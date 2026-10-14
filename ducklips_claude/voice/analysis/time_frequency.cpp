#include "time_frequency.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace tf {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

std::uint32_t read_u32(const unsigned char* bytes) {
  return static_cast<std::uint32_t>(bytes[0]) | (static_cast<std::uint32_t>(bytes[1]) << 8U) |
         (static_cast<std::uint32_t>(bytes[2]) << 16U) | (static_cast<std::uint32_t>(bytes[3]) << 24U);
}

std::uint16_t read_u16(const unsigned char* bytes) {
  return static_cast<std::uint16_t>(bytes[0] | (bytes[1] << 8U));
}

// Fine log-frequency row of a frequency, or -1 outside the grid.
std::int32_t fine_row(const Grid& grid, double hz, std::int32_t rows) {
  if (hz <= 0.0) {
    return -1;
  }
  double position = std::log2(hz / grid.low_hz) * grid.fine_per_octave;
  std::int32_t row = static_cast<std::int32_t>(std::lround(position));
  return row >= 0 && row < rows ? row : -1;
}

std::int32_t fine_row_count(const Grid& grid) {
  return static_cast<std::int32_t>(std::floor(std::log2(grid.high_hz / grid.low_hz) * grid.fine_per_octave)) + 1;
}

std::vector<float> fine_hz(const Grid& grid) {
  std::int32_t rows = fine_row_count(grid);
  std::vector<float> hz(static_cast<std::size_t>(rows), 0.0F);
  for (std::int32_t row = 0; row < rows; ++row) {
    hz[static_cast<std::size_t>(row)] =
        static_cast<float>(grid.low_hz * std::pow(2.0, static_cast<double>(row) / grid.fine_per_octave));
  }
  return hz;
}

Field make_field(std::int32_t rows, std::int32_t columns, const std::vector<float>& hz,
                 const std::vector<float>& seconds) {
  Field field;
  field.rows = rows;
  field.columns = columns;
  field.data.assign(static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns), 0.0F);
  field.hz = hz;
  field.seconds = seconds;
  return field;
}

// Places a band field (rows at arbitrary centre frequencies) onto the fine
// log grid by linear interpolation in log frequency.
Field band_to_fine(const Field& bands, const Grid& grid) {
  std::vector<float> hz = fine_hz(grid);
  std::int32_t rows = static_cast<std::int32_t>(hz.size());
  Field fine = make_field(rows, bands.columns, hz, bands.seconds);
  std::int32_t lower = 0;
  for (std::int32_t row = 0; row < rows; ++row) {
    double target = std::log2(static_cast<double>(hz[static_cast<std::size_t>(row)]));
    while (lower + 1 < bands.rows - 1 &&
           std::log2(static_cast<double>(bands.hz[static_cast<std::size_t>(lower + 1)])) < target) {
      ++lower;
    }
    double low = std::log2(static_cast<double>(bands.hz[static_cast<std::size_t>(lower)]));
    double high = std::log2(static_cast<double>(bands.hz[static_cast<std::size_t>(lower + 1)]));
    double weight = std::clamp((target - low) / (high - low), 0.0, 1.0);
    for (std::int32_t column = 0; column < bands.columns; ++column) {
      std::size_t below = static_cast<std::size_t>(lower) * static_cast<std::size_t>(bands.columns) +
                          static_cast<std::size_t>(column);
      std::size_t above = below + static_cast<std::size_t>(bands.columns);
      fine.data[static_cast<std::size_t>(row) * static_cast<std::size_t>(bands.columns) +
                static_cast<std::size_t>(column)] =
          static_cast<float>((1.0 - weight) * bands.data[below] + weight * bands.data[above]);
    }
  }
  return fine;
}

double quantile(std::vector<float> values, double fraction) {
  if (values.empty()) {
    return 0.0;
  }
  std::size_t index = static_cast<std::size_t>(fraction * static_cast<double>(values.size() - 1));
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
  return values[index];
}

// The enhanced-STFT recipe carried over: cartoon = 3x3 average of the fused
// plain field; texture = fused reassigned amplitude; scaled to the same
// 99.5th percentile and added.
Field cartoon_plus_texture(const Field& fused_plain, const Field& fused_squeezed) {
  Field result = fused_plain;
  std::int32_t rows = fused_plain.rows;
  std::int32_t columns = fused_plain.columns;
  for (std::int32_t row = 0; row < rows; ++row) {
    for (std::int32_t column = 0; column < columns; ++column) {
      double sum = 0.0;
      std::int32_t count = 0;
      for (std::int32_t dr = -1; dr <= 1; ++dr) {
        for (std::int32_t dc = -1; dc <= 1; ++dc) {
          std::int32_t r = row + dr;
          std::int32_t c = column + dc;
          if (r >= 0 && r < rows && c >= 0 && c < columns) {
            sum += fused_plain.data[static_cast<std::size_t>(r) * static_cast<std::size_t>(columns) +
                                    static_cast<std::size_t>(c)];
            ++count;
          }
        }
      }
      result.data[static_cast<std::size_t>(row) * static_cast<std::size_t>(columns) +
                  static_cast<std::size_t>(column)] = static_cast<float>(sum / count);
    }
  }
  double scale = quantile(result.data, 0.995) / std::max(quantile(fused_squeezed.data, 0.995), 1.0e-30);
  for (std::size_t index = 0; index < result.data.size(); ++index) {
    result.data[index] += static_cast<float>(scale * fused_squeezed.data[index]);
  }
  return result;
}

void maximum_into(Field& target, const Field& source) {
  for (std::size_t index = 0; index < target.data.size(); ++index) {
    target.data[index] = std::max(target.data[index], source.data[index]);
  }
}

// Synchrosqueezed power accumulator on the fine grid, converted to amplitude.
void finish_squeeze(Field& squeezed) {
  for (float& value : squeezed.data) {
    value = std::sqrt(std::max(value, 0.0F));
  }
}

// Linear interpolation of a complex sequence sampled every `spacing` samples.
Complex sample_at(const std::vector<Complex>& values, double position) {
  if (position <= 0.0) {
    return values.front();
  }
  double last = static_cast<double>(values.size() - 1);
  if (position >= last) {
    return values.back();
  }
  std::size_t index = static_cast<std::size_t>(position);
  double weight = position - static_cast<double>(index);
  return (1.0 - weight) * values[index] + weight * values[index + 1];
}

// Constant-Q order N quality and atom length (Garces 2020, eqs. 11 and A37).
double order_quality(double order) {
  return 1.0 / (std::pow(2.0, 1.0 / (2.0 * order)) - std::pow(2.0, -1.0 / (2.0 * order)));
}

double atom_cycles(double order) {
  return 2.0 * std::sqrt(std::log(2.0)) * order_quality(order);
}

}  // namespace

bool read_wave(const std::string& path, Signal& signal) {
  std::FILE* file = std::fopen(path.c_str(), "rb");
  if (file == nullptr) {
    return false;
  }
  std::vector<unsigned char> bytes;
  unsigned char buffer[65536];
  std::size_t count = 0;
  while ((count = std::fread(buffer, 1, sizeof(buffer), file)) > 0) {
    bytes.insert(bytes.end(), buffer, buffer + count);
  }
  std::fclose(file);
  if (bytes.size() < 44 || std::memcmp(bytes.data(), "RIFF", 4) != 0 || std::memcmp(bytes.data() + 8, "WAVE", 4) != 0) {
    return false;
  }
  std::uint16_t format = 0;
  std::uint16_t channels = 0;
  std::uint16_t bits = 0;
  std::uint32_t rate = 0;
  std::size_t position = 12;
  while (position + 8 <= bytes.size()) {
    std::uint32_t chunk = read_u32(bytes.data() + position + 4);
    if (std::memcmp(bytes.data() + position, "fmt ", 4) == 0) {
      format = read_u16(bytes.data() + position + 8);
      channels = read_u16(bytes.data() + position + 10);
      rate = read_u32(bytes.data() + position + 12);
      bits = read_u16(bytes.data() + position + 22);
    } else if (std::memcmp(bytes.data() + position, "data", 4) == 0) {
      std::size_t end = std::min(bytes.size(), position + 8 + chunk);
      std::size_t frame = static_cast<std::size_t>(channels) * (bits / 8U);
      if (frame == 0) {
        return false;
      }
      signal.samples.clear();
      for (std::size_t offset = position + 8; offset + frame <= end; offset += frame) {
        double sum = 0.0;
        for (std::uint16_t channel = 0; channel < channels; ++channel) {
          const unsigned char* sample = bytes.data() + offset + channel * (bits / 8U);
          if (format == 1 && bits == 16) {
            sum += static_cast<std::int16_t>(read_u16(sample)) / 32768.0;
          } else if (format == 3 && bits == 32) {
            float value = 0.0F;
            std::memcpy(&value, sample, 4);
            sum += value;
          } else if (format == 1 && bits == 32) {
            sum += static_cast<std::int32_t>(read_u32(sample)) / 2147483648.0;
          } else {
            return false;
          }
        }
        signal.samples.push_back(sum / channels);
      }
      signal.rate = rate;
      return true;
    }
    position += 8 + chunk + (chunk & 1U);
  }
  return false;
}

bool write_npy(const std::string& path, const std::vector<float>& data, std::int32_t rows, std::int32_t columns) {
  std::FILE* file = std::fopen(path.c_str(), "wb");
  if (file == nullptr) {
    return false;
  }
  char header[128];
  int length = std::snprintf(header, sizeof(header), "{'descr': '<f4', 'fortran_order': False, 'shape': (%d, %d), }",
                             rows, columns);
  int total = 10 + length + 1;
  int padding = (64 - total % 64) % 64;
  std::string text(header, static_cast<std::size_t>(length));
  text.append(static_cast<std::size_t>(padding), ' ');
  text.push_back('\n');
  std::uint16_t header_length = static_cast<std::uint16_t>(text.size());
  std::fwrite("\x93NUMPY\x01\x00", 1, 8, file);
  unsigned char size_bytes[2] = {static_cast<unsigned char>(header_length & 0xFFU),
                                 static_cast<unsigned char>(header_length >> 8U)};
  std::fwrite(size_bytes, 1, 2, file);
  std::fwrite(text.data(), 1, text.size(), file);
  std::fwrite(data.data(), sizeof(float), data.size(), file);
  return std::fclose(file) == 0;
}

std::int32_t next_power_of_two(std::int32_t value) {
  std::int32_t power = 1;
  while (power < value) {
    power <<= 1;
  }
  return power;
}

Fft::Fft(std::int32_t size) : size_(size) {
  twiddles_.resize(static_cast<std::size_t>(size / 2));
  for (std::int32_t index = 0; index < size / 2; ++index) {
    twiddles_[static_cast<std::size_t>(index)] = std::polar(1.0, -kTwoPi * index / size);
  }
  reversed_.resize(static_cast<std::size_t>(size));
  std::int32_t bits = 0;
  while ((1 << bits) < size) {
    ++bits;
  }
  for (std::int32_t index = 0; index < size; ++index) {
    std::int32_t reversed = 0;
    for (std::int32_t bit = 0; bit < bits; ++bit) {
      if ((index >> bit) & 1) {
        reversed |= 1 << (bits - 1 - bit);
      }
    }
    reversed_[static_cast<std::size_t>(index)] = reversed;
  }
}

void Fft::forward(std::vector<Complex>& values) const { run(values, false); }

void Fft::inverse(std::vector<Complex>& values) const { run(values, true); }

void Fft::run(std::vector<Complex>& values, bool inverse) const {
  for (std::int32_t index = 0; index < size_; ++index) {
    std::int32_t other = reversed_[static_cast<std::size_t>(index)];
    if (other > index) {
      std::swap(values[static_cast<std::size_t>(index)], values[static_cast<std::size_t>(other)]);
    }
  }
  for (std::int32_t length = 2; length <= size_; length <<= 1) {
    std::int32_t half = length / 2;
    std::int32_t stride = size_ / length;
    for (std::int32_t start = 0; start < size_; start += length) {
      for (std::int32_t offset = 0; offset < half; ++offset) {
        Complex twiddle = twiddles_[static_cast<std::size_t>(offset * stride)];
        if (inverse) {
          twiddle = std::conj(twiddle);
        }
        Complex& even = values[static_cast<std::size_t>(start + offset)];
        Complex& odd = values[static_cast<std::size_t>(start + offset + half)];
        Complex rotated = odd * twiddle;
        odd = even - rotated;
        even += rotated;
      }
    }
  }
}

std::vector<Complex> analytic_signal(const Signal& signal) {
  std::int32_t size = next_power_of_two(static_cast<std::int32_t>(signal.samples.size()));
  std::vector<Complex> values(static_cast<std::size_t>(size), Complex(0.0, 0.0));
  for (std::size_t index = 0; index < signal.samples.size(); ++index) {
    values[index] = Complex(signal.samples[index], 0.0);
  }
  Fft fft(size);
  fft.forward(values);
  for (std::int32_t index = 1; index < size / 2; ++index) {
    values[static_cast<std::size_t>(index)] *= 2.0;
  }
  for (std::int32_t index = size / 2 + 1; index < size; ++index) {
    values[static_cast<std::size_t>(index)] = Complex(0.0, 0.0);
  }
  fft.inverse(values);
  for (Complex& value : values) {
    value /= static_cast<double>(size);
  }
  values.resize(signal.samples.size());
  return values;
}

std::vector<double> band_centres(const Grid& grid) {
  std::vector<double> centres;
  for (std::int32_t index = 0;; ++index) {
    double hz = grid.low_hz * std::pow(2.0, static_cast<double>(index) / grid.bands_per_octave);
    if (hz > grid.high_hz) {
      break;
    }
    centres.push_back(hz);
  }
  return centres;
}

std::vector<float> column_times(const Signal& signal, const Grid& grid) {
  double duration = static_cast<double>(signal.samples.size()) / signal.rate;
  std::vector<float> seconds;
  for (double time = 0.0; time < duration; time += grid.hop_seconds) {
    seconds.push_back(static_cast<float>(time));
  }
  return seconds;
}

// ---------------------------------------------------------------------------
// Short-time Fourier transform, enhanced by phase-advance reassignment and
// fusion across three windows.

Analysis stft(const Signal& signal, const Grid& grid, double window_seconds) {
  std::vector<float> seconds = column_times(signal, grid);
  std::int32_t columns = static_cast<std::int32_t>(seconds.size());
  std::vector<double> centres = band_centres(grid);
  std::int32_t bands = static_cast<std::int32_t>(centres.size());
  std::vector<float> band_hz(centres.begin(), centres.end());
  Field fused_plain;
  Field fused_squeezed;
  Analysis result;
  const double scales[3] = {0.5, 1.0, 2.0};
  std::int32_t length = static_cast<std::int32_t>(signal.samples.size());
  for (std::int32_t variant = 0; variant < 3; ++variant) {
    std::int32_t window = static_cast<std::int32_t>(window_seconds * scales[variant] * signal.rate);
    std::int32_t size = next_power_of_two(window);
    Fft fft(size);
    std::vector<double> taper(static_cast<std::size_t>(window), 0.0);
    double taper_sum = 0.0;
    for (std::int32_t index = 0; index < window; ++index) {
      taper[static_cast<std::size_t>(index)] = 0.5 - 0.5 * std::cos(kTwoPi * (index + 0.5) / window);
      taper_sum += taper[static_cast<std::size_t>(index)];
    }
    Field plain = make_field(bands, columns, band_hz, seconds);
    Field squeezed = make_field(fine_row_count(grid), columns, fine_hz(grid), seconds);
    std::vector<Complex> now(static_cast<std::size_t>(size));
    std::vector<Complex> later(static_cast<std::size_t>(size));
    double bin_hz = signal.rate / size;
    for (std::int32_t column = 0; column < columns; ++column) {
      std::int32_t centre = static_cast<std::int32_t>(std::lround(seconds[static_cast<std::size_t>(column)] * signal.rate));
      std::fill(now.begin(), now.end(), Complex(0.0, 0.0));
      std::fill(later.begin(), later.end(), Complex(0.0, 0.0));
      for (std::int32_t index = 0; index < window; ++index) {
        std::int32_t sample = centre - window / 2 + index;
        double value = sample >= 0 && sample < length ? signal.samples[static_cast<std::size_t>(sample)] : 0.0;
        double next = sample + 1 >= 0 && sample + 1 < length ? signal.samples[static_cast<std::size_t>(sample + 1)] : 0.0;
        now[static_cast<std::size_t>(index)] = Complex(value * taper[static_cast<std::size_t>(index)], 0.0);
        later[static_cast<std::size_t>(index)] = Complex(next * taper[static_cast<std::size_t>(index)], 0.0);
      }
      fft.forward(now);
      fft.forward(later);
      double normal = 2.0 / taper_sum;
      for (std::int32_t band = 0; band < bands; ++band) {
        double position = centres[static_cast<std::size_t>(band)] / bin_hz;
        std::int32_t bin = static_cast<std::int32_t>(position);
        double weight = position - bin;
        double magnitude = (1.0 - weight) * std::abs(now[static_cast<std::size_t>(bin)]) +
                           weight * std::abs(now[static_cast<std::size_t>(bin + 1)]);
        plain.data[static_cast<std::size_t>(band) * static_cast<std::size_t>(columns) +
                   static_cast<std::size_t>(column)] = static_cast<float>(magnitude * normal);
      }
      std::int32_t top = std::min(size / 2 - 1, static_cast<std::int32_t>(grid.high_hz * 1.05 / bin_hz));
      for (std::int32_t bin = 1; bin <= top; ++bin) {
        Complex a = now[static_cast<std::size_t>(bin)];
        double power = std::norm(a) * normal * normal;
        if (power <= 0.0) {
          continue;
        }
        double advance = std::arg(later[static_cast<std::size_t>(bin)] * std::conj(a));
        double hz = advance * signal.rate / kTwoPi;
        std::int32_t row = fine_row(grid, hz, squeezed.rows);
        if (row >= 0) {
          squeezed.data[static_cast<std::size_t>(row) * static_cast<std::size_t>(columns) +
                        static_cast<std::size_t>(column)] += static_cast<float>(power);
        }
      }
    }
    finish_squeeze(squeezed);
    Field plain_fine = band_to_fine(plain, grid);
    if (variant == 0) {
      fused_plain = plain_fine;
      fused_squeezed = squeezed;
    } else {
      maximum_into(fused_plain, plain_fine);
      maximum_into(fused_squeezed, squeezed);
    }
    if (variant == 1) {
      result.plain = plain;
    }
  }
  result.enhanced = cartoon_plus_texture(fused_plain, fused_squeezed);
  return result;
}

// ---------------------------------------------------------------------------
// Constant-Q Gabor atoms (Garces 2020). Each band is computed by a small
// inverse FFT of only the spectrum under its Gaussian, shifted to baseband
// and decimated, so the cost per band follows its bandwidth.

Analysis constant_q(const Signal& signal, const Grid& grid) {
  std::vector<float> seconds = column_times(signal, grid);
  std::int32_t columns = static_cast<std::int32_t>(seconds.size());
  std::vector<double> centres = band_centres(grid);
  std::int32_t bands = static_cast<std::int32_t>(centres.size());
  std::vector<float> band_hz(centres.begin(), centres.end());
  std::int32_t length = static_cast<std::int32_t>(signal.samples.size());
  std::int32_t size = next_power_of_two(length);
  std::vector<Complex> spectrum(static_cast<std::size_t>(size), Complex(0.0, 0.0));
  for (std::int32_t index = 0; index < length; ++index) {
    spectrum[static_cast<std::size_t>(index)] = Complex(signal.samples[static_cast<std::size_t>(index)], 0.0);
  }
  Fft full(size);
  full.forward(spectrum);
  double bin_hz = signal.rate / size;
  double hop_samples = grid.hop_seconds * signal.rate;

  const double orders[3] = {3.0, 6.0, 12.0};
  Field fused_plain;
  Field fused_squeezed;
  Analysis result;
  std::vector<Complex> baseband;
  for (std::int32_t variant = 0; variant < 3; ++variant) {
    double cycles = atom_cycles(orders[variant]);
    Field plain = make_field(bands, columns, band_hz, seconds);
    Field squeezed = make_field(fine_row_count(grid), columns, fine_hz(grid), seconds);
    for (std::int32_t band = 0; band < bands; ++band) {
      double centre = centres[static_cast<std::size_t>(band)];
      double spread = centre / cycles;  // Gaussian standard deviation in Hz
      std::int32_t centre_bin = static_cast<std::int32_t>(std::lround(centre / bin_hz));
      std::int32_t half_span = static_cast<std::int32_t>(std::ceil(4.0 * spread / bin_hz)) + 1;
      std::int32_t minimum = next_power_of_two(static_cast<std::int32_t>(std::ceil(2.0 * size / hop_samples)));
      std::int32_t small = std::min(size, std::max(next_power_of_two(4 * half_span), minimum));
      double decimation = static_cast<double>(size) / small;
      baseband.assign(static_cast<std::size_t>(small), Complex(0.0, 0.0));
      for (std::int32_t offset = -half_span; offset <= half_span; ++offset) {
        std::int32_t bin = centre_bin + offset;
        if (bin <= 0 || bin >= size / 2) {
          continue;
        }
        double ratio = bin * bin_hz / centre - 1.0;
        double gain = std::exp(-0.5 * cycles * cycles * ratio * ratio);
        std::int32_t slot = ((offset % small) + small) % small;
        baseband[static_cast<std::size_t>(slot)] += 2.0 * gain * spectrum[static_cast<std::size_t>(bin)];
      }
      Fft small_fft(small);
      small_fft.inverse(baseband);
      for (Complex& value : baseband) {
        value /= static_cast<double>(size);
      }
      for (std::int32_t column = 0; column < columns; ++column) {
        double sample = seconds[static_cast<std::size_t>(column)] * signal.rate;
        double position = sample / decimation;
        Complex value = sample_at(baseband, position);
        Complex next = sample_at(baseband, position + 1.0);
        double magnitude = std::abs(value);
        plain.data[static_cast<std::size_t>(band) * static_cast<std::size_t>(columns) +
                   static_cast<std::size_t>(column)] = static_cast<float>(magnitude);
        double advance = std::arg(next * std::conj(value));
        double hz = centre_bin * bin_hz + advance * signal.rate / (kTwoPi * decimation);
        std::int32_t row = fine_row(grid, hz, squeezed.rows);
        if (row >= 0) {
          squeezed.data[static_cast<std::size_t>(row) * static_cast<std::size_t>(columns) +
                        static_cast<std::size_t>(column)] += static_cast<float>(magnitude * magnitude);
        }
      }
    }
    finish_squeeze(squeezed);
    Field plain_fine = band_to_fine(plain, grid);
    if (variant == 0) {
      fused_plain = plain_fine;
      fused_squeezed = squeezed;
    } else {
      maximum_into(fused_plain, plain_fine);
      maximum_into(fused_squeezed, squeezed);
    }
    if (variant == 1) {
      result.plain = plain;
    }
  }
  result.enhanced = cartoon_plus_texture(fused_plain, fused_squeezed);
  return result;
}

// ---------------------------------------------------------------------------
// Time-causal Gabor analogue (Lindeberg 2025): demodulate by each band's
// frequency, then smooth with the time-causal limit kernel realised as K
// first-order recursive filters in cascade (c = sqrt 2). The output is read
// back by the kernel's mean delay, and reassigned by its phase advance.

Analysis causal_gabor(const Signal& signal, const Grid& grid) {
  std::vector<float> seconds = column_times(signal, grid);
  std::int32_t columns = static_cast<std::int32_t>(seconds.size());
  std::vector<double> centres = band_centres(grid);
  std::int32_t bands = static_cast<std::int32_t>(centres.size());
  std::vector<float> band_hz(centres.begin(), centres.end());
  std::vector<Complex> analytic = analytic_signal(signal);
  std::int32_t length = static_cast<std::int32_t>(analytic.size());
  constexpr std::int32_t kStages = 8;
  const double kRatio = std::sqrt(2.0);

  const double orders[3] = {3.0, 6.0, 12.0};
  Field fused_plain;
  Field fused_squeezed;
  Analysis result;
  std::vector<Complex> smoothed(static_cast<std::size_t>(length));
  for (std::int32_t variant = 0; variant < 3; ++variant) {
    double cycles = atom_cycles(orders[variant]);
    Field plain = make_field(bands, columns, band_hz, seconds);
    Field squeezed = make_field(fine_row_count(grid), columns, fine_hz(grid), seconds);
    for (std::int32_t band = 0; band < bands; ++band) {
      double centre = centres[static_cast<std::size_t>(band)];
      // Same temporal spread as the Gabor atom of this order: sigma (samples).
      double sigma = cycles * signal.rate / (kTwoPi * centre);
      double variance = sigma * sigma;
      double time_constant[kStages];
      double previous = 0.0;
      double delay = 0.0;
      for (std::int32_t stage = 1; stage <= kStages; ++stage) {
        double level = std::pow(kRatio, 2.0 * (stage - kStages)) * variance;
        double increment = level - previous;
        previous = level;
        time_constant[stage - 1] = (std::sqrt(1.0 + 4.0 * increment) - 1.0) / 2.0;
        delay += time_constant[stage - 1];
      }
      double omega = kTwoPi * centre / signal.rate;
      Complex rotation = std::polar(1.0, -omega);
      Complex phase(1.0, 0.0);
      Complex states[kStages];
      for (std::int32_t stage = 0; stage < kStages; ++stage) {
        states[stage] = Complex(0.0, 0.0);
      }
      for (std::int32_t index = 0; index < length; ++index) {
        Complex value = analytic[static_cast<std::size_t>(index)] * phase;
        phase *= rotation;
        if ((index & 1023) == 0) {
          phase /= std::abs(phase);
        }
        for (std::int32_t stage = 0; stage < kStages; ++stage) {
          states[stage] += (value - states[stage]) / (1.0 + time_constant[stage]);
          value = states[stage];
        }
        smoothed[static_cast<std::size_t>(index)] = value;
      }
      for (std::int32_t column = 0; column < columns; ++column) {
        double position = seconds[static_cast<std::size_t>(column)] * signal.rate + delay;
        Complex value = sample_at(smoothed, position);
        Complex next = sample_at(smoothed, position + 1.0);
        double magnitude = std::abs(value);
        plain.data[static_cast<std::size_t>(band) * static_cast<std::size_t>(columns) +
                   static_cast<std::size_t>(column)] = static_cast<float>(magnitude);
        double hz = centre + std::arg(next * std::conj(value)) * signal.rate / kTwoPi;
        std::int32_t row = fine_row(grid, hz, squeezed.rows);
        if (row >= 0) {
          squeezed.data[static_cast<std::size_t>(row) * static_cast<std::size_t>(columns) +
                        static_cast<std::size_t>(column)] += static_cast<float>(magnitude * magnitude);
        }
      }
    }
    finish_squeeze(squeezed);
    Field plain_fine = band_to_fine(plain, grid);
    if (variant == 0) {
      fused_plain = plain_fine;
      fused_squeezed = squeezed;
    } else {
      maximum_into(fused_plain, plain_fine);
      maximum_into(fused_squeezed, squeezed);
    }
    if (variant == 1) {
      result.plain = plain;
    }
  }
  result.enhanced = cartoon_plus_texture(fused_plain, fused_squeezed);
  return result;
}

// ---------------------------------------------------------------------------
// Smoothed pseudo-Wigner-Ville: lag window h, triangular time smoothing g.
// Enhanced: the minimum of the positive parts over three (h, g) shapes.
// Auto-terms sit in the same place in every shape; interference terms move
// and change sign, so the consensus keeps components and drops cross-terms.

Analysis pseudo_wigner(const Signal& signal, const Grid& grid) {
  std::vector<float> seconds = column_times(signal, grid);
  std::int32_t columns = static_cast<std::int32_t>(seconds.size());
  std::vector<double> centres = band_centres(grid);
  std::int32_t bands = static_cast<std::int32_t>(centres.size());
  std::vector<float> band_hz(centres.begin(), centres.end());
  std::vector<Complex> analytic = analytic_signal(signal);
  std::int32_t length = static_cast<std::int32_t>(analytic.size());
  const double lag_seconds[3] = {0.012, 0.020, 0.030};
  const double smooth_seconds[3] = {0.002, 0.004, 0.006};
  Field consensus;
  Analysis result;
  for (std::int32_t variant = 0; variant < 3; ++variant) {
    std::int32_t lag_half = static_cast<std::int32_t>(lag_seconds[variant] * signal.rate / 2.0);
    std::int32_t smooth_half = std::max(1, static_cast<std::int32_t>(smooth_seconds[variant] * signal.rate / 2.0));
    std::int32_t size = next_power_of_two(2 * lag_half + 2);
    Fft fft(size);
    std::vector<double> lag_window(static_cast<std::size_t>(lag_half + 1));
    double lag_sum = 0.0;
    for (std::int32_t lag = 0; lag <= lag_half; ++lag) {
      double value = 0.5 + 0.5 * std::cos(kPi * lag / (lag_half + 1));
      lag_window[static_cast<std::size_t>(lag)] = value;
      lag_sum += lag == 0 ? value : 2.0 * value;
    }
    std::vector<double> smoothing(static_cast<std::size_t>(2 * smooth_half + 1));
    double smoothing_sum = 0.0;
    for (std::int32_t offset = -smooth_half; offset <= smooth_half; ++offset) {
      double value = 1.0 - std::abs(static_cast<double>(offset)) / (smooth_half + 1);
      smoothing[static_cast<std::size_t>(offset + smooth_half)] = value;
      smoothing_sum += value;
    }
    Field plain = make_field(bands, columns, band_hz, seconds);
    std::vector<Complex> kernel(static_cast<std::size_t>(size));
    double bin_hz = signal.rate / (2.0 * size);
    for (std::int32_t column = 0; column < columns; ++column) {
      std::int32_t centre = static_cast<std::int32_t>(std::lround(seconds[static_cast<std::size_t>(column)] * signal.rate));
      std::fill(kernel.begin(), kernel.end(), Complex(0.0, 0.0));
      for (std::int32_t lag = 0; lag <= lag_half; ++lag) {
        Complex sum(0.0, 0.0);
        for (std::int32_t offset = -smooth_half; offset <= smooth_half; ++offset) {
          std::int32_t ahead = centre + offset + lag;
          std::int32_t behind = centre + offset - lag;
          if (ahead < 0 || ahead >= length || behind < 0 || behind >= length) {
            continue;
          }
          sum += smoothing[static_cast<std::size_t>(offset + smooth_half)] *
                 analytic[static_cast<std::size_t>(ahead)] * std::conj(analytic[static_cast<std::size_t>(behind)]);
        }
        sum *= lag_window[static_cast<std::size_t>(lag)] / smoothing_sum;
        kernel[static_cast<std::size_t>(lag)] = sum;
        if (lag > 0) {
          kernel[static_cast<std::size_t>(size - lag)] = std::conj(sum);
        }
      }
      fft.forward(kernel);
      for (std::int32_t band = 0; band < bands; ++band) {
        double position = centres[static_cast<std::size_t>(band)] / bin_hz;
        std::int32_t bin = static_cast<std::int32_t>(position);
        double weight = position - bin;
        double power = (1.0 - weight) * kernel[static_cast<std::size_t>(bin)].real() +
                       weight * kernel[static_cast<std::size_t>(bin + 1)].real();
        plain.data[static_cast<std::size_t>(band) * static_cast<std::size_t>(columns) +
                   static_cast<std::size_t>(column)] =
            static_cast<float>(std::sqrt(std::max(power, 0.0) / lag_sum));
      }
    }
    Field plain_fine = band_to_fine(plain, grid);
    if (variant == 0) {
      consensus = plain_fine;
    } else {
      for (std::size_t index = 0; index < consensus.data.size(); ++index) {
        consensus.data[index] = std::min(consensus.data[index], plain_fine.data[index]);
      }
    }
    if (variant == 1) {
      result.plain = plain;
    }
  }
  result.enhanced = consensus;
  return result;
}

}  // namespace tf
