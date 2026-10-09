// Measurement tool for DuckLips output and control recordings.
//   spectrogram:  log-magnitude STFT as a PNG (no external libraries)
//   pitch:        autocorrelation F0 over voiced frames
//   formants:     LPC envelope peaks on a 12 kHz copy, median over voiced frames
// Works on any PCM16/24/float WAV; channels are averaged.
#include "ducklips/wav.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

struct Arguments final {
    std::string input{};
    std::string spectrogram{};
    double from{0.0};
    double to{-1.0};
    double maximum_hz{8000.0};
    std::size_t lpc_order{12};
    bool clicks{false};
    bool track{false};
};

// ---------------------------------------------------------------- PNG output

std::uint32_t crc32(const unsigned char* data, std::size_t size, std::uint32_t crc) noexcept {
    crc = ~crc;
    for (std::size_t index = 0; index < size; ++index) {
        crc ^= data[index];
        for (std::size_t bit = 0; bit < 8; ++bit) {
            std::uint32_t mask = (crc & 1u) != 0 ? 0xEDB88320u : 0u;
            crc = (crc >> 1) ^ mask;
        }
    }
    std::uint32_t result = ~crc;
    return result;
}

void put_be32(std::vector<unsigned char>& bytes, std::uint32_t value) {
    bytes.push_back(static_cast<unsigned char>(value >> 24));
    bytes.push_back(static_cast<unsigned char>(value >> 16));
    bytes.push_back(static_cast<unsigned char>(value >> 8));
    bytes.push_back(static_cast<unsigned char>(value));
}

void put_chunk(std::vector<unsigned char>& file, const char* type, const std::vector<unsigned char>& data) {
    put_be32(file, static_cast<std::uint32_t>(data.size()));
    std::vector<unsigned char> body(type, type + 4);
    body.insert(body.end(), data.begin(), data.end());
    file.insert(file.end(), body.begin(), body.end());
    put_be32(file, crc32(body.data(), body.size(), 0));
}

// RGB rows, stored (uncompressed) deflate blocks inside a zlib stream.
bool write_png(const std::string& path, std::size_t width, std::size_t height, const std::vector<unsigned char>& rgb) {
    std::vector<unsigned char> raw{};
    raw.reserve(height * (1 + 3 * width));
    for (std::size_t row = 0; row < height; ++row) {
        raw.push_back(0); // filter: none
        raw.insert(raw.end(), rgb.begin() + static_cast<std::ptrdiff_t>(row * 3 * width),
                   rgb.begin() + static_cast<std::ptrdiff_t>((row + 1) * 3 * width));
    }
    std::vector<unsigned char> zlib{0x78, 0x01};
    std::size_t offset = 0;
    while (offset < raw.size() || raw.empty()) {
        std::size_t block = std::min<std::size_t>(65535, raw.size() - offset);
        bool last = offset + block >= raw.size();
        zlib.push_back(last ? 1 : 0);
        zlib.push_back(static_cast<unsigned char>(block & 0xFF));
        zlib.push_back(static_cast<unsigned char>(block >> 8));
        zlib.push_back(static_cast<unsigned char>(~block & 0xFF));
        zlib.push_back(static_cast<unsigned char>((~block >> 8) & 0xFF));
        zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(offset),
                    raw.begin() + static_cast<std::ptrdiff_t>(offset + block));
        offset += block;
        if (raw.empty()) {
            break;
        }
    }
    std::uint32_t a = 1;
    std::uint32_t b = 0;
    for (std::size_t index = 0; index < raw.size(); ++index) {
        a = (a + raw[index]) % 65521u;
        b = (b + a) % 65521u;
    }
    put_be32(zlib, (b << 16) | a);
    std::vector<unsigned char> file{0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
    std::vector<unsigned char> header{};
    put_be32(header, static_cast<std::uint32_t>(width));
    put_be32(header, static_cast<std::uint32_t>(height));
    header.push_back(8);
    header.push_back(2);
    header.push_back(0);
    header.push_back(0);
    header.push_back(0);
    put_chunk(file, "IHDR", header);
    put_chunk(file, "IDAT", zlib);
    put_chunk(file, "IEND", std::vector<unsigned char>{});
    std::ofstream stream(path, std::ios::binary);
    stream.write(reinterpret_cast<const char*>(file.data()), static_cast<std::streamsize>(file.size()));
    bool written = static_cast<bool>(stream);
    return written;
}

// ---------------------------------------------------------------- FFT

void fft(std::vector<std::complex<double>>& data) noexcept {
    const std::size_t size = data.size();
    for (std::size_t index = 1, reversed = 0; index < size; ++index) {
        std::size_t bit = size >> 1;
        for (; (reversed & bit) != 0; bit >>= 1) {
            reversed ^= bit;
        }
        reversed ^= bit;
        if (index < reversed) {
            std::swap(data[index], data[reversed]);
        }
    }
    for (std::size_t length = 2; length <= size; length <<= 1) {
        double angle = -2.0 * pi / static_cast<double>(length);
        std::complex<double> unit(std::cos(angle), std::sin(angle));
        for (std::size_t start = 0; start < size; start += length) {
            std::complex<double> twiddle(1.0, 0.0);
            for (std::size_t offset = 0; offset < length / 2; ++offset) {
                std::complex<double> even = data[start + offset];
                std::complex<double> odd = data[start + offset + length / 2] * twiddle;
                data[start + offset] = even + odd;
                data[start + offset + length / 2] = even - odd;
                twiddle *= unit;
            }
        }
    }
}

// Colour map from dark blue through red to pale yellow; position in [0, 1].
std::array<unsigned char, 3> heat(double position) noexcept {
    double t = std::clamp(position, 0.0, 1.0);
    double red = std::clamp(1.6 * t - 0.1, 0.0, 1.0);
    double green = std::clamp(2.0 * t - 1.0, 0.0, 1.0);
    double blue = std::clamp(0.35 + 0.9 * t - 1.6 * t * t + 1.0 * std::max(t - 0.8, 0.0) * 4.0, 0.0, 1.0);
    std::array<unsigned char, 3> colour{static_cast<unsigned char>(255.0 * red), static_cast<unsigned char>(255.0 * green),
                                        static_cast<unsigned char>(255.0 * blue)};
    return colour;
}

bool render_spectrogram(const std::vector<double>& signal, double rate, double maximum_hz, const std::string& path) {
    const std::size_t window = 1024;
    const std::size_t height = 400;
    std::size_t hop = std::max<std::size_t>(64, signal.size() / 1400);
    if (signal.size() < window) {
        return false;
    }
    std::size_t frames = (signal.size() - window) / hop + 1;
    std::vector<double> magnitude(frames * height, 0.0);
    std::vector<std::complex<double>> buffer(window);
    double loudest = -300.0;
    for (std::size_t frame = 0; frame < frames; ++frame) {
        for (std::size_t index = 0; index < window; ++index) {
            double hann = 0.5 - 0.5 * std::cos(2.0 * pi * static_cast<double>(index) / static_cast<double>(window));
            buffer[index] = std::complex<double>(signal[frame * hop + index] * hann, 0.0);
        }
        fft(buffer);
        for (std::size_t row = 0; row < height; ++row) {
            double frequency = maximum_hz * (static_cast<double>(row) + 0.5) / static_cast<double>(height);
            std::size_t bin = static_cast<std::size_t>(frequency / rate * static_cast<double>(window));
            double level = 10.0 * std::log10(std::norm(buffer[bin]) + 1.0e-20);
            magnitude[frame * height + row] = level;
            loudest = std::max(loudest, level);
        }
    }
    std::vector<unsigned char> rgb(frames * height * 3, 0);
    for (std::size_t frame = 0; frame < frames; ++frame) {
        for (std::size_t row = 0; row < height; ++row) {
            double level = magnitude[frame * height + row];
            double position = (level - (loudest - 75.0)) / 75.0;
            std::array<unsigned char, 3> colour = heat(position);
            std::size_t y = height - 1 - row;
            std::size_t pixel = (y * frames + frame) * 3;
            rgb[pixel] = colour[0];
            rgb[pixel + 1] = colour[1];
            rgb[pixel + 2] = colour[2];
        }
    }
    // Faint gridlines every kilohertz.
    for (double grid = 1000.0; grid < maximum_hz; grid += 1000.0) {
        std::size_t row = static_cast<std::size_t>(grid / maximum_hz * static_cast<double>(height));
        std::size_t y = height - 1 - row;
        for (std::size_t frame = 0; frame < frames; frame += 3) {
            std::size_t pixel = (y * frames + frame) * 3;
            rgb[pixel] = 255;
            rgb[pixel + 1] = 255;
            rgb[pixel + 2] = 255;
        }
    }
    bool written = write_png(path, frames, height, rgb);
    return written;
}

// ---------------------------------------------------------------- pitch

double frame_pitch(const std::vector<double>& signal, std::size_t start, std::size_t length, double rate,
                   double& strength) noexcept {
    const std::size_t shortest = static_cast<std::size_t>(rate / 500.0);
    const std::size_t longest = static_cast<std::size_t>(rate / 60.0);
    double energy = 0.0;
    for (std::size_t index = 0; index < length; ++index) {
        energy += signal[start + index] * signal[start + index];
    }
    strength = 0.0;
    if (energy <= 1.0e-12) {
        return 0.0;
    }
    double best = 0.0;
    std::size_t best_lag = 0;
    for (std::size_t lag = shortest; lag <= longest && lag < length; ++lag) {
        double sum = 0.0;
        double lagged_energy = 0.0;
        for (std::size_t index = 0; index + lag < length; ++index) {
            sum += signal[start + index] * signal[start + index + lag];
            lagged_energy += signal[start + index + lag] * signal[start + index + lag];
        }
        double normalized = sum / std::sqrt(energy * lagged_energy + 1.0e-30);
        if (normalized > best) {
            best = normalized;
            best_lag = lag;
        }
    }
    strength = best;
    if (best_lag == 0) {
        return 0.0;
    }
    double pitch = rate / static_cast<double>(best_lag);
    return pitch;
}

// ---------------------------------------------------------------- formants

std::vector<double> downsample_by_four(const std::vector<double>& signal) {
    // Windowed-sinc low-pass at 5.5 kHz for 48 kHz input, then decimate.
    const std::size_t taps = 63;
    const double cutoff = 5500.0 / 48000.0;
    std::vector<double> kernel(taps, 0.0);
    double sum = 0.0;
    for (std::size_t index = 0; index < taps; ++index) {
        double centre = static_cast<double>(index) - static_cast<double>(taps - 1) / 2.0;
        double sinc = centre == 0.0 ? 2.0 * cutoff : std::sin(2.0 * pi * cutoff * centre) / (pi * centre);
        double window = 0.54 - 0.46 * std::cos(2.0 * pi * static_cast<double>(index) / static_cast<double>(taps - 1));
        kernel[index] = sinc * window;
        sum += kernel[index];
    }
    std::vector<double> output{};
    for (std::size_t position = taps; position < signal.size(); position += 4) {
        double value = 0.0;
        for (std::size_t index = 0; index < taps; ++index) {
            value += kernel[index] * signal[position - index];
        }
        output.push_back(value / sum);
    }
    return output;
}

// LPC envelope peaks (Hz) of one frame: autocorrelation method, Levinson-Durbin.
std::vector<double> frame_formants(const std::vector<double>& signal, std::size_t start, std::size_t length,
                                   std::size_t order, double rate) {
    std::vector<double> frame(length, 0.0);
    for (std::size_t index = 0; index < length; ++index) {
        double previous = start + index > 0 ? signal[start + index - 1] : 0.0;
        double emphasized = signal[start + index] - 0.97 * previous;
        double window = 0.54 - 0.46 * std::cos(2.0 * pi * static_cast<double>(index) / static_cast<double>(length - 1));
        frame[index] = emphasized * window;
    }
    std::vector<double> correlation(order + 1, 0.0);
    for (std::size_t lag = 0; lag <= order; ++lag) {
        double sum = 0.0;
        for (std::size_t index = lag; index < length; ++index) {
            sum += frame[index] * frame[index - lag];
        }
        correlation[lag] = sum;
    }
    std::vector<double> peaks{};
    if (correlation[0] <= 1.0e-14) {
        return peaks;
    }
    correlation[0] *= 1.0001;
    std::vector<double> coefficients(order + 1, 0.0);
    std::vector<double> scratch(order + 1, 0.0);
    coefficients[0] = 1.0;
    double error = correlation[0];
    for (std::size_t step = 1; step <= order; ++step) {
        double accumulation = correlation[step];
        for (std::size_t index = 1; index < step; ++index) {
            accumulation += coefficients[index] * correlation[step - index];
        }
        double reflection = -accumulation / error;
        scratch = coefficients;
        for (std::size_t index = 1; index < step; ++index) {
            coefficients[index] = scratch[index] + reflection * scratch[step - index];
        }
        coefficients[step] = reflection;
        error *= 1.0 - reflection * reflection;
    }
    const std::size_t grid = 1100;
    std::vector<double> envelope(grid, 0.0);
    for (std::size_t point = 0; point < grid; ++point) {
        double frequency = 5500.0 * static_cast<double>(point) / static_cast<double>(grid);
        double angle = 2.0 * pi * frequency / rate;
        std::complex<double> response(0.0, 0.0);
        for (std::size_t index = 0; index <= order; ++index) {
            response += coefficients[index] * std::polar(1.0, -angle * static_cast<double>(index));
        }
        envelope[point] = -10.0 * std::log10(std::norm(response) + 1.0e-30);
    }
    for (std::size_t point = 1; point + 1 < grid; ++point) {
        double frequency = 5500.0 * static_cast<double>(point) / static_cast<double>(grid);
        if (frequency > 180.0 && envelope[point] > envelope[point - 1] && envelope[point] >= envelope[point + 1]) {
            peaks.push_back(frequency);
        }
    }
    return peaks;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    double middle = values[values.size() / 2];
    return middle;
}

bool parse(int count, char** values, Arguments& arguments) {
    for (int index = 1; index < count; ++index) {
        std::string flag = values[index];
        bool has_value = index + 1 < count;
        if (flag == "--spectrogram" && has_value) {
            arguments.spectrogram = values[++index];
        } else if (flag == "--from" && has_value) {
            arguments.from = std::atof(values[++index]);
        } else if (flag == "--to" && has_value) {
            arguments.to = std::atof(values[++index]);
        } else if (flag == "--max-hz" && has_value) {
            arguments.maximum_hz = std::atof(values[++index]);
        } else if (flag == "--track") {
            arguments.track = true;
        } else if (flag == "--clicks") {
            arguments.clicks = true;
        } else if (flag == "--lpc-order" && has_value) {
            arguments.lpc_order = static_cast<std::size_t>(std::atoi(values[++index]));
        } else if (arguments.input.empty()) {
            arguments.input = flag;
        } else {
            return false;
        }
    }
    bool valid = !arguments.input.empty();
    return valid;
}

}

int main(int argc, char** argv) {
    Arguments arguments{};
    if (!parse(argc, argv, arguments)) {
        std::printf("usage: ducklips_analyze FILE.wav [--spectrogram OUT.png] [--from S] [--to S] [--max-hz HZ] [--clicks] [--track]\n");
        return 2;
    }
    ducklips::AudioBuffer audio{};
    if (!ducklips::read_wav(arguments.input, audio)) {
        std::fprintf(stderr, "cannot read %s\n", arguments.input.c_str());
        return 1;
    }
    std::size_t frames = audio.samples.size() / audio.channels;
    std::vector<double> mono(frames, 0.0);
    for (std::size_t frame = 0; frame < frames; ++frame) {
        double sum = 0.0;
        for (std::size_t channel = 0; channel < audio.channels; ++channel) {
            sum += static_cast<double>(audio.samples[frame * audio.channels + channel]);
        }
        mono[frame] = sum / static_cast<double>(audio.channels);
    }
    double rate = audio.sample_rate;
    std::size_t first = static_cast<std::size_t>(std::max(0.0, arguments.from) * rate);
    std::size_t last = frames;
    if (arguments.to > 0.0) {
        last = std::min(frames, static_cast<std::size_t>(arguments.to * rate));
    }
    first = std::min(first, last);
    std::vector<double> span(mono.begin() + static_cast<std::ptrdiff_t>(first), mono.begin() + static_cast<std::ptrdiff_t>(last));

    double peak = 0.0;
    double energy = 0.0;
    for (std::size_t index = 0; index < span.size(); ++index) {
        peak = std::max(peak, std::fabs(span[index]));
        energy += span[index] * span[index];
    }
    double rms = std::sqrt(energy / std::max<double>(1.0, static_cast<double>(span.size())));
    std::printf("file %s: %.3f s, %.0f Hz, %zu channels; range %.3f-%.3f s peak %.3f rms %.4f\n",
                arguments.input.c_str(), static_cast<double>(frames) / rate, rate, audio.channels,
                static_cast<double>(first) / rate, static_cast<double>(last) / rate, peak, rms);

    if (!arguments.spectrogram.empty()) {
        if (render_spectrogram(span, rate, arguments.maximum_hz, arguments.spectrogram)) {
            std::printf("spectrogram %s (0-%.0f Hz, 75 dB range)\n", arguments.spectrogram.c_str(), arguments.maximum_hz);
        }
    }

    if (arguments.clicks) {
        // Impulsive events: the second difference of the signal against its
        // local 20 ms rms. Speech has no sample-to-sample corners that sharp.
        const std::size_t half = static_cast<std::size_t>(0.010 * rate);
        std::vector<double> curvature(span.size(), 0.0);
        for (std::size_t index = 1; index + 1 < span.size(); ++index) {
            curvature[index] = std::fabs(span[index + 1] - 2.0 * span[index] + span[index - 1]);
        }
        std::vector<std::pair<double, std::size_t>> events{};
        for (std::size_t index = half; index + half < span.size(); index += 1) {
            double local = 0.0;
            for (std::size_t offset = index - half; offset < index + half; offset += 8) {
                local += curvature[offset] * curvature[offset];
            }
            double local_rms = std::sqrt(local / static_cast<double>(2 * half / 8)) + 1.0e-9;
            double ratio = curvature[index] / local_rms;
            if (ratio > 8.0) {
                events.push_back(std::pair<double, std::size_t>(ratio, index));
                index += half;
            }
        }
        std::sort(events.rbegin(), events.rend());
        std::printf("impulsive events (second difference / local rms > 8): %zu\n", events.size());
        for (std::size_t index = 0; index < events.size() && index < 12; ++index) {
            std::printf("  %.4f s  ratio %.1f  sample %.3f\n", static_cast<double>(first + events[index].second) / rate,
                        events[index].first, span[events[index].second]);
        }
    }

    if (arguments.track && rate == 48000.0) {
        // Formant track every 5 ms (LPC peaks on a 12 kHz copy, 20 ms frames),
        // plus the spectral centroid above 1 kHz of a 10 ms frame (burst cue).
        std::vector<double> low = downsample_by_four(span);
        const std::size_t frame = static_cast<std::size_t>(0.020 * 12000.0);
        const std::size_t step = static_cast<std::size_t>(0.005 * 12000.0);
        const std::size_t burst = 512;
        std::vector<std::complex<double>> buffer(burst);
        std::printf("time      rms     F1    F2    F3   centroid>1k\n");
        for (std::size_t start = 0; start + frame < low.size(); start += step) {
            double energy_sum = 0.0;
            for (std::size_t index = 0; index < frame; ++index) {
                energy_sum += low[start + index] * low[start + index];
            }
            double frame_rms = std::sqrt(energy_sum / static_cast<double>(frame));
            std::vector<double> peaks = frame_formants(low, start, frame, arguments.lpc_order, 12000.0);
            std::size_t centre = (start + frame / 2) * 4;
            double centroid = 0.0;
            if (centre + burst / 2 < span.size() && centre >= burst / 2) {
                for (std::size_t index = 0; index < burst; ++index) {
                    double hann = 0.5 - 0.5 * std::cos(2.0 * pi * static_cast<double>(index) / static_cast<double>(burst));
                    buffer[index] = std::complex<double>(span[centre - burst / 2 + index] * hann, 0.0);
                }
                fft(buffer);
                double weighted = 0.0;
                double total = 0.0;
                for (std::size_t bin = 1; bin < burst / 2; ++bin) {
                    double frequency = rate * static_cast<double>(bin) / static_cast<double>(burst);
                    if (frequency < 1000.0) {
                        continue;
                    }
                    double power = std::norm(buffer[bin]);
                    weighted += frequency * power;
                    total += power;
                }
                centroid = total > 0.0 ? weighted / total : 0.0;
            }
            double f1 = peaks.size() > 0 ? peaks[0] : 0.0;
            double f2 = peaks.size() > 1 ? peaks[1] : 0.0;
            double f3 = peaks.size() > 2 ? peaks[2] : 0.0;
            double time = static_cast<double>(first) / rate + static_cast<double>(start + frame / 2) / 12000.0;
            std::printf("%.3f  %.4f  %4.0f  %4.0f  %4.0f   %5.0f\n", time, frame_rms, f1, f2, f3, centroid);
        }
    }

    // Pitch and formants over 10 ms hops, voiced and loud frames only.
    std::vector<double> pitches{};
    std::vector<std::vector<double>> formant_tracks(4);
    std::vector<double> decimated = rate == 48000.0 ? downsample_by_four(span) : std::vector<double>{};
    const std::size_t pitch_window = static_cast<std::size_t>(0.040 * rate);
    const std::size_t hop = static_cast<std::size_t>(0.010 * rate);
    std::size_t voiced_frames = 0;
    for (std::size_t start = 0; start + pitch_window < span.size(); start += hop) {
        double frame_energy = 0.0;
        for (std::size_t index = 0; index < pitch_window; ++index) {
            frame_energy += span[start + index] * span[start + index];
        }
        double frame_rms = std::sqrt(frame_energy / static_cast<double>(pitch_window));
        if (frame_rms < 0.1 * rms) {
            continue;
        }
        double strength = 0.0;
        double pitch = frame_pitch(span, start, pitch_window, rate, strength);
        if (strength < 0.5 || pitch <= 0.0) {
            continue;
        }
        voiced_frames += 1;
        pitches.push_back(pitch);
        if (!decimated.empty()) {
            std::size_t low_start = start / 4 + static_cast<std::size_t>(0.0075 * 12000.0);
            std::size_t low_length = static_cast<std::size_t>(0.025 * 12000.0);
            if (low_start + low_length < decimated.size()) {
                std::vector<double> peaks = frame_formants(decimated, low_start, low_length, arguments.lpc_order, 12000.0);
                for (std::size_t index = 0; index < peaks.size() && index < formant_tracks.size(); ++index) {
                    formant_tracks[index].push_back(peaks[index]);
                }
            }
        }
    }
    std::printf("voiced frames %zu; median F0 %.1f Hz", voiced_frames, median(pitches));
    if (!pitches.empty()) {
        std::sort(pitches.begin(), pitches.end());
        std::printf(" (10th %.1f, 90th %.1f)", pitches[pitches.size() / 10], pitches[(pitches.size() * 9) / 10]);
    }
    std::printf("\n");
    if (!decimated.empty()) {
        std::printf("median LPC peaks: F1 %.0f  F2 %.0f  F3 %.0f  F4 %.0f Hz\n", median(formant_tracks[0]),
                    median(formant_tracks[1]), median(formant_tracks[2]), median(formant_tracks[3]));
    }
    return 0;
}
