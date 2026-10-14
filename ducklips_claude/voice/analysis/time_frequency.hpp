#pragma once

// Four time-frequency views of a sound, built for comparing synthetic speech
// with human recordings. Each view sits at a different corner of the Gabor
// uncertainty box:
//
//   stft            - fixed window, linear frequency (the reference view)
//   constant_q      - Garces (2020) quantized constant-Q Gabor atoms:
//                     Gaussian atoms whose duration scales with period
//   causal_gabor    - Lindeberg (2025) time-causal analogue of the Gabor
//                     transform: a cascade of first-order recursive filters
//                     (the time-causal limit kernel) on demodulated signal
//   pseudo_wigner   - smoothed pseudo-Wigner-Ville distribution: quadratic,
//                     the sharpest joint resolution, with cross-terms
//
// Each transform has a plain field (rows = analysis bands) and an enhanced
// field on a fine grid that reaches pixel resolution: synchrosqueezed
// (energy moved to its instantaneous frequency) and fused across several
// resolutions, as the enhanced STFT fuses its apertures.

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace tf {

using Complex = std::complex<double>;

struct Signal {
  std::vector<double> samples;
  double rate = 0.0;
};

// A row-major field: value(row, column) = data[row * columns + column].
struct Field {
  std::int32_t rows = 0;
  std::int32_t columns = 0;
  std::vector<float> data;
  std::vector<float> hz;       // centre frequency of each row
  std::vector<float> seconds;  // centre time of each column
};

struct Grid {
  double hop_seconds = 0.002;   // column spacing
  double low_hz = 60.0;         // lowest band centre
  double high_hz = 7600.0;      // highest band centre
  std::int32_t bands_per_octave = 24;
  std::int32_t fine_per_octave = 96;  // rows of the enhanced (reassigned) grid
};

[[nodiscard]] bool read_wave(const std::string& path, Signal& signal);
[[nodiscard]] bool write_npy(const std::string& path, const std::vector<float>& data, std::int32_t rows,
                             std::int32_t columns);

// Discrete Fourier transform of a power-of-two length, in place.
class Fft {
 public:
  explicit Fft(std::int32_t size);
  void forward(std::vector<Complex>& values) const;
  void inverse(std::vector<Complex>& values) const;
  [[nodiscard]] std::int32_t size() const { return size_; }

 private:
  void run(std::vector<Complex>& values, bool inverse) const;
  std::int32_t size_ = 0;
  std::vector<Complex> twiddles_;
  std::vector<std::int32_t> reversed_;
};

std::int32_t next_power_of_two(std::int32_t value);
std::vector<Complex> analytic_signal(const Signal& signal);
std::vector<double> band_centres(const Grid& grid);
std::vector<float> column_times(const Signal& signal, const Grid& grid);

struct Analysis {
  Field plain;
  Field enhanced;
};

Analysis stft(const Signal& signal, const Grid& grid, double window_seconds);
Analysis constant_q(const Signal& signal, const Grid& grid);
Analysis causal_gabor(const Signal& signal, const Grid& grid);
Analysis pseudo_wigner(const Signal& signal, const Grid& grid);

}  // namespace tf
