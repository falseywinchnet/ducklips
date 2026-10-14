#pragma once
// Signals and measurements shared by both models: the jet source, the
// decimator to 48 kHz, spectra and the summary measures.
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace boundary {

// Jet noise at the internal rate: DuckLips TurbulenceNoise (Gaussian, 500 Hz
// one-pole high-pass) through a 12 kHz Butterworth low-pass so the loop is
// not driven above the band of interest. Fills output (resized once).
void make_jet_noise(std::uint64_t seed, double rate, std::size_t count, std::vector<double>& output);
// A unit impulse through the same 12 kHz low-pass.
void make_band_impulse(double rate, std::size_t count, std::vector<double>& output);

// Windowed-sinc (Blackman) decimator by an integer factor.
class Decimator final {
public:
    void prepare(std::size_t factor, double input_rate, double cutoff_hz);
    void reset() noexcept;
    // Returns true when an output sample was produced into output.
    bool push(double input, double& output) noexcept;

private:
    std::size_t factor_{1};
    std::size_t phase_{0};
    std::size_t write_{0};
    std::vector<double> taps_{};
    std::vector<double> history_{}; // ring of taps_.size()
};

void fft_in_place(std::vector<std::complex<double>>& data);

// Welch power spectrum (Hann, 50% overlap) of signal[begin, end). psd has
// fft_size/2 + 1 bins; bin spacing is rate / fft_size.
void welch_power(const std::vector<double>& signal, std::size_t begin, std::size_t end, std::size_t fft_size,
                 std::vector<double>& psd);
// Moving average over +- half_width bins.
void smooth_bins(const std::vector<double>& input, std::size_t half_width, std::vector<double>& output);

struct SpectrumSummary final {
    double rms_db{0.0};
    double peak_hz{0.0};          // maximum of the 300 Hz-smoothed spectrum, 1.5-11 kHz
    double band_2500_4000_db{0.0};
    double band_4000_8000_db{0.0};
    double centroid_hz{0.0};      // 1-11 kHz
};

[[nodiscard]] SpectrumSummary summarize(const std::vector<double>& signal, double rate, double from_s, double to_s);

struct Peak final {
    double hz{0.0};
    double db{0.0};
    double bandwidth_hz{0.0};     // -3 dB width, 0 when not resolved
};

// Peaks of the impulse response's magnitude spectrum in [low, high] Hz within
// floor_db of the maximum, smoothed over +- smooth_hz. Clears and fills peaks.
void transfer_peaks(const std::vector<double>& impulse_response, double rate, double low, double high,
                    double floor_db, double smooth_hz, std::vector<Peak>& peaks);

[[nodiscard]] double rms(const std::vector<double>& signal, std::size_t begin, std::size_t end);

}
