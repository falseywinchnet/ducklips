#include "analysis.hpp"
#include "ducklips/physics.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <cmath>

namespace boundary {

namespace {

// RBJ low-pass biquad, direct form I.
struct Biquad final {
    double b0{1.0};
    double b1{0.0};
    double b2{0.0};
    double a1{0.0};
    double a2{0.0};
    double x1{0.0};
    double x2{0.0};
    double y1{0.0};
    double y2{0.0};
};

Biquad butterworth_lowpass(double cutoff, double rate) {
    double omega = 2.0 * ducklips::pi * cutoff / rate;
    double alpha = std::sin(omega) / (2.0 * std::sqrt(0.5));
    double cosine = std::cos(omega);
    double a0 = 1.0 + alpha;
    Biquad filter{};
    filter.b0 = (1.0 - cosine) * 0.5 / a0;
    filter.b1 = (1.0 - cosine) / a0;
    filter.b2 = filter.b0;
    filter.a1 = -2.0 * cosine / a0;
    filter.a2 = (1.0 - alpha) / a0;
    return filter;
}

double run_biquad(Biquad& filter, double x) {
    double y = filter.b0 * x + filter.b1 * filter.x1 + filter.b2 * filter.x2 - filter.a1 * filter.y1 - filter.a2 * filter.y2;
    filter.x2 = filter.x1;
    filter.x1 = x;
    filter.y2 = filter.y1;
    filter.y1 = y;
    return y;
}

double decibels(double power) {
    double result = 10.0 * std::log10(std::max(power, 1.0e-300));
    return result;
}

}

void make_jet_noise(std::uint64_t seed, double rate, std::size_t count, std::vector<double>& output) {
    ducklips::TurbulenceNoise noise{};
    noise.prepare(seed, 500.0, rate);
    Biquad first = butterworth_lowpass(12000.0, rate);
    Biquad second = butterworth_lowpass(12000.0, rate);
    output.assign(count, 0.0);
    for (std::size_t i = 0; i < count; ++i) {
        double x = noise.next();
        double y = run_biquad(first, x);
        output[i] = run_biquad(second, y);
    }
}

void make_band_impulse(double rate, std::size_t count, std::vector<double>& output) {
    Biquad first = butterworth_lowpass(12000.0, rate);
    Biquad second = butterworth_lowpass(12000.0, rate);
    output.assign(count, 0.0);
    for (std::size_t i = 0; i < count; ++i) {
        double x = 0.0;
        if (i == 0) {
            x = 1.0;
        }
        double y = run_biquad(first, x);
        output[i] = run_biquad(second, y);
    }
}

void Decimator::prepare(std::size_t factor, double input_rate, double cutoff_hz) {
    factor_ = factor;
    std::size_t length = 16 * factor + 1;
    taps_.assign(length, 0.0);
    double fc = cutoff_hz / input_rate;
    double middle = 0.5 * static_cast<double>(length - 1);
    double sum = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        double x = static_cast<double>(i) - middle;
        double sinc = 2.0 * fc;
        if (std::fabs(x) > 1.0e-12) {
            sinc = std::sin(2.0 * ducklips::pi * fc * x) / (ducklips::pi * x);
        }
        double phase = 2.0 * ducklips::pi * static_cast<double>(i) / static_cast<double>(length - 1);
        double window = 0.42 - 0.5 * std::cos(phase) + 0.08 * std::cos(2.0 * phase);
        taps_[i] = sinc * window;
        sum += taps_[i];
    }
    for (std::size_t i = 0; i < length; ++i) {
        taps_[i] /= sum;
    }
    history_.assign(length, 0.0);
    reset();
}

void Decimator::reset() noexcept {
    std::fill(history_.begin(), history_.end(), 0.0);
    phase_ = 0;
    write_ = 0;
}

bool Decimator::push(double input, double& output) noexcept {
    const std::size_t length = history_.size();
    history_[write_] = input;
    write_ += 1;
    if (write_ == length) {
        write_ = 0;
    }
    phase_ += 1;
    if (phase_ < factor_) {
        return false;
    }
    phase_ = 0;
    double sum = 0.0;
    std::size_t index = write_; // oldest sample
    for (std::size_t k = 0; k < length; ++k) {
        sum += taps_[k] * history_[index];
        index += 1;
        if (index == length) {
            index = 0;
        }
    }
    output = sum;
    return true;
}

void fft_in_place(std::vector<std::complex<double>>& data) {
    const std::size_t n = data.size();
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t bit = n >> 1;
        while ((j & bit) != 0) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }
    for (std::size_t length = 2; length <= n; length <<= 1) {
        double angle = -2.0 * ducklips::pi / static_cast<double>(length);
        std::complex<double> step(std::cos(angle), std::sin(angle));
        for (std::size_t start = 0; start < n; start += length) {
            std::complex<double> w(1.0, 0.0);
            std::size_t half = length / 2;
            for (std::size_t k = 0; k < half; ++k) {
                std::complex<double> even = data[start + k];
                std::complex<double> odd = w * data[start + k + half];
                data[start + k] = even + odd;
                data[start + k + half] = even - odd;
                w *= step;
            }
        }
    }
}

void welch_power(const std::vector<double>& signal, std::size_t begin, std::size_t end, std::size_t fft_size,
                 std::vector<double>& psd) {
    psd.assign(fft_size / 2 + 1, 0.0);
    std::vector<std::complex<double>> frame(fft_size);
    std::vector<double> window(fft_size, 0.0);
    double window_power = 0.0;
    for (std::size_t i = 0; i < fft_size; ++i) {
        window[i] = 0.5 - 0.5 * std::cos(2.0 * ducklips::pi * static_cast<double>(i) / static_cast<double>(fft_size));
        window_power += window[i] * window[i];
    }
    std::size_t hop = fft_size / 2;
    std::size_t frames = 0;
    for (std::size_t start = begin; start + fft_size <= end; start += hop) {
        for (std::size_t i = 0; i < fft_size; ++i) {
            frame[i] = std::complex<double>(signal[start + i] * window[i], 0.0);
        }
        fft_in_place(frame);
        for (std::size_t k = 0; k < psd.size(); ++k) {
            psd[k] += std::norm(frame[k]);
        }
        frames += 1;
    }
    if (frames == 0) {
        return;
    }
    // Normalized so that the sum over bins (one side, doubled) equals the mean square.
    double scale = 1.0 / (static_cast<double>(frames) * window_power * static_cast<double>(fft_size));
    for (std::size_t k = 0; k < psd.size(); ++k) {
        psd[k] *= scale * static_cast<double>(fft_size);
    }
}

void smooth_bins(const std::vector<double>& input, std::size_t half_width, std::vector<double>& output) {
    const std::size_t n = input.size();
    output.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t low = 0;
        if (i > half_width) {
            low = i - half_width;
        }
        std::size_t high = std::min(n - 1, i + half_width);
        double sum = 0.0;
        for (std::size_t k = low; k <= high; ++k) {
            sum += input[k];
        }
        output[i] = sum / static_cast<double>(high - low + 1);
    }
}

double rms(const std::vector<double>& signal, std::size_t begin, std::size_t end) {
    double sum = 0.0;
    std::size_t last = std::min(end, signal.size());
    for (std::size_t i = begin; i < last; ++i) {
        sum += signal[i] * signal[i];
    }
    double count = static_cast<double>(std::max<std::size_t>(last - std::min(begin, last), 1));
    double result = std::sqrt(sum / count);
    return result;
}

SpectrumSummary summarize(const std::vector<double>& signal, double rate, double from_s, double to_s) {
    SpectrumSummary summary{};
    std::size_t begin = static_cast<std::size_t>(std::llround(from_s * rate));
    std::size_t end = std::min(signal.size(), static_cast<std::size_t>(std::llround(to_s * rate)));
    summary.rms_db = 20.0 * std::log10(std::max(rms(signal, begin, end), 1.0e-300));
    const std::size_t fft_size = 1024;
    std::vector<double> psd{};
    welch_power(signal, begin, end, fft_size, psd);
    double bin_hz = rate / static_cast<double>(fft_size);
    std::vector<double> smooth{};
    smooth_bins(psd, static_cast<std::size_t>(std::llround(150.0 / bin_hz)), smooth);
    double best = -1.0;
    double low_band = 0.0;
    double high_band = 0.0;
    double weighted = 0.0;
    double total = 0.0;
    for (std::size_t k = 0; k < psd.size(); ++k) {
        double hz = static_cast<double>(k) * bin_hz;
        if (hz >= 1500.0 && hz <= 11000.0 && smooth[k] > best) {
            best = smooth[k];
            summary.peak_hz = hz;
        }
        if (hz >= 2500.0 && hz < 4000.0) {
            low_band += psd[k];
        }
        if (hz >= 4000.0 && hz < 8000.0) {
            high_band += psd[k];
        }
        if (hz >= 1000.0 && hz <= 11000.0) {
            weighted += hz * psd[k];
            total += psd[k];
        }
    }
    summary.band_2500_4000_db = decibels(low_band);
    summary.band_4000_8000_db = decibels(high_band);
    if (total > 0.0) {
        summary.centroid_hz = weighted / total;
    }
    return summary;
}

void transfer_peaks(const std::vector<double>& impulse_response, double rate, double low, double high,
                    double floor_db, double smooth_hz, std::vector<Peak>& peaks) {
    peaks.clear();
    std::size_t size = 1;
    while (size < impulse_response.size()) {
        size <<= 1;
    }
    std::vector<std::complex<double>> spectrum(size);
    for (std::size_t i = 0; i < size; ++i) {
        double x = 0.0;
        if (i < impulse_response.size()) {
            x = impulse_response[i];
        }
        spectrum[i] = std::complex<double>(x, 0.0);
    }
    fft_in_place(spectrum);
    std::size_t bins = size / 2 + 1;
    std::vector<double> power(bins, 0.0);
    for (std::size_t k = 0; k < bins; ++k) {
        power[k] = std::norm(spectrum[k]);
    }
    double bin_hz = rate / static_cast<double>(size);
    std::vector<double> smooth{};
    smooth_bins(power, static_cast<std::size_t>(std::llround(smooth_hz / bin_hz)), smooth);
    std::size_t first = static_cast<std::size_t>(std::ceil(low / bin_hz));
    std::size_t last = std::min(bins - 2, static_cast<std::size_t>(std::floor(high / bin_hz)));
    double maximum = 0.0;
    for (std::size_t k = first; k <= last; ++k) {
        maximum = std::max(maximum, smooth[k]);
    }
    double maximum_db = decibels(maximum);
    std::size_t neighbourhood = static_cast<std::size_t>(std::llround(250.0 / bin_hz));
    for (std::size_t k = first; k <= last; ++k) {
        double level = smooth[k];
        if (decibels(level) < maximum_db - floor_db) {
            continue;
        }
        bool is_peak = true;
        std::size_t from = 0;
        if (k > neighbourhood) {
            from = k - neighbourhood;
        }
        std::size_t to = std::min(bins - 1, k + neighbourhood);
        for (std::size_t m = from; m <= to; ++m) {
            if (smooth[m] > level) {
                is_peak = false;
                break;
            }
        }
        if (!is_peak) {
            continue;
        }
        Peak peak{};
        peak.hz = static_cast<double>(k) * bin_hz;
        peak.db = decibels(level) - maximum_db;
        double half = 0.5 * level;
        std::size_t left = k;
        while (left > 0 && smooth[left] > half) {
            left -= 1;
        }
        std::size_t right = k;
        while (right + 1 < bins && smooth[right] > half) {
            right += 1;
        }
        if (left > 0 && right + 1 < bins) {
            peak.bandwidth_hz = static_cast<double>(right - left) * bin_hz;
        }
        peaks.push_back(peak);
    }
}

}
