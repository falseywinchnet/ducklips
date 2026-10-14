// Runs the four time-frequency views on one WAV and saves every field.
// usage: tf_analyze IN.wav OUT_PREFIX [--hop s] [--low hz] [--high hz]
// Writes OUT_PREFIX_<view>_{plain,enhanced}.npy, *_hz.npy and OUT_PREFIX_seconds.npy.

#include "time_frequency.hpp"

#include <chrono>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <string>

namespace {

bool save_field(const std::string& prefix, const std::string& name, const tf::Field& field) {
  bool good = tf::write_npy(prefix + "_" + name + ".npy", field.data, field.rows, field.columns);
  good = tf::write_npy(prefix + "_" + name + "_hz.npy", field.hz, 1, static_cast<std::int32_t>(field.hz.size())) && good;
  return good;
}

double seconds_since(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

}  // namespace

int main(int count, char** values) {
  if (count < 3) {
    std::fprintf(stderr, "usage: tf_analyze IN.wav OUT_PREFIX [--hop s] [--low hz] [--high hz]\n");
    return 2;
  }
  tf::Grid grid;
  std::string views = "stft,cq,causal,wigner";
  for (int index = 3; index + 1 < count; index += 2) {
    std::string flag = values[index];
    if (flag == "--views") {
      views = values[index + 1];
      continue;
    }
    double value = std::strtod(values[index + 1], nullptr);
    if (flag == "--hop") {
      grid.hop_seconds = value;
    } else if (flag == "--low") {
      grid.low_hz = value;
    } else if (flag == "--high") {
      grid.high_hz = value;
    }
  }
  tf::Signal signal;
  if (!tf::read_wave(values[1], signal)) {
    std::fprintf(stderr, "cannot read %s\n", values[1]);
    return 1;
  }
  grid.high_hz = std::min(grid.high_hz, 0.47 * signal.rate);
  std::string prefix = values[2];
  std::vector<float> seconds = tf::column_times(signal, grid);
  bool good = tf::write_npy(prefix + "_seconds.npy", seconds, 1, static_cast<std::int32_t>(seconds.size()));

  // Only the requested views are computed (a tuner needs two of four).
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  tf::Analysis view;
  if (views.find("stft") != std::string::npos) {
    view = tf::stft(signal, grid, 0.025);
    good = save_field(prefix, "stft_plain", view.plain) && save_field(prefix, "stft_enhanced", view.enhanced) && good;
  }
  if (views.find("cq") != std::string::npos) {
    view = tf::constant_q(signal, grid);
    good = save_field(prefix, "cq_plain", view.plain) && save_field(prefix, "cq_enhanced", view.enhanced) && good;
  }
  if (views.find("causal") != std::string::npos) {
    view = tf::causal_gabor(signal, grid);
    good = save_field(prefix, "causal_plain", view.plain) && save_field(prefix, "causal_enhanced", view.enhanced) && good;
  }
  if (views.find("wigner") != std::string::npos) {
    view = tf::pseudo_wigner(signal, grid);
    good = save_field(prefix, "wigner_plain", view.plain) && save_field(prefix, "wigner_enhanced", view.enhanced) && good;
  }
  std::fprintf(stderr, "%s: %.3f s for %.2f s of %.0f Hz audio\n", views.c_str(), seconds_since(start),
               signal.samples.size() / signal.rate, signal.rate);
  return good ? 0 : 1;
}
