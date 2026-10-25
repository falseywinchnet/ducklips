// Several voices at once, each on its own thread: every voice speaks its
// own score (or the same score with its own seed and pitch), and the run
// reports whether all of them keep up with real time together.
// usage: voices OUT.wav SCORE [SCORE...]
//   PITCH_STEP=semitones moves each further voice's pitch (default 3);
//   COPIES=n speaks the first score n times when only one is given.

#include "tract/speaker.hpp"
#include "voice/score.hpp"
#include "wave_out.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

namespace {

bool read_score(const char* path, voice::Score& score) {
  std::ifstream input(path);
  if (!input) {
    return false;
  }
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!voice::parse_score_line(line, score)) {
      return false;
    }
  }
  voice::sort_score(score);
  return true;
}

// Moves every pitch target of a score by a ratio (another speaker's range).
void transpose(voice::Score& score, double ratio) {
  for (voice::PartEvent& event : score.parts) {
    if (event.part == voice::Part::pitch) {
      event.a = static_cast<float>(event.a * ratio);
    }
  }
}

struct Job {
  voice::Score score;
  std::uint32_t seed = 1U;
  std::vector<double> samples;  // 48 kHz
  double seconds = 0.0;         // wall time spent
};

// One voice from start to end; owns nothing shared.
void run_voice(Job* job) {
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  tract::BodyLayout layout = tract::speech_layout();
  layout.seed = job->seed;
  tract::Speaker speaker(layout, job->score, job->seed);
  Biquad first = lowpass(18000.0, tract::kStepRate, 0.5412);
  Biquad second = lowpass(18000.0, tract::kStepRate, 1.3066);
  long total = static_cast<long>(std::ceil(job->score.duration * tract::kStepRate));
  job->samples.reserve(static_cast<std::size_t>(total / 2 + 1));
  for (long step = 0; step < total; ++step) {
    double filtered = run(second, run(first, speaker.step()));
    if (step % 2 == 0) {
      job->samples.push_back(std::isfinite(filtered) ? filtered : 0.0);
    }
  }
  job->seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
}

void flush_subnormals() {
#if defined(__aarch64__)
  std::uint64_t control = 0;
  __asm__ volatile("mrs %0, fpcr" : "=r"(control));
  control |= (1ULL << 24);
  __asm__ volatile("msr fpcr, %0" : : "r"(control));
#endif
}

void voice_thread(Job* job) {
  flush_subnormals();
  run_voice(job);
}

}  // namespace

int main(int count, char** values) {
  if (count < 3) {
    std::fprintf(stderr, "usage: voices OUT.wav SCORE [SCORE...]\n");
    return 2;
  }
  double step_semitones = std::getenv("PITCH_STEP") != nullptr ? std::strtod(std::getenv("PITCH_STEP"), nullptr) : 3.0;
  int copies = std::getenv("COPIES") != nullptr ? std::atoi(std::getenv("COPIES")) : 1;
  std::vector<std::string> paths;
  for (int index = 2; index < count; ++index) {
    paths.push_back(values[index]);
  }
  if (paths.size() == 1) {
    for (int copy = 1; copy < copies; ++copy) {
      paths.push_back(paths[0]);
    }
  }
  std::vector<Job> jobs(paths.size());
  for (std::size_t index = 0; index < paths.size(); ++index) {
    if (!read_score(paths[index].c_str(), jobs[index].score)) {
      std::fprintf(stderr, "cannot read %s\n", paths[index].c_str());
      return 1;
    }
    transpose(jobs[index].score, std::pow(2.0, step_semitones * static_cast<double>(index) / 12.0));
    jobs[index].seed = 73021U + 7919U * static_cast<std::uint32_t>(index);
  }
  std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
  std::vector<std::thread> threads;
  for (Job& job : jobs) {
    threads.emplace_back(voice_thread, &job);
  }
  for (std::thread& thread : threads) {
    thread.join();
  }
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - began).count();
  // Mix, each voice levelled by its own peak.
  std::size_t length = 0;
  for (const Job& job : jobs) {
    length = std::max(length, job.samples.size());
  }
  std::vector<double> mix(length, 0.0);
  double longest = 0.0;
  for (const Job& job : jobs) {
    double peak = 1e-12;
    for (double value : job.samples) {
      peak = std::fmax(peak, std::fabs(value));
    }
    for (std::size_t index = 0; index < job.samples.size(); ++index) {
      mix[index] += job.samples[index] / peak / static_cast<double>(jobs.size());
    }
    longest = std::fmax(longest, job.score.duration);
    std::fprintf(stderr, "voice seed %u: %.2f s of speech in %.2f s (%.1fx)\n", job.seed, job.score.duration,
                 job.seconds, job.score.duration / job.seconds);
  }
  std::vector<float> samples;
  samples.reserve(mix.size());
  for (double value : mix) {
    samples.push_back(static_cast<float>(0.9 * value));
  }
  write_wave(values[1], samples, static_cast<int>(tract::kAudioRate));
  std::fprintf(stderr, "%zu voices together: %.2f s of speech in %.2f s wall (%.1fx real time)\n", jobs.size(),
               longest, wall, longest / wall);
  return 0;
}
