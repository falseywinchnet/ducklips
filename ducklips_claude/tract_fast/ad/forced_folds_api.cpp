// C entry points (ctypes) for the forced folds and their forward-mode
// tangents: analysis/gradients.py calls them in-process, in place of a
// sustain FORCEWAVE render (a subprocess, temporary files, and the whole
// body's vocal tube, nose and branches, none of which reach forced folds).
//
// Arrays are row-major float64. Tangent arrays are [step][direction].
//
// Tangent blocks start late: a block whose inputs' tangents are all zero
// before step s0 has zero tangents there too (the folds are causal), so it
// starts from the primal state at the last checkpoint at or before s0,
// promoted (zero tangents). A phone's log areas reach the sound only from
// its keyframes on: in a three-vowel word the third phone's blocks start
// about 60% of the way in.

#include "forced_folds.hpp"

#include <cmath>
#include <cstring>
#include <thread>
#include <vector>

namespace {

using tract_ad::Constant;
using tract_ad::Dual;
using tract_ad::ForcedFolds;
using tract_ad::PostureCopy;

constexpr long kCheckpointSteps = 1024;  // primal checkpoints this far apart (each about 8 kB)

// The posture target of sustain FOLDKEYS at a 1 ms tick (96 steps): the
// targets' row for this tick (the last row past the end).
template <class S>
PostureCopy<S> posture_at(const double* targets, long tick, long ticks) {
  long row = tick < ticks ? tick : ticks - 1;
  PostureCopy<S> posture;
  posture.lower_rest = Constant<S>::make(targets[3 * row]);
  posture.upper_rest = Constant<S>::make(targets[3 * row + 1]);
  posture.chink = Constant<S>::make(targets[3 * row + 2]);
  return posture;
}

// A block's lanes: directions first .. first + width - 1 (lanes past width
// stay zero).
struct Lanes {
  int directions = 0;
  int first = 0;
  int width = 0;
};

// The same with tangents: targets_dot is [tick][3][directions].
template <int N>
PostureCopy<Dual<N>> posture_at_dual(const double* targets, const double* targets_dot, long tick, long ticks,
                                     const Lanes& lanes) {
  PostureCopy<Dual<N>> posture = posture_at<Dual<N>>(targets, tick, ticks);
  if (targets_dot == nullptr) {
    return posture;
  }
  long row = tick < ticks ? tick : ticks - 1;
  Dual<N>* fields[3] = {&posture.lower_rest, &posture.upper_rest, &posture.chink};
  for (int field = 0; field < 3; ++field) {
    const double* source = targets_dot + (3 * row + field) * static_cast<long>(lanes.directions) + lanes.first;
    for (int k = 0; k < lanes.width; ++k) {
      fields[field]->d[k] = source[k];
    }
  }
  return posture;
}

template <int N>
Dual<N> load_dual(double value, const double* row, const Lanes& lanes) {
  Dual<N> out = Constant<Dual<N>>::make(value);
  const double* source = row + lanes.first;
  for (int k = 0; k < lanes.width; ++k) {
    out.d[k] = source[k];
  }
  return out;
}

// Everything a block reads.
struct FoldsWork {
  long steps = 0;
  int vocal_count = 0;
  double pitch = 120.0;
  double lung = 800.0;
  const double* waves = nullptr;
  const double* waves_dot = nullptr;
  const double* area0 = nullptr;
  const double* area0_dot = nullptr;
  const double* targets = nullptr;
  const double* targets_dot = nullptr;
  long ticks = 0;
  int directions = 0;
  const int* block_starts = nullptr;  // block b: directions block_starts[b] .. block_starts[b + 1] - 1
  int blocks = 0;
  const std::vector<ForcedFolds<double>>* checkpoints = nullptr;  // the primal state at every kCheckpointSteps
  int offset = 0;  // this thread's blocks: every `stride`-th from `offset`
  int stride = 1;
  double* flow_dot = nullptr;
};

// The first step where any of a block's input tangents is nonzero (steps if
// none).
long first_moved_step(const FoldsWork& work, const Lanes& lanes) {
  for (long step = 0; step < work.steps; ++step) {
    const double* wave_row = work.waves_dot + step * static_cast<long>(lanes.directions) + lanes.first;
    const double* area_row = work.area0_dot + step * static_cast<long>(lanes.directions) + lanes.first;
    for (int k = 0; k < lanes.width; ++k) {
      if (wave_row[k] != 0.0 || area_row[k] != 0.0) {
        return step;
      }
    }
    if (work.targets_dot != nullptr && step % 96 == 0) {
      long tick = step / 96 < work.ticks ? step / 96 : work.ticks - 1;
      for (int field = 0; field < 3; ++field) {
        const double* source = work.targets_dot + (3 * tick + field) * static_cast<long>(lanes.directions) + lanes.first;
        for (int k = 0; k < lanes.width; ++k) {
          if (source[k] != 0.0) {
            return step;
          }
        }
      }
    }
  }
  return work.steps;
}

// One block of up to N tangent directions, from its checkpoint to the end.
template <int N>
void run_block(const FoldsWork& work, const Lanes& lanes) {
  long moved = first_moved_step(work, lanes);
  if (moved >= work.steps) {
    return;  // no input tangent: the flow's tangents stay zero
  }
  long checkpoint = moved / kCheckpointSteps;
  long start = checkpoint * kCheckpointSteps;
  ForcedFolds<Dual<N>> folds(work.vocal_count, work.pitch, work.lung);
  folds.promote_from((*work.checkpoints)[static_cast<std::size_t>(checkpoint)]);
  for (long step = start; step < work.steps; ++step) {
    if (work.targets != nullptr && step % 96 == 0) {
      folds.set_fold_target(posture_at_dual<N>(work.targets, work.targets_dot, step / 96, work.ticks, lanes));
    }
    const double* wave_row = work.waves_dot + step * static_cast<long>(lanes.directions);
    const double* area_row = work.area0_dot + step * static_cast<long>(lanes.directions);
    Dual<N> wave = load_dual<N>(work.waves[step], wave_row, lanes);
    Dual<N> area = load_dual<N>(work.area0[step], area_row, lanes);
    Dual<N> out = folds.step(wave, area);
    double* out_row = work.flow_dot + step * static_cast<long>(lanes.directions) + lanes.first;
    for (int k = 0; k < lanes.width; ++k) {
      out_row[k] = out.d[k];
    }
  }
}

// A thread's blocks, each at the narrowest instantiated width that holds it.
void run_blocks(const FoldsWork* work) {
  for (int block = work->offset; block < work->blocks; block += work->stride) {
    Lanes lanes;
    lanes.directions = work->directions;
    lanes.first = work->block_starts[block];
    lanes.width = work->block_starts[block + 1] - lanes.first;
    if (lanes.width <= 0) {
      continue;
    }
    if (lanes.width <= 8) {
      run_block<8>(*work, lanes);
    } else if (lanes.width <= 16) {
      run_block<16>(*work, lanes);
    } else if (lanes.width <= 32) {
      run_block<32>(*work, lanes);
    } else {
      run_block<48>(*work, lanes);
    }
  }
}

void run_threads(const FoldsWork& shared, int threads) {
  if (threads <= 1) {
    run_blocks(&shared);
    return;
  }
  std::vector<FoldsWork> shares(static_cast<std::size_t>(threads), shared);
  std::vector<std::thread> workers;
  for (int thread = 0; thread < threads; ++thread) {
    shares[static_cast<std::size_t>(thread)].offset = thread;
    shares[static_cast<std::size_t>(thread)].stride = threads;
    workers.emplace_back(run_blocks, &shares[static_cast<std::size_t>(thread)]);
  }
  for (std::thread& worker : workers) {
    worker.join();
  }
}

}  // namespace

extern "C" {

// The forced folds' flow per step (sustain RAW's first column under
// FORCEWAVE): waves and area0 per step; targets (ticks x 3: lower rest,
// upper rest, chink) for sustain FOLDKEYS or null for the modal posture.
int ff_primal(long steps, int vocal_count, double pitch, double lung, const double* waves, const double* area0,
              const double* targets, long ticks, double* flow) {
  ForcedFolds<double> folds(vocal_count, pitch, lung);
  for (long step = 0; step < steps; ++step) {
    if (targets != nullptr && step % 96 == 0) {
      folds.set_fold_target(posture_at<double>(targets, step / 96, ticks));
    }
    flow[step] = folds.step(waves[step], area0[step]);
  }
  return 0;
}

// ff_primal with the ear's gain set (diagnosis: 0 holds the tension open
// loop, so whole-step closure timing cannot reach the flow).
int ff_primal_feedback(long steps, int vocal_count, double pitch, double lung, const double* waves, const double* area0,
                       double pitch_feedback, double* flow) {
  ForcedFolds<double> folds(vocal_count, pitch, lung);
  folds.set_pitch_feedback(pitch_feedback);
  for (long step = 0; step < steps; ++step) {
    flow[step] = folds.step(waves[step], area0[step]);
  }
  return 0;
}

// Diagnosis: ff_primal_feedback with each step's ForcedFolds::trace (13
// values) into trace [step][13].
int ff_trace(long steps, int vocal_count, double pitch, double lung, const double* waves, const double* area0,
             double pitch_feedback, double* flow, double* trace) {
  ForcedFolds<double> folds(vocal_count, pitch, lung);
  folds.set_pitch_feedback(pitch_feedback);
  for (long step = 0; step < steps; ++step) {
    flow[step] = folds.step(waves[step], area0[step]);
    folds.trace(trace + 13 * step);
  }
  return 0;
}

// The flow and its derivatives along `directions` tangents of the inputs
// (waves_dot, area0_dot: [step][direction]; targets_dot [tick][3][direction]
// or null), the directions in `blocks` blocks (block_starts: blocks + 1
// ascending starts, each block at most 48 wide) shared among `threads`
// threads. flow_dot is zeroed here first.
int ff_tangent(long steps, int vocal_count, double pitch, double lung, const double* waves, const double* waves_dot,
               const double* area0, const double* area0_dot, const double* targets, const double* targets_dot,
               long ticks, int directions, const int* block_starts, int blocks, int threads, double* flow,
               double* flow_dot) {
  for (int block = 0; block < blocks; ++block) {
    if (block_starts[block + 1] - block_starts[block] > 48 || block_starts[block + 1] < block_starts[block]) {
      return 1;
    }
  }
  // The primal pass, with its checkpoints.
  std::vector<ForcedFolds<double>> checkpoints;
  ForcedFolds<double> folds(vocal_count, pitch, lung);
  for (long step = 0; step < steps; ++step) {
    if (step % kCheckpointSteps == 0) {
      checkpoints.push_back(folds);
    }
    if (targets != nullptr && step % 96 == 0) {
      folds.set_fold_target(posture_at<double>(targets, step / 96, ticks));
    }
    flow[step] = folds.step(waves[step], area0[step]);
  }
  std::memset(flow_dot, 0, sizeof(double) * static_cast<std::size_t>(steps) * static_cast<std::size_t>(directions));
  FoldsWork work;
  work.steps = steps;
  work.vocal_count = vocal_count;
  work.pitch = pitch;
  work.lung = lung;
  work.waves = waves;
  work.waves_dot = waves_dot;
  work.area0 = area0;
  work.area0_dot = area0_dot;
  work.targets = targets;
  work.targets_dot = targets_dot;
  work.ticks = ticks;
  work.directions = directions;
  work.block_starts = block_starts;
  work.blocks = blocks;
  work.checkpoints = &checkpoints;
  work.flow_dot = flow_dot;
  run_threads(work, threads);
  return 0;
}

}  // extern "C"
