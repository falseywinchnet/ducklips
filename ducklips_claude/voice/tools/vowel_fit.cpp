// Finds tongue and lip postures whose closed-glottis tract resonances match
// measured vowel formants (Hillenbrand et al. 1995, adult male means). The
// body itself is the model: each candidate posture is rendered and its
// impulse response read for peaks.
//
// usage: vowel_fit [VOWEL] | bounded [VOWEL] | front [VOWEL [F1 F2 F3]] | rhotic

#include "voice/body.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

struct VowelTarget {
  const char* name;
  double f1;
  double f2;
  double f3;
  double protrusion;
  double place;
  double area;
  double lips;
};

// Starting postures are rough guesses; the fit moves them.
constexpr VowelTarget kTargets[] = {
    {"IY", 329, 2387, 3125, 0.0, 13.0, 0.4, 3.0}, {"IH", 434, 2128, 2854, 0.0, 12.5, 0.9, 3.0},
    {"EY", 471, 2168, 2848, 0.0, 12.5, 0.8, 3.0}, {"EH", 594, 1785, 2883, 0.0, 12.0, 1.6, 4.0},
    {"AE", 581, 1885, 2759, 0.0, 12.0, 2.0, 5.0}, {"AA", 761, 1361, 3000, 0.0, 5.5, 0.7, 5.0},
    {"AO", 600, 900, 2450, 0.6, 6.0, 0.6, 2.0},   {"OW", 515, 1500, 2600, 0.8, 8.5, 0.7, 0.5},
    {"UH", 469, 1122, 2434, 0.3, 9.5, 0.9, 1.5},  {"UW", 347, 1800, 2600, 0.8, 10.5, 0.5, 0.6},
    {"AH", 586, 1588, 2700, 0.0, 7.0, 1.0, 3.0},  {"AX", 500, 1600, 2600, 0.0, 9.5, 2.5, 3.0},
    {"ER", 474, 1379, 1710, 0.3, 11.5, 0.7, 1.5},
};

struct Resonances {
  std::array<double, 4> frequency{};
  int count = 0;
};

struct Extra {
  bool active = false;
  double dorsum_place = 13.0;
  double dorsum_area = 0.5;
  double floor = 0.5;
};

Extra g_extra;

Resonances measure(double place, double area, double lips, double protrusion) {
  voice::BodyOptions options;
  options.hold_noise = 0.0;
  options.aspiration = 0.0;
  options.frication = 0.0;
  options.skin = 0.0;
  voice::Body body(options);
  voice::PartEvent event;
  event.part = voice::Part::tongue;
  event.a = static_cast<float>(place);
  event.b = static_cast<float>(area);
  body.apply(event);
  event.part = voice::Part::lips;
  event.a = static_cast<float>(lips);
  event.b = static_cast<float>(protrusion);
  body.apply(event);
  // The glottis slightly apart, as averaged over a voiced cycle: fully
  // closed, the subglottal tract couples through the posterior chink and
  // splits F1 into a peak and a broad bump that the peak picker mistakes
  // for F1 (EH read 580 Hz where voiced speech gives 440-460).
  event.part = voice::Part::glottis;
  event.a = 0.01F;
  event.b = 1.0F;
  body.apply(event);
  if (g_extra.active) {
    event.part = voice::Part::dorsum;
    event.a = static_cast<float>(g_extra.dorsum_place);
    event.b = static_cast<float>(g_extra.dorsum_area);
    body.apply(event);
    event.part = voice::Part::floor;
    event.a = static_cast<float>(g_extra.floor);
    body.apply(event);
  }
  body.settle();
  std::vector<float> response(6000, 0.0F);
  body.render(response.data(), 480);
  body.kick(1000.0);
  body.render(response.data(), 6000);
  // Spectrum from 150 Hz to 4 kHz on a 10 Hz grid, skipping the direct pulse,
  // with the radiation slope (+6 dB/octave) removed: volume-velocity peaks.
  constexpr int kBins = 386;
  std::array<double, kBins> level{};
  for (int bin = 0; bin < kBins; ++bin) {
    double frequency = 150.0 + 10.0 * bin;
    double omega = 2.0 * voice::kPi * frequency / voice::kAudioRate;
    double real = 0.0;
    double imaginary = 0.0;
    for (int index = 72; index < 6000; ++index) {
      double value = response[static_cast<std::size_t>(index)];
      real += value * std::cos(omega * index);
      imaginary -= value * std::sin(omega * index);
    }
    level[bin] = (real * real + imaginary * imaginary) / (frequency * frequency);
  }
  Resonances result;
  for (int bin = 2; bin < kBins - 2 && result.count < 4; ++bin) {
    if (level[bin] > level[bin - 1] && level[bin] >= level[bin + 1] && level[bin] > level[bin - 2] &&
        level[bin] >= level[bin + 2]) {
      result.frequency[static_cast<std::size_t>(result.count)] = 150.0 + 10.0 * bin;
      ++result.count;
    }
  }
  return result;
}

double cost(const VowelTarget& target, double place, double area, double lips) {
  Resonances found = measure(place, area, lips, target.protrusion);
  if (found.count < 3) {
    return 100.0;
  }
  double e1 = std::log(found.frequency[0] / target.f1);
  double e2 = std::log(found.frequency[1] / target.f2);
  double e3 = std::log(found.frequency[2] / target.f3);
  return e1 * e1 + e2 * e2 + 0.5 * e3 * e3;
}

void fit(const VowelTarget& target) {
  double place = target.place;
  double log_area = std::log(target.area);
  double log_lips = std::log(target.lips);
  double best = cost(target, place, std::exp(log_area), std::exp(log_lips));
  // Coarse multi-start: keep the best of a grid if it beats the guess.
  const double start_places[5] = {3.5, 6.0, 8.5, 11.0, 13.0};
  const double start_areas[3] = {0.4, 1.0, 2.5};
  const double start_lips[3] = {0.3, 1.5, 5.0};
  for (double start_place : start_places) {
    for (double start_area : start_areas) {
      for (double start_lip : start_lips) {
        double trial = cost(target, start_place, start_area, start_lip);
        if (trial < best) {
          best = trial;
          place = start_place;
          log_area = std::log(start_area);
          log_lips = std::log(start_lip);
        }
      }
    }
  }
  double steps[3] = {1.0, 0.3, 0.3};
  for (int round = 0; round < 40; ++round) {
    bool improved = false;
    for (int axis = 0; axis < 3; ++axis) {
      for (int sign = -1; sign <= 1; sign += 2) {
        double trial_place = place;
        double trial_area = log_area;
        double trial_lips = log_lips;
        if (axis == 0) {
          trial_place = std::clamp(place + sign * steps[0], 3.0, 15.0);
        } else if (axis == 1) {
          trial_area = std::clamp(log_area + sign * steps[1], std::log(0.15), std::log(6.0));
        } else {
          trial_lips = std::clamp(log_lips + sign * steps[2], std::log(0.15), std::log(10.0));
        }
        double trial = cost(target, trial_place, std::exp(trial_area), std::exp(trial_lips));
        if (trial < best) {
          best = trial;
          place = trial_place;
          log_area = trial_area;
          log_lips = trial_lips;
          improved = true;
        }
      }
    }
    if (!improved) {
      steps[0] *= 0.5;
      steps[1] *= 0.5;
      steps[2] *= 0.5;
      if (steps[0] < 0.05) {
        break;
      }
    }
  }
  Resonances found = measure(place, std::exp(log_area), std::exp(log_lips), target.protrusion);
  std::printf("{\"%s\", %.2f, %.3f, %.3f, %.2f},  // %.0f %.0f %.0f -> %.0f %.0f %.0f %.0f  cost %.4f\n",
              target.name, place, std::exp(log_area), std::exp(log_lips), target.protrusion, target.f1,
              target.f2, target.f3, found.frequency[0], found.frequency[1], found.frequency[2],
              found.frequency[3], best);
  std::fflush(stdout);
}

}  // namespace

// Rhotic: a bunched dorsum near the palate, the tongue root drawn back into
// the pharynx, rounded lips and a deep space under the raised blade.
void fit_rhotic() {
  VowelTarget target = {"ER", 474, 1379, 1710, 0.3, 5.0, 1.0, 1.0};
  g_extra.active = true;
  // The tuned inventory's ER, for reference.
  g_extra.dorsum_place = 14.0;
  g_extra.dorsum_area = 0.15;
  g_extra.floor = 4.0;
  target.protrusion = 0.3;
  Resonances current = measure(3.5, 1.2, 1.0, 0.3);
  std::printf("current ER -> %.0f %.0f %.0f %.0f\n", current.frequency[0], current.frequency[1],
              current.frequency[2], current.frequency[3]);
  double best = 1.0e9;
  double best_values[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double dorsum_places[6] = {12.0, 13.0, 14.0, 14.5, 15.0, 15.5};
  const double dorsum_areas[5] = {0.08, 0.15, 0.25, 0.4, 0.7};
  const double floors[3] = {1.0, 3.0, 6.0};
  const double root_places[3] = {3.5, 5.0, 6.5};
  const double root_areas[3] = {0.4, 0.9, 2.0};
  const double lip_areas[4] = {0.3, 0.6, 1.2, 2.4};
  const double protrusions[2] = {0.3, 0.7};
  for (double dorsum_place : dorsum_places) {
    for (double dorsum_area : dorsum_areas) {
      for (double floor : floors) {
        for (double protrusion : protrusions) {
          g_extra.dorsum_place = dorsum_place;
          g_extra.dorsum_area = dorsum_area;
          g_extra.floor = floor;
          target.protrusion = protrusion;
          for (double root_place : root_places) {
            for (double root_area : root_areas) {
              for (double lip_area : lip_areas) {
                double trial = cost(target, root_place, root_area, lip_area);
                if (trial < best) {
                  best = trial;
                  best_values[0] = dorsum_place;
                  best_values[1] = dorsum_area;
                  best_values[2] = root_place;
                  best_values[3] = root_area;
                  best_values[4] = lip_area;
                  best_values[5] = floor;
                  best_values[6] = protrusion;
                  std::printf("  dorsum %.1f %.2f floor %.1f root %.1f %.2f lips %.2f protrusion %.1f cost %.4f\n",
                              dorsum_place, dorsum_area, floor, root_place, root_area, lip_area, protrusion, best);
                  std::fflush(stdout);
                }
              }
            }
          }
        }
      }
    }
  }
  g_extra.dorsum_place = best_values[0];
  g_extra.dorsum_area = best_values[1];
  g_extra.floor = best_values[5];
  Resonances found = measure(best_values[2], best_values[3], best_values[4], best_values[6]);
  std::printf("ER dorsum %.1f %.2f floor %.1f root %.1f %.2f lips %.2f protrusion %.1f -> %.0f %.0f %.0f %.0f cost %.4f\n",
              best_values[0], best_values[1], best_values[5], best_values[2], best_values[3], best_values[4],
              best_values[6], found.frequency[0], found.frequency[1], found.frequency[2], found.frequency[3], best);
}

// Lip openings held to anthropometric ranges (male, citation vowels): wide
// openings radiate F3 and F4 away (bandwidths near 250 Hz at 9 cm^2).
// Targets: Hillenbrand 1995 men; UW toward the fronted American /u/.
struct BoundedTarget {
  VowelTarget target;
  double lip_low;
  double lip_high;
};

constexpr BoundedTarget kBounded[] = {
    {{"IY", 342, 2322, 3000, 0.0, 13.0, 0.4, 2.5}, 2.0, 3.5},
    {{"IH", 427, 2034, 2684, 0.0, 12.5, 0.9, 3.0}, 2.5, 4.5},
    {{"EY", 476, 2089, 2691, 0.0, 12.5, 0.8, 3.0}, 2.5, 5.0},
    {{"EH", 580, 1799, 2605, 0.0, 11.0, 1.4, 4.0}, 3.0, 6.0},
    {{"AE", 588, 1952, 2601, 0.0, 11.0, 1.8, 5.0}, 4.0, 7.5},
    {{"AA", 768, 1333, 2522, 0.0, 4.5, 1.0, 5.5}, 4.5, 7.0},
    {{"AO", 652, 997, 2538, 0.3, 5.5, 0.5, 2.0}, 1.5, 3.0},
    {{"OW", 497, 910, 2459, 0.6, 6.5, 0.7, 1.2}, 0.8, 1.6},
    {{"UH", 469, 1122, 2434, 0.4, 8.5, 1.0, 1.5}, 1.0, 2.0},
    {{"UW", 360, 1500, 2300, 0.6, 13.0, 0.4, 0.3}, 0.15, 0.6},
    {{"AH", 623, 1200, 2550, 0.0, 7.0, 0.8, 4.0}, 3.0, 5.0},
    {{"AX", 500, 1500, 2500, 0.0, 9.5, 1.5, 3.0}, 2.0, 4.0},
};

void fit_bounded(const BoundedTarget& bounded) {
  const VowelTarget& target = bounded.target;
  double low_lips = std::log(bounded.lip_low);
  double high_lips = std::log(bounded.lip_high);
  double place = target.place;
  double log_area = std::log(target.area);
  double log_lips = std::clamp(std::log(target.lips), low_lips, high_lips);
  double best = cost(target, place, std::exp(log_area), std::exp(log_lips));
  const double start_places[9] = {3.5, 5.0, 6.5, 8.0, 9.5, 11.0, 12.0, 13.0, 14.0};
  const double start_areas[4] = {0.25, 0.6, 1.2, 2.5};
  const double start_lips[2] = {bounded.lip_low, bounded.lip_high};
  for (double start_place : start_places) {
    for (double start_area : start_areas) {
      for (double start_lip : start_lips) {
        double trial = cost(target, start_place, start_area, start_lip);
        if (trial < best) {
          best = trial;
          place = start_place;
          log_area = std::log(start_area);
          log_lips = std::log(start_lip);
        }
      }
    }
  }
  double steps[3] = {1.0, 0.3, 0.2};
  for (int round = 0; round < 60; ++round) {
    bool improved = false;
    for (int axis = 0; axis < 3; ++axis) {
      for (int sign = -1; sign <= 1; sign += 2) {
        double trial_place = place;
        double trial_area = log_area;
        double trial_lips = log_lips;
        if (axis == 0) {
          trial_place = std::clamp(place + sign * steps[0], 3.0, 15.0);
        } else if (axis == 1) {
          trial_area = std::clamp(log_area + sign * steps[1], std::log(0.15), std::log(6.0));
        } else {
          trial_lips = std::clamp(log_lips + sign * steps[2], low_lips, high_lips);
        }
        double trial = cost(target, trial_place, std::exp(trial_area), std::exp(trial_lips));
        if (trial < best) {
          best = trial;
          place = trial_place;
          log_area = trial_area;
          log_lips = trial_lips;
          improved = true;
        }
      }
    }
    if (!improved) {
      steps[0] *= 0.5;
      steps[1] *= 0.5;
      steps[2] *= 0.5;
      if (steps[0] < 0.03) {
        break;
      }
    }
  }
  Resonances found = measure(place, std::exp(log_area), std::exp(log_lips), target.protrusion);
  std::printf("constexpr Posture k%s = {%.2f, %.3f, %.3f, %.1f};  // %.0f %.0f %.0f -> %.0f %.0f %.0f  cost %.4f\n",
              target.name, place, std::exp(log_area), std::exp(log_lips), target.protrusion, target.f1, target.f2,
              target.f3, found.frequency[0], found.frequency[1], found.frequency[2], best);
  std::fflush(stdout);
}

// Front vowels with a raised front dorsum as a second shaping point, lips
// within their anthropometric range.
void fit_front(const BoundedTarget& bounded) {
  const VowelTarget& target = bounded.target;
  g_extra.active = true;
  g_extra.floor = 0.5;
  double best = 1.0e9;
  double values[5] = {};
  const double dorsum_places[8] = {11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0};
  const double dorsum_areas[7] = {0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2};
  const double places[5] = {5.0, 6.5, 8.0, 9.5, 11.0};
  const double areas[4] = {0.8, 1.4, 2.2, 3.5};
  const double lip_areas[3] = {bounded.lip_low, std::sqrt(bounded.lip_low * bounded.lip_high), bounded.lip_high};
  for (double dorsum_place : dorsum_places) {
    for (double dorsum_area : dorsum_areas) {
      for (double place : places) {
        for (double area : areas) {
          for (double lip : lip_areas) {
            g_extra.dorsum_place = dorsum_place;
            g_extra.dorsum_area = dorsum_area;
            double trial = cost(target, place, area, lip);
            if (trial < best) {
              best = trial;
              values[0] = dorsum_place;
              values[1] = dorsum_area;
              values[2] = place;
              values[3] = area;
              values[4] = lip;
            }
          }
        }
      }
    }
  }
  g_extra.dorsum_place = values[0];
  g_extra.dorsum_area = values[1];
  Resonances found = measure(values[2], values[3], values[4], 0.0);
  std::printf("%s dorsum %.1f %.2f tongue %.1f %.2f lips %.2f -> %.0f %.0f %.0f cost %.4f\n", target.name, values[0],
              values[1], values[2], values[3], values[4], found.frequency[0], found.frequency[1], found.frequency[2],
              best);
}

int main(int argument_count, char** arguments) {
  if (argument_count > 1 && std::strcmp(arguments[1], "front") == 0) {
    for (const BoundedTarget& bounded : kBounded) {
      if (argument_count > 2 && std::strcmp(arguments[2], bounded.target.name) != 0) {
        continue;
      }
      // vowel_fit front NAME F1 F2 F3: other formant targets for that vowel.
      BoundedTarget chosen = bounded;
      if (argument_count > 5) {
        chosen.target.f1 = std::strtod(arguments[3], nullptr);
        chosen.target.f2 = std::strtod(arguments[4], nullptr);
        chosen.target.f3 = std::strtod(arguments[5], nullptr);
      }
      fit_front(chosen);
    }
    return 0;
  }
  if (argument_count > 1 && std::strcmp(arguments[1], "bounded") == 0) {
    for (const BoundedTarget& bounded : kBounded) {
      if (argument_count > 2 && std::strcmp(arguments[2], bounded.target.name) != 0) {
        continue;
      }
      fit_bounded(bounded);
    }
    return 0;
  }
  if (argument_count > 1 && std::strcmp(arguments[1], "rhotic") == 0) {
    fit_rhotic();
    return 0;
  }
  for (const VowelTarget& target : kTargets) {
    if (argument_count > 1 && std::strcmp(arguments[1], target.name) != 0) {
      continue;
    }
    fit(target);
  }
  return 0;
}
