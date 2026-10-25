// The physical tongue alone: does it hold still at rest, and does it settle
// on a held target? Prints the distance from the target (rms and largest,
// cm), the feedforward's own reach estimate and the area every 20 ms.
// usage: tongue_probe [key=value ...]
//   rx ry rtipx rtipy: the first posture held; tx ty tipx tipy jaw: the
//   target posture (default /aa/-like); every: print interval (ms);
//   gain integral force act pliancy: TongueSettings; seconds.

#include "tract/articulation.hpp"
#include "tract/tongue.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {

double value_of(int count, char** values, const char* key, double fallback) {
  std::size_t length = std::strlen(key);
  for (int index = 1; index < count; ++index) {
    if (std::strncmp(values[index], key, length) == 0 && values[index][length] == '=') {
      return std::strtod(values[index] + length + 1, nullptr);
    }
  }
  return fallback;
}

}  // namespace

int main(int count, char** values) {
  tract::MidsagittalTract anatomy;
  tract::TractProfile profile;
  tract::Articulators rest;
  rest.jaw_angle = value_of(count, values, "jaw", 6.0);
  rest.tongue_x = value_of(count, values, "rx", 4.0);
  rest.tongue_y = value_of(count, values, "ry", 6.6);
  rest.tip_x = value_of(count, values, "rtipx", 6.9);
  rest.tip_y = value_of(count, values, "rtipy", 7.4);
  anatomy.compute(rest, profile);
  tract::TongueContour relaxed;
  anatomy.tongue_outline(rest, relaxed);
  tract::TongueSettings settings;
  settings.gain = value_of(count, values, "gain", 3.0);
  settings.integral_rate = value_of(count, values, "integral", 10.0);
  settings.muscle_force = value_of(count, values, "force", 4000.0);
  settings.activation_time = value_of(count, values, "act", 0.015);
  settings.tip_pliancy = value_of(count, values, "pliancy", 0.2);
  settings.support = value_of(count, values, "support", settings.support);
  settings.elasticity = value_of(count, values, "elasticity", settings.elasticity);
  settings.bending = value_of(count, values, "bending", settings.bending);
  settings.section_mass = value_of(count, values, "mass", settings.section_mass);
  settings.tip_mass = value_of(count, values, "tipmass", settings.tip_mass);
  settings.damping_gain = value_of(count, values, "damping", settings.damping_gain);
  settings.tip_stretch = value_of(count, values, "stretch", settings.tip_stretch);
  settings.statics = value_of(count, values, "statics", settings.statics);
  settings.contact_stiffness = value_of(count, values, "soft", settings.contact_stiffness);
  settings.contact_loss = value_of(count, values, "softloss", settings.contact_loss);
  tract::Tongue tongue;
  tongue.reset(relaxed, settings);
  tongue.place(relaxed);
  tract::Articulators aimed = rest;
  aimed.tongue_x = value_of(count, values, "tx", 3.37);
  aimed.tongue_y = value_of(count, values, "ty", 6.66);
  aimed.tip_x = value_of(count, values, "tipx", 6.71);
  aimed.tip_y = value_of(count, values, "tipy", 6.63);
  tract::TongueContour target;
  anatomy.tongue_outline(aimed, target);
  double seconds = value_of(count, values, "seconds", 1.0);
  int frames = static_cast<int>(seconds / 1.0e-3);
  int every = static_cast<int>(value_of(count, values, "every", 20.0));
  for (int phase = 0; phase < 2; ++phase) {
    std::printf("%s\n", phase == 0 ? "== holding the rest" : "== holding the target");
    for (int frame = 0; frame < frames; ++frame) {
      tongue.set_target(phase == 0 ? relaxed : target);
      tract::Articulators held = rest;
      held.contour = &tongue.contour();
      anatomy.compute(held, profile);
      tongue.step(1.0e-3, rest, anatomy);
      if (frame % every == 0) {
        double rms = 0.0;
        double largest = 0.0;
        tongue.target_error(rms, largest);
        double narrowest = 1.0e9;
        double where = 0.0;
        for (int cut = 0; cut < profile.count; ++cut) {
          std::size_t c = static_cast<std::size_t>(cut);
          if (profile.arc[c] > 10.5 && profile.arc[c] < profile.teeth_arc + 0.5 && profile.area[c] < narrowest) {
            narrowest = profile.area[c];
            where = profile.arc[c];
          }
        }
        std::size_t tip = static_cast<std::size_t>(tongue.contour().count - 1);
        std::printf("t %.3f err %.3f largest %.3f reach %.3f  narrowest %.4f cm2 at %.1f cm  tip %.2f,%.2f (aim %.2f,%.2f)\n",
                    frame * 1.0e-3, rms, largest, tongue.reach_error(), narrowest, where, tongue.contour().x[tip],
                    tongue.contour().y[tip], target.x[tip], target.y[tip]);
        std::printf("    superior longitudinal %.3f (feedforward %.3f, integral %.3f)  inferior %.3f (%.3f, %.3f)\n",
                    tongue.activation(7), tongue.feedforward(7), tongue.integral(7), tongue.activation(8),
                    tongue.feedforward(8), tongue.integral(8));
      }
    }
  }
  return 0;
}
