// The midsagittal articulatory model, one shape per input line.
// usage: geometry_probe [outline]   (reads lines of key=value from stdin)
//   keys: jaw tx ty tipx tipy ap pro lw drop lat groove, and area rules as
//   rule.REGION=alpha,beta (REGION: larynx pharynx velar palatal alveolar
//   teeth), anat.scale=s (the speaker's size), anat.root=k (root follows body).
// For each line it prints "length L teeth T lips P velum V tip Q" and then
// one "cut arc area width region x y" line per cut; with `outline`, also the
// hull and organ outlines ("hull x y", "organ x y"). Lines end with "end".

#include "tract/articulation.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

namespace {

const char* region_name(tract::CutRegion region) {
  switch (region) {
    case tract::CutRegion::larynx:
      return "larynx";
    case tract::CutRegion::pharynx:
      return "pharynx";
    case tract::CutRegion::velar:
      return "velar";
    case tract::CutRegion::palatal:
      return "palatal";
    case tract::CutRegion::alveolar:
      return "alveolar";
    case tract::CutRegion::teeth:
      return "teeth";
    case tract::CutRegion::lips:
      return "lips";
  }
  return "?";
}

bool region_of(const std::string& name, tract::CutRegion& region) {
  if (name == "larynx") {
    region = tract::CutRegion::larynx;
  } else if (name == "pharynx") {
    region = tract::CutRegion::pharynx;
  } else if (name == "velar") {
    region = tract::CutRegion::velar;
  } else if (name == "palatal") {
    region = tract::CutRegion::palatal;
  } else if (name == "alveolar") {
    region = tract::CutRegion::alveolar;
  } else if (name == "teeth") {
    region = tract::CutRegion::teeth;
  } else {
    return false;
  }
  return true;
}

void apply(const std::string& key, const std::string& value, tract::Articulators& organs,
           tract::MidsagittalTract& model) {
  double number = std::strtod(value.c_str(), nullptr);
  if (key == "jaw") {
    organs.jaw_angle = number;
  } else if (key == "tx") {
    organs.tongue_x = number;
  } else if (key == "ty") {
    organs.tongue_y = number;
  } else if (key == "tipx") {
    organs.tip_x = number;
  } else if (key == "tipy") {
    organs.tip_y = number;
  } else if (key == "ap") {
    organs.lip_aperture = number;
  } else if (key == "pro") {
    organs.lip_protrusion = number;
  } else if (key == "lw") {
    organs.lip_width = number;
  } else if (key == "drop") {
    organs.larynx_drop = number;
  } else if (key == "lat") {
    organs.lateral = number;
  } else if (key == "groove") {
    organs.groove = number;
  } else if (key == "anat.scale") {
    model.set_scale(number);
  } else if (key == "anat.root") {
    model.set_root_follow(number);
  } else if (key.rfind("rule.", 0) == 0) {
    tract::CutRegion region = tract::CutRegion::pharynx;
    std::size_t comma = value.find(',');
    if (region_of(key.substr(5), region) && comma != std::string::npos) {
      tract::AreaRule rule;
      rule.alpha = std::strtod(value.c_str(), nullptr);
      rule.beta = std::strtod(value.c_str() + comma + 1, nullptr);
      model.set_area_rule(region, rule);
    }
  }
}

}  // namespace

int main(int count, char** values) {
  bool outline = count > 1 && std::strcmp(values[1], "outline") == 0;
  tract::MidsagittalTract model;
  tract::TractProfile profile;
  std::string line;
  while (std::getline(std::cin, line)) {
    tract::Articulators organs;
    std::istringstream fields(line);
    std::string field;
    while (fields >> field) {
      std::size_t equals = field.find('=');
      if (equals != std::string::npos) {
        apply(field.substr(0, equals), field.substr(equals + 1), organs, model);
      }
    }
    model.compute(organs, profile);
    std::printf("length %.3f teeth %.3f lips %.3f velum %.3f tip %.3f\n", profile.length, profile.teeth_arc,
                profile.lip_arc, profile.velum_arc, profile.tip_arc);
    for (int cut = 0; cut < profile.count; ++cut) {
      std::size_t index = static_cast<std::size_t>(cut);
      std::printf("cut %.3f %.4f %.4f %s %.3f %.3f\n", profile.arc[index], profile.area[index], profile.width[index],
                  region_name(profile.region[index]), profile.mid_x[index], profile.mid_y[index]);
    }
    if (outline) {
      for (int point = 0; point < model.hull_points(); ++point) {
        std::printf("hull %.3f %.3f\n", model.hull_x(point), model.hull_y(point));
      }
      for (int point = 0; point < model.organ_points(); ++point) {
        std::printf("organ %.3f %.3f\n", model.organ_x(point), model.organ_y(point));
      }
    }
    std::printf("end\n");
    std::fflush(stdout);
  }
  return 0;
}
