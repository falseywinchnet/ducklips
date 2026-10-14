#pragma once

// The vocal tract's shape from its articulators: one area per oral section
// (glottis to lips). Shared by the body and by offline tools (codebooks,
// inversion) so they all use the same mouth.

#include "voice/body.hpp"

#include <array>

namespace voice {

struct TractShape {
  double tongue_place = 10.75;  // body constriction, cm from glottis
  double tongue_area = 2.5;     // cm^2
  double tip_place = 16.0;
  double tip_root = 2.83;       // signed root area (negative presses)
  double dorsum_place = 12.0;
  double dorsum_root = 2.83;
  double lip_root = 1.27;       // signed root area of the lip opening
  double protrusion = 0.0;      // cm
  double teeth_gap = 8.0;       // cm^2
  double jaw_closure = 0.0;     // 0 open .. 1 raised
  double free_front = 1.0;      // 0..1: how far a raised tip frees the space ahead of it from the body
};

struct TractLayout {
  double lip_start = 0.0;    // cm from glottis where the lips begin
  double teeth_start = 0.0;
  double teeth_end = 0.0;
};

TractLayout tract_areas(const TractShape& shape, std::array<double, kOralSections>& areas);

}  // namespace voice
