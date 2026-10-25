#pragma once

// The tongue's roles in the mouth as an instrument. Each sound asks a
// different job of each piece of the tongue (root, dorsum, blade, tip):
//   shaper - the resonator: it sets the cavities' areas, slowly and stably
//            (vowels, glides, /r/, /h/, rest);
//   valve  - the occluder: it strikes the hull, seals and releases
//            (/t d n l/, /th dh/, the flap with the tip and blade; /k g ng/
//            with the dorsum);
//   nozzle - the jet's aperture: it holds a precise, still channel against
//            the ridge and teeth (/s z sh zh/, /ch jh/).
// Each role is represented either by the body (the physical tongue's
// sections, moved by its muscles) or pinned to the geometric tongue's outline
// (the planner's shape as the motors hold it), as the tongue's settings
// choose (TongueSettings shaper_pinned, valve_pinned, nozzle_pinned). Around
// each phone the role comes into force over role_ramp, role_lead before the
// phone, and leaves role_hold after it.

#include "tract/tongue.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace tract {

enum class TongueRole : std::uint8_t { shaper, valve, nozzle };

constexpr int kTonguePieces = 4;  // root, dorsum, blade, tip

struct PhoneSpan {
  double start = 0.0;
  double duration = 0.0;
  std::array<TongueRole, kTonguePieces> roles{};
};

// The roles a phone (ARPAbet; stress digits and a "-1" suffix are ignored)
// gives each piece.
void tongue_roles(const std::string& phone, std::array<TongueRole, kTonguePieces>& roles);

// A score's phone comment ("# NAME start duration pitch", the stress in the name).
bool phone_span_from_comment(const std::string& line, PhoneSpan& span);

// Each piece's representation at `time`: how strongly it is pinned to the
// geometric outline (0 the body, 1 pinned), and the body's stiffness toward
// its target and damping ratio, blended between the roles in force.
struct RoleState {
  std::array<double, kTonguePieces> pins{};
  std::array<double, kTonguePieces> stiffness{};
  std::array<double, kTonguePieces> damping{};
  std::array<double, kTonguePieces> shaping{};  // the share of the resonator's role
};

// `cursor` advances through the time-ordered phones.
void role_state(const std::vector<PhoneSpan>& phones, std::size_t& cursor, double time,
                const TongueSettings& settings, RoleState& state);

}  // namespace tract
