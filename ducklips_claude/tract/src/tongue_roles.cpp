#include "tract/tongue_roles.hpp"

#include <algorithm>
#include <cstdio>

namespace tract {

namespace {

constexpr int kDorsum = 1;
constexpr int kBlade = 2;
constexpr int kTip = 3;

std::string base_phone(const std::string& phone) {
  std::string base = phone;
  std::size_t dash = base.find('-');
  if (dash != std::string::npos) {
    base.erase(dash);
  }
  while (!base.empty() && base.back() >= '0' && base.back() <= '9') {
    base.pop_back();
  }
  return base;
}

bool one_of(const std::string& phone, const char* const* names, int count) {
  for (int index = 0; index < count; ++index) {
    if (phone == names[index]) {
      return true;
    }
  }
  return false;
}

// A role's representation: pinned share, stiffness and damping ratio.
void represent(TongueRole role, const TongueSettings& settings, double& pinned, double& stiffness, double& damping) {
  if (role == TongueRole::valve) {
    pinned = settings.valve_pinned;
    stiffness = settings.valve_stiffness;
    damping = settings.valve_damping;
  } else if (role == TongueRole::nozzle) {
    pinned = settings.nozzle_pinned;
    stiffness = settings.nozzle_stiffness;
    damping = settings.nozzle_damping;
  } else {
    pinned = settings.shaper_pinned;
    stiffness = settings.shaper_stiffness;
    damping = settings.shaper_damping;
  }
}

// 0..1: how far a phone's role is in force at `time` (smoothstep ramps).
double role_presence(const PhoneSpan& phone, double time, const TongueSettings& settings) {
  double ramp = std::max(settings.role_ramp, 1.0e-4);
  double on = phone.start - settings.role_lead;
  double off = phone.start + phone.duration + settings.role_hold;
  double rise = std::clamp((time - (on - ramp)) / ramp, 0.0, 1.0);
  double fall = std::clamp((off + ramp - time) / ramp, 0.0, 1.0);
  double share = std::min(rise, fall);
  return share * share * (3.0 - 2.0 * share);
}

}  // namespace

void tongue_roles(const std::string& phone, std::array<TongueRole, kTonguePieces>& roles) {
  static const char* const kTipValves[] = {"T", "D", "N", "L", "DH", "TH", "DX"};
  static const char* const kDorsumValves[] = {"K", "G", "NG"};
  static const char* const kNozzles[] = {"S", "Z", "SH", "ZH", "CH", "JH"};
  roles.fill(TongueRole::shaper);
  std::string base = base_phone(phone);
  if (one_of(base, kTipValves, 7)) {
    roles[static_cast<std::size_t>(kBlade)] = TongueRole::valve;
    roles[static_cast<std::size_t>(kTip)] = TongueRole::valve;
  } else if (one_of(base, kDorsumValves, 3)) {
    roles[static_cast<std::size_t>(kDorsum)] = TongueRole::valve;
  } else if (one_of(base, kNozzles, 6)) {
    roles[static_cast<std::size_t>(kBlade)] = TongueRole::nozzle;
    roles[static_cast<std::size_t>(kTip)] = TongueRole::nozzle;
  }
}

bool phone_span_from_comment(const std::string& line, PhoneSpan& span) {
  char name[32] = {};
  double start = 0.0;
  double duration = 0.0;
  double pitch = 0.0;
  // The name carries the stress ("EH0", "DH-1"); the header line
  // ("# phones: name stress start ...") fails the numbers.
  if (std::sscanf(line.c_str(), "# %31s %lf %lf %lf", name, &start, &duration, &pitch) != 4) {
    return false;
  }
  span.start = start;
  span.duration = duration;
  tongue_roles(name, span.roles);
  return true;
}

void role_state(const std::vector<PhoneSpan>& phones, std::size_t& cursor, double time,
                const TongueSettings& settings, RoleState& state) {
  double ramp = std::max(settings.role_ramp, 1.0e-4);
  while (cursor < phones.size() &&
         phones[cursor].start + phones[cursor].duration + settings.role_hold + ramp < time) {
    ++cursor;
  }
  std::array<double, kTonguePieces> strongest{};
  std::array<TongueRole, kTonguePieces> role{};
  role.fill(TongueRole::shaper);
  for (std::size_t look = cursor; look < phones.size() && phones[look].start - settings.role_lead - ramp <= time;
       ++look) {
    double presence = role_presence(phones[look], time, settings);
    for (int piece = 0; piece < kTonguePieces; ++piece) {
      std::size_t p = static_cast<std::size_t>(piece);
      if (phones[look].roles[p] != TongueRole::shaper && presence > strongest[p]) {
        strongest[p] = presence;
        role[p] = phones[look].roles[p];
      }
    }
  }
  double base_pinned = 0.0;
  double base_stiffness = 0.0;
  double base_damping = 0.0;
  represent(TongueRole::shaper, settings, base_pinned, base_stiffness, base_damping);
  for (int piece = 0; piece < kTonguePieces; ++piece) {
    std::size_t p = static_cast<std::size_t>(piece);
    double pinned = 0.0;
    double stiffness = 0.0;
    double damping = 0.0;
    represent(role[p], settings, pinned, stiffness, damping);
    double share = strongest[p];
    state.pins[p] = base_pinned + share * (pinned - base_pinned);
    state.stiffness[p] = base_stiffness + share * (stiffness - base_stiffness);
    state.damping[p] = base_damping + share * (damping - base_damping);
    state.shaping[p] = 1.0 - share;
  }
}

}  // namespace tract
