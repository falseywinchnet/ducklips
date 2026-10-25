#pragma once

// The planner's parameter tree: the timing, glottal, breath, release and
// articulation settings of PlanOptions, the consonant shapes
// (ConsonantShapes) and a duration factor per phone, as a block file
// (voice/block_file.hpp):
//
//   timing {
//     vowel_hold 0.6
//   }
//   articulation {
//     sibilant_x 7.05
//   }
//   duration {
//     CH 1
//   }
//
// `say --plan-params FILE` reads one over the defaults (and the flags before
// it); `say --plan-params-out FILE` writes the tree as planned. Every path is
// declared once in plan_parameters.cpp; an undeclared one is an error.

#include "voice/block_file.hpp"
#include "voice/planner.hpp"

#include <string>
#include <vector>

namespace voice {

// Sets every path in `entries`; false, naming it in `error`, at the first
// path the planner does not declare.
bool apply_plan_parameters(const std::vector<ParameterEntry>& entries, PlanOptions& options, std::string& error);
// The whole declared tree with the values the options hold now.
void collect_plan_parameters(const PlanOptions& options, std::vector<ParameterEntry>& entries);

}  // namespace voice
