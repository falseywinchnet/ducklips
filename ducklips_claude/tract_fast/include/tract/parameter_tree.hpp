#pragma once

// The engine's parameters as a tree of named blocks, read from and written
// to a block file:
//
//   # comment
//   fold {
//     cover_mass 0.09
//     support_stiffness 150000
//   }
//   glottal_noise {
//     aspiration 8
//   }
//
// The file parses to a tree (ParameterNode); its leaves are addressed by
// dotted paths ("fold.cover_mass"). Every path the engine accepts is
// declared once in parameter_tree.cpp, bound to a field of the body's layout
// (applied before the body is built) or to a setting of the speaker
// (applied after). An undeclared path is an error, never ignored: whatever
// writes a file (the behavioural fit, tract/analysis/fit_blocks.py) learns at
// once that a name does not exist.

#include "tract/speaker.hpp"
#include "tract/voice_body.hpp"
#include "voice/block_file.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace tract {

// The block file's tree and leaves (voice/block_file.hpp).
using ParameterNode = voice::ParameterNode;
using ParameterEntry = voice::ParameterEntry;
using voice::flatten_parameters;
using voice::parse_parameters;
using voice::write_parameters;

// Sets the declared layout fields named in `entries` (the speaker's are
// skipped). Returns false, naming the path in `error`, at the first path the
// engine does not declare.
bool apply_layout_parameters(const std::vector<ParameterEntry>& entries, BodyLayout& layout, std::string& error);
// Sets the declared speaker settings named in `entries` (the layout's are
// skipped).
bool apply_speaker_parameters(const std::vector<ParameterEntry>& entries, Speaker& speaker, std::string& error);
// The whole declared tree with the values the engine holds now.
void collect_parameters(const BodyLayout& layout, Speaker& speaker, std::vector<ParameterEntry>& entries);

}  // namespace tract
