#pragma once

// Block files: named blocks of numbers, the parameter trees of the engine
// (tract/parameter_tree.hpp) and the planner (voice/plan_parameters.hpp).
//
//   # comment
//   fold {
//     cover_mass 0.09
//   }
//
// A file parses to a tree (ParameterNode); its leaves are addressed by
// dotted paths ("fold.cover_mass").

#include <iosfwd>
#include <string>
#include <vector>

namespace voice {

// One node of a block file: a block of named children, or a number.
struct ParameterNode {
  std::string name;
  bool block = false;
  double value = 0.0;
  std::vector<ParameterNode> children;
};

// A leaf by its dotted path.
struct ParameterEntry {
  std::string path;
  double value = 0.0;
};

// Parses a block file into `root` (a block named ""). On failure returns
// false with the line and reason in `error`.
bool parse_parameters(std::istream& input, ParameterNode& root, std::string& error);
// Every leaf under `root`, in file order.
void flatten_parameters(const ParameterNode& root, std::vector<ParameterEntry>& entries);
// Reads and flattens a block file; false (with the reason) when it cannot.
bool read_parameter_file(const std::string& path, std::vector<ParameterEntry>& entries, std::string& error);
// Writes entries (grouped by their paths' blocks, in the given order) as a
// block file.
void write_parameters(std::ostream& output, const std::vector<ParameterEntry>& entries);

}  // namespace voice
