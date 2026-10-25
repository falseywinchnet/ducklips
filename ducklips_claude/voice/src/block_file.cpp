#include "voice/block_file.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <ostream>

namespace voice {

namespace {

// Splits a line into words, braces always words of their own; a '#' ends it.
void split_line(const std::string& line, std::vector<std::string>& words) {
  std::string word;
  for (char letter : line) {
    if (letter == '#') {
      break;
    }
    bool space = letter == ' ' || letter == '\t' || letter == '\r';
    bool brace = letter == '{' || letter == '}';
    if ((space || brace) && !word.empty()) {
      words.push_back(word);
      word.clear();
    }
    if (brace) {
      words.push_back(std::string(1, letter));
    } else if (!space) {
      word.push_back(letter);
    }
  }
  if (!word.empty()) {
    words.push_back(word);
  }
}

struct Word {
  std::string text;
  int line = 0;
};

// items := (NAME '{' items '}' | NAME NUMBER)*
bool parse_items(const std::vector<Word>& words, std::size_t& at, ParameterNode& block, std::string& error) {
  while (at < words.size()) {
    const Word& name = words[at];
    if (name.text == "}") {
      return true;
    }
    if (name.text == "{") {
      error = "line " + std::to_string(name.line) + ": a block needs a name";
      return false;
    }
    if (at + 1 >= words.size()) {
      error = "line " + std::to_string(name.line) + ": " + name.text + " has no value";
      return false;
    }
    const Word& next = words[at + 1];
    ParameterNode child;
    child.name = name.text;
    if (next.text == "{") {
      child.block = true;
      at += 2;
      if (!parse_items(words, at, child, error)) {
        return false;
      }
      if (at >= words.size()) {
        error = "block " + child.name + " is not closed";
        return false;
      }
      at += 1;  // its '}'
    } else {
      char* end = nullptr;
      child.value = std::strtod(next.text.c_str(), &end);
      if (end == next.text.c_str() || *end != '\0') {
        error = "line " + std::to_string(next.line) + ": " + name.text + ": not a number: " + next.text;
        return false;
      }
      at += 2;
    }
    block.children.push_back(child);
  }
  return true;
}

void flatten_into(const ParameterNode& node, const std::string& prefix, std::vector<ParameterEntry>& entries) {
  for (const ParameterNode& child : node.children) {
    std::string path = prefix.empty() ? child.name : prefix + "." + child.name;
    if (child.block) {
      flatten_into(child, path, entries);
    } else {
      ParameterEntry entry;
      entry.path = path;
      entry.value = child.value;
      entries.push_back(entry);
    }
  }
}

void split_path(const std::string& path, std::vector<std::string>& parts) {
  parts.clear();
  std::size_t begin = 0;
  while (true) {
    std::size_t dot = path.find('.', begin);
    parts.push_back(path.substr(begin, dot == std::string::npos ? std::string::npos : dot - begin));
    if (dot == std::string::npos) {
      return;
    }
    begin = dot + 1;
  }
}

}  // namespace

bool parse_parameters(std::istream& input, ParameterNode& root, std::string& error) {
  root = ParameterNode();
  root.block = true;
  std::vector<Word> words;
  std::vector<std::string> line_words;
  std::string line;
  int number = 0;
  while (std::getline(input, line)) {
    number += 1;
    line_words.clear();
    split_line(line, line_words);
    for (const std::string& text : line_words) {
      Word word;
      word.text = text;
      word.line = number;
      words.push_back(word);
    }
  }
  std::size_t at = 0;
  if (!parse_items(words, at, root, error)) {
    return false;
  }
  if (at < words.size()) {
    error = "line " + std::to_string(words[at].line) + ": unmatched }";
    return false;
  }
  return true;
}

void flatten_parameters(const ParameterNode& root, std::vector<ParameterEntry>& entries) {
  flatten_into(root, "", entries);
}

bool read_parameter_file(const std::string& path, std::vector<ParameterEntry>& entries, std::string& error) {
  std::ifstream file(path);
  if (!file) {
    error = "cannot read " + path;
    return false;
  }
  ParameterNode tree;
  if (!parse_parameters(file, tree, error)) {
    error = path + ": " + error;
    return false;
  }
  entries.clear();
  flatten_parameters(tree, entries);
  return true;
}

void write_parameters(std::ostream& output, const std::vector<ParameterEntry>& entries) {
  std::vector<std::string> open;
  std::vector<std::string> parts;
  char number[64];
  for (const ParameterEntry& entry : entries) {
    split_path(entry.path, parts);
    std::size_t shared = 0;
    while (shared < open.size() && shared + 1 < parts.size() && open[shared] == parts[shared]) {
      shared += 1;
    }
    while (open.size() > shared) {
      open.pop_back();
      output << std::string(2 * open.size(), ' ') << "}\n";
    }
    while (open.size() + 1 < parts.size()) {
      output << std::string(2 * open.size(), ' ') << parts[open.size()] << " {\n";
      open.push_back(parts[open.size()]);
    }
    std::snprintf(number, sizeof(number), "%.10g", entry.value);
    output << std::string(2 * open.size(), ' ') << parts.back() << " " << number << "\n";
  }
  while (!open.empty()) {
    open.pop_back();
    output << std::string(2 * open.size(), ' ') << "}\n";
  }
}

}  // namespace voice
