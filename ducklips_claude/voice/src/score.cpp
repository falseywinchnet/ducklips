#include "voice/score.hpp"

#include <algorithm>
#include <cstdio>
#include <sstream>

namespace voice {

namespace {

struct PartName {
  Part part;
  const char* name;
};

constexpr PartName kPartNames[] = {
    {Part::lungs, "lungs"},   {Part::pitch, "pitch"}, {Part::glottis, "glottis"},
    {Part::tongue, "tongue"}, {Part::tip, "tip"},     {Part::dorsum, "dorsum"},
    {Part::lips, "lips"},     {Part::velum, "velum"}, {Part::floor, "floor"},  {Part::jaw, "jaw"},
    {Part::lateral, "lateral"},
};

struct MarkName {
  Mark mark;
  const char* name;
};

constexpr MarkName kMarkNames[] = {
    {Mark::breath, "breath"}, {Mark::phrase, "phrase"}, {Mark::word, "word"}, {Mark::stop, "stop"}};

bool part_event_earlier(const PartEvent& left, const PartEvent& right) {
  return left.time < right.time;
}

bool mark_event_earlier(const MarkEvent& left, const MarkEvent& right) {
  return left.time < right.time;
}

}  // namespace

const char* part_name(Part part) {
  for (const PartName& entry : kPartNames) {
    if (entry.part == part) {
      return entry.name;
    }
  }
  return "unknown";
}

const char* mark_name(Mark mark) {
  for (const MarkName& entry : kMarkNames) {
    if (entry.mark == mark) {
      return entry.name;
    }
  }
  return "unknown";
}

bool part_from_name(const std::string& name, Part& part) {
  for (const PartName& entry : kPartNames) {
    if (name == entry.name) {
      part = entry.part;
      return true;
    }
  }
  return false;
}

bool mark_from_name(const std::string& name, Mark& mark) {
  for (const MarkName& entry : kMarkNames) {
    if (name == entry.name) {
      mark = entry.mark;
      return true;
    }
  }
  return false;
}

void sort_score(Score& score) {
  std::stable_sort(score.parts.begin(), score.parts.end(), part_event_earlier);
  std::stable_sort(score.marks.begin(), score.marks.end(), mark_event_earlier);
}

std::string score_text(const Score& score) {
  std::string text;
  std::size_t part_index = 0;
  std::size_t mark_index = 0;
  char line[256];
  while (part_index < score.parts.size() || mark_index < score.marks.size()) {
    bool take_mark = mark_index < score.marks.size() &&
                     (part_index >= score.parts.size() ||
                      score.marks[mark_index].time <= score.parts[part_index].time);
    if (take_mark) {
      const MarkEvent& event = score.marks[mark_index];
      std::snprintf(line, sizeof(line), "%.4f mark %s %.3f", event.time, mark_name(event.mark),
                    static_cast<double>(event.a));
      text += line;
      if (!event.text.empty()) {
        text += " ";
        text += event.text;
      }
      text += "\n";
      ++mark_index;
    } else {
      const PartEvent& event = score.parts[part_index];
      std::snprintf(line, sizeof(line), "%.4f %s %.3f %.3f %.3f\n", event.time, part_name(event.part),
                    static_cast<double>(event.a), static_cast<double>(event.b),
                    static_cast<double>(event.move));
      text += line;
      ++part_index;
    }
  }
  return text;
}

bool parse_score_line(const std::string& line, Score& score) {
  std::istringstream stream(line);
  double time = 0.0;
  std::string word;
  if (!(stream >> time >> word)) {
    return false;
  }
  if (word == "mark") {
    MarkEvent event;
    event.time = time;
    std::string name;
    if (!(stream >> name) || !mark_from_name(name, event.mark)) {
      return false;
    }
    stream >> event.a;
    std::string rest;
    std::getline(stream, rest);
    std::size_t start = rest.find_first_not_of(' ');
    if (start != std::string::npos) {
      event.text = rest.substr(start);
    }
    score.marks.push_back(event);
    score.duration = std::max(score.duration, time + static_cast<double>(event.a));
    return true;
  }
  PartEvent event;
  event.time = time;
  if (!part_from_name(word, event.part)) {
    return false;
  }
  if (!(stream >> event.a >> event.b >> event.move)) {
    return false;
  }
  score.parts.push_back(event);
  score.duration = std::max(score.duration, time);
  return true;
}

}  // namespace voice
