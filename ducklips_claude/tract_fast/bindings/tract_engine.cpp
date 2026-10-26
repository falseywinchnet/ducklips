// The in-process engine's C interface (tract_engine.h) over the sustain
// tool's own run (tract/sustain_run.hpp).

#include "tract_engine.h"

#include "tract/sustain_run.hpp"

#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

struct tract_settings {
  tract::SustainSettings settings;
};

struct tract_run {
  explicit tract_run(const tract::SustainSettings& settings) : run(settings) {}
  tract_run(const tract_run& other) : run(other.run), files_opened(true) {}
  tract::SustainRun run;
  // The per-step files named by the settings (RAW, PORTWAVES, PORTDUMP,
  // GEODUMP), opened at this handle's first advance; a copy streams none.
  std::unique_ptr<tract::StreamFiles> files;
  bool files_opened = false;
};

namespace {

thread_local std::string last_error;

int fail(const std::string& text) {
  last_error = text;
  return 1;
}

bool same(const char* name, const char* expected) { return std::strcmp(name, expected) == 0; }

// The settings field a numeric variable sets (nullptr: none of them).
std::optional<double>* optional_field(tract::SustainSettings& settings, const char* name) {
  struct Entry {
    const char* name;
    std::optional<double> tract::SustainSettings::*field;
  };
  static const Entry kEntries[] = {
      {"PITCHEND", &tract::SustainSettings::pitch_end},
      {"ASPIRATION", &tract::SustainSettings::aspiration},
      {"FRICATION", &tract::SustainSettings::frication},
      {"FLUX", &tract::SustainSettings::flux},
      {"CLOSING", &tract::SustainSettings::closing},
      {"TREMOR", &tract::SustainSettings::tremor},
      {"ADDTREMOR", &tract::SustainSettings::adductor_tremor},
      {"BREATHTREMOR", &tract::SustainSettings::breath_tremor},
      {"PITCHFB", &tract::SustainSettings::pitch_feedback},
      {"PORTRESIST", &tract::SustainSettings::port_resistance},
      {"REVJET", &tract::SustainSettings::reversible_jet},
      {"REVJETCLOSE", &tract::SustainSettings::reversible_jet_closing},
      {"FOLDTURB", &tract::SustainSettings::fold_turbulence},
      {"SKIN", &tract::SustainSettings::skin},
      {"WALLVAR", &tract::SustainSettings::wall_variation},
      {"VISCOUS", &tract::SustainSettings::viscous},
      {"ACTIVELENGTH", &tract::SustainSettings::active_length},
      {"STATICAREA", &tract::SustainSettings::static_area},
      {"REVNOSE", &tract::SustainSettings::reversible_nose},
      {"SEAL", &tract::SustainSettings::seal},
      {"SEALBLEND", &tract::SustainSettings::seal_blend},
      {"JETBOUND", &tract::SustainSettings::jet_bound},
      {"ENTERINGONLY", &tract::SustainSettings::entering_from},
      {"NOSEVISCOUS", &tract::SustainSettings::nose_viscous},
      {"NOSESHELF", &tract::SustainSettings::nose_shelf},
      {"NOSEHOP", &tract::SustainSettings::nose_hop},
      {"HOPLOSS", &tract::SustainSettings::hop_loss},
      {"CHINK", &tract::SustainSettings::chink},
      {"VELUM", &tract::SustainSettings::velum},
  };
  for (const Entry& entry : kEntries) {
    if (same(name, entry.name)) {
      return &(settings.*(entry.field));
    }
  }
  return nullptr;
}

double* plain_field(tract::SustainSettings& settings, const char* name) {
  struct Entry {
    const char* name;
    double tract::SustainSettings::*field;
  };
  static const Entry kEntries[] = {
      {"PITCH", &tract::SustainSettings::pitch},
      {"LUNG", &tract::SustainSettings::lung},
      {"SECONDS", &tract::SustainSettings::seconds},
      {"RELEASE", &tract::SustainSettings::release},
      {"TRACE", &tract::SustainSettings::trace_ms},
      {"WALLS", &tract::SustainSettings::wall_scale},
      {"WALLMASS", &tract::SustainSettings::mass_scale},
      {"AXIAL", &tract::SustainSettings::axial_scale},
      {"OBSERVE", &tract::SustainSettings::observe_gain},
  };
  for (const Entry& entry : kEntries) {
    if (same(name, entry.name)) {
      return &(settings.*(entry.field));
    }
  }
  return nullptr;
}

std::string* path_field(tract::SustainSettings& settings, const char* name) {
  struct Entry {
    const char* name;
    std::string tract::SustainSettings::*field;
  };
  static const Entry kEntries[] = {
      {"WAVE", &tract::SustainSettings::wave_path},
      {"SOURCE", &tract::SustainSettings::source_path},
      {"RAW", &tract::SustainSettings::raw_path},
      {"PORTWAVES", &tract::SustainSettings::port_waves_path},
      {"PORTDUMP", &tract::SustainSettings::port_dump_path},
      {"GEODUMP", &tract::SustainSettings::geometry_path},
      {"AREAS", &tract::SustainSettings::areas_path},
      {"WALLDUMP", &tract::SustainSettings::wall_dump_path},
      {"NOSEDUMP", &tract::SustainSettings::nose_dump_path},
  };
  for (const Entry& entry : kEntries) {
    if (same(name, entry.name)) {
      return &(settings.*(entry.field));
    }
  }
  return nullptr;
}

// A caller's list of variables.
class ListedEnvironment : public tract::EnvironmentLookup {
 public:
  ListedEnvironment(int count, const char* const* names, const char* const* values)
      : count_(count), names_(names), values_(values) {}
  const char* find(const char* name) const override {
    const char* found = nullptr;
    for (int index = 0; index < count_; ++index) {
      if (same(names_[index], name)) {
        found = values_[index];
      }
    }
    return found;
  }

 private:
  int count_;
  const char* const* names_;
  const char* const* values_;
};

// Writes each step's rows into the caller's buffers (from their starts at
// the advance's first step) and forwards it to the run's files.
class BufferSink : public tract::SustainSink {
 public:
  BufferSink(long first, double* raw, double* port_waves, double* geometry, tract::StreamFiles* files)
      : first_(first), raw_(raw), port_waves_(port_waves), geometry_(geometry), files_(files) {}

  void record(long step, tract::VoiceBody& body, const tract::StepFlows& flows) override {
    std::size_t offset = static_cast<std::size_t>(step - first_);
    if (raw_ != nullptr) {
      raw_[2 * offset] = flows.entering;
      raw_[2 * offset + 1] = flows.radiated;
    }
    if (port_waves_ != nullptr) {
      tract::port_wave_row(body, port_waves_ + 6 * offset);
    }
    if (geometry_ != nullptr) {
      tract::Tube& vocal = body.vocal();
      std::size_t count = static_cast<std::size_t>(vocal.count());
      for (tract::Index node = 0; node < vocal.count(); ++node) {
        tract::geometry_row(vocal, node, geometry_ + 12 * (offset * count + static_cast<std::size_t>(node)));
      }
    }
    if (files_ != nullptr) {
      files_->record(step, body, flows);
    }
    written_ = static_cast<long>(offset) + 1;
  }

  long written() const { return written_; }

 private:
  long first_;
  double* raw_;
  double* port_waves_;
  double* geometry_;
  tract::StreamFiles* files_;
  long written_ = 0;
};

bool any_stream(const tract::SustainSettings& settings) {
  return !settings.raw_path.empty() || !settings.port_waves_path.empty() || !settings.port_dump_path.empty() ||
         !settings.geometry_path.empty();
}

}  // namespace

extern "C" {

const char* tract_error(void) { return last_error.c_str(); }

tract_settings* tract_settings_new(void) { return new tract_settings(); }

tract_settings* tract_settings_copy(const tract_settings* settings) { return new tract_settings(*settings); }

void tract_settings_free(tract_settings* settings) { delete settings; }

int tract_settings_number(tract_settings* handle, const char* name, double value) {
  tract::SustainSettings& settings = handle->settings;
  if (std::optional<double>* field = optional_field(settings, name)) {
    *field = value;
    return 0;
  }
  if (double* field = plain_field(settings, name)) {
    *field = value;
    return 0;
  }
  if (same(name, "RAWDOUBLE")) {
    settings.raw_double = value != 0.0;
    return 0;
  }
  if (same(name, "RAWCORD")) {
    settings.raw_cord = true;  // the tool reads its presence alone
    return 0;
  }
  return fail(std::string("unknown number ") + name);
}

int tract_settings_pair(tract_settings* handle, const char* name, double first, double second) {
  tract::SustainSettings& settings = handle->settings;
  if (same(name, "NOSERAD")) {
    settings.nose_radiation = true;
    settings.nose_radiation_mass = first;
    settings.nose_radiation_resistance = second;
  } else if (same(name, "REST")) {
    settings.rest = true;
    settings.rest_lower = first;
    settings.rest_upper = second;
  } else if (same(name, "FOLDDAMP")) {
    settings.fold_damping = true;
    settings.fold_damping_lower = first;
    settings.fold_damping_upper = second;
  } else if (same(name, "GLIDE")) {
    settings.glide_start = first;
    settings.glide_duration = second;
  } else {
    return fail(std::string("unknown pair ") + name);
  }
  return 0;
}

int tract_settings_text(tract_settings* handle, const char* name, const char* text) {
  tract::SustainSettings& settings = handle->settings;
  if (same(name, "BRANCHES")) {
    settings.branches = std::string(text);
    return 0;
  }
  if (std::string* field = path_field(settings, name)) {
    *field = text;
    return 0;
  }
  return fail(std::string("unknown text ") + name);
}

int tract_settings_array(tract_settings* handle, const char* name, const double* values, int64_t rows,
                         int64_t columns) {
  tract::SustainSettings& settings = handle->settings;
  std::size_t row_count = static_cast<std::size_t>(rows);
  std::size_t column_count = static_cast<std::size_t>(columns);
  std::size_t count = row_count * column_count;
  if (rows < 0 || columns < 0) {
    return fail(std::string("negative shape for ") + name);
  }
  if (same(name, "AREAFILE")) {
    if (rows != 1 || columns < 2) {
      return fail("AREAFILE is one row: the section length, then the areas");
    }
    settings.section_length = values[0];
    settings.areas.assign(values + 1, values + count);
  } else if (same(name, "AREAFILE2")) {
    if (rows != 1 || columns < 2) {
      return fail("AREAFILE2 is one row: a length, then the areas");
    }
    settings.areas_to.assign(values + 1, values + count);
    settings.glide = true;
  } else if (same(name, "AREASEQ")) {
    if (rows > 0 && column_count != settings.areas.size() + 1) {
      return fail("AREASEQ rows are a time and the tube's count of areas (set AREAFILE first)");
    }
    settings.key_times.clear();
    settings.key_areas.clear();
    for (std::size_t row = 0; row < row_count; ++row) {
      const double* line = values + row * column_count;
      settings.key_times.push_back(line[0]);
      settings.key_areas.insert(settings.key_areas.end(), line + 1, line + column_count);
    }
  } else if (same(name, "VELUMSEQ") || same(name, "JETSEQ")) {
    if (rows > 0 && columns != 2) {
      return fail(std::string(name) + " rows are a time and a value");
    }
    std::vector<double>& times = same(name, "VELUMSEQ") ? settings.velum_times : settings.jet_times;
    std::vector<double>& levels = same(name, "VELUMSEQ") ? settings.velum_areas : settings.jet_gates;
    times.clear();
    levels.clear();
    for (std::size_t row = 0; row < row_count; ++row) {
      times.push_back(values[2 * row]);
      levels.push_back(values[2 * row + 1]);
    }
  } else if (same(name, "FOLDKEYS")) {
    if (rows > 0 && columns != 4) {
      return fail("FOLDKEYS rows are a time, the lower and upper rests and the chink");
    }
    settings.fold_key_times.clear();
    settings.fold_key_postures.clear();
    for (std::size_t row = 0; row < row_count; ++row) {
      const double* line = values + 4 * row;
      tract::CordPosture posture;
      posture.lower_rest = line[1];
      posture.upper_rest = line[2];
      posture.chink = line[3];
      settings.fold_key_times.push_back(line[0]);
      settings.fold_key_postures.push_back(posture);
    }
  } else if (same(name, "FOLDSEQ")) {
    if (rows > 0 && columns != 5) {
      return fail("FOLDSEQ rows are a time, the lower and upper rests, the chink and a move time");
    }
    settings.fold_times.clear();
    settings.fold_postures.clear();
    settings.fold_moves.clear();
    for (std::size_t row = 0; row < row_count; ++row) {
      const double* line = values + 5 * row;
      tract::CordPosture posture;
      posture.lower_rest = line[1];
      posture.upper_rest = line[2];
      posture.chink = line[3];
      settings.fold_times.push_back(line[0]);
      settings.fold_postures.push_back(posture);
      settings.fold_moves.push_back(line[4]);
    }
  } else if (same(name, "AREATICKS")) {
    if (rows > 0 && column_count != settings.areas.size()) {
      return fail("AREATICKS rows are the tube's count of areas (set AREAFILE first)");
    }
    settings.tick_areas.assign(values, values + count);
  } else if (same(name, "VELUMTICKS")) {
    if (rows > 0 && columns != 1) {
      return fail("VELUMTICKS rows are one area");
    }
    settings.tick_velum.assign(values, values + count);
  } else if (same(name, "FOLDTICKS")) {
    if (rows > 0 && columns != 3) {
      return fail("FOLDTICKS rows are the lower and upper rests and the chink");
    }
    settings.tick_folds.clear();
    for (std::size_t row = 0; row < row_count; ++row) {
      const double* line = values + 3 * row;
      tract::CordPosture posture;
      posture.lower_rest = line[0];
      posture.upper_rest = line[1];
      posture.chink = line[2];
      settings.tick_folds.push_back(posture);
    }
  } else if (same(name, "FORCEWAVE")) {
    settings.forced_waves.assign(values, values + count);
  } else if (same(name, "FORCEFLOW")) {
    settings.observed_flows.assign(values, values + count);
  } else {
    return fail(std::string("unknown array ") + name);
  }
  return 0;
}

int tract_settings_environment(tract_settings* handle, int count, const char* const* names,
                               const char* const* values) {
  ListedEnvironment environment(count, names, values);
  std::string error;
  if (environment.find("AREAFILE") != nullptr) {
    // The tool's AREAFILE reading, through its own parser.
    std::vector<std::string> arguments = {"aa", ""};
    int status = 0;
    tract::SustainSettings parsed;
    if (!tract::settings_from_environment(arguments, environment, parsed, error, status)) {
      return fail(error);
    }
    parsed.wave_path = handle->settings.wave_path;
    parsed.pitch = handle->settings.pitch;
    parsed.lung = handle->settings.lung;
    parsed.seconds = handle->settings.seconds;
    handle->settings = parsed;
    return 0;
  }
  if (!tract::options_from_environment(environment, handle->settings, error)) {
    return fail(error);
  }
  return 0;
}

int64_t tract_settings_steps(const tract_settings* settings) {
  return static_cast<int64_t>(static_cast<long>(settings->settings.seconds * tract::kStepRate));
}

int64_t tract_settings_vocal_nodes(const tract_settings* settings) {
  tract::BodyLayout layout = tract::vowel_layout(settings->settings.areas, settings->settings.section_length);
  return static_cast<int64_t>(layout.vocal_areas.size());
}

int64_t tract_settings_divergence(const tract_settings* first, const tract_settings* second) {
  return static_cast<int64_t>(tract::first_divergent_step(first->settings, second->settings));
}

tract_run* tract_run_new(const tract_settings* settings) {
  std::string error;
  if (!tract::validate_settings(settings->settings, error)) {
    fail(error);
    return nullptr;
  }
  return new tract_run(settings->settings);
}

tract_run* tract_run_copy(const tract_run* run) { return new tract_run(*run); }

void tract_run_free(tract_run* run) { delete run; }

int64_t tract_run_step(const tract_run* run) { return run->run.step(); }

int64_t tract_run_total(const tract_run* run) { return run->run.total(); }

int64_t tract_run_advance(tract_run* run, const tract_settings* settings, int64_t end, double* raw,
                          double* port_waves, double* geometry, int32_t* status) {
  if (!run->files_opened) {
    run->files_opened = true;
    if (any_stream(settings->settings)) {
      run->files = std::make_unique<tract::StreamFiles>(settings->settings, false);
    }
  }
  BufferSink sink(run->run.step(), raw, port_waves, geometry, run->files.get());
  bool finite = run->run.advance(settings->settings, static_cast<long>(end), sink);
  if (status != nullptr) {
    *status = finite ? 0 : 1;
  }
  return sink.written();
}

void tract_run_finish(const tract_run* run, const tract_settings* settings) { run->run.finish(settings->settings); }

int64_t tract_run_nose(const tract_run* handle, double* header, double* rows, int64_t capacity) {
  tract::VoiceBody& body = const_cast<tract::VoiceBody&>(handle->run.body());
  tract::Tube& nose = body.nose();
  if (header != nullptr) {
    header[0] = static_cast<double>(body.nose_port_node());
    header[1] = static_cast<double>(nose.count());
    header[2] = nose.wave_smoothing();
    header[3] = nose.viscous_scale();
  }
  for (tract::Index node = 0; rows != nullptr && node < nose.count() && node < capacity; ++node) {
    double* row = rows + 12 * static_cast<std::size_t>(node);
    for (int column = 0; column < 12; ++column) {
      row[column] = 0.0;
    }
    row[0] = nose.area(node);
    row[2] = node + 1 < nose.count() ? nose.link_resistance(node) : 0.0;
    nose.lining(node, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11]);
  }
  return static_cast<int64_t>(nose.count());
}

int64_t tract_render(const tract_settings* settings, double* raw, double* port_waves, double* geometry,
                     int32_t* status) {
  tract_run* run = tract_run_new(settings);
  if (run == nullptr) {
    if (status != nullptr) {
      *status = 2;
    }
    return -1;
  }
  int64_t written = tract_run_advance(run, settings, run->run.total(), raw, port_waves, geometry, status);
  run->files.reset();
  tract_run_finish(run, settings);
  tract_run_free(run);
  return written;
}

}  // extern "C"
