#pragma once

// The sustain tool's run, as a reusable object: the command line
// (tools/sustain.cpp) and the in-process engine (bindings/tract_engine.cpp)
// both fill a SustainSettings and step a SustainRun, so the two cannot
// drift. The settings carry every input the command line takes from its
// arguments, its environment and the files the environment names, already
// parsed (doubles, never text); settings_from_environment is the one
// parser of the environment's text, used by the command line and by the
// engine for the variables a caller's environment passes through.
//
// A SustainRun is a value: copying it snapshots the whole run (the voice
// body and everything it owns, the loop's filters, step and audio), and a
// copy resumes exactly where the original stood.

#include "tract/cord.hpp"
#include "tract/voice_body.hpp"

#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace tract {

// A two-pole section of the 96 -> 48 kHz decimation filter (the arithmetic
// of tools/wave_out.hpp's Biquad, here in the library's namespace).
struct DecimationSection {
  double b0 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
  double z1 = 0.0;
  double z2 = 0.0;
};

DecimationSection decimation_lowpass(double corner, double rate, double q);
double run_section(DecimationSection& section, double input);
// A 16-bit mono WAV, as tools/wave_out.hpp writes it.
void write_wave_file(const char* path, const std::vector<float>& samples, int rate);

struct SustainSettings {
  // Command-line arguments: VOWEL OUT.wav [pitch_hz] [lung_pa] [seconds].
  std::string vowel = "aa";
  std::string wave_path;  // empty: no WAV is written
  double pitch = 120.0;
  double lung = 800.0;
  double seconds = 1.5;
  // The tube: AREAFILE's (or the Story vowel's) section length and areas.
  double section_length = 0.0;
  std::vector<double> areas;
  // AREAFILE2 with GLIDE=start,duration (both present).
  bool glide = false;
  std::vector<double> areas_to;
  double glide_start = 0.0;
  double glide_duration = 0.0;
  // AREASEQ: keyframe times and their areas (row-major, areas.size() a row).
  std::vector<double> key_times;
  std::vector<double> key_areas;
  // VELUMSEQ.
  std::vector<double> velum_times;
  std::vector<double> velum_areas;
  // FOLDSEQ: trigger times (s), postures and move times.
  std::vector<double> fold_times;
  std::vector<CordPosture> fold_postures;
  std::vector<double> fold_moves;
  // JETSEQ.
  std::vector<double> jet_times;
  std::vector<double> jet_gates;
  // FOLDKEYS.
  std::vector<double> fold_key_times;
  std::vector<CordPosture> fold_key_postures;
  // AREATICKS, VELUMTICKS, FOLDTICKS: a control field sampled at the 1 ms
  // ticks the run reads its controls at (row k at tick k, the last row held
  // after). A sampled field enters exactly: the engine reads controls only at
  // ticks. They take the place of AREASEQ, VELUMSEQ and FOLDKEYS when set.
  std::vector<double> tick_areas;
  std::vector<double> tick_velum;
  std::vector<CordPosture> tick_folds;
  // Scalars of the environment (absent: the variable is unset).
  std::optional<double> pitch_end;
  double release = 0.0;
  double trace_ms = 0.0;
  std::optional<std::string> branches;
  std::optional<double> aspiration;
  std::optional<double> frication;
  std::optional<double> flux;
  std::optional<double> closing;
  std::optional<double> tremor;
  std::optional<double> adductor_tremor;
  std::optional<double> breath_tremor;
  std::optional<double> pitch_feedback;
  std::optional<double> port_resistance;
  std::optional<double> reversible_jet;
  std::optional<double> reversible_jet_closing;
  std::optional<double> fold_turbulence;
  std::optional<double> skin;
  double wall_scale = 1.0;
  double mass_scale = 1.0;
  double axial_scale = 1.0;
  std::optional<double> wall_variation;
  std::optional<double> viscous;
  std::optional<double> active_length;
  std::optional<double> static_area;
  std::optional<double> reversible_nose;
  std::optional<double> seal;
  std::optional<double> seal_blend;
  std::optional<double> jet_bound;   // JETBOUND=1 (VoiceBody::set_jet_bound, from the lungs' pressure)  // SEALBLEND=1 (Tube::set_seal_blend)
  bool nose_radiation = false;
  double nose_radiation_mass = 1.0;
  double nose_radiation_resistance = 1.0;
  std::optional<double> nose_viscous;
  std::optional<double> nose_shelf;
  std::optional<double> nose_hop;
  std::optional<double> hop_loss;
  bool rest = false;
  double rest_lower = 0.0;
  double rest_upper = 0.0;
  std::optional<double> chink;
  bool fold_damping = false;
  double fold_damping_lower = 0.0;
  double fold_damping_upper = 0.0;
  std::optional<double> velum;
  // FORCEWAVE (per step, NaN: free) and FORCEFLOW with OBSERVE.
  std::vector<double> forced_waves;
  std::vector<double> observed_flows;
  double observe_gain = 0.0;
  bool raw_cord = false;    // RAWCORD (present)
  bool raw_double = false;  // RAWDOUBLE (nonzero)
  // Output paths (empty: not written).
  std::string source_path;      // SOURCE
  std::string raw_path;         // RAW
  std::string port_waves_path;  // PORTWAVES
  std::string port_dump_path;   // PORTDUMP
  std::string geometry_path;    // GEODUMP
  std::string areas_path;       // AREAS
  std::string wall_dump_path;   // WALLDUMP
  std::string nose_dump_path;   // NOSEDUMP
  // (Every field is either part of the body's construction or output,
  // compared whole by same_construction, or an input over time, compared
  // step by step by first_divergent_step: a new field must join one of
  // them, or a resumed run could reuse a prefix it does not share.)
};

// The inputs a run applies at a control tick (every 96 steps, 1 ms): the
// command line's keyframe interpolations, computed once for applying and
// for comparing two settings.
struct TickControls {
  double pitch = 0.0;
  bool fold_key = false;
  CordPosture fold_posture;
  bool jet = false;
  double jet_gate = 0.0;
  bool velum = false;
  double velum_area = 0.0;
  int area_source = 0;  // 0 none, 1 AREASEQ keyframes, 2 AREAFILE2 glide, 3 AREATICKS field
  std::vector<double> rest_areas;
};

// The tick's controls at `step` (a multiple of 96) of a run of `total` steps.
void tick_controls(const SustainSettings& settings, long step, long total, TickControls& controls);
// The step a FOLDSEQ trigger at `time` (s) fires at.
long trigger_step(double time);
// The steps a run of these settings takes and the step its RELEASE starts.
long total_steps(const SustainSettings& settings);
long release_step(const SustainSettings& settings);

// Whether two settings build the same body and record the same outputs
// (every field but the inputs over time), bit for bit.
bool same_construction(const SustainSettings& first, const SustainSettings& second);
// The first step whose inputs differ between two settings (the controls
// applied before it steps, its forced wave and observed flow): a run of
// `second` may resume from any state of a run of `first` at a step at or
// before it. -1 when they build different bodies (no state of the one is
// the other's); the shorter run's total when nothing differs within it.
long first_divergent_step(const SustainSettings& first, const SustainSettings& second);

// Where a variable's value is looked up (the process environment, or a
// caller's list). Returns nullptr for an unset variable.
class EnvironmentLookup {
 public:
  virtual ~EnvironmentLookup() = default;
  virtual const char* find(const char* name) const = 0;
};

// The process environment (std::getenv).
class ProcessEnvironment : public EnvironmentLookup {
 public:
  const char* find(const char* name) const override;
};

// Fills `settings` from the arguments and the environment exactly as the
// sustain tool reads them (files named by the environment are read).
// `arguments` are argv[1..]; returns false with `error` set (and the
// command line's exit status in `status`) when the tool would stop.
bool settings_from_environment(const std::vector<std::string>& arguments, const EnvironmentLookup& environment,
                               SustainSettings& settings, std::string& error, int& status);

// Every variable but AREAFILE, read as the command line reads it (the tube's
// areas must already be in `settings`: AREAFILE2 and AREASEQ are checked
// against their count). Unset variables leave their defaults.
bool options_from_environment(const EnvironmentLookup& environment, SustainSettings& settings, std::string& error);

// Checks what the command line checks after reading (counts of keyframe
// rows, the glide); returns false with `error` set.
bool validate_settings(const SustainSettings& settings, std::string& error);

// Steps since the run started: the RAW pair's flows for one step.
struct StepFlows {
  double entering = 0.0;
  double radiated = 0.0;
};

// Receives each step's outputs (the per-step streams). The body is the
// run's own, as it stands after the step; read it, never step it.
class SustainSink {
 public:
  virtual ~SustainSink() = default;
  virtual void record(long step, VoiceBody& body, const StepFlows& flows) = 0;
  // After the step's output passed the finiteness check (the command
  // line's TRACE lines come here).
  virtual void checked(long step, VoiceBody& body);
};

// PORTWAVES' six columns for the body as it stands.
void port_wave_row(VoiceBody& body, double row[6]);
// GEODUMP's twelve columns for one vocal node.
void geometry_row(Tube& vocal, Index node, double row[12]);

// The per-step files the command line writes (RAW, PORTWAVES, PORTDUMP,
// GEODUMP) and its TRACE lines, each when its path (or TRACE) is set.
class StreamFiles : public SustainSink {
 public:
  StreamFiles(const SustainSettings& settings, bool trace);
  ~StreamFiles() override;
  StreamFiles(const StreamFiles&) = delete;
  StreamFiles& operator=(const StreamFiles&) = delete;
  void record(long step, VoiceBody& body, const StepFlows& flows) override;
  void checked(long step, VoiceBody& body) override;

 private:
  std::FILE* raw_ = nullptr;
  std::FILE* port_waves_ = nullptr;
  std::FILE* port_dump_ = nullptr;
  std::FILE* geometry_ = nullptr;
  bool raw_double_ = false;
  long trace_every_ = 0;
};

class SustainRun {
 public:
  // Builds and configures the body as the command line does before its
  // first step.
  explicit SustainRun(const SustainSettings& settings);

  // Steps from step() up to `end` (exclusive; at most total()), reading the
  // time-varying inputs from `settings` (the run's own or a changed copy:
  // inputs acting from the current step on may differ from the ones the
  // run began with). Returns false at a non-finite output, after recording
  // that step (as the command line stops after writing its RAW pair).
  bool advance(const SustainSettings& settings, long end, SustainSink& sink);

  // The command line's ending: the WAV (and SOURCE) normalized to 0.7 of
  // their peaks, written when their paths are set.
  void finish(const SustainSettings& settings) const;

  long step() const { return step_; }
  long total() const { return total_; }
  bool failed() const { return failed_; }
  double peak() const { return peak_; }
  VoiceBody& body() { return body_; }
  const VoiceBody& body() const { return body_; }

 private:
  void control(const SustainSettings& settings, long step);
  void first_step_dumps(const SustainSettings& settings);

  VoiceBody body_;
  DecimationSection first_;
  DecimationSection second_;
  DecimationSection source_first_;
  DecimationSection source_second_;
  long total_ = 0;
  long release_step_ = 0;
  long step_ = 0;
  bool failed_ = false;
  bool keep_audio_ = false;
  bool keep_source_ = false;
  double peak_ = 0.0;
  double source_peak_ = 0.0;
  double previous_flow_ = 0.0;
  std::vector<double> audio_;
  std::vector<double> source_;
  TickControls tick_;  // the current tick's controls (kept to reuse its buffer)
};

}  // namespace tract
