#pragma once

// The connected system: lungs -> trachea (tube of ring modules) -> cord
// module -> vocal tube (modules from the epilarynx to the lips) -> radiation.
// The lungs are the only sustained energy source.

#include "tract/cord.hpp"
#include "tract/tube.hpp"
#include "tract/variation.hpp"

#include <algorithm>
#include <vector>

namespace tract {

struct BodyLayout {
  std::vector<double> trachea_areas;   // cm^2 per node, lungs first
  std::vector<double> vocal_areas;     // cm^2 per node, glottis first
  std::vector<int> vocal_materials;    // index into the vocal materials, per node
  std::vector<SegmentMaterial> vocal_material_list;
  double skin_gain = 1.0;              // share of the lining's motion heard through the skin
  double skin_corner = 250.0;          // Hz: two poles; neck and cheek tissue pass only the lows
  // Bounded depths of the slow variation of causes (0 turns each off).
  double breath_variation = 0.03;      // share of lung pressure
  double tension_variation = 0.015;    // share of cord tension
  double posture_variation = 0.0015;   // cm of fold rest gap
  double wall_variation = 0.10;        // share of wall support stiffness
  double wall_shared = 0.3;            // part of the wall variation common to the whole body
  // Glottal turbulence gain (VocalTractLab = 1). With the round-jet corner,
  // 8 puts /h/ above 2 kHz at -29 dB re the loudest vowel, bdl's median -27
  // (51 /hh/ in arctic). (On the slit-height corner the gain was 2.)
  double aspiration = 8.0;
  // Decelerating jets are unstable: while the folds cut the flow off, the
  // glottal jet's turbulence grows by this factor times the deceleration
  // (over 1 ms, relative to the flow). A steady jet (/h/) is unaffected.
  double closing_turbulence = 0.0;
  // Glottal noise corner from the jet's hydraulic diameter. On the slit height
  // the corner sat at 3-12 kHz: /hh/ came out 13 dB short at 2-3 kHz and
  // 5 dB strong above 5 kHz against bdl (267 tokens); round, with gain 8,
  // /hh/ is within 2 dB of bdl's shape except 2-3 kHz (-6).
  bool glottal_round_jet = true;
  // Turbulent flux through the glottis: rms share of the flow, times the
  // turbulence onset (Reynolds number), shaped like the jet noise.
  double jet_flux = 0.0;
  // Its spectrum's corner (Hz; 0 leaves it shaped like the jet noise, which
  // at the glottis's narrow slit is nearly white to 8 kHz: vowels' 6-8 kHz
  // stood 11-16 dB over bdl's). Two low-pass stages of white noise,
  // unity gain below the corner.
  double flux_corner = 0.0;
  // Physiological tremor, 2-15 Hz, rms: of the stretching muscle (share of
  // cord tension), the adducting muscles (cm of fold rest gap) and the
  // breathing muscles (share of lung pressure).
  double tremor = 0.005;
  double adductor_tremor = 0.0007;
  double breath_tremor = 0.015;
  // Auditory pitch feedback: share of the log pitch error corrected per
  // glottal cycle (0 leaves tension open-loop), and its bound.
  double pitch_feedback = 0.08;
  double pitch_correction_limit = 0.3;
  double fold_turbulence = 1.0;        // share of the jet's turbulent pressure felt by the folds
  // Turbulence of jets from oral constrictions (the front end's tuned jet
  // model, ported): overall gain, Strouhal number of the noise corner on the
  // slit height, the band below corner / frication_band rising 12 dB/oct,
  // distances over which a jet becomes turbulent and decays (cm), how far it
  // stays coherent enough to strike the incisors, and the extra source where
  // a lip channel meets the upper teeth.
  // 1.5: /s/ /sh/ at -14 dB re the loudest vowel (bdl -13; at 0.6, -22),
  // /f/ -34 (bdl -32), /z/ -16 (-13), /th/ -34 (-29); 40 bdl sentences.
  // 0.8 (2026-10-13; 1.5 before): with the folds damped and the tube's
  // highs lossier, /s/ sat 7 dB over bdl at 4.5-7 kHz and /sh/ 3-13 dB
  // over at 1.5-3 kHz; at 0.8 both sit within about 3 dB of bdl.
  double frication = 0.8;
  // s, the jet noise's build-up time (0: follows the velocity at once).
  // 0.008: /p/ /k/ releases rise 31 dB per 5 ms (bdl 29-34; 41 at 0) and
  // peak about 16 dB under the loudest vowel (bdl 20-21; 13-15 at 0).
  double jet_attack = 0.008;
  double strouhal = 0.2;
  double frication_band = 3.0;
  double jet_growth = 0.4;
  double jet_spread = 0.0;
  double jet_corner_ceiling = 12000.0; // Hz: the highest oral jet-noise corner
  double lip_slit_width = 2.5;         // cm: the width of a lip (or labiodental) slit             // cm added to the slit height where the turbulence forms
  double jet_decay = 1.2;
  double obstacle_reach = 0.8;
  double labiodental = 1.2;            // (0.6 in the front end's body; with /f/ at 0.06 cm^2, 1.2 gives bdl's /f/ level)
  bool piriform = true;                // side branches (diagnostics can switch them off)
  double piriform_loss = 0.5;          // the piriform lining's loss (fraction of critical)
  double piriform_port = 0.6;          // cm^2, both openings into the hypopharynx together
  double piriform_scale = 1.0;         // the recesses' cross-sections (and volume) times this
  // The recesses' frequency-growing loss (their hops' low-pass k, 0 off):
  // narrow, mucosa-lined pockets lose to their boundary layers, so their
  // zero is a dip, not a cut (lossless it took 17-23 dB at 4.6-5 kHz).
  double piriform_smoothing = 0.0;
  bool under_tongue = true;
  bool sinuses = true;
  // Scales the sinus ostia's areas (Dang & Honda's 0.11-0.185 cm^2). At
  // full size their side resonances cut nasal murmurs 4-7 dB under bdl at
  // 0.5-2 kHz (and 12 dB more than without sinuses); at 0.25 (0.03-0.05
  // cm^2, small but measured ostia) with the nose's hop loss 0.4 the murmur
  // sits within 2 dB of bdl in five bands (111 nasals, a0001-a0030).
  double sinus_neck_scale = 0.25;
  std::uint32_t seed = 73021U;
  // The source: false, the self-oscillating triangular glottis; true, the
  // kinematic folds (Cord::set_kinematic) with these settings.
  bool kinematic_glottis = false;
  // The fold tube (Cord::set_fold_tube): the folds as a short tube of the
  // body's three pieces, blown by the lungs alone. Takes precedence over the
  // kinematic folds.
  bool fold_tube = false;
  FoldTubeSettings fold;
  double kinematic_amplitude = 0.0007;  // cm per sqrt(dyn/cm^2): 0.06 cm at 750 Pa
  double kinematic_phase_lag = 1.22;    // rad, upper edge behind the lower
  double kinematic_softness = 0.002;    // cm, width of the soft contact
  double kinematic_jitter = 0.17;       // rad per sqrt(s) of phase random walk
};

// A rest layout for a vowel from an MRI area profile (sections of equal
// length, glottis first), resampled onto the node grid. Materials: firm
// larynx walls for the first 2 cm, soft pharynx and mouth, lips last 1 cm.
BodyLayout vowel_layout(const std::vector<double>& section_areas, double section_length);

class VoiceBody {
 public:
  explicit VoiceBody(const BodyLayout& layout);

  void set_lung_pressure(double pascal, double move_seconds);
  // Velum: area of the velopharyngeal port, cm^2 (0 closes the nose).
  void set_velum(double area) { velum_area_ = std::max(area, 0.0); }
  // Opening of the space under the tongue blade, cm^2: only a raised tip
  // (/s sh r/) opens it; for vowels the tip rests behind the lower teeth.
  void set_under_tongue(double area) { under_tongue_area_ = std::max(area, 0.0); }
  // Where the nasal port and the pocket under the tongue open (cm from the
  // glottis), for a body whose geometry moves them.
  void set_velum_place(double place);
  // The node where the lip opening is (a body whose open end flares past the
  // lips); jets there are lip jets. -1: the tube's last active node.
  void set_lip_end(Index node) { lip_end_ = node; }
  Index lip_end() const { return lip_end_; }
  void set_under_tongue_place(double place);
  // Release a share of the walls' support stiffness over nodes first..last of
  // the vocal tube (muscles letting go), applied at the next control tick.
  void set_wall_release(Index first, Index last, double share) {
    release_first_ = first;
    release_last_ = last;
    release_share_ = std::clamp(share, 0.0, 0.95);
  }
  // Muscles pull the walls of nodes first..last outward by `offset` cm (the
  // support frame's anchor moves; the tissue follows through its own mass
  // against the air pressure). Applied at the next control tick.
  void set_wall_pull(Index first, Index last, double offset) {
    pull_first_ = first;
    pull_last_ = last;
    pull_offset_ = offset;
  }
  // Where the incisors are (cm from the glottis): the obstacle a jet strikes.
  void set_teeth(double start, double end) {
    teeth_start_ = start;
    teeth_end_ = end;
  }
  // Intended pitch, Hz: sets the stretching muscle's open-loop tension; the
  // heard period then trims it.
  void set_pitch(double hz) {
    intended_pitch_ = hz;
    cord_.set_pitch(hz);
  }
  double intended_pitch() const { return intended_pitch_; }
  double pitch_correction() const { return pitch_correction_; }
  double lung_pressure() const { return lung_; }  // dyn/cm^2
  Cord& cord() { return cord_; }
  Tube& vocal() { return vocal_; }
  Tube& trachea() { return trachea_; }
  Tube& nose() { return branch_tubes_[static_cast<std::size_t>(nose_)]; }
  double velum_port() const { return velum_now_; }  // cm^2, the nasal port as it is now

  // One state step; returns the radiated sound (arbitrary units) at the
  // step rate.
  double step();

 private:
  // A closed or radiating side tube joined to a node of a parent tube.
  struct Branch {
    int parent = 0;          // 0 vocal tube, 1 nose
    int port = 0;            // side-port index on the parent
    Index node = 0;          // parent node the port opens from
    double port_area = 0.0;  // cm^2; -1: the velum sets it, -2: the tongue tip sets it
    int tube = 0;            // index into branches_
    int group = 0;           // 0 always stepped, 1 behind the velum, 2 behind the tongue tip
  };

  // State of one jet-noise shaping filter.
  struct JetNoise {
    double low = 0.0;
    double lower = 0.0;
    double band = 0.0;
    double second = 0.0;
  };
  void control_tick();
  void hear_pitch();
  double aspiration_sample();
  double tremor_white();
  void move_branch_group(int group, double place);
  void place_jet();
  void sound_jet();
  double jet_gauss();
  double shaped_jet_noise(JetNoise& state, double velocity, double width);
  double tremor_band(int channel);
  void add_branch(int parent, Index node, double port_area, Tube tube);
  void exchange_branches();

  Tube trachea_;
  Tube vocal_;
  Cord cord_;
  std::vector<Tube> branch_tubes_;
  std::vector<Branch> branches_;
  int nose_ = -1;               // index of the nose in branch_tubes_
  double velum_area_ = 0.0;
  double under_tongue_area_ = 0.0;
  double velum_now_ = 0.0;
  double velum_guide_ = 0.0;  // first of the velum glide's two stages
  // Steps since the velum / the space under the tongue last opened: the
  // cavities behind a port shut for 50 ms hold no sound worth computing.
  long velum_quiet_ = 0;
  long under_quiet_ = 0;
  std::vector<int> tube_group_;  // per branch tube: its group
  double under_tongue_now_ = 0.0;
  double under_tongue_guide_ = 0.0;
  double skin_gain_ = 1.0;
  double skin_pole_gain_ = 1.0;  // one skin low-pass pole's gain per step
  BodyLayout settings_;
  Variation variation_;
  long tick_ = 0;
  // Glottal turbulence (dipole 1.5 cm above the glottis).
  Index aspiration_node_ = 4;
  double aspiration_flow_ = 0.0;      // source amplitude low-passed at 500 Hz
  double aspiration_root_ = 0.0;      // sqrt of the source amplitude, 3 ms rise
  // Tremor channels: 0 stretching muscle, 1 adductors, 2 breathing muscles.
  static constexpr int kTremorChannels = 3;
  // White noise through two one-pole low-passes at 2 Hz and two at 15 Hz;
  // the band is the difference of the two second stages.
  double tremor_slow_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_slower_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_fast_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_faster_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_value_[kTremorChannels] = {0.0, 0.0, 0.0};
  double tremor_scale_ = 1.0;  // makes the band unit rms
  std::uint32_t tremor_random_ = 0x2545F491U;
  double glottal_flow_previous_ = 0.0;
  double deceleration_ = 0.0;  // smoothed -dU/dt of the glottal flow, cm^3/s^2
  double flow_peak_ = 0.0;     // recent peak glottal flow, cm^3/s
  double noise_low_ = 0.0;     // aspiration shaping: the band below corner / 3 removed
  double noise_lower_ = 0.0;
  double noise_z1_ = 0.0;     // and two low-pass stages at the corner
  double noise_z2_ = 0.0;
  double flux_z1_ = 0.0;      // the flux noise's own two low-pass stages
  double flux_z2_ = 0.0;
  double area_floor_ = 1.0;   // running minimum of the glottal area (15 ms release), cm^2
  std::uint32_t noise_random_ = 0x9E3779B9U;
  // Pitch heard from the closure instants (contact onset).
  double intended_pitch_ = 120.0;
  double pitch_correction_ = 0.0;  // log share added to tension
  bool was_closed_ = false;
  long last_closure_ = -1;
  // Oral jet turbulence: where the jet forms (narrowest rest section past the
  // pharynx, or -1), its source weight on each link downstream, and the
  // shaping filter of the shared turbulence signal.
  Index release_first_ = 0;
  Index release_last_ = -1;
  double release_share_ = 0.0;
  Index pull_first_ = 0;
  Index pull_last_ = -1;
  double pull_offset_ = 0.0;
  double teeth_start_ = 16.1;
  double teeth_end_ = 16.6;
  Index jet_narrowest_ = -1;
  std::vector<double> jet_weights_;           // free-jet turbulence per link
  std::vector<double> jet_obstacle_weights_;  // the jet striking the incisors per link
  JetNoise jet_free_;
  JetNoise jet_obstacle_;
  std::uint32_t jet_random_ = 0x6C078965U;
  double lung_flow_slow_ = 0.0;  // flow from the lungs below 20 Hz, cm^3/s
  // Muscle causes interpolated across each control tick: 0 cord tension
  // factor, 1 fold rest offset (cm), 2 breath factor.
  static constexpr int kCauses = 3;
  double cause_from_[kCauses] = {1.0, 0.0, 1.0};
  double cause_to_[kCauses] = {1.0, 0.0, 1.0};
  double breath_factor_ = 1.0;   // variation and tremor of breath, gliding (2 ms) to its 1 kHz target
  double lung_ = 0.0;          // dyn/cm^2
  double lung_target_ = 0.0;
  double lung_ramp_ = 0.0;       // the command ramped across the tick
  double lung_ramp_from_ = 0.0;
  int lung_ramp_left_ = 0;
  double lung_gain_ = 1.0;
  double jet_strength_ = 0.0;  // the jet noise's strength after its attack
  double from_cord_below_ = 0.0;
  double from_cord_above_ = 0.0;
  double previous_flow_ = 0.0;
  Index lip_end_ = -1;
  double skin_low_ = 0.0;    // the lining's flow through the first of two skin low-pass poles
  double skin_lower_ = 0.0;  // and through the second
};

}  // namespace tract
