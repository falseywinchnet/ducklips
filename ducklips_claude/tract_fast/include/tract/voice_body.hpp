#pragma once

// The connected system: lungs -> trachea (tube of ring modules) -> cord
// module -> vocal tube (modules from the epilarynx to the lips) -> radiation.
// The lungs are the only sustained energy source.

#include "tract/cord.hpp"
#include "tract/tongue.hpp"
#include "tract/tube.hpp"
#include "tract/variation.hpp"

#include <algorithm>
#include <limits>
#include <vector>

namespace tract {

struct BodyLayout {
  std::vector<double> trachea_areas;   // cm^2 per node, lungs first
  std::vector<double> vocal_areas;     // cm^2 per node, glottis first
  std::vector<int> vocal_materials;    // index into the vocal materials, per node
  std::vector<SegmentMaterial> vocal_material_list;
  double wall_stiffness_scale = 1.0;   // every vocal wall's support stiffness times this
  double wall_loss_scale = 1.0;        // ... and its support loss
  double skin_gain = 1.0;              // share of the lining's motion heard through the skin
  double skin_corner = 250.0;          // Hz: two poles; neck and cheek tissue pass only the lows
  // Hz: the recording chain's low cut, a two-pole Butterworth high-pass after
  // the room (speak); 0: none. (bdl's ARCTIC readings sit 14-19 dB under
  // their 150-300 Hz band from 10 to 90 Hz, ours 0-5: our voiceless sounds'
  // slow flow swings, 2026-10-18.)
  double recording_highpass = 0.0;
  // The recording microphone's presence rise (speak's room; dB, Hz, Q; 0 dB:
  // none). (Every vowel of ours sat 7-16 dB under bdl at 3.0-3.4 kHz whatever
  // the folds, losses or epilarynx; speech microphones lift about there.)
  double recording_presence_db = 0.0;
  double recording_presence_hz = 3500.0;
  double recording_presence_q = 0.7;
  double frame_corner = 250.0;          // Hz, the articulators' pushed air (Tube::set_frame_corner)
  double support_log_glide = 0.0;      // s; walls' support glides in its logarithm (Tube::set_support_log_glide); 0: linear, 2 ms
  // A side port that moves along the vocal tube with the articulation (the
  // velum, the space under the tongue tip, the cheeks) hands its area from
  // the node it leaves to the next over about this many seconds; 0: it jumps
  // (as before). Measured 2026-10-18: the velum's port stepped from node 21
  // to 20 as it lowered for "uniform"'s /n/, the nasal flow (about 25
  // cm^3/s) landed on another node in one step and the click the owner
  // heard spread from there (NATO 22.323 s).
  double branch_glide = 0.0;
  // Warmth, after the Canopy Falls harp's W2 (playsuite-games 912b87f; the
  // owner: "w2 is nice. really nice", and here: "its the 1/n^2"): heard
  // mostly as the voice's own swing, the flow leaving the lips (harmonics
  // falling as the flow's, about 1/n^2), rather than its rate of change
  // (radiation's, one power of n brighter). The flow is level-matched to
  // the rate of change at the intended pitch (times 2 pi f) and freed of
  // the breath's steady stream by two high-pass poles at warmth_corner.
  // 0: the rate of change alone (as before); the harp keeps 0.85.
  double warmth = 0.0;
  double warmth_corner = 50.0;         // Hz
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
  // A listening ear (ear_latency > 0): auditory feedback as people have it.
  // Each heard period's log error is smoothed over ear_smoothing seconds,
  // reaches the muscle ear_latency later, and moves the correction at
  // ear_rate per second (an integrator). 0 latency: the cycle-by-cycle ear
  // above (pitch_feedback per closure, no delay), which chased every period
  // the moving tongue perturbed: the owner heard a warble, gone with it off.
  double ear_latency = 0.0;
  double ear_smoothing = 0.1;
  double ear_rate = 3.0;
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
  // The reversible jet (0: off): each link downstream of the narrowest node
  // is a turbulence source driven only by its own flow of the step before,
  // the narrowest area and a counter-hashed noise, so a peel from the lips
  // reaches every source knowing what drives it (VoiceBody::reversible_jet).
  // The value is its gain, as `frication`.
  double reversible_jet = 0.0;
  // With a seal its jet fades out as the narrowest node closes: none at a
  // quarter of this area (cm^2), full from it (0: no fade). Each link's
  // velocity is its own flow over the narrowest node's area, so with that
  // node sealed (1e-4 cm^2) the cavity in front's small flows read as huge
  // jets and the render ran away (/aa k! aa s! aa).
  double reversible_jet_closing = 0.0;
  // A factor on the velum port's viscous resistance (1: the viscous law's).
  double port_resistance_scale = 1.0;
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
  // The nose's hop loss (Speaker: Tube::set_wave_smoothing on the nose; 0.4
  // was fixed there) and the vocal tube's (below 0: the Speaker's own, 0.02
  // for the articulatory body). (bdl4: nasal murmurs too bright, the vowel-
  // to-murmur drop at 0.5-2 kHz 5 dB against bdl's 10.)
  double nose_smoothing = 0.4;
  double vocal_hop_loss = -1.0;
  // Experiments on the tube (2026-10-16; 0 or 1 leaves each off):
  // the hop loss follows each link's cross-section, as boundary layers do
  // (Tube::set_smoothing_by_area; 1 on);
  double hop_loss_by_area = 0.0;
  // the cheek pockets (buccal vestibules) beside the teeth, both together,
  // as closed side branches moving with the teeth: their cross-section,
  // port and volume times this (0 off; 1: 3 cm deep, 0.8 cm^2, port 0.4);
  double buccal_scale = 0.0;
  // the velum's least opening while the mouth is open (cm^2): the
  // velopharynx seldom seals fully in vowels; it closes under the oral
  // pressure of an obstruent (fading as a constriction narrows to 0.05 cm^2);
  double velum_leak = 0.0;
  // Below the trachea, the bronchial tree (subglottal.*): bronchial_length
  // cm of tube whose total cross-section grows (each generation branches)
  // from the trachea's to bronchial_area at the lungs; the lungs' end is a
  // port lung_port_scale times wider than the last section (a large, soft
  // volume that reflects sound back up, as an open end does). 0 / 1: the
  // trachea ends in a matched (reflectionless) port, as before.
  // A factor on the vocal tube's viscous resistance: the law (A^-2.6 when
  // narrow) is a steady-flow one; at formant frequencies the boundary
  // layer's loss in a 0.2 cm^2 channel is about a fifth of it.
  double viscous_scale = 1.0;
  PinchSettings pinch;  // the pinch tube in the cord's place (pinch.on)
  double bronchial_length = 0.0;   // cm
  double bronchial_area = 20.0;    // cm^2, total, at the lungs
  double lung_port_scale = 1.0;
  // the lips' radiation load: the end correction's air mass and the
  // resistance, times these (1, 1: a piston in an infinite baffle; a head
  // radiates more like an unflanged end, about 0.72 and 1.04);
  double radiation_mass_scale = 1.0;
  double radiation_resistance_scale = 1.0;
  // The tongue as a body moved by its muscles (tract/tongue.hpp; 1 on, 0 the
  // geometric tongue the planner's targets place directly).
  double physical_tongue = 0.0;
  TongueSettings tongue;
  // With the physical tongue off: 1 draws the tongue as VocalTractLab does
  // (tract/articulation.hpp VtlTongueShape), 0 the old body circle and tip.
  // Every length of the midsagittal anatomy times this (a smaller speaker:
  // a shorter tract; MidsagittalTract::set_scale). 1: as drawn.
  double anatomy_scale = 1.0;
  // The epilarynx: its height (cm; below 0 the anatomy's own 2.3) and its
  // cross-section's scale (MidsagittalTract::set_larynx_height / _area_scale).
  double epilarynx_height = -1.0;
  double epilarynx_area = 1.0;
  double lateral_area = 0.12;          // cm^2, an /l/'s side channels (MidsagittalTract::set_lateral_area)
  double lateral_reach = 2.0;          // cm behind the tip they run (MidsagittalTract::set_lateral_reach)
  // An /l/'s closed midline space over the blade, behind the tip contact: a
  // side branch whose quarter wave puts the lateral's zero near 2 kHz. Its
  // port (cm^2, at full lateral; 0: no pocket), its length and cross-section
  // (cm, cm^2), and how far behind the tip it joins (cm). (bdl's /l/ dips to
  // -39 to -49 dB re his peak at 1.4-2.5 kHz; ours, one narrowed tube with
  // no branch, -28 to -35, 2026-10-19.)
  double lateral_pocket = 0.0;
  double lateral_pocket_length = 4.0;
  double lateral_pocket_area = 0.6;
  double lateral_pocket_back = 3.0;
  double vtl_tongue = 0.0;
  double vtl_body_radius = 1.8;   // cm
  double vtl_tip_radius = 0.2;    // cm
  double vtl_blade_along = 0.61;
  double vtl_blade_lift = 1.15;   // cm
  double vtl_root_bulge = 0.3;    // cm
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
  // VocalTractLab 2019's skew, double pulsing, flutter and length scaling
  // (Cord::set_kinematic_shape); all 0: as before. Their most natural male:
  // skew 0, double pulsing 0.05.
  double kinematic_skew = 0.0;
  double kinematic_double_pulsing = 0.0;
  double kinematic_flutter = 0.0;
  double kinematic_length_scaling = 0.0;
  // VocalTractLab 2.4's "Geometric glottis 2025" on the kinematic folds
  // (Cord::set_geometric_2025): 1 on; its pulse shape and relative amplitude.
  double geometric_2025 = 0.0;
  double pulse_shape = 0.5;
  double relative_amplitude = 1.0;
  double contact_softness_2025 = 0.0;   // cm; 0: cut at the midline (theirs)
  // Cycle-to-cycle variation of the 2025 folds (Cord::set_cycle_variation): sd
  // of the pulse shape, the lag (rad), the swing (share), the rests (cm).
  double cycle_shape = 0.0;
  double cycle_lag = 0.0;
  double cycle_amplitude = 0.0;
  double cycle_rest = 0.0;
  // The 2025 folds' size (Cord::set_geometry_2025): JD2's male defaults.
  // The folds stop swinging as the arytenoids abduct: full swing at the
  // planner's gap up to abduction_fade_start (cm), none from
  // abduction_fade_end (0: no fade, as before). Measured 2026-10-18 (bdl
  // profile): with the 2025 folds swinging at every gap, /f/ carried +28 dB
  // and /hh/ +17 dB of low bands against bdl (the swing follows the
  // commanded breath, which an open glottis does not lower).
  double abduction_fade_start = 0.0;
  double abduction_fade_end = 0.0;
  double fold_length_2025 = 1.6;
  double fold_thickness_2025 = 0.45;
  double fold_rest_pitch_2025 = 120.0;
  // VocalTractLab 2.4's glottal flow (Cord::set_vtl_flow) and aspiration
  // (VoiceBody::vtl_aspiration_sample, in place of ours, the flux noise and
  // the turbulence on the folds). 0: off.
  double vtl_glottal_flow = 0.0;
  double vtl_inertance_floor = 1.0e-4;  // cm^2
  double vtl_aspiration = 0.0;
  double vtl_aspiration_gain = 1.0;     // times their 0.2e-6
  double vtl_aspiration_corner = 1000.0; // Hz, their noise filter's corner (and gain (1000 / f_c)^1.5)
  double vtl_area_limit = 0.0;          // 1: VocalTractLab's area-rate limit on the vocal tube and nose (Tube::set_area_limit)
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
  // The reversible jet's gate (a factor on its gain, 1 by default): frication
  // as a gesture, on for a fricative and off for the vowels around it.
  void set_jet_gate(double gate) { jet_gate_ = gate; }
  // JETBOUND: the reversible jet's speed bounded by the fastest the lungs'
  // pressure can drive air, sqrt(2 p / rho) (cm/s; 0: unbounded). The jet's
  // speed at a link is its flow over the narrowest area; as the narrowest
  // closes under the jet's gate that ratio, cubed in the strength, ran away
  // (Oscar's smooth /s k/: non-finite at 0.2573 s).
  void set_jet_bound(double speed) { jet_bound_ = speed; }
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
  // The lateral pocket's port (cm^2) and where it joins (cm from the glottis).
  void set_lateral_pocket(double area) { lateral_pocket_area_ = std::max(area, 0.0); }
  void set_lateral_pocket_place(double place);
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
  void set_teeth(double start, double end);
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
  Index nose_port_node() const;  // the vocal node the nose opens from
  double under_tongue_port() const { return under_tongue_now_; }  // cm^2, the sublingual port as it is now

  // One state step; returns the radiated sound (arbitrary units) at the
  // step rate.
  double step();

  // The folds driven by a given tract (diagnostic, for the reversible
  // engine's equation-error test): with forcing on, the folds read `wave`
  // as the wave the vocal tube returns toward the glottis this step instead
  // of the tube's own, and forced_flow() is the flow they then pass into
  // the tube, (sent - wave) A0 / rho c. The tube still steps; nothing it
  // does reaches the folds. Off by default.
  void force_return(bool on, double wave) {
    forcing_ = on;
    forced_wave_ = wave;
  }
  double forced_flow() const { return forced_flow_; }
  // An observer on the forced folds: with `gain` above 0 and a finite
  // `implied` flow (the flow the sound implies this step), the folds feel
  // the forced wave plus gain (rho c / A0) (forced flow - implied flow): a
  // mismatch pulls them back toward the flow the sound implies, so a fold
  // phase left wrong by a stretch they ran free through decays instead of
  // persisting. At the truth the two flows are one and nothing changes;
  // forced_flow() stays the flow the forced wave alone implies.
  void observe_flow(double gain, double implied) {
    observe_gain_ = gain;
    observed_flow_ = implied;
  }

 private:
  // A closed or radiating side tube joined to a node of a parent tube.
  struct Branch {
    int parent = 0;          // 0 vocal tube, 1 nose
    int port = 0;            // side-port index on the parent
    Index node = 0;          // parent node the port opens from
    double port_area = 0.0;  // cm^2; -1: the velum sets it, -2: the tongue tip sets it
    int tube = 0;            // index into branches_
    int group = 0;           // 0 always stepped, 1 behind the velum, 2 behind the tongue tip
    // A moving port's crossfade (BodyLayout::branch_glide): the node it
    // left keeps a shadow port with this share of the area, fading out.
    int shadow_port = -1;
    Index shadow_node = 0;
    double shadow_share = 0.0;
  };

  // State of one jet-noise shaping filter.
  struct JetNoise {
    double low = 0.0;
    double lower = 0.0;
    double band = 0.0;
    double second = 0.0;
    // The gains of the last corner this filter was given (a memo: the
    // corner sits clamped at 1500 Hz or the ceiling most of the time, and
    // the per-link exps were about a tenth of a jet word's run, sample
    // 2026-10-20; the same corner gives the same gains, so the filter
    // evolves exactly as without it: 9% fewer instructions, MEASURED).
    double corner = std::numeric_limits<double>::quiet_NaN();
    double corner_gain = 0.0;
    double floor_gain = 0.0;
    double passed = 0.0;
  };
  // A corner's gains into `state`'s memo (the shaping filters' arithmetic).
  void jet_gains(JetNoise& state, double corner) const;
  void control_tick();
  void hear_pitch();
  double aspiration_sample();
  double vtl_aspiration_sample();
  double tremor_white();
  void move_branch_group(int group, double place);
  void place_jet();
  void sound_jet();
  double jet_gauss();
  double shaped_jet_noise(JetNoise& state, double velocity, double width);
  void reversible_jet();
  double tremor_band(int channel);
  void add_branch(int parent, Index node, double port_area, Tube tube);
  void exchange_branches();

  Tube trachea_;
  Tube vocal_;
  Cord cord_;
  std::vector<Tube> branch_tubes_;
  std::vector<Branch> branches_;
  int nose_ = -1;               // index of the nose in branch_tubes_
  bool buccal_ = false;         // the cheek pockets are present (branch group 3, at the teeth)
  double oral_open_ = 1.0;      // 0..1: how open the mouth is (the velum's leak follows it)
  double velum_area_ = 0.0;
  double under_tongue_area_ = 0.0;
  double velum_now_ = 0.0;
  double jet_gate_ = 1.0;
  double jet_bound_ = 0.0;
  double velum_guide_ = 0.0;  // first of the velum glide's two stages
  // Steps since the velum / the space under the tongue last opened: the
  // cavities behind a port shut for 50 ms hold no sound worth computing.
  long velum_quiet_ = 0;
  long under_quiet_ = 0;
  double lateral_pocket_area_ = 0.0;
  double lateral_pocket_guide_ = 0.0;
  double lateral_pocket_now_ = 0.0;
  long lateral_quiet_ = 0;
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
  double vtl_noise_flow_ = 0.0;       // the glottal outflow through 500 Hz (VocalTractLab's)
  double vtl_noise_amp_ = 0.0;        // its source amplitude, dyn/cm^2
  double vtl_noise_state_[2] = {0.0, 0.0};  // the 1 kHz low-pass
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
  double heard_error_ = 0.0;        // the listening ear's smoothed log error
  std::vector<double> ear_delay_;   // its errors on their way to the muscle, one per control tick
  std::size_t ear_slot_ = 0;
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
  std::vector<JetNoise> link_jets_;    // the reversible jet's per-link shaping filters
  long jet_step_ = 0;                  // its step counter (the noise's hash key)
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
  bool forcing_ = false;
  double forced_wave_ = 0.0;
  double forced_flow_ = 0.0;
  double observe_gain_ = 0.0;
  double observed_flow_ = 0.0;
  double previous_flow_ = 0.0;
  double warm_low_ = 0.0;     // the flow through the first warmth pole (the part removed)
  double warm_lower_ = 0.0;   // and the second
  double warm_pitch_ = 120.0; // Hz, the intended pitch followed over 20 ms
  Index lip_end_ = -1;
  double skin_low_ = 0.0;    // the lining's flow through the first of two skin low-pass poles
  double skin_lower_ = 0.0;  // and through the second
};

}  // namespace tract
