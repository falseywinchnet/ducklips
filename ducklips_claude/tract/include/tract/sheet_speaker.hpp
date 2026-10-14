#pragma once

// The sheet speaker: reads the body language (voice::Score, written by the
// text front end: phrases, breath marks, phonemes, timing, pitch) and
// deforms the sheet tube to satisfy it.
//
// Every part of the score is a held coordinate of a SomaticHold (opposing
// motor pools, sensing delay, efference copy, relaxing strain memory): its
// posture is reached and held with effort, and wanders with the effort. Once
// per millisecond the held posture becomes the vocal tract's area function
// (the front end's geometry), and every tract section's ring support moves
// to it; the glottis part sets the slit cord's adduction, the pitch intent
// its stretch, the lungs the breath. Turbulence forms where a jet leaves a
// narrowing (the glottis, an oral constriction) and is placed there as a
// pressure source. Sound itself comes only from the sheet tube.

#include "tract/sheet_tube.hpp"
#include "tract/somatic_hold.hpp"
#include "voice/score.hpp"

#include <array>
#include <cstdio>
#include <cstdint>
#include <vector>

namespace tract {

// The speaking tube: Astra's subglottal sections, a slit cord, then the tract
// (17.5 cm, the front end's oral grid) as `modules` modules of three sections
// (elastic, support, mass) at the neutral 3 cm^2, walls carrying the tissue
// behind them (three times the vowel tubes' first masses: /a/ 720/1020/2960,
// /i/ 310/2140/2930, /u/ 300/770/2360 Hz).
struct SheetVoice;
SheetDesign speech_sheet_design(const SheetVoice& voice);
// The nose: velum to nostrils, 11.7 cm (the front end's nasal tract:
// nasopharynx 1.6 cm^2, main cavity 2.6, vestibule narrowing to 1.0), as
// three-section modules; it radiates at the nostrils.
// area_share and length_share scale one side's cross-sections and length.
SheetDesign nasal_sheet_design(const SheetVoice& voice, double area_share, double length_share);
// Astra's maxillary sinus (a neck and a chamber, firm) and the space under
// the tongue (an entry and a floor pocket, soft): closed side tubes.
SheetDesign sinus_sheet_design(const SheetVoice& voice, double radius_share);
SheetDesign under_tongue_sheet_design();
// The piriform recesses (the front end's): blind pockets beside the larynx,
// `depth` m deep, tapering from 0.9 to 0.4 cm^2; their quarter-wave zeros
// carve the voiced spectrum near 4-5 kHz.
SheetDesign piriform_sheet_design(double depth);
// The lateral pocket over the blade behind an /l/ contact (the front end's):
// a 0.5 cm neck into 2.5 cm^3.
SheetDesign lateral_sheet_design();

struct SheetVoice {
  // The cord's pitch at its base stretch, and F0 ~ stretch^exponent
  // (measured on /a/ and /i/ at 800 Pa: stretch 0.9-2.2 gave 98-145 Hz).
  double base_pitch = 101.5;
  double pitch_exponent = 0.42;
  double stretch_low = 0.9;   // below it the cord stops voicing steadily
  double stretch_high = 2.3;
  double pressed = -0.0001;   // m, the cord's voicing adduction (rest half-gap 0.1 mm + this)
  // Motor-pool fluctuation and common drive (held /a/ with the triangular
  // glottis: 0.1 and 0.015 gave jitter 0.25%, shimmer 1.0%, HNR 18 dB, F0
  // sd 7 cents; 0.2 and 0.03 gave HNR 8 and 22 cents of fast wander).
  double hold_noise = 0.1;
  double hold_common = 0.02;
  double frication = 1.0;     // oral jet turbulence gain (with the monopole source; series: 3)
  // The jet's turbulence along its path (the front end's): it grows over
  // jet_growth and decays over jet_decay (cm) past the constriction's exit,
  // four times stronger where a jet still coherent (within obstacle_reach)
  // strikes the incisors; a lip channel adds the labiodental weight.
  double jet_growth = 0.4;
  double jet_decay = 1.2;
  double obstacle_reach = 0.8;
  double obstacle_gain = 3.0;  // a coherent jet striking the incisors: weight x (1 + this)
  double labiodental = 0.6;
  // The oral jet noise's corner stays within these (Hz; the front end's,
  // fitted to bdl's fricatives). Up to 12 kHz, /sh ch/ rolled off below
  // about 4 kHz: weak at 1-3 kHz, too strong at 3-7 kHz.
  // The oral jet's turbulence as a local flow (monopole, m^3/s from the
  // pressure strength over the local impedance) instead of a pressure in
  // series: a series source drives the cavity behind the constriction as
  // strongly as the one ahead (ideal model: /sh/'s 1.1 kHz back resonance
  // dominated). false: series pressure (the front end's).
  // Adopted with frication 1: /sh/ 33 -> 50%, /ch/ ~20 -> 44%, /s/ 74%,
  // /s/ within 2-3 dB of bdl at 3-7 kHz; overall unchanged (29.8%, three
  // seeds). /sh/ remains 14 dB loud at 3-4.5 kHz.
  bool jet_monopole = true;
  double jet_corner_floor = 1500.0;
  double jet_corner_ceiling = 5500.0;
  double aspiration = 0.3;    // glottal jet turbulence gain
  // 96 kHz and 16 modules: the tract's 48 sections are the front end's
  // oral grid (3.6 mm). At 48 kHz and 8 modules (7.3 mm sections) speech
  // was 79-85% voiced and 2.3 dB from the owner's spectrum; here 92-96%
  // and 1.6 dB (front end engine: 89%, 2.0 dB). 1.8x real time.
  double rate = 96000.0;
  int modules = 16;
  // The cord's source settings (see SheetDesign). A 0.3 mm taper with
  // elastic contact: closure zips from the front, and the voice's long-term
  // spectrum came within 3.1 dB of the owner's (4.4 with parallel folds
  // stopped at contact: 9.5 dB too strong at 1.6-2.5 kHz).
  double cord_taper = 0.0003;
  double cord_contact = 1.0;
  double cord_loss = -1.0;
  double lung_scale = 0.9;  // with the wider subglottal airway (1.0 before it)
  // The cord section: thickness (m), mass and elasticity factors, support.
  double cord_thickness = 0.004;
  double cord_mass = 3.0;
  double cord_elasticity = 0.4;
  double cord_support = 0.2;
  // Scales the subglottal bands' support offsets (1: Astra's, airway
  // 0.25-0.39 cm^2). 0.7 with lung_scale 0.9: H1-H2 +4.2 -> +2.0 dB (the
  // owner -3.9), spectrum 1.8 -> 1.6 dB, identification 33.6% (34.3).
  double subglottal_narrowing = 0.7;
  double lung_resistance = 1.5e5;  // Pa s/m^3; 0 leaves Astra's matched inlet alone
  double lung_corner = 20.0;  // Hz, the lung branch's R/L corner (SheetDesign::lung_corner)
  bool side_tubes = true;   // the nose, sinus and under-tongue joints
  bool lateral_pocket = false;  // the lateral pocket alone (without the piriform recesses)
  bool recesses = false;    // the piriform recesses and the lateral pocket (REJECTED as default: spectrum 1.7 -> 2.7 dB, -10 dB at 2.5-4 kHz; /l/ 16% -> 7%)
  double tissue_relaxation = 0.1;
  bool jet_loss = true;
  // The folds: the triangular glottis (true), or the sheet's own slit cord.
  bool triangular_glottis = true;
  // The folds' elongation by deflection (Cord::set_elongation, 1/cm^2).
  // 200: a held /a/ rises 121 -> 140 Hz over 390 -> 1040 Pa (about 3 Hz
  // per cmH2O; people 2-6), where it fell 150 -> 122 without. Over 29 ARCTIC
  // sentences: vowel pitch error 77 -> 69 cents, swing 120 -> 113, onset
  // error 138 -> 111, voicing 88 -> 89%. (300: the same pitch figures,
  // identification 31.4% against 32.4.)
  double fold_elongation = 200.0;
  double fold_elongation_reference = 1.6e-3;  // cm^2 (2e-3 before the wider airway and softer breath)
  double fold_elongation_release = 0.01;      // s (30 ms: pitch error 80 cents)
  double fold_elongation_floor = 0.9;        // the stretch factor's floor
  // The nose's air-path taps: spread (a gentle low-pass) but summing to
  // about 1, so the steady flow passes. Flat losses (0.96, 0.005, 0.002)
  // bled it: with the mouth sealed, pressure built behind the seal, the folds
  // stopped in /n m/, and the murmur read 15 dB under bdl's at 200-400 Hz.
  // These (with velum_speed 3, nasal_wall_loss 8): /n/ 26%, /m/ 9% (0
  // before), 200-400 Hz at bdl's level.
  double nasal_taps[3] = {0.4, 0.35, 0.248};
  // Hz, the nose's air-path low-pass (0: off). /n/ re bdl at 400-1500 Hz was
  // 11-17 dB too bright without it and 3-14 dB too dark at 2000; the lossy
  // taps replaced it.
  double nasal_lowpass = 0.0;
  // The nose's wall masses, scaled (bone under mucosa; x4-x15 changed
  // nothing: the murmur's peak and identification stayed).
  double nasal_mass = 1.0;
  double nasal_wall_loss = 8.0;  // the nose walls' damping, scaled (widens its resonances, keeps the flow)
  // Two passages behind the velum (the septum), 55% and 45% of the
  // single nose's cross-sections, the second 6% longer, each with its own
  // sinus (the second's chamber 1.25 times as wide): their unequal
  // resonances and zeros damp the murmur above 400 Hz.
  bool split_nose = false;
  double nasal_length = 0.117;  // m, velum to nostrils (the front end's 11.7 cm)
  // The oral tract's air-path taps (Astra's). Lossier taps bled the steady
  // flow: voicing fell to 58-89% and the spectrum went 3.5-10 dB off.
  double oral_taps[3] = {0.985, 0.010, 0.003};
  // The airways below the folds (SheetDesign::subglottal_taps; a first value
  // below 0: the oral taps).
  double subglottal_taps[3] = {-1.0, 0.0, 0.0};
  // The tract walls' damping, scaled (x4: voiced 82-99%, spectrum 1.6 dB;
  // B1 stays 85-90 Hz).
  double wall_loss = 4.0;
  double tract_mass = 1.0;  // the oral tract walls' masses, scaled
  // cm^3/s, the throat's and cheeks' expansion behind a voiced oral closure
  // (people enlarge the cavity by a few cm^3 over a /b/; 0: off). 60 keeps
  // the folds vibrating through a closure, but weakens the bursts:
  // identification 29.8% (30: 31.3).
  double closure_expansion = 0.0;
  double closure_volume = 1.0e9;  // cm^3, the walls' give per closure (huge: unlimited)
  // The share of the tract walls' volume velocity radiated through the
  // skin of the throat and cheeks (0: the walls are silent). 1: /b/ re bdl
  // -23.6 -> -4.6 dB, /g/ -13 -> -9, spectrum 2.2 -> 2.1 dB, identification
  // 31.5% (31.8 without).
  double skin_radiation = 1.0;
  double skin_corner = 400.0;  // Hz, the skin's two-pole low-pass
  double velum_gain = 1.0;     // the velum port's area, scaled
  // Where the space under the tongue opens into the tract, cm from the
  // glottis. Ahead of /sh/'s tip (15.2 or 15.8 cm) the pocket resonated
  // near 1.2 kHz and /sh/ fell from 33% to 22-28%, /ch/ from 19% to 0-6%.
  double under_tongue_at = 13.5;
  double sinus_radius = 0.008; // m, the sinus chamber (Astra's: about 7 cm^3; a maxillary sinus is ~15)
  double control_rate = 8000.0;  // Hz, the sheet's pieces follow the moving supports
  // The body is shaped this often (Hz), its posture interpolated between the
  // motor pools' 1 ms ticks (stepped at each tick, a held /f/ rang at 1 kHz
  // and its harmonics).
  double shape_rate = 8000.0;
  int listen = 0;  // 0: mouth and nostrils; 1: the mouth only; 2: the nostrils only (diagnosis)
  // Its voicing posture: rest displacements of the folds' lower and upper
  // edges at the vocal processes (cm), as the segment-module speaker's.
  double voicing_lower = 0.004;
  double voicing_upper = -0.008;
  double voiced_open_gain = 0.7;  // the opening per cm of glottis setting below 0.03 (voiced fricatives)
  // The ear: at each closure the heard period corrects the folds' tension
  // by this share of log(intended / heard), within the limit (the
  // segment-module body's; the glottis's open-loop pitch was calibrated in
  // another tube: 141 Hz for an intended 120 here).
  // (0.08 and 0.3 left our vowels a median 157 cents off the plan, with
  // 113 cents of swing inside a vowel: tone-like, "vietnamese"; the front
  // end's 79 and 80. 0.2 and 1.2: 89 and 79 cents on the two sentences.)
  // 0.12 once the ear was gated to phonation: at 0.2 its corrections
  // chased the nose's pull on the folds next to nasals (/ow/ in "hello"
  // swung 498 cents, 144 with the ear off); 0.12: swing 142 -> 106 cents
  // on sentence 1 at the same 57-cent accuracy.
  double pitch_feedback = 0.12;
  double pitch_correction_limit = 1.2;
  // Stiffening for a voiceless consonant (glottis wider than voicing,
  // narrower than devoicing_widest, below which breaths are excluded):
  // tension x exp(raise), relaxing with devoicing_relax once the glottis
  // closes. bdl's vowels after voiceless sounds begin about +125 cents over
  // their median in the first 15 ms and +20 by 15-30 ms.
  // 0.25 with 40 ms: identification 31.8% (31.5 without; 0.15: 31.7).
  double devoicing_raise = 0.25;
  double devoicing_relax = 0.04;   // s
  double devoicing_widest = 0.3;
  // Motor commands lead their acoustic events (EMG precedes movement): each
  // event acts this long before its time.
  double lead = 0.0;  // four-view identification: 23.4% at 0, 18.8% at 0.02, 8.4% at 0.04
  double fastest = 400.0;  // 1/s, the fastest held movement
  // The oral articulators (tongue, tip, dorsum, lips, jaw) move this many
  // times the score's rate: a constriction's area goes as the square of its
  // root, so the last percent of travel sets a channel, and at the score's
  // rate /s/ reached its 0.08 cm^2 channel only as it ended.
  double articulator_speed = 1.5;
  double glottis_speed = 1.0;  // the same for the arytenoids (the glottis part)
  // The velum: at the score's rate (100 ms moves, starting ~27 ms ahead) it
  // was barely open (0-0.1 cm^2) as each nasal began and reached 0.3-0.4
  // only at its end: a sealed mouth and a shut nose, heard as /d b/.
  double velum_speed = 3.0;
  // A falling breath (the lungs' target dropping) moves this many times the
  // score's rate: words ended a median 49 ms after their last phone, the
  // front end's engine at it (-30 dB re peak).
  double breath_release_speed = 3.0;  // owner NATO 10.44 -> 10.20; x6 cut voicing to 80%
  // The breath is commanded this much earlier than the articulators: our
  // words began a median 65 ms after their first phone, the front end's 30
  // (the folds need pressure before they sound). Owner NATO 10.20 -> 9.30,
  // identification 26.1 -> 25.5% (seed range about 1).
  double laryngeal_lead = 0.03;
  // A voicing onset (the glottis closing toward phonation) is commanded this
  // much earlier, and moves voicing_speed times faster. The score's onsets
  // were timed for the old engine's quick folds: after /k/ in "speak with"
  // the glottis stayed 0.14 cm open through all of /w/ (whispered, nearly
  // silent: "with" garbled). Openings for voiceless sounds keep their time.
  double voicing_lead = 0.03;  // identification 31.3 -> 32.7% (three seeds); 0.05: 31.1
  double voicing_speed = 2.0;
  // A closure is aimed past contact, pressing (this root, sqrt cm^2, below
  // zero): contact is then made on the way, not approached forever (the
  // score's area 0 left the tip of /n/ sealing only as /n/ ended: heard /d/).
  double closure_press = 0.6;
  // Voiced closures (the segment-module speaker's): with the mouth sealed
  // above the pharynx and the folds set to vibrate, the pharyngeal walls
  // (3.5-10 cm) are drawn outward by up to this (m), over ~50 ms, and back
  // over ~40 ms once the mouth opens: the growing cavity keeps air flowing
  // across the folds (bdl voices /b d g/ about 61% of their closures; ours
  // /b/ 12%, /g/ 38%, /d/ 60%).
  double closure_pull = 0.0;
  // ... and their ring support is released toward (1 - this) (the
  // segment-module speaker's wall_release).
  double closure_release = 0.0;
  // An opening constriction (from closed or pressed toward open) moves this
  // many times faster: tissue under pressure parts abruptly. Pressed past
  // contact and released at the score's rate, /ch/ in "chair" opened 30 ms
  // after the score's release and its frication was cut in half.
  double release_speed = 3.0;  // identification 28.9 -> 29.4% at x2 (/ch/ 6 -> 25%), 29.9% at x4
};

class SheetSpeaker {
 public:
  SheetSpeaker(const voice::Score& score, std::uint32_t seed, const SheetVoice& voice);

  // Streaming: a further event, placed in time order among those not yet
  // reached (one already past acts at once).
  void push(const voice::PartEvent& event);
  void set_duration(double seconds) { score_.duration = seconds; }
  double time() const { return time_; }
  const SheetTube& tube() const { return tube_; }
  const SheetDesign& design() const { return design_; }
  bool finished() const { return time_ >= score_.duration; }
  // One step; returns the radiated sound (the rate of change of the
  // radiated flow, m^3/s^2).
  double step();
  bool fault() const {
    return tube_.fault() || nose_.fault() || sinus_.fault() || nose_right_.fault() || sinus_right_.fault() || under_.fault() || piriform_left_.fault() ||
           piriform_right_.fault() || lateral_.fault();
  }
  SheetTube& tube() { return tube_; }
  // One line of state: time, held lungs (Pa), glottis (cm), pitch (Hz),
  // glottal flow (cm^3/s) and area (cm^2), pressure below the cord (Pa),
  // narrowest tract area (cm^2) and its node.
  void trace(std::FILE* out) const;

 private:
  struct Turbulence {
    double low = 0.0;
    double lower = 0.0;
    double band = 0.0;
    double second = 0.0;
  };

  void add_coordinate(voice::Part part, int which, double origin, double scale, double low, double high,
                      int group);
  int tract_node_at(double cm) const;
  int coordinate(voice::Part part, int which) const {
    return coordinates_[static_cast<std::size_t>(2 * static_cast<int>(part) + which)];
  }
  double held(voice::Part part, int which) const;
  void apply(const voice::PartEvent& event);
  void shape_body();
  void sound_glottis();
  void sound_constriction();
  double turbulence(Turbulence& state, double velocity, double width, double lowest, double highest);
  double gauss();

  voice::Score score_;
  SheetVoice voice_;
  SheetDesign design_;
  SheetTube tube_;
  SheetTube nose_;
  SheetTube sinus_;
  SheetTube nose_right_;   // the second passage of a split nose
  SheetTube sinus_right_;
  SheetTube under_;
  SheetTube piriform_left_;
  SheetTube piriform_right_;
  SheetTube lateral_;
  SheetJoint piriform_left_joint_;
  SheetJoint piriform_right_joint_;
  SheetJoint lateral_joint_;
  Cord glottis_;
  SheetJoint velum_joint_;
  SheetJoint sinus_joint_;
  SheetJoint velum_right_joint_;
  SheetJoint sinus_right_joint_;
  SheetJoint under_joint_;
  double previous_nasal_flow_ = 0.0;
  double closure_relaxation_ = 0.0;  // 0..1, how far the pharyngeal walls have let go
  double expansion_ = 0.0;           // m^3/s, the walls' expansion behind a voiced closure
  double expanded_volume_ = 0.0;    // m^3, how far the walls have given in this closure
  int pharynx_node_ = 0;             // where the expansion is drawn (5 and 11 cm)
  int mouth_node_ = 0;
  double previous_volume_ = 0.0;     // m^3, the tract's volume (skin radiation)
  double previous_wall_flow_ = 0.0;  // m^3/s
  double skin_low_[2] = {0.0, 0.0};  // m^3/s, the walls' flow through the skin's low-pass
  double skin_gain_ = 0.0;
  void hear_pitch();
  double intended_pitch_ = 120.0;
  double pitch_correction_ = 0.0;  // log share of tension
  double devoicing_ = 0.0;         // log share of tension, the stiffening after a voiceless consonant
  long last_closure_ = -1;
  bool was_closed_ = false;
  double target_narrowest_ = 0.0;  // cm^2, the area function's narrowest (trace)
  double jet_velocity_ = 0.0;      // m/s, the oral jet (trace)
  double teeth_start_ = 0.0;       // cm from the glottis
  double teeth_end_ = 0.0;
  std::vector<double> jet_weights_;  // per tube node, the jet's turbulence weight
  std::vector<double> node_center_;  // per tube node, cm from the glottis (tract nodes)
  double jet_strength_ = 0.0;      // Pa, its turbulence's strength (trace)
  SomaticHold hold_;
  std::array<int, 2 * static_cast<std::size_t>(voice::kPartCount)> coordinates_{};
  std::array<double, SomaticHold::kMaxCoordinates> previous_{};  // held values at the last tick
  double blend_ = 1.0;  // 0..1 from the last tick's values to this tick's
  int per_shape_ = 12;
  int cord_band_ = -1;
  int first_tract_band_ = -1;
  int first_tract_node_ = 0;
  std::vector<double> band_center_;  // cm from the glottis, per band (tract bands)
  std::size_t next_event_ = 0;
  std::vector<char> applied_;  // per score event, applied early (laryngeal lead)
  long tick_ = 0;
  int per_control_ = 48;
  double time_ = 0.0;
  double previous_flow_ = 0.0;
  std::uint32_t random_ = 1U;
  Turbulence glottal_noise_;
  Turbulence oral_noise_;
};

}  // namespace tract
