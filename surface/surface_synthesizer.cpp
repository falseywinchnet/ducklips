#include "surface_synthesizer.hpp"
#include <algorithm>
#include <cmath>

namespace ducklips {
namespace {

constexpr double critical_reynolds = 1800.0;
constexpr double minimum_cycle_amplitude = 2.0e-6; // m^2 of glottal area swing for a counted cycle

// Dynamic pressure of a turbulent jet that has crossed the critical Reynolds
// number; zero below it.
double turbulence_pressure(double jet_velocity, double hydraulic_diameter, const AirProperties& air) noexcept {
    if (jet_velocity <= 0.0 || hydraulic_diameter <= 0.0) {
        return 0.0;
    }
    double kinematic_viscosity = air.viscosity / air.density;
    double reynolds = jet_velocity * hydraulic_diameter / kinematic_viscosity;
    if (reynolds <= critical_reynolds) {
        return 0.0;
    }
    double ratio = critical_reynolds / reynolds;
    double gate = 1.0 - ratio * ratio;
    double pressure = 0.5 * air.density * jet_velocity * jet_velocity * gate;
    return pressure;
}

// Three cascaded one-pole low-passes (unity gain at low frequency): turbulent
// noise falls steeply above its peak (MEASURED against a human recording:
// two poles left 16-22 kHz 19 dB too strong).
double shape_noise(std::array<double, 3>& stages, double coefficient, double input) noexcept {
    stages[0] += coefficient * (input - stages[0]);
    stages[1] += coefficient * (stages[0] - stages[1]);
    stages[2] += coefficient * (stages[1] - stages[2]);
    return stages[2];
}

float finite_sample(double value) noexcept {
    if (!std::isfinite(value)) {
        return 0.0f;
    }
    double clamped = std::clamp(value, -1.0, 1.0);
    float sample = static_cast<float>(clamped);
    return sample;
}

double sum_of(std::span<const double> values) noexcept {
    double total = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
        total += values[index];
    }
    return total;
}

}

void SurfaceSynthesizer::prepare(const VoiceProfile& voice, const SurfaceSynthesisOptions& options) {
    voice_ = voice;
    voice_.folds.medial_roundness = options.medial_roundness;
    voice_.folds.contact_softness = options.contact_softness;
    options_ = options;
    air_ = options.air;
    sample_period_ = 1.0 / options.rate.samples_per_second;
    internal_rate_ = options.rate.samples_per_second * static_cast<double>(options.oversample);
    internal_period_ = 1.0 / internal_rate_;
    block_internal_ = options.control_block * options.oversample;
    vowels_.build(voice_, options.oral_sections, options.nasal_sections, air_);
    model_.prepare(voice_, options.oral_sections, options.nasal_sections);
    model_.set_floor(options.floor);
    model_.allocate(shape_);
    model_.build_vowel(vowels_.shape(Phone::AX), shape_);

    // The cavity surface: the airway, the nose from the velar port, and the
    // side cavities from the stations they open beside.
    CavityGeometry geometry{};
    TubeGeometry airway{};
    airway.stations = options.oral_sections;
    airway.length = sum_of(std::span<const double>(shape_.oral_length));
    airway.ring_pieces = options.ring_pieces;
    airway.end = TubeEnd::radiating;
    if (options.solid_airway) {
        airway.solid = true;
        airway.bent = true;
        airway.cross_paths = true;
        airway.axial_share = options.axial_share;
        airway.cross_timing = options.cross_timing;
        airway.edge_tension = options.edge_tension;
        // The ellipses' own perimeter replaces the flattening factor.
        airway.perimeter_factor = 1.0;
    }
    geometry.tubes.push_back(airway);
    TubeGeometry nose{};
    nose.stations = 1 + options.nasal_sections;
    nose.length = shape_.port_length + sum_of(std::span<const double>(shape_.nasal_length));
    nose.ring_pieces = options.branch_ring_pieces;
    nose.parent = 0;
    nose.parent_station = model_.port_section();
    nose.end = TubeEnd::radiating;
    // Turbinates multiply the lining; most of it is mucosa on bone.
    nose.perimeter_factor = 3.0;
    nose.yielding_fraction = 0.3;
    nose_tube_ = geometry.tubes.size();
    geometry.tubes.push_back(nose);
    std::span<const SideCavity> sides = model_.side_cavities();
    side_tube_.clear();
    std::size_t side_section = 0;
    for (std::size_t index = 0; index < sides.size(); ++index) {
        TubeGeometry side{};
        side.stations = sides[index].sections;
        side.length = sum_of(std::span<const double>(shape_.side_length).subspan(side_section, sides[index].sections));
        side.ring_pieces = options.branch_ring_pieces;
        side.perimeter_factor = sides[index].perimeter_factor;
        if (sides[index].attach_node < options.oral_sections) {
            // Beside the larynx: soft-walled pockets.
            side.parent = 0;
            side.parent_station = sides[index].attach_node;
        } else {
            // From the nose: bony sinuses.
            side.parent = nose_tube_;
            side.parent_station = 1 + (sides[index].attach_node - options.oral_sections);
            side.yielding = false;
        }
        side_tube_.push_back(geometry.tubes.size());
        geometry.tubes.push_back(side);
        side_section += sides[index].sections;
    }
    // The space under the tongue: one station beside each oral section that
    // can hold it, closed at the incisors.
    floor_count_ = 0;
    if (options.floor) {
        model_.floor_sections(floor_first_, floor_count_);
    }
    if (floor_count_ > 1) {
        TubeGeometry under{};
        under.stations = floor_count_;
        under.length = 0.0;
        for (std::size_t section = floor_first_; section < floor_first_ + floor_count_; ++section) {
            under.length += shape_.oral_length[section];
        }
        under.ring_pieces = 1;
        under.parent = 0;
        under.parent_station = floor_first_;
        under.parallel = true;
        under.end = TubeEnd::closed;
        // The underside of the tongue and the floor of the mouth: all soft.
        under.yielding_fraction = 1.0;
        geometry.tubes.push_back(under);
    } else {
        floor_count_ = 0;
    }
    std::size_t stations = 0;
    for (std::size_t tube = 0; tube < geometry.tubes.size(); ++tube) {
        stations += geometry.tubes[tube].stations;
    }
    opening_.assign(stations, 0.0);
    station_area_.assign(stations, 0.0);
    write_station_areas();
    geometry.station_area = station_area_;
    AirConstants constants{};
    constants.density = air_.density;
    constants.sound_speed = air_.sound_speed;
    constants.viscosity = air_.viscosity;
    network_.build(geometry, internal_rate_, constants);

    queue_.prepare(options.cue_capacity);
    air_path_.prepare(options.air_path, air_, options.rate);
    double block_rate = options.rate.samples_per_second / static_cast<double>(options.control_block);
    jitter_noise_.prepare(voice.seed * 7 + 1, 12.0, block_rate);
    shimmer_noise_.prepare(voice.seed * 7 + 2, 8.0, block_rate);
    tremor_noise_.prepare(voice.seed * 7 + 6, 40.0, internal_rate_);
    adduction_noise_.prepare(voice.seed * 7 + 7, 20.0, internal_rate_);
    roughness_noise_.prepare(voice.seed * 7 + 8, 150.0, internal_rate_);
    roughness_wander_.prepare(voice.seed * 7 + 9, 150.0, internal_rate_);
    aspiration_noise_.prepare(voice.seed * 7 + 3, 300.0, internal_rate_);
    frication_noise_.prepare(voice.seed * 7 + 4, 300.0, internal_rate_);
    buffet_noise_.seed(voice.seed * 7 + 5);
    // A whispering glottis is held wide open: its jet is broad and its
    // turbulence brighter than the breath in a voiced sound (the human whisper
    // carries noise to about 6 kHz).
    double aspiration_corner = voice_.phonation ? options.aspiration_corner : options.whisper_aspiration_corner;
    aspiration_coefficient_ = 1.0 - std::exp(-2.0 * pi * aspiration_corner / internal_rate_);
    frication_coefficient_ = 1.0 - std::exp(-2.0 * pi * options.highest_frication_corner / internal_rate_);
    exit_smoothing_ = 1.0 - std::exp(-2.0 * pi * 1000.0 / internal_rate_);

    // Decimation kernel: pass to 0.47 of the output rate, Blackman window.
    const std::size_t taps = 64 * options.oversample + 1;
    kernel_.assign(taps, 0.0);
    double kernel_sum = 0.0;
    for (std::size_t index = 0; index < taps; ++index) {
        double centre = static_cast<double>(index) - static_cast<double>(taps - 1) / 2.0;
        double cutoff = 0.47 / static_cast<double>(options.oversample);
        double sinc = centre == 0.0 ? 2.0 * cutoff : std::sin(2.0 * pi * cutoff * centre) / (pi * centre);
        double phase = 2.0 * pi * static_cast<double>(index) / static_cast<double>(taps - 1);
        double window = 0.42 - 0.5 * std::cos(phase) + 0.08 * std::cos(2.0 * phase);
        kernel_[index] = sinc * window;
        kernel_sum += kernel_[index];
    }
    for (std::size_t index = 0; index < taps; ++index) {
        kernel_[index] /= kernel_sum;
    }
    history_size_ = taps;
    radiator_history_.assign(radiator_count * 2 * history_size_, 0.0);
    // Oldest sample first: tap t of the reversed kernel meets the sample
    // taps - 1 - t steps back.
    reversed_kernel_.assign(taps, 0.0);
    for (std::size_t index = 0; index < taps; ++index) {
        reversed_kernel_[index] = kernel_[taps - 1 - index];
    }
    unit_pitch_ = calibrate_folds();
    trace_.clear();
    trace_.reserve(options.trace_samples);
    trace_begin_ = -1;
    if (options.trace_start >= 0.0) {
        trace_begin_ = static_cast<SampleCount>(options.trace_start * internal_rate_);
    }
    reset();
}

double SurfaceSynthesizer::calibrate_folds() const noexcept {
    // Run the folds alone at Q = 1, conversational lung pressure and voicing
    // adduction, and measure their cycle rate after the onset transient. The
    // pitch controller corrects for the tract's loading while speaking.
    TwoMassGlottis folds{};
    folds.reset(voice_.folds);
    NoiseGenerator turbulence{};
    turbulence.seed(voice_.seed * 7 + 5);
    const double turbulence_scale = 0.01 * voice_.lung_pressure;
    const std::size_t total = static_cast<std::size_t>(0.5 * internal_rate_);
    const std::size_t settle = static_cast<std::size_t>(0.2 * internal_rate_);
    double last = 0.0;
    double trough = 1.0;
    bool rising = false;
    std::size_t first_peak = 0;
    std::size_t last_peak = 0;
    bool started = false;
    std::size_t cycles = 0;
    for (std::size_t sample = 0; sample < total; ++sample) {
        double surface = turbulence_scale * turbulence.gaussian();
        folds.advance(voice_.lung_pressure, voice_.adducted_area, 1.0, internal_period_, options_.glottis_substeps,
                      surface);
        double area = folds.minimum_area();
        if (area > last) {
            rising = true;
        } else if (rising && area < last) {
            rising = false;
            if (last - trough > minimum_cycle_amplitude && sample > settle) {
                if (!started) {
                    first_peak = sample;
                    started = true;
                } else {
                    cycles += 1;
                }
                last_peak = sample;
            }
            trough = last;
        }
        trough = std::min(trough, area);
        last = area;
    }
    if (cycles == 0) {
        return 120.0;
    }
    double seconds = static_cast<double>(last_peak - first_peak) / internal_rate_;
    double pitch = static_cast<double>(cycles) / seconds;
    return pitch;
}

void SurfaceSynthesizer::reset() noexcept {
    articulators_ = ArticulatorState{};
    articulators_.apply_speeds(voice_.speeds);
    articulators_.place_vowel(vowels_.shape(Phone::AX));
    articulators_.glottal_rest_area.place(voice_.abducted_area);
    articulators_.lung_pressure.place(0.0);
    articulators_.pitch.place(voice_.base_pitch);
    articulators_.effort.place(1.0);
    articulators_.loudness.place(1.0);
    articulators_.port_area.place(0.0);
    for (std::size_t slot = 0; slot < constriction_slot_count; ++slot) {
        articulators_.slot[slot].activation.place(0.0);
    }
    glottis_.reset(voice_.folds);
    network_.reset();
    air_path_.reset();
    aspiration_shaper_ = std::array<double, 3>{};
    frication_shaper_ = std::array<double, 3>{};
    for (std::size_t index = 0; index < radiator_history_.size(); ++index) {
        radiator_history_[index] = 0.0;
    }
    history_index_ = 0;
    position_ = 0;
    finish_time_ = -1;
    block_remaining_ = 0;
    internal_in_block_ = 0;
    lung_pressure_ = 0.0;
    previous_lung_pressure_ = 0.0;
    rest_area_ = voice_.abducted_area;
    previous_rest_area_ = voice_.abducted_area;
    tension_ = 1.0;
    previous_tension_ = 1.0;
    output_gain_ = 1.0;
    previous_output_gain_ = 1.0;
    jet_area_ = 0.0;
    exit_flow_smoothed_ = 0.0;
    exit_speed_smoothed_ = 0.0;
    previous_frication_gain_ = 0.0;
    previous_jet_area_ = 0.0;
    jet_at_lips_ = false;
    free_jet_flow_ = 0.0;
    frication_gain_ = 0.0;
    tremor_phase_ = 0.0;
    transglottal_ = 0.0;
    glottal_flow_ = 0.0;
    glottal_area_ = 0.0;
    previous_lip_flow_ = 0.0;
    previous_nose_flow_ = 0.0;
    previous_wall_flow_ = 0.0;
    pitch_integral_ = 0.0;
    last_area_ = 0.0;
    trough_area_ = 0.0;
    rising_ = false;
    last_peak_ = -1;
    last_voiced_ = 0;
    internal_position_ = 0;
    measured_pitch_ = 0.0;
    probe_ = SurfaceProbe{};
    trace_.clear();
}

bool SurfaceSynthesizer::enqueue(const Cue& cue) noexcept {
    bool accepted = queue_.push(cue);
    return accepted;
}

std::size_t SurfaceSynthesizer::queue_space() const noexcept {
    return queue_.free_space();
}

void SurfaceSynthesizer::apply_cue(const Cue& cue) noexcept {
    std::size_t slot = static_cast<std::size_t>(cue.slot);
    CriticallyDampedTracker* tracker = nullptr;
    switch (cue.kind) {
    case CueKind::vowel:
        if (cue.angular_rate <= 0.0) {
            articulators_.place_vowel(cue.vowel);
        } else {
            articulators_.target_vowel(cue.vowel, cue.angular_rate);
        }
        return;
    case CueKind::constriction:
        articulators_.slot[slot].begin(cue.gesture, cue.value, cue.angular_rate);
        return;
    case CueKind::release:
        tracker = &articulators_.slot[slot].activation;
        break;
    case CueKind::velum:
        tracker = &articulators_.port_area;
        break;
    case CueKind::glottis:
        tracker = &articulators_.glottal_rest_area;
        break;
    case CueKind::lungs:
        tracker = &articulators_.lung_pressure;
        break;
    case CueKind::pitch:
        tracker = &articulators_.pitch;
        break;
    case CueKind::effort:
        tracker = &articulators_.effort;
        break;
    case CueKind::loudness:
        tracker = &articulators_.loudness;
        break;
    case CueKind::finish:
        finish_time_ = cue.time;
        return;
    }
    if (cue.angular_rate <= 0.0) {
        (*tracker).place(cue.value);
    } else {
        (*tracker).target = cue.value;
        (*tracker).angular_rate = cue.angular_rate;
    }
}

void SurfaceSynthesizer::write_station_areas() noexcept {
    std::size_t flat = 0;
    for (std::size_t section = 0; section < shape_.oral_area.size(); ++section) {
        station_area_[flat] = shape_.oral_area[section];
        flat += 1;
    }
    station_area_[flat] = shape_.port_area;
    flat += 1;
    for (std::size_t section = 0; section < shape_.nasal_area.size(); ++section) {
        station_area_[flat] = shape_.nasal_area[section];
        flat += 1;
    }
    for (std::size_t section = 0; section < shape_.side_area.size(); ++section) {
        station_area_[flat] = shape_.side_area[section];
        flat += 1;
    }
    for (std::size_t index = 0; index < floor_count_; ++index) {
        std::size_t section = floor_first_ + index;
        station_area_[flat] = std::max(shape_.sublingual_area[section], options_.floor_minimum_area);
        opening_[flat] = shape_.sublingual_opening[section] * options_.floor_opening_scale;
        flat += 1;
    }
}

void SurfaceSynthesizer::begin_block() noexcept {
    const SampleCount block_end = position_ + static_cast<SampleCount>(options_.control_block);
    Cue cue{};
    while (queue_.peek(cue) && cue.time < block_end) {
        apply_cue(cue);
        queue_.pop();
    }
    const double block_seconds = static_cast<double>(options_.control_block) * sample_period_;
    articulators_.advance(block_seconds);
    model_.build(articulators_, shape_);
    // The tongue body follows the blade's pull with its own inertia.
    Point2 pull = model_.required_drag();
    articulators_.body_offset_x.target = pull.x;
    articulators_.body_offset_y.target = pull.y;
    write_station_areas();
    network_.set_openings(std::span<const double>(opening_));
    network_.reshape(std::span<const double>(station_area_), block_internal_);

    // The jet leaves the narrowest section; a sibilant jet striking the
    // teeth is louder.
    std::size_t constriction = model_.narrowest_section();
    previous_jet_area_ = jet_area_;
    jet_area_ = std::max(model_.narrowest_area(), 0.0);
    jet_station_ = network_.tube_station(0, constriction);
    // Pressures two sections either side of the narrowest one drive the jet.
    jet_behind_ = network_.tube_station(0, constriction >= 2 ? constriction - 2 : 0);
    jet_beyond_ = constriction + 2 < options_.oral_sections ? network_.tube_station(0, constriction + 2) : no_parent;
    // While a gesture with an obstacle is under way the source sits at the
    // teeth and the obstacle's gain blends in with the gesture's progress, so
    // the noise rises and falls with the airflow. Switching place and gain at
    // one activation made fricatives start and stop like rectangles.
    std::size_t source_station = jet_station_;
    double frication_multiplier = 1.0;
    for (std::size_t slot = 0; slot < constriction_slot_count; ++slot) {
        ConstrictionGesture gesture = articulators_.slot[slot].current();
        double activation = articulators_.slot[slot].activation.value;
        if (activation > 0.05 && gesture.obstacle_from_lips >= 0.0) {
            std::size_t obstacle = model_.section_at_distance(gesture.obstacle_from_lips);
            if (obstacle > constriction) {
                double progress = std::clamp((activation - 0.05) / 0.6, 0.0, 1.0);
                double blend = progress * progress * (3.0 - 2.0 * progress);
                frication_multiplier = 1.0 + blend * (voice_.obstacle_gain * options_.obstacle_boost - 1.0);
                if (gesture.target_area < 1.0e-6) {
                    // A stop's release burst (MEASURED against the human reading:
                    // 6-10 dB under the vowel; at the sibilant's gain it peaked
                    // above the vowel and was heard as a click).
                    frication_multiplier *= options_.burst_scale;
                }
                source_station = network_.tube_station(0, obstacle);
            }
        }
    }
    network_.set_turbulence_exit(source_station, block_internal_);
    jet_at_lips_ = source_station == network_.tube_station(0, options_.oral_sections - 1);
    // Noise corner from the jet's speed and width.
    if (jet_area_ > 1.0e-10) {
        double speed = std::sqrt(2.0 * std::fabs(exit_flow_smoothed_) / air_.density);
        double width = 2.0 * std::sqrt(jet_area_ / pi);
        double corner = std::clamp(options_.frication_strouhal * speed / width, options_.lowest_frication_corner,
                                   options_.highest_frication_corner);
        frication_coefficient_ = 1.0 - std::exp(-2.0 * pi * corner / internal_rate_);
    }
    previous_frication_gain_ = frication_gain_;
    frication_gain_ = voice_.frication_gain * options_.frication_scale * frication_multiplier;

    // Larynx: breathing pressure with shimmer; tension from intended pitch,
    // the pitch controller's correction, jitter and tremor.
    previous_lung_pressure_ = lung_pressure_;
    previous_rest_area_ = rest_area_;
    previous_tension_ = tension_;
    previous_output_gain_ = output_gain_;
    double shimmer = shimmer_noise_.next();
    double effort = std::clamp(articulators_.effort.value, 0.3, 2.0);
    lung_pressure_ = articulators_.lung_pressure.value * effort * (1.0 + voice_.shimmer * shimmer);
    // Effort above 1 presses the folds together; below 1 lets them part.
    double adduction = std::clamp(1.0 - 1.5 * (effort - 1.0), 0.4, 1.6);
    rest_area_ = std::max(articulators_.glottal_rest_area.value * adduction, -2.0e-6);
    output_gain_ = std::max(articulators_.loudness.value, 0.0);
    tremor_phase_ += 2.0 * pi * voice_.tremor_rate * block_seconds;
    if (tremor_phase_ > 2.0 * pi) {
        tremor_phase_ -= 2.0 * pi;
    }
    double wander = options_.slow_wander_scale * voice_.jitter * jitter_noise_.next()
        + voice_.tremor_depth * std::sin(tremor_phase_);
    double intended = std::max(articulators_.pitch.value, 30.0);
    // Without voicing for 30 ms the correction relaxes toward none (time
    // constant about 50 ms), so the next onset starts from the feedforward.
    if (internal_position_ - last_voiced_ > static_cast<SampleCount>(0.03 * internal_rate_)) {
        pitch_integral_ *= std::exp(-block_seconds / 0.05);
    }
    double feedforward = intended / unit_pitch_;
    tension_ = std::clamp(feedforward * std::exp(pitch_integral_) * (1.0 + wander), 0.3, 6.0);

    probe_.lung_pressure = lung_pressure_;
    probe_.intended_pitch = intended;
    probe_.measured_pitch = measured_pitch_;
    probe_.narrowest_area = jet_area_;
    probe_.jet_section = constriction;
    probe_.lip_area = shape_.oral_area[shape_.oral_area.size() - 1];
    probe_.port_area = shape_.port_area;
    probe_.intraoral_pressure = network_.station_pressure(network_.tube_station(0, constriction > 0 ? constriction - 1 : 0));
    probe_.jet_flow = network_.exit_flow(jet_station_);
    probe_.glottal_flow = glottal_flow_;
    probe_.glottal_area = glottal_area_;
    block_remaining_ = options_.control_block;
    internal_in_block_ = 0;
}

void SurfaceSynthesizer::track_pitch(double glottal_area) noexcept {
    if (glottal_area > last_area_) {
        rising_ = true;
    } else if (rising_ && glottal_area < last_area_) {
        rising_ = false;
        if (last_area_ - trough_area_ > minimum_cycle_amplitude) {
            SampleCount peak = internal_position_ - 1;
            if (last_peak_ >= 0) {
                double period = static_cast<double>(peak - last_peak_) * internal_period_;
                if (period > 1.0 / 900.0 && period < 1.0 / 35.0) {
                    measured_pitch_ = 1.0 / period;
                    double intended = std::max(articulators_.pitch.value, 30.0);
                    double ratio = measured_pitch_ / intended;
                    // Only a cycle near the intended pitch is voicing; slow
                    // wobbles of an open glottis while breathing, and octave
                    // jumps, must not wind the controller up.
                    if (ratio > 0.6 && ratio < 1.6) {
                        // Integral control at a fixed rate per second, not per
                        // cycle: per cycle, a high voice corrects so often that
                        // it chases its own jitter (period alternation).
                        double error = std::log(intended / measured_pitch_);
                        pitch_integral_ = std::clamp(pitch_integral_ + options_.pitch_correction_rate * period * error, -1.2, 1.2);
                        last_voiced_ = internal_position_;
                    }
                }
            }
            last_peak_ = peak;
        }
        trough_area_ = last_area_;
    }
    trough_area_ = std::min(trough_area_, glottal_area);
    last_area_ = glottal_area;
}

void SurfaceSynthesizer::advance_internal(double fraction) noexcept {
    // Block values ramp across the block so nothing steps at the block rate.
    double lungs = previous_lung_pressure_ + fraction * (lung_pressure_ - previous_lung_pressure_);
    double rest = previous_rest_area_ + fraction * (rest_area_ - previous_rest_area_);
    double tension = previous_tension_ + fraction * (tension_ - previous_tension_);
    double unsteadiness = voice_.jitter / 0.006;
    double rough = voice_.folds.roughness * roughness_noise_.next();
    tension *= 1.0 + unsteadiness * options_.tremor * tremor_noise_.next() + rough;
    double opening = rest * std::max(0.0, 1.0 + unsteadiness * options_.adduction_wander * adduction_noise_.next()
                                              + 3.0 * voice_.folds.roughness * roughness_wander_.next());

    // Turbulence of the glottal jet of the previous sample: on the fold
    // surface and, downstream, as aspiration with a falling spectrum.
    double aspiration = 0.0;
    if (glottal_area_ > 1.0e-10) {
        double jet = std::fabs(glottal_flow_) / glottal_area_;
        double slit = 2.0 * glottal_area_ / voice_.folds.length;
        double direction = glottal_flow_ < 0.0 ? options_.inhaled_aspiration : 1.0;
        aspiration = direction * voice_.aspiration_gain * options_.aspiration_scale
            * turbulence_pressure(jet, slit, air_)
            * shape_noise(aspiration_shaper_, aspiration_coefficient_, aspiration_noise_.next());
    }
    double area = std::max(opening, 0.0);
    if (voice_.phonation) {
        double buffet = 0.002 * lungs * buffet_noise_.gaussian();
        glottis_.advance(transglottal_, opening, tension, internal_period_, options_.glottis_substeps,
                         buffet + aspiration);
        area = std::max(glottis_.minimum_area(), 0.0);
    }
    glottal_area_ = area;
    track_pitch(area);

    // Frication at the jet exit, from the exit's own flow.
    double frication = 0.0;
    // The jet's area moves through the block like the geometry does.
    const double jet_area = previous_jet_area_ + fraction * (jet_area_ - previous_jet_area_);
    if (jet_area > 1.0e-10) {
        // The jet's speed from the pressure drop across the constriction
        // (Bernoulli), smoothed: flow over area understates a jet through a
        // gap that has just opened, whose flow is still building against the
        // air's inertia (MEASURED: a /d/ release read 5 m/s with 280 Pa
        // behind it, and made no burst). Never beyond what the lungs drive.
        double behind = network_.station_pressure(jet_behind_);
        double beyond = jet_beyond_ != no_parent ? network_.station_pressure(jet_beyond_) : 0.0;
        exit_flow_smoothed_ += exit_smoothing_ * (behind - beyond - exit_flow_smoothed_);
        double limit = std::sqrt(2.0 * std::max(lungs, 0.0) / air_.density);
        double velocity = std::min(std::sqrt(2.0 * std::fabs(exit_flow_smoothed_) / air_.density), limit);
        // A closing gap: the drop across it is large but its flow is dying
        // away, so the jet is no faster than its flow over its area (with the
        // drop alone, every closure flashed like a release).
        exit_speed_smoothed_ += exit_smoothing_ * (std::fabs(network_.exit_flow(jet_station_)) / jet_area
                                                   - exit_speed_smoothed_);
        if (jet_area_ < previous_jet_area_) {
            velocity = std::min(velocity, exit_speed_smoothed_);
        }
        double diameter = 2.0 * std::sqrt(jet_area / pi);
        // Relative to a mid-sized fricative constriction of this body: the
        // same gesture in a smaller body is the same relative source.
        double area_strength = std::sqrt(jet_area / (0.1 * square_centimetre * model_.area_scale()));
        double gain = previous_frication_gain_ + fraction * (frication_gain_ - previous_frication_gain_);
        frication = gain * turbulence_pressure(velocity, diameter, air_) * area_strength
            * shape_noise(frication_shaper_, frication_coefficient_, frication_noise_.next());
    }

    network_.gather();
    double above = network_.glottis_free_pressure();
    double impedance = network_.glottis_impedance();
    double flow = orifice_flow(area, lungs - above, impedance, air_);
    transglottal_ = lungs - (above + impedance * flow);
    glottal_flow_ = flow;
    network_.scatter(flow, aspiration, frication);
    if (trace_begin_ >= 0 && internal_position_ >= trace_begin_ && trace_.size() < trace_.capacity()) {
        trace_.push_back(SurfaceTraceSample{.lip_area = network_.opening_area(0),
                                            .lip_resistance = network_.opening_resistance(0),
                                            .lip_flow = network_.radiated_flow(0), .glottal_flow = flow,
                                            .frication = frication, .aspiration = aspiration,
                                            .entry_last = network_.entry_resistance(options_.oral_sections - 1),
                                            .entry_before = network_.entry_resistance(options_.oral_sections - 2),
                                            .lip_pressure = network_.station_pressure(options_.oral_sections - 1)});
    }
    internal_position_ += 1;

    // Radiators: rate of change of the flow through the lips, the nostrils
    // and (through the neck and cheeks) the yielding walls.
    double free_jet = jet_at_lips_ ? options_.free_jet_radiation * frication * jet_area_
        / (air_.density * air_.sound_speed) : 0.0;
    // Fade the free-jet term like any other source: it must not step.
    free_jet_flow_ += 0.02 * (free_jet - free_jet_flow_);
    double lip_flow = network_.radiated_flow(0) + free_jet_flow_;
    double nose_flow = network_.radiated_flow(1);
    double wall_flow = network_.wall_flow();
    std::array<double, radiator_count> values{(lip_flow - previous_lip_flow_) * internal_rate_,
                                             (nose_flow - previous_nose_flow_) * internal_rate_,
                                             voice_.tissue_transmission * (wall_flow - previous_wall_flow_) * internal_rate_};
    for (std::size_t radiator = 0; radiator < radiator_count; ++radiator) {
        double* history = radiator_history_.data() + radiator * 2 * history_size_;
        history[history_index_] = values[radiator];
        history[history_index_ + history_size_] = values[radiator];
    }
    history_index_ = (history_index_ + 1) % history_size_;
    previous_lip_flow_ = lip_flow;
    previous_nose_flow_ = nose_flow;
    previous_wall_flow_ = wall_flow;
}

double SurfaceSynthesizer::decimated(std::size_t radiator) const noexcept {
    // The history is stored twice end to end, so the newest history_size_
    // samples are contiguous: history[index .. index + size) oldest first,
    // against the time-reversed kernel. Four partial sums let the compiler
    // keep several multiply-adds in flight.
    const double* history = radiator_history_.data() + radiator * 2 * history_size_ + history_index_;
    const double* kernel = reversed_kernel_.data();
    std::array<double, 4> partial{};
    std::size_t tap = 0;
    for (; tap + 4 <= history_size_; tap += 4) {
        partial[0] += kernel[tap] * history[tap];
        partial[1] += kernel[tap + 1] * history[tap + 1];
        partial[2] += kernel[tap + 2] * history[tap + 2];
        partial[3] += kernel[tap + 3] * history[tap + 3];
    }
    double value = partial[0] + partial[1] + partial[2] + partial[3];
    for (; tap < history_size_; ++tap) {
        value += kernel[tap] * history[tap];
    }
    return value;
}

void SurfaceSynthesizer::render(std::span<float> stereo) noexcept {
    const std::size_t frames = stereo.size() / 2;
    std::array<double, radiator_count> radiated{};
    std::array<double, microphone_count> microphone{};
    for (std::size_t frame = 0; frame < frames; ++frame) {
        if (block_remaining_ == 0) {
            begin_block();
        }
        for (std::size_t step = 0; step < options_.oversample; ++step) {
            internal_in_block_ += 1;
            double fraction = static_cast<double>(internal_in_block_) / static_cast<double>(block_internal_);
            advance_internal(fraction);
        }
        for (std::size_t radiator = 0; radiator < radiator_count; ++radiator) {
            bool heard = (options_.radiator_mask & (1u << radiator)) != 0u;
            radiated[radiator] = heard ? decimated(radiator) : 0.0;
        }
        air_path_.process(radiated, microphone);
        double ramp = static_cast<double>(internal_in_block_) / static_cast<double>(block_internal_);
        double gain = previous_output_gain_ + ramp * (output_gain_ - previous_output_gain_);
        stereo[2 * frame] = finite_sample(gain * microphone[0]);
        stereo[2 * frame + 1] = finite_sample(gain * microphone[1]);
        block_remaining_ -= 1;
        position_ += 1;
    }
}

bool SurfaceSynthesizer::finished() const noexcept {
    bool done = finish_time_ >= 0 && position_ >= finish_time_;
    return done;
}

SampleCount SurfaceSynthesizer::position() const noexcept {
    return position_;
}

const VowelLibrary& SurfaceSynthesizer::vowels() const noexcept {
    return vowels_;
}

const SurfaceProbe& SurfaceSynthesizer::probe() const noexcept {
    return probe_;
}

double SurfaceSynthesizer::unit_tension_pitch() const noexcept {
    return unit_pitch_;
}

std::size_t SurfaceSynthesizer::piece_count() const noexcept {
    return network_.piece_count();
}

const std::vector<SurfaceTraceSample>& SurfaceSynthesizer::trace() const noexcept {
    return trace_;
}

std::size_t SurfaceSynthesizer::link_count() const noexcept {
    return network_.link_count();
}

}

namespace ducklips {

std::span<const double> SurfaceSynthesizer::station_areas() const noexcept {
    return std::span<const double>(station_area_);
}

}
