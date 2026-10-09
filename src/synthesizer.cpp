#include "ducklips/synthesizer.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {

constexpr double critical_reynolds = 1800.0;
constexpr double minimum_cycle_amplitude = 2.0e-6; // m^2 of glottal area swing for a counted cycle

// Dynamic pressure of a turbulent jet that has crossed the critical Reynolds
// number; zero below it. The returned value is scaled by a source gain and
// multiplied by unit noise to form a series pressure source.
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

float finite_sample(double value, SynthesisStatistics& statistics) noexcept {
    if (!std::isfinite(value)) {
        ++statistics.nonfinite_samples;
        return 0.0f;
    }
    if (std::fabs(value) > 1.0) { ++statistics.clipped_samples; }
    double clamped = std::clamp(value, -1.0, 1.0);
    float sample = static_cast<float>(clamped);
    return sample;
}

}

void Synthesizer::prepare(const VoiceProfile& voice, const SynthesisOptions& options) {
    if (options.canonical_substeps<1 || options.canonical_substeps>64) { throw std::invalid_argument("canonical substeps outside [1,64]"); }
    voice_ = voice;
    options_ = options;
    air_ = options.air;
    sample_period_ = 1.0 / options.rate.samples_per_second;
    vowels_.build(voice, options.oral_sections, options.nasal_sections, air_);
    model_.prepare(voice, options.oral_sections, options.nasal_sections);
    model_.allocate(shape_);
    model_.allocate(start_shape_);
    model_.allocate(geometry_shape_);
    for (std::size_t index = 0; index < networks_.size(); ++index) {
        networks_[index].prepare(options.oral_sections, options.nasal_sections, model_.port_section(),
                                 model_.side_cavities());
        bases_[index].prepare(networks_[index].node_count(), options.mode_capacity);
    }
    matrix_.prepare(networks_[0].node_count(), networks_[0].tree_edges(), model_.port_section());
    solver_.prepare(networks_[0].node_count());
    resonator_.prepare(options.mode_capacity, networks_[0].node_count(), networks_[0].branch_count());
    if (options.canonical_motion) {
        model_.build(articulators_,geometry_shape_);
        networks_[0].assemble(geometry_shape_,air_,voice_.losses);
        canonical_.prepare(networks_[0],air_);
    }
    air_path_.prepare(options.air_path, air_, options.rate);
    queue_.prepare(options.cue_capacity);
    double block_rate = options.rate.samples_per_second / static_cast<double>(options.control_block);
    jitter_noise_.prepare(voice.seed * 7 + 1, 12.0, block_rate);
    shimmer_noise_.prepare(voice.seed * 7 + 2, 8.0, block_rate);
    aspiration_noise_.prepare(voice.seed * 7 + 3, 300.0, options.rate.samples_per_second);
    frication_noise_.prepare(voice.seed * 7 + 4, 300.0, options.rate.samples_per_second);
    // Two one-pole stages at 20 Hz separate the mean glottal flow from its
    // acoustic fluctuation; at 100 Hz the residual ripple is below 4 %.
    mean_coefficient_ = 1.0 - std::exp(-2.0 * pi * 20.0 / options.rate.samples_per_second);
    unit_pitch_ = calibrate_folds();
    trace_.clear();
    trace_.reserve(options.trace_samples);
    trace_begin_ = -1;
    if (options.trace_start >= 0.0) {
        trace_begin_ = samples_from_seconds(options.rate, options.trace_start);
    }
    reset();
}

double Synthesizer::calibrate_folds() const noexcept {
    // Run the folds alone at Q = 1, conversational lung pressure and voicing
    // adduction, and measure their cycle rate after the onset transient.
    TwoMassGlottis folds{};
    folds.reset(voice_.folds);
    NoiseGenerator turbulence{};
    turbulence.seed(voice_.seed * 7 + 5);
    const double turbulence_scale = 0.01 * voice_.lung_pressure;
    const SampleCount total = samples_from_seconds(options_.rate, 0.5);
    const SampleCount settle = samples_from_seconds(options_.rate, 0.2);
    double last = 0.0;
    double trough = 1.0;
    bool rising = false;
    SampleCount first_peak = -1;
    SampleCount last_peak = -1;
    std::size_t cycles = 0;
    for (SampleCount sample = 0; sample < total; ++sample) {
        double surface = turbulence_scale * turbulence.gaussian();
        folds.advance(voice_.lung_pressure, voice_.adducted_area, 1.0, sample_period_, options_.glottis_substeps, surface);
        double area = folds.minimum_area();
        if (area > last) {
            rising = true;
        } else if (rising && area < last) {
            rising = false;
            if (last - trough > minimum_cycle_amplitude && sample > settle) {
                if (first_peak < 0) {
                    first_peak = sample;
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
    double seconds = seconds_from_samples(options_.rate, last_peak - first_peak);
    double pitch = static_cast<double>(cycles) / seconds;
    return pitch;
}

void Synthesizer::reset() noexcept {
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
    intraoral_.reset();
    resonator_.clear();
    canonical_ready_=false;
    if (options_.canonical_motion) { canonical_.reset(); }
    bases_[0].mode_count = 0;
    bases_[1].mode_count = 0;
    current_ = 0;
    air_path_.reset();
    position_ = 0;
    finish_time_ = -1;
    block_remaining_ = 0;
    geometry_remaining_ = 0;
    geometry_steps_ = 1;
    geometry_step_index_ = 0;
    geometry_ready_ = false;
    geometry_lip_metric_jump_ = 0.0;
    geometry_lip_projection_jump_ = 0.0;
    geometry_pressure_jump_ = 0.0;
    geometry_lip_slope_jump_ = 0.0;
    geometry_energy_ratio_ = 0.0;
    geometry_parameter_change_ = 0.0;
    transglottal_ = 0.0;
    glottal_pressure_ = 0.0;
    glottal_flow_ = 0.0;
    glottal_area_ = 0.0;
    mean_stage_ = 0.0;
    mean_flow_ = 0.0;
    oral_outlet_ = 0.0;
    nasal_outlet_ = 0.0;
    previous_oral_outlet_ = 0.0;
    previous_nasal_outlet_ = 0.0;
    previous_oral_outflow_ = 0.0;
    previous_nasal_outflow_ = 0.0;
    oral_stage_ = 0.0;
    nasal_stage_ = 0.0;
    pitch_integral_ = 0.0;
    last_area_ = 0.0;
    trough_area_ = 0.0;
    rising_ = false;
    last_peak_ = -1;
    measured_pitch_ = 0.0;
    tremor_phase_ = 0.0;
    output_gain_ = 1.0;
    lung_pressure_ = 0.0;
    shimmer_multiplier_ = 1.0;
    pitch_multiplier_ = 1.0;
    pitch_correction_ = 1.0;
    tap_previous_flow_ = 0.0;
    flow_sum_ = 0.0;
    area_sum_ = 0.0;
    probe_samples_ = 0;
    statistics_ = SynthesisStatistics{};
    probe_ = SynthesisProbe{};
    trace_.clear();
}

bool Synthesizer::enqueue(const Cue& cue) noexcept {
    bool accepted = queue_.push(cue);
    return accepted;
}

std::size_t Synthesizer::queue_space() const noexcept {
    return queue_.free_space();
}

void Synthesizer::apply_cue(const Cue& cue) noexcept {
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

void Synthesizer::begin_block() noexcept {
    const SampleCount block_end = position_ + static_cast<SampleCount>(options_.control_block);
    Cue cue{};
    while (queue_.peek(cue) && cue.time < block_end) {
        apply_cue(cue);
        queue_.pop();
    }
    const double block_seconds = static_cast<double>(options_.control_block) * sample_period_;
    // Geometry solves remain at control rate; physiological forces follow
    // their continuous trajectories in render_samples, not a block hold.
    articulators_.advance_geometry(block_seconds);
    model_.build(articulators_, shape_);
    // The tongue body follows the blade's pull with its own inertia.
    Point2 pull = model_.required_drag();
    articulators_.body_offset_x.target = pull.x;
    articulators_.body_offset_y.target = pull.y;

    // Adaptive geometry stepping. A closure changes the mode shapes so fast
    // that one projection per block scatters energy out of the modes it
    // belongs to; the modes of a smoothly closing tube follow adiabatically.
    // The step count grows with the largest log-area change in the block.
    std::size_t steps = 1;
    if (geometry_ready_) {
        const double floor = 0.01 * square_centimetre * model_.area_scale();
        double largest = 0.0;
        for (std::size_t section = 0; section < shape_.oral_area.size(); ++section) {
            double ratio = (shape_.oral_area[section] + floor) / (geometry_shape_.oral_area[section] + floor);
            largest = std::max(largest, std::fabs(std::log(ratio)));
        }
        double port_ratio = (shape_.port_area + floor) / (geometry_shape_.port_area + floor);
        largest = std::max(largest, std::fabs(std::log(port_ratio)));
        double wanted = std::ceil(largest / options_.geometry_step_limit);
        steps = static_cast<std::size_t>(std::clamp(wanted, 1.0, static_cast<double>(options_.maximum_geometry_steps)));
    }
    start_shape_.port_area = geometry_shape_.port_area;
    for (std::size_t section = 0; section < shape_.oral_area.size(); ++section) {
        start_shape_.oral_area[section] = geometry_shape_.oral_area[section];
        start_shape_.oral_length[section] = geometry_shape_.oral_length[section];
    }
    if (!geometry_ready_) {
        // First block: start from the target itself.
        start_shape_.port_area = shape_.port_area;
        for (std::size_t section = 0; section < shape_.oral_area.size(); ++section) {
            start_shape_.oral_area[section] = shape_.oral_area[section];
            start_shape_.oral_length[section] = shape_.oral_length[section];
        }
        geometry_ready_ = true;
    }
    if (options_.canonical_motion) { steps=1; }
    geometry_steps_ = steps;
    geometry_step_index_ = 0;
    statistics_.geometry_updates += steps;
    statistics_.most_geometry_steps = std::max(statistics_.most_geometry_steps, steps);

    // Source placement: the jet leaves the narrowest section; a sibilant jet
    // strikes the teeth further downstream.
    const AcousticNetwork& layout = networks_[current_];
    const std::size_t oral = model_.oral_sections();
    std::size_t constriction = model_.narrowest_section();
    frication_branch_ = layout.lip_branch();
    if (constriction + 1 < oral) {
        frication_branch_ = layout.oral_branch(constriction);
    }
    double frication_multiplier = 1.0;
    for (std::size_t slot = 0; slot < constriction_slot_count; ++slot) {
        ConstrictionGesture gesture = articulators_.slot[slot].current();
        if (articulators_.slot[slot].activation.value > 0.3 && gesture.obstacle_from_lips >= 0.0) {
            std::size_t obstacle = model_.section_at_distance(gesture.obstacle_from_lips);
            if (obstacle > constriction) {
                frication_branch_ = layout.lip_branch();
                if (obstacle + 1 < oral) {
                    frication_branch_ = layout.oral_branch(obstacle);
                }
                frication_multiplier = voice_.obstacle_gain;
            }
        }
    }
    aspiration_branch_ = layout.oral_branch(0);

    // Aerodynamics: the oral constriction and the nasal path vent one cavity.
    double nasal_outlet = 0.0;
    double nostrils = shape_.nasal_area[shape_.nasal_area.size() - 1]; // nostril entrance
    if (shape_.port_area > 0.0 && nostrils > 0.0) {
        double series = 1.0 / (shape_.port_area * shape_.port_area) + 1.0 / (nostrils * nostrils);
        nasal_outlet = 1.0 / std::sqrt(series);
    }
    previous_oral_outlet_ = oral_outlet_;
    previous_nasal_outlet_ = nasal_outlet_;
    constriction_area_ = std::max(model_.narrowest_area(), 0.0);
    oral_outlet_ = constriction_area_;
    nasal_outlet_ = nasal_outlet;
    const double bulk_stiffness = air_.density * air_.sound_speed * air_.sound_speed;
    previous_cavity_compliance_=cavity_compliance_;
    previous_cavity_volume_=cavity_volume_;
    cavity_volume_=model_.oral_volume();
    cavity_compliance_=cavity_volume_/bulk_stiffness+voice_.static_wall_compliance;
    if (options_.canonical_motion && !canonical_ready_) {
        previous_cavity_compliance_=cavity_compliance_;
        previous_cavity_volume_=cavity_volume_;
    }
    frication_scale_ = voice_.frication_gain * frication_multiplier;

    // Larynx: breathing pressure with shimmer; tension from intended pitch,
    // the pitch controller's correction, jitter and tremor.
    double shimmer = shimmer_noise_.next();
    shimmer_multiplier_ = 1.0 + voice_.shimmer * shimmer;
    double effort = std::clamp(articulators_.effort.value, 0.3, 2.0);
    lung_pressure_ = articulators_.lung_pressure.value * effort * shimmer_multiplier_;
    // Effort above 1 presses the folds together; below 1 lets them part.
    double adduction = std::clamp(1.0 - 1.5 * (effort - 1.0), 0.4, 1.6);
    rest_area_ = std::max(articulators_.glottal_rest_area.value * adduction, -2.0e-6);
    output_gain_ = std::max(articulators_.loudness.value, 0.0);
    tremor_phase_ += 2.0 * pi * voice_.tremor_rate * block_seconds;
    if (tremor_phase_ > 2.0 * pi) {
        tremor_phase_ -= 2.0 * pi;
    }
    double wander = voice_.jitter * jitter_noise_.next() + voice_.tremor_depth * std::sin(tremor_phase_);
    double intended = std::max(articulators_.pitch.value, 30.0);
    double feedforward = intended / unit_pitch_;
    pitch_multiplier_ = 1.0 + wander;
    tension_ = std::clamp(feedforward * pitch_correction_ * pitch_multiplier_, 0.3, 6.0);

    // Probe of the block that just finished.
    if (probe_samples_ > 0) {
        probe_.glottal_flow = flow_sum_ / static_cast<double>(probe_samples_);
        probe_.glottal_area = area_sum_ / static_cast<double>(probe_samples_);
    }
    flow_sum_ = 0.0;
    area_sum_ = 0.0;
    probe_samples_ = 0;
    probe_.intraoral_pressure = intraoral_.pressure();
    probe_.lung_pressure = lung_pressure_;
    probe_.intended_pitch = intended;
    probe_.measured_pitch = measured_pitch_;
    probe_.tension = tension_;
    probe_.narrowest_area = constriction_area_;
    probe_.acoustic_energy = options_.canonical_motion ? canonical_.energy() : resonator_.energy();
    probe_.transglottal = transglottal_;
    probe_.glottal_pressure = glottal_pressure_;
    probe_.modes = bases_[current_].mode_count;
    for (std::size_t index = 0; index < probe_.lowest_mode_hz.size(); ++index) {
        double frequency = 0.0;
        if (index < bases_[current_].mode_count) {
            frequency = bases_[current_].angular_frequency[index] / (2.0 * pi);
        }
        probe_.lowest_mode_hz[index] = frequency;
    }
    block_remaining_ = options_.control_block;
}

void Synthesizer::begin_geometry_step() noexcept {
    geometry_step_index_ += 1;
    const double fraction = static_cast<double>(geometry_step_index_) / static_cast<double>(geometry_steps_);
    geometry_shape_.port_area = start_shape_.port_area + fraction * (shape_.port_area - start_shape_.port_area);
    for (std::size_t section = 0; section < shape_.oral_area.size(); ++section) {
        double start_area = start_shape_.oral_area[section];
        double start_length = start_shape_.oral_length[section];
        geometry_shape_.oral_area[section] = start_area + fraction * (shape_.oral_area[section] - start_area);
        geometry_shape_.oral_length[section] = start_length + fraction * (shape_.oral_length[section] - start_length);
    }
    for (std::size_t section = 0; section < shape_.nasal_area.size(); ++section) {
        geometry_shape_.nasal_area[section] = shape_.nasal_area[section];
        geometry_shape_.nasal_length[section] = shape_.nasal_length[section];
    }

    const std::size_t next = 1 - current_;
    AcousticNetwork& network = networks_[next];
    network.assemble(geometry_shape_, air_, voice_.losses);
    if (options_.canonical_motion) {
        if (!canonical_ready_) {
            networks_[current_].assemble(geometry_shape_,air_,voice_.losses);
            canonical_ready_=true;
        }
        current_=next;geometry_remaining_=options_.control_block;
        ++statistics_.blocks;statistics_.fewest_modes=0;
        return;
    }
    network.write_operator(matrix_);
    const double maximum = 2.0 * pi * options_.maximum_mode_frequency;
    ModalSolveReport report = solver_.solve(network, matrix_, maximum, bases_[current_], bases_[next]);
    statistics_.blocks += 1;
    statistics_.warm_modes += report.warm_modes;
    statistics_.bisected_modes += report.bisected_modes;
    statistics_.factorizations += report.factorizations;
    statistics_.fewest_modes = std::min(statistics_.fewest_modes, bases_[next].mode_count);
    statistics_.most_modes = std::max(statistics_.most_modes, bases_[next].mode_count);
    if (options_.diagnostics) {
        const ModalBasis& basis = bases_[next];
        double defect = 0.0;
        for (std::size_t first = 0; first < basis.mode_count; ++first) {
            for (std::size_t second = first; second < basis.mode_count; ++second) {
                double sum = 0.0;
                for (std::size_t node = 0; node < network.node_count(); ++node) {
                    sum += network.nodes[node].compliance * basis.pressure(first, node) * basis.pressure(second, node);
                }
                double expected = 0.0;
                if (first == second) {
                    expected = 1.0;
                }
                defect = std::max(defect, std::fabs(sum - expected));
            }
        }
        if (defect > statistics_.worst_orthogonality) {
            statistics_.worst_orthogonality = defect;
            statistics_.worst_orthogonality_time = seconds_from_samples(options_.rate, position_);
        }
    }
    const bool trace_geometry = options_.trace_samples > 0 && bases_[current_].mode_count > 0;
    double previous_pressure = 0.0, previous_lips = 0.0, expected_lips = 0.0;
    double previous_slope = 0.0, previous_energy = 0.0;
    if (trace_geometry) {
        const AcousticNetwork& old_network = networks_[current_];
        const std::size_t lip_node = old_network.oral_sections()-1;
        previous_energy = resonator_.energy();
        geometry_parameter_change_ = 0.0;
        for (std::size_t node=0;node<network.node_count();++node) {
            geometry_parameter_change_ = std::max(geometry_parameter_change_,
                std::abs(network.nodes[node].compliance/old_network.nodes[node].compliance-1));
        }
        for (std::size_t branch=0;branch<network.branch_count();++branch) {
            double a=old_network.branches[branch].conductance,b=network.branches[branch].conductance;
            if (a+b>0) { geometry_parameter_change_ = std::max(geometry_parameter_change_,std::abs(a-b)/(a+b)); }
        }
        double lip_pressure = 0.0;
        for (std::size_t mode=0;mode<bases_[current_].mode_count;++mode) {
            const double amplitude = resonator_.state[mode].real();
            previous_pressure += bases_[current_].pressure(mode,0)*amplitude;
            lip_pressure += bases_[current_].pressure(mode,lip_node)*amplitude;
            previous_slope += resonator_.lip_output[mode]*(-bases_[current_].damping[mode]*amplitude
                -bases_[current_].angular_frequency[mode]*resonator_.state[mode].imag());
        }
        previous_lips = lip_pressure*old_network.branches[old_network.lip_branch()].conductance;
        double ratio = 1.0;
        if (options_.invariant == ReprojectionInvariant::passive) {
            ratio = std::min(old_network.nodes[lip_node].compliance,network.nodes[lip_node].compliance)
                    /network.nodes[lip_node].compliance;
        } else if (options_.invariant == ReprojectionInvariant::charge_momentum) {
            ratio = old_network.nodes[lip_node].compliance/network.nodes[lip_node].compliance;
        }
        expected_lips = ratio*lip_pressure*network.branches[network.lip_branch()].conductance;
    }
    if (bases_[current_].mode_count > 0) {
        resonator_.reproject(networks_[current_], bases_[current_], network, bases_[next], options_.invariant);
    }
    resonator_.configure(network, bases_[next], sample_period_, aspiration_branch_, frication_branch_);
    if (trace_geometry) {
        double pressure=0.0,lips=0.0,slope=0.0;
        for (std::size_t mode=0;mode<bases_[next].mode_count;++mode) {
            pressure += resonator_.glottal_weight[mode]*resonator_.state[mode].real();
            lips += resonator_.lip_output[mode]*resonator_.state[mode].real();
            slope += resonator_.lip_output[mode]*(-bases_[next].damping[mode]*resonator_.state[mode].real()
                -bases_[next].angular_frequency[mode]*resonator_.state[mode].imag());
        }
        geometry_pressure_jump_ = pressure-previous_pressure;
        geometry_lip_slope_jump_ = slope-previous_slope;
        geometry_energy_ratio_ = previous_energy>0 ? resonator_.energy()/previous_energy : 1.0;
        geometry_lip_metric_jump_ = expected_lips-previous_lips;
        geometry_lip_projection_jump_ = lips-expected_lips;
    }
    current_ = next;

    // Split the block's samples as evenly as possible among its steps.
    std::size_t base = options_.control_block / geometry_steps_;
    std::size_t extra = options_.control_block % geometry_steps_;
    geometry_remaining_ = base;
    if (geometry_step_index_ <= extra) {
        geometry_remaining_ += 1;
    }
}

void Synthesizer::track_pitch(double glottal_area) noexcept {
    if (glottal_area > last_area_) {
        rising_ = true;
    } else if (rising_ && glottal_area < last_area_) {
        rising_ = false;
        if (last_area_ - trough_area_ > minimum_cycle_amplitude) {
            SampleCount peak = position_ - 1;
            if (last_peak_ >= 0) {
                double period = seconds_from_samples(options_.rate, peak - last_peak_);
                if (period > 1.0 / 900.0 && period < 1.0 / 35.0) {
                    measured_pitch_ = 1.0 / period;
                    double intended = std::max(articulators_.pitch.value, 30.0);
                    double error = std::log(intended / measured_pitch_);
                    pitch_integral_ = std::clamp(pitch_integral_ + 0.25 * error, -1.2, 1.2);
                    pitch_correction_ = std::exp(pitch_integral_);
                    statistics_.pitch_cycles += 1;
                }
            }
            last_peak_ = peak;
        }
        trough_area_ = last_area_;
    }
    trough_area_ = std::min(trough_area_, glottal_area);
    last_area_ = glottal_area;
}

void Synthesizer::render_samples(std::span<float> stereo) noexcept {
    const std::size_t frames = stereo.size() / 2;
    const std::size_t modes = resonator_.mode_count;
    const double period = sample_period_;
    double half_impedance = 0.5 * resonator_.glottal_impedance;
    const double constriction_diameter = 2.0 * std::sqrt(constriction_area_ / pi);
    // Turbulent source strength grows with the jet's cross-section (Stevens):
    // sqrt(A / 0.1 cm^2) relative to a mid-sized fricative constriction.
    const double area_strength = std::sqrt(constriction_area_ / (0.1 * square_centimetre));
    const double block_length = static_cast<double>(options_.control_block);
    const double root_two_over_density = std::sqrt(2.0 / air_.density);
    std::array<double, radiator_count> radiated{};
    std::array<double, microphone_count> microphone{};
    for (std::size_t frame = 0; frame < frames; ++frame) {
        articulators_.advance_physiology(period);
        const double effort = std::clamp(articulators_.effort.value, 0.3, 2.0);
        lung_pressure_ = articulators_.lung_pressure.value * effort * shimmer_multiplier_;
        const double adduction = std::clamp(1.0 - 1.5 * (effort - 1.0), 0.4, 1.6);
        rest_area_ = std::max(articulators_.glottal_rest_area.value * adduction, -2.0e-6);
        tension_ = std::clamp(std::max(articulators_.pitch.value, 30.0)/unit_pitch_
                             * pitch_correction_ * pitch_multiplier_, 0.3, 6.0);
        output_gain_ = std::max(articulators_.loudness.value, 0.0);
        // Turbulence of the glottal jet of the previous sample. The same
        // fluctuating pressure acts on the fold surface and, downstream, as
        // the aspiration source in the tract.
        double aspiration = 0.0;
        if (glottal_area_ > 1.0e-10) {
            double jet = std::fabs(glottal_flow_) / glottal_area_;
            double slit = 2.0 * glottal_area_ / voice_.folds.length;
            aspiration = voice_.aspiration_gain * turbulence_pressure(jet, slit, air_) * aspiration_noise_.next();
        }

        // Folds move under the pressure drop of the previous sample.
        double area = std::max(rest_area_, 0.0);
        if (voice_.phonation) {
            glottis_.advance(transglottal_, rest_area_, tension_, period, options_.glottis_substeps, aspiration);
            area = std::max(glottis_.minimum_area(), 0.0);
        }
        glottal_area_ = area;
        track_pitch(area);
        double cavity = intraoral_.pressure();
        double constriction_jet = std::sqrt(2.0 * std::fabs(cavity) / air_.density);
        double frication = frication_scale_ * turbulence_pressure(constriction_jet, constriction_diameter, air_)
            * area_strength * frication_noise_.next();

        // Free acoustic evolution with series pressure sources.
        double free_pressure = 0.0;
        if (options_.canonical_motion) {
            const double fraction=(block_length-double(block_remaining_-frame))/block_length;
            canonical_.free_step(networks_[1-current_],networks_[current_],fraction,1/block_length,period,
                                 options_.canonical_substeps,aspiration_branch_,frication_branch_,aspiration,frication);
        }
        for (std::size_t mode = 0; mode < modes; ++mode) {
            double series = resonator_.aspiration_weight[mode] * aspiration + resonator_.frication_weight[mode] * frication;
            std::complex<double> drive(0.0, series);
            std::complex<double> value = resonator_.propagator[mode] * resonator_.state[mode]
                + resonator_.input_gain[mode] * drive;
            resonator_.state[mode] = value;
            free_pressure += resonator_.glottal_weight[mode] * value.real();
        }

        // Glottal flow against the midpoint supraglottal pressure, solved
        // implicitly with the tract's instantaneous response.
        // Only the fluctuation of the glottal flow about its slow mean drives
        // the acoustic modes; the mean flow belongs to the aerodynamic cavity.
        double free_drop = lung_pressure_ - cavity - 0.5 * (glottal_pressure_ + free_pressure)
            + half_impedance * mean_flow_;
        if (options_.canonical_motion) {
            // The full network includes modes above output Nyquist. Use the
            // work-conjugate pressure integral, not an endpoint trapezoid.
            half_impedance=canonical_.load_impedance();
            free_drop=lung_pressure_-cavity-canonical_.load_pressure()+half_impedance*mean_flow_;
        }
        double flow = orifice_flow(area, free_drop, half_impedance, air_);
        double fluctuation = flow - mean_flow_;

        double next_pressure = 0.0;
        double lips = 0.0;
        double nostrils = 0.0;
        double throat = 0.0;
        for (std::size_t mode = 0; mode < modes; ++mode) {
            std::complex<double> injected = resonator_.input_gain[mode] * (resonator_.glottal_weight[mode] * fluctuation);
            std::complex<double> value = resonator_.state[mode] + injected;
            resonator_.state[mode] = value;
            double amplitude = value.real();
            next_pressure += resonator_.glottal_weight[mode] * amplitude;
            lips += resonator_.lip_output[mode] * amplitude;
            nostrils += resonator_.nostril_output[mode] * amplitude;
            throat += resonator_.wall_output[mode] * amplitude;
        }
        if (options_.canonical_motion) {
            canonical_.inject(fluctuation);next_pressure=canonical_.pressure();
            const AcousticNetwork& layout=networks_[current_];
            lips=canonical_.flow_derivative(layout.lip_branch(),period);
            nostrils=canonical_.flow_derivative(layout.nostril_branch(),period);
            for (std::size_t node=0;node<layout.oral_sections();++node) { throat+=canonical_.flow_derivative(layout.wall_branch(node),period); }
        }
        transglottal_ = lung_pressure_ - cavity - 0.5 * (glottal_pressure_ + next_pressure);
        if (options_.canonical_motion) {
            transglottal_=lung_pressure_-cavity-canonical_.load_pressure()-canonical_.load_impedance()*fluctuation;
        }
        glottal_pressure_ = next_pressure;
        glottal_flow_ = flow;
        mean_stage_ += mean_coefficient_ * (flow - mean_stage_);
        mean_flow_ += mean_coefficient_ * (mean_stage_ - mean_flow_);

        // Aerodynamics on the slow flow. Outlet areas ramp across the block so
        // the radiated mean-flow derivative has no block-rate steps.
        double ramp = (static_cast<double>(block_length) - static_cast<double>(block_remaining_ - frame) + 1.0) / block_length;
        double oral_area = previous_oral_outlet_ + ramp * (oral_outlet_ - previous_oral_outlet_);
        double nasal_area = previous_nasal_outlet_ + ramp * (nasal_outlet_ - previous_nasal_outlet_);
        // The cavity conserves mass, so it receives the true glottal flow; a
        // lagging mean would overfill it past lung pressure during closures.
        double vented,displacement_flow=0,mass_residual=0;
        if (options_.canonical_motion) {
            const double change=cavity_compliance_-previous_cavity_compliance_;
            const double first=previous_cavity_compliance_+(ramp-1/block_length)*change;
            const double last=previous_cavity_compliance_+ramp*change;
            const double displaced=(cavity_volume_-previous_cavity_volume_)/block_length;
            vented=intraoral_.advance_moving(flow,oral_area+nasal_area,first,last,displaced,period,air_);
            displacement_flow=displaced/period;
            const double outlet=std::copysign((oral_area+nasal_area)*std::sqrt(2*std::abs(vented)/air_.density),vented);
            mass_residual=(last*vented-first*cavity+displaced)/period-flow+outlet;
            statistics_.maximum_cavity_mass_residual=std::max(statistics_.maximum_cavity_mass_residual,std::abs(mass_residual));
        } else {
            vented=intraoral_.advance(flow,oral_area+nasal_area,cavity_compliance_,period,air_);
        }
        double jet = root_two_over_density * std::sqrt(std::fabs(vented));
        if (vented < 0.0) {
            jet = -jet;
        }
        // Mean outflow radiates through its rate of change: the release of a
        // stop vents stored air and is heard as a low-frequency pop. The
        // outflow follows the glottal pulses when the tract is open; those
        // already reach the air through the modes, so only its slow part
        // (the same two-pole 20 Hz smoothing) is radiated here.
        oral_stage_ += mean_coefficient_ * (oral_area * jet - oral_stage_);
        double oral_outflow = previous_oral_outflow_ + mean_coefficient_ * (oral_stage_ - previous_oral_outflow_);
        nasal_stage_ += mean_coefficient_ * (nasal_area * jet - nasal_stage_);
        double nasal_outflow = previous_nasal_outflow_ + mean_coefficient_ * (nasal_stage_ - previous_nasal_outflow_);
        if (options_.radiate_mean_flow) {
            lips += (oral_outflow - previous_oral_outflow_) / period;
            nostrils += (nasal_outflow - previous_nasal_outflow_) / period;
        }
        previous_oral_outflow_ = oral_outflow;
        previous_nasal_outflow_ = nasal_outflow;

        if (resonator_.frication_at_lips) {
            lips += resonator_.lip_conductance * frication;
        }
        radiated[static_cast<std::size_t>(Radiator::lips)] = lips;
        radiated[static_cast<std::size_t>(Radiator::nostrils)] = nostrils;
        radiated[static_cast<std::size_t>(Radiator::throat)] = voice_.tissue_transmission * throat;
        air_path_.process(radiated, microphone);
        if (trace_begin_ >= 0 && position_ >= trace_begin_ && trace_.size() < trace_.capacity()) {
            trace_.push_back(TraceSample{.lips = lips, .nostrils = nostrils, .throat = throat,
                                         .glottal_pressure = next_pressure, .energy = options_.canonical_motion ? canonical_.energy() : resonator_.energy(),
                                         .modes = modes, .glottal_area = area, .glottal_flow = flow,
                                         .tension = tension_, .geometry_lip_metric_jump = geometry_lip_metric_jump_,
                                         .geometry_lip_projection_jump = geometry_lip_projection_jump_,
                                         .geometry_pressure_jump = geometry_pressure_jump_,
                                         .geometry_lip_slope_jump = geometry_lip_slope_jump_,
                                         .geometry_energy_ratio = geometry_energy_ratio_,
                                         .geometry_parameter_change = geometry_parameter_change_,
                                         .frication_pressure = frication, .frication_branch = frication_branch_,
                                         .constriction_area = constriction_area_, .cavity_pressure = cavity,
                                         .cavity_displacement_flow = displacement_flow, .cavity_mass_residual = mass_residual});
        }
        geometry_lip_metric_jump_ = 0.0;
        geometry_lip_projection_jump_ = 0.0;
        geometry_pressure_jump_ = 0.0;
        geometry_lip_slope_jump_ = 0.0;
        geometry_energy_ratio_ = 0.0;
        geometry_parameter_change_ = 0.0;
        double gain = output_gain_;
        double left = microphone[0];
        double right = microphone[1];
        if (options_.tap == OutputTap::lips) {
            // Lip radiation as pressure at 0.4 m, without delays or the second ear.
            double direct = air_.density / (4.0 * pi * 0.4) * lips / options_.air_path.full_scale_pressure;
            left = direct;
            right = direct;
        } else if (options_.tap == OutputTap::glottal_source) {
            double derivative = (flow - tap_previous_flow_) / period;
            double direct = air_.density / (4.0 * pi * 0.4) * derivative / options_.air_path.full_scale_pressure;
            left = direct;
            right = direct;
        }
        tap_previous_flow_ = flow;
        stereo[2 * frame] = finite_sample(gain * left, statistics_);
        stereo[2 * frame + 1] = finite_sample(gain * right, statistics_);

        flow_sum_ += flow;
        area_sum_ += area;
        probe_samples_ += 1;
        position_ += 1;
    }
}

void Synthesizer::render(std::span<float> stereo) noexcept {
    const std::size_t frames = stereo.size() / 2;
    std::size_t frame = 0;
    while (frame < frames) {
        if (block_remaining_ == 0) {
            begin_block();
        }
        if (geometry_remaining_ == 0) {
            begin_geometry_step();
        }
        std::size_t run = std::min(geometry_remaining_, frames - frame);
        render_samples(stereo.subspan(2 * frame, 2 * run));
        frame += run;
        block_remaining_ -= run;
        geometry_remaining_ -= run;
    }
}

bool Synthesizer::finished() const noexcept {
    bool done = finish_time_ >= 0 && position_ >= finish_time_;
    return done;
}

SampleCount Synthesizer::position() const noexcept {
    return position_;
}

const SynthesisStatistics& Synthesizer::statistics() const noexcept {
    return statistics_;
}

const SynthesisProbe& Synthesizer::probe() const noexcept {
    return probe_;
}

const VowelLibrary& Synthesizer::vowels() const noexcept {
    return vowels_;
}

const VoiceProfile& Synthesizer::voice() const noexcept {
    return voice_;
}

const SynthesisOptions& Synthesizer::options() const noexcept {
    return options_;
}

const std::vector<TraceSample>& Synthesizer::trace() const noexcept {
    return trace_;
}

double Synthesizer::unit_tension_pitch() const noexcept {
    return unit_pitch_;
}

}
