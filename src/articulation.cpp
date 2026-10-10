#include "ducklips/articulation.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {

constexpr double reference_tract_length = 0.175; // m
constexpr double epilarynx_fraction = 0.12;
constexpr double lip_fraction = 0.06;
constexpr double minimum_vowel_area = 0.15; // cm^2, reference scale
constexpr double constriction_core = 0.0015; // m, flat-topped half width at reference scale

double lip_storage_area(const LipProfile& profile,double first,double last) noexcept {
    // Fixed eight-point volume quadrature on each smooth piece. Contact at
    // the exterior does not empty the partially open section behind it.
    constexpr std::array<double,4> x{.1834346424956498,.525532409916329,.7966664774136267,.9602898564975363};
    constexpr std::array<double,4> w{.362683783378362,.3137066458778873,.2223810344533745,.1012285362903763};
    std::array<double,4> cuts{first,last,std::clamp(profile.place-profile.core,first,last),
                            std::clamp(profile.place+profile.core,first,last)};
    std::sort(cuts.begin(),cuts.end());
    double volume=0;
    for(std::size_t piece=1;piece<cuts.size();++piece){
        const double mid=.5*(cuts[piece]+cuts[piece-1]),half=.5*(cuts[piece]-cuts[piece-1]);
        for(std::size_t q=0;q<x.size();++q){volume+=half*w[q]*(profile.area(mid-half*x[q])+profile.area(mid+half*x[q]));}
    }
    return volume/(last-first);
}

struct Knot final {
    double position{}; // fraction of length
    double value{};
};

// Neutral reference area function, cm^2 against fraction of length from the
// glottis. A modelling choice shaped after published adult male area
// functions; not a measurement of a particular speaker.
constexpr std::array<Knot, 7> neutral_knots{{
    {0.00, 0.6}, {0.10, 0.8}, {0.16, 2.6}, {0.45, 3.0}, {0.60, 3.4}, {0.88, 3.2}, {1.00, 3.0},
}};

// Nasal reference area function, cm^2 from the velar port to the nostril
// entrance; the last value is replaced by the voice's nostril area.
constexpr std::array<Knot, 6> nasal_knots{{
    {0.00, 1.2}, {0.15, 2.0}, {0.40, 4.0}, {0.70, 4.0}, {0.90, 2.0}, {1.00, 1.4},
}};

template <std::size_t count>
double interpolate_log(const std::array<Knot, count>& knots, double position) noexcept {
    if (position <= knots[0].position) {
        return std::log(knots[0].value);
    }
    for (std::size_t index = 1; index < count; ++index) {
        if (position <= knots[index].position) {
            const Knot& left = knots[index - 1];
            const Knot& right = knots[index];
            double t = (position - left.position) / (right.position - left.position);
            double value = (1.0 - t) * std::log(left.value) + t * std::log(right.value);
            return value;
        }
    }
    double last = std::log(knots[count - 1].value);
    return last;
}

double smoothstep(double t) noexcept {
    double clamped = std::clamp(t, 0.0, 1.0);
    double value = clamped * clamped * (3.0 - 2.0 * clamped);
    return value;
}

}

void CriticallyDampedTracker::advance(double step) noexcept {
    // Three equal poles settle later than two at the same rate; the 1.4
    // factor keeps the 90 % settling time close to the second-order tracker
    // the cue rates were chosen for.
    double rate = angular_rate;
    if (speed_limit > 0.0) {
        rate = std::min(rate, speed_limit);
    }
    rate *= 1.4;
    // Exact motion of the three repeated poles under a held target. The old
    // Euler steps gave geometry endpoints a velocity inconsistent with the
    // tracker's stored velocity. A radiating moving wall exposes that error.
    const double displacement=value-target;
    const double command_error=command-target;
    const double slope=velocity+rate*displacement;
    const double curvature=rate*rate*command_error;
    const double polynomial=displacement+step*(slope+.5*step*curvature);
    const double decay=std::exp(-rate*step);
    value=target+decay*polynomial;
    velocity=decay*(slope+step*curvature-rate*polynomial);
    command=target+decay*command_error;
}

void CriticallyDampedTracker::place(double position) noexcept {
    value = position;
    target = position;
    command = position;
    velocity = 0.0;
}

void ConstrictionSlotState::begin(const ConstrictionGesture& gesture, double activation_target,
                                  double angular_rate) noexcept {
    // Rate at which an active articulator glides to a new place, and the much
    // faster rate at which it changes degree in place: going from a fricative
    // gap to a seal is a millimetre of motion (/z/ to /d/ in "is DuckLips").
    constexpr double glide_rate = 160.0;  // rad/s
    constexpr double degree_rate = 450.0; // rad/s
    constexpr double idle_activation = 0.05;
    if (activation.value < idle_activation) {
        place.place(gesture.place_from_lips);
        target_area.place(gesture.target_area);
        width.place(gesture.width);
        tip_down.place(gesture.tip_down);
    } else {
        place.target = gesture.place_from_lips;
        place.angular_rate = glide_rate;
        target_area.target = gesture.target_area;
        target_area.angular_rate = degree_rate;
        width.target = gesture.width;
        width.angular_rate = glide_rate;
        tip_down.target = gesture.tip_down;
        tip_down.angular_rate = glide_rate;
    }
    obstacle_from_lips = gesture.obstacle_from_lips;
    if (angular_rate<=0) {
        activation.place(activation_target);
    } else {
        activation.target = activation_target;
        activation.angular_rate = angular_rate;
    }
}

ConstrictionGesture ConstrictionSlotState::current() const noexcept {
    ConstrictionGesture gesture{};
    gesture.place_from_lips = place.value;
    gesture.width = std::max(width.value, 1.0e-4);
    gesture.target_area = std::max(target_area.value, 0.0);
    gesture.obstacle_from_lips = obstacle_from_lips;
    gesture.tip_down = std::clamp(tip_down.value, 0.0, 1.0);
    return gesture;
}

void ConstrictionSlotState::advance(double step) noexcept {
    activation.advance(step);
    place.advance(step);
    target_area.advance(step);
    width.advance(step);
    tip_down.advance(step);
}

void ArticulatorState::advance_geometry(double step) noexcept {
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        shape[parameter].advance(step);
    }
    lip_area.advance(step);
    lip_protrusion.advance(step);
    for (std::size_t index = 0; index < constriction_slot_count; ++index) {
        slot[index].advance(step);
    }
    port_area.advance(step);
    body_offset_x.advance(step);
    body_offset_y.advance(step);
}

void ArticulatorState::advance_physiology(double step) noexcept {
    glottal_rest_area.advance(step);
    lung_pressure.advance(step);
    pitch.advance(step);
    effort.advance(step);
    loudness.advance(step);
}

void ArticulatorState::advance(double step) noexcept {
    advance_geometry(step);
    advance_physiology(step);
}

void ArticulatorState::apply_speeds(const ArticulatorSpeeds& speeds) noexcept {
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        shape[parameter].speed_limit = speeds.tongue_body;
    }
    shape[jaw_parameter].speed_limit = speeds.jaw;
    shape[blade_parameter].speed_limit = speeds.tongue_tip;
    lip_area.speed_limit = speeds.lips;
    lip_protrusion.speed_limit = speeds.lips;
    const std::array<double, constriction_slot_count> slot_speed{speeds.lips, speeds.tongue_tip, speeds.tongue_body};
    for (std::size_t index = 0; index < constriction_slot_count; ++index) {
        slot[index].activation.speed_limit = slot_speed[index];
        slot[index].place.speed_limit = slot_speed[index];
        slot[index].target_area.speed_limit = slot_speed[index];
        slot[index].width.speed_limit = slot_speed[index];
    }
    body_offset_x.speed_limit = speeds.tongue_body;
    body_offset_y.speed_limit = speeds.tongue_body;
    body_offset_x.angular_rate = speeds.tongue_body;
    body_offset_y.angular_rate = speeds.tongue_body;
    port_area.speed_limit = speeds.velum;
    glottal_rest_area.speed_limit = speeds.glottis;
    lung_pressure.speed_limit = speeds.lungs;
    pitch.speed_limit = speeds.larynx_pitch;
}

void ArticulatorState::place_vowel(const VowelShape& vowel) noexcept {
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        shape[parameter].place(vowel.coefficients[parameter]);
    }
    lip_area.place(vowel.lip_area);
    lip_protrusion.place(vowel.lip_protrusion);
}

void ArticulatorState::target_vowel(const VowelShape& vowel, double angular_rate) noexcept {
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        shape[parameter].target = vowel.coefficients[parameter];
        shape[parameter].angular_rate = angular_rate;
    }
    lip_area.target = vowel.lip_area;
    lip_area.angular_rate = angular_rate * 1.3;
    lip_protrusion.target = vowel.lip_protrusion;
    lip_protrusion.angular_rate = angular_rate * 1.3;
}

void TractModel::prepare(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections) {
    if (oral_sections < 24 || nasal_sections < 6) {
        throw std::invalid_argument("TractModel needs at least 24 oral and 6 nasal sections");
    }
    oral_sections_ = oral_sections;
    nasal_sections_ = nasal_sections;
    length_scale_ = voice.tract_length / reference_tract_length;
    area_scale_ = length_scale_ * length_scale_;
    section_length_ = voice.tract_length / static_cast<double>(oral_sections);
    double count = static_cast<double>(oral_sections);
    lip_sections_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::lround(lip_fraction * count)));
    epilarynx_sections_ = static_cast<std::size_t>(std::lround(epilarynx_fraction * count));
    std::size_t port = static_cast<std::size_t>(std::lround(voice.port_position * count));
    port_section_ = std::clamp(port, epilarynx_sections_ + 1, oral_sections - lip_sections_ - 2);

    neutral_log_area_.assign(oral_sections, 0.0);
    arc_.assign(oral_sections, 0.0);
    gap_.assign(oral_sections, 0.0);
    lower_gap_.assign(oral_sections, 0.0);
    open_above_.assign(oral_sections, 0.0);
    residual_basis_.assign(oral_sections * anatomy_parameter_count, 0.0);
    double shaped_start = static_cast<double>(epilarynx_sections_) / count;
    double shaped_end = 1.0 - static_cast<double>(lip_sections_) / count;
    for (std::size_t section = 0; section < oral_sections; ++section) {
        double position = (static_cast<double>(section) + 0.5) / count;
        neutral_log_area_[section] = interpolate_log(neutral_knots, position);
        arc_[section] = position * MidsagittalTongue::length;
        if (position < shaped_start || position >= shaped_end) {
            continue;
        }
        double local = (position - shaped_start) / (shaped_end - shaped_start);
        double taper = smoothstep(std::min(local, 1.0 - local) / 0.1);
        for (std::size_t mode = 0; mode < anatomy_parameter_count; ++mode) {
            double wave = std::cos(static_cast<double>(mode + 1) * pi * local);
            residual_basis_[section * anatomy_parameter_count + mode] = wave * taper;
        }
    }

    nasal_area_.assign(nasal_sections, 0.0);
    nasal_length_.assign(nasal_sections, voice.nasal_length / static_cast<double>(nasal_sections));
    for (std::size_t section = 0; section < nasal_sections; ++section) {
        double position = (static_cast<double>(section) + 0.5) / static_cast<double>(nasal_sections);
        double area = std::exp(interpolate_log(nasal_knots, position)) * square_centimetre * area_scale_;
        nasal_area_[section] = area;
    }
    nasal_area_[nasal_sections - 1] = voice.nostril_area;

    // Side cavities, reference geometry scaled with the voice. Values are
    // modelling choices shaped after published adult measurements, not a
    // particular speaker.
    side_cavities_.clear();
    side_area_.clear();
    side_length_.clear();
    if (voice.side_cavities) {
        // Piriform fossae (both sides as one branch): open beside the
        // epilarynx, about 2 cm deep; a quarter-wave antiresonance near 4-5 kHz.
        std::size_t laryngeal = std::max<std::size_t>(1, epilarynx_sections_);
        side_cavities_.push_back(SideCavity{.attach_node = laryngeal, .sections = 2, .perimeter_factor = 1.3});
        for (std::size_t section = 0; section < 2; ++section) {
            side_area_.push_back(1.0 * square_centimetre * area_scale_);
            side_length_.push_back(1.0 * centimetre * length_scale_);
        }
        // Maxillary sinuses (both as one Helmholtz resonator): a narrow ostium
        // into a large bony cavity, opening from the middle of the nasal tract.
        std::size_t nasal_middle = oral_sections + static_cast<std::size_t>(std::lround(0.4 * static_cast<double>(nasal_sections)));
        side_cavities_.push_back(SideCavity{.attach_node = nasal_middle, .sections = 3, .perimeter_factor = 1.5});
        side_area_.push_back(0.14 * square_centimetre * area_scale_);
        side_length_.push_back(0.5 * centimetre * length_scale_);
        for (std::size_t section = 0; section < 2; ++section) {
            side_area_.push_back(10.0 * square_centimetre * area_scale_);
            side_length_.push_back(1.5 * centimetre * length_scale_);
        }
    }
    distance_.assign(oral_sections, 0.0);
    allocate(scratch_);
    narrowest_ = 0;
    narrowest_area_ = 0.0;
    oral_volume_ = 0.0;
}

void TractModel::allocate(TractShape& shape) const {
    shape.oral_area.assign(oral_sections_, 0.0);
    shape.oral_length.assign(oral_sections_, section_length_);
    shape.nasal_area.assign(nasal_sections_, 0.0);
    shape.nasal_length.assign(nasal_sections_, 0.0);
    shape.port_area = 0.0;
    shape.lip_aperture = -1.0;
    shape.port_length = 0.01 * length_scale_;
    shape.side_area = side_area_;
    shape.side_length = side_length_;
    shape.sublingual_area.assign(oral_sections_, 0.0);
    shape.sublingual_opening.assign(oral_sections_, 0.0);
}

TongueShape TractModel::pose(const TonguePosture& posture,const TongueGesture& body,const TongueGesture& tip) noexcept {
    if (!conservative_body_) { return MidsagittalTongue::pose(posture,body,tip); }
    const Point2 centre{posture.body_x+posture.offset_x,posture.body_y-.9*std::clamp(posture.jaw,0.0,1.0)+posture.offset_y};
    if (!body_surface_.configure(centre,body)) {
        geometry_error_=body_surface_.error();TongueShape failed{};failed.valid=false;return failed;
    }
    const TongueShape result=MidsagittalTongue::pose(posture,body,tip,&body_surface_);
    if (!result.valid) { geometry_error_=body_surface_.error()?body_surface_.error():"blade attachment material mapping failed"; }
    return result;
}

void TractModel::write_posture(const TonguePosture& posture, const TongueGesture& body, const TongueGesture& tip,
                               std::span<const double> residual, double lip_area, double lip_protrusion,
                               bool vowel_floor, TractShape& shape) noexcept {
    std::size_t lip_start = oral_sections_ - lip_sections_;
    tongue_ = pose(posture, body, tip);
    if (!tongue_.valid) { return; }
    MidsagittalTongue::gaps(tongue_, posture, std::span<const double>(arc_.data(), lip_start),
                            std::span<double>(gap_.data(), lip_start));
    for (std::size_t i=0;i<lip_start;++i) {
        if (!std::isfinite(gap_[i])) { geometry_error_="mapped airway intersection failed";return; }
    }
    // The space under the tongue: a channel about as wide as the arch of the
    // lower teeth (2.2 cm in the reference speaker) between the underside
    // and the floor (a modelling choice).
    const double sublingual_width = 2.2; // cm
    bool with_floor = floor_ && shape.sublingual_area.size() == oral_sections_;
    if (with_floor) {
        MidsagittalTongue::floor_space(tongue_, posture, std::span<const double>(arc_.data(), lip_start),
                                       std::span<double>(gap_.data(), lip_start),
                                       std::span<double>(lower_gap_.data(), lip_start),
                                       std::span<double>(open_above_.data(), lip_start));
        for (std::size_t section = 0; section < oral_sections_; ++section) {
            double lower = section < lip_start && section >= epilarynx_sections_ ? lower_gap_[section] : 0.0;
            double face = section < lip_start && section >= epilarynx_sections_ ? open_above_[section] : 0.0;
            shape.sublingual_area[section] = sublingual_width * lower * square_centimetre * area_scale_;
            shape.sublingual_opening[section] = sublingual_width * face * square_centimetre * area_scale_;
        }
    }
    for (std::size_t section = 0; section < lip_start; ++section) {
        double area = 0.0;
        if (section < epilarynx_sections_) {
            area = std::exp(neutral_log_area_[section]);
        } else {
            area = MidsagittalTongue::area_from_gap(gap_[section], arc_[section]);
            double log_residual = 0.0;
            for (std::size_t mode = 0; mode < anatomy_parameter_count; ++mode) {
                log_residual += residual[mode] * residual_basis_[section * anatomy_parameter_count + mode];
            }
            area *= std::exp(log_residual);
            if (vowel_floor) {
                // A vowel never seals the tract; contact belongs to consonants.
                area = std::max(area, minimum_vowel_area);
            }
        }
        shape.oral_area[section] = area * square_centimetre * area_scale_;
        shape.oral_length[section] = section_length_;
    }
    double protrusion = std::max(lip_protrusion, 0.0) * length_scale_;
    double lip_length = section_length_ + protrusion / static_cast<double>(lip_sections_);
    double opening = std::max(lip_area, 0.0) * area_scale_;
    for (std::size_t section = lip_start; section < oral_sections_; ++section) {
        shape.oral_area[section] = opening;
        shape.oral_length[section] = lip_length;
    }
    double running = 0.5 * shape.oral_length[oral_sections_ - 1];
    distance_[oral_sections_ - 1] = running;
    for (std::size_t section = oral_sections_ - 1; section > 0; --section) {
        running += 0.5 * (shape.oral_length[section] + shape.oral_length[section - 1]);
        distance_[section - 1] = running;
    }
}

void TractModel::apply_constriction(const ConstrictionGesture& gesture, double activation, TractShape& shape) noexcept {
    // Tissue contact: the narrowing follows a smoothstep of activation, so the
    // area meets its target (and a closure meets zero) with zero rate of
    // change. A linear map closes with a kink that splatters energy into the
    // highest modes, heard as a click.
    double amount = smoothstep(activation);
    if (amount <= 0.0) {
        return;
    }
    double place = gesture.place_from_lips * length_scale_;
    double width = gesture.width * length_scale_;
    double core = constriction_core * length_scale_;
    double target = gesture.target_area * area_scale_;
    for (std::size_t section = epilarynx_sections_; section < oral_sections_; ++section) {
        double separation = std::fabs(distance_[section] - place);
        double weight = 1.0;
        if (separation > core) {
            double normalized = (separation - core) / width;
            weight = std::exp(-normalized * normalized);
        }
        if (weight < 1.0e-4) {
            continue;
        }
        double area = shape.oral_area[section];
        if (target < area) {
            shape.oral_area[section] = area + amount * weight * (target - area);
        }
    }
}

void TractModel::finish(TractShape& shape) noexcept {
    for (std::size_t section = 0; section < nasal_sections_; ++section) {
        shape.nasal_area[section] = nasal_area_[section];
        shape.nasal_length[section] = nasal_length_[section];
    }
    narrowest_ = epilarynx_sections_;
    narrowest_area_ = shape.oral_area[epilarynx_sections_];
    double volume = 0.0;
    for (std::size_t section = 0; section < oral_sections_; ++section) {
        volume += shape.oral_area[section] * shape.oral_length[section];
        if (section >= epilarynx_sections_ && shape.oral_area[section] < narrowest_area_) {
            narrowest_area_ = shape.oral_area[section];
            narrowest_ = section;
        }
    }
    oral_volume_ = volume;
    if (shape.lip_aperture >= 0 && shape.lip_aperture < narrowest_area_) {
        narrowest_area_ = shape.lip_aperture;
        narrowest_ = oral_sections_ - 1;
    }
}

namespace {

TongueGesture tongue_gesture(const ConstrictionSlotState& slot) noexcept {
    ConstrictionGesture gesture = slot.current();
    TongueGesture result{};
    result.activation = std::max(slot.activation.value, 0.0);
    result.place_from_lips = gesture.place_from_lips / centimetre;
    double arc = MidsagittalTongue::length - result.place_from_lips;
    result.gap = MidsagittalTongue::gap_from_area(gesture.target_area / square_centimetre, arc);
    result.tip_down = gesture.tip_down;
    result.width = gesture.width / centimetre;
    return result;
}

TonguePosture posture_from(const std::array<double, shape_parameter_count>& coefficients) noexcept {
    TonguePosture posture{};
    posture.body_x = coefficients[body_x_parameter];
    posture.body_y = coefficients[body_y_parameter];
    posture.jaw = coefficients[jaw_parameter];
    posture.blade = coefficients[blade_parameter];
    posture.root = coefficients[root_parameter];
    return posture;
}

}

void TractModel::build(const ArticulatorState& state, TractShape& shape) noexcept {
    geometry_error_=nullptr;
    std::array<double, shape_parameter_count> coefficients{};
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        coefficients[parameter] = state.shape[parameter].value;
    }
    const ConstrictionSlotState& tip_slot = state.slot[static_cast<std::size_t>(ConstrictionSlot::tongue_tip)];
    const ConstrictionSlotState& body_slot = state.slot[static_cast<std::size_t>(ConstrictionSlot::tongue_body)];
    const ConstrictionSlotState& lip_slot = state.slot[static_cast<std::size_t>(ConstrictionSlot::lips)];
    std::span<const double> residual(coefficients.data() + tongue_parameter_count, anatomy_parameter_count);
    TonguePosture posture = posture_from(coefficients);
    TongueGesture body_gesture = tongue_gesture(body_slot);
    TongueGesture tip_gesture = tongue_gesture(tip_slot);
    // What the tip demands of the posture alone, and the posture as the body
    // actually is: displaced by its lingering offset. The blade constraint
    // still holds exactly while the tip is down.
    required_drag_ = pose(posture, body_gesture, tip_gesture).drag_vector;
    if (geometry_error_) { return; }
    posture.offset_x = state.body_offset_x.value;
    posture.offset_y = state.body_offset_y.value;
    write_posture(posture, body_gesture, tip_gesture, residual, state.lip_area.value, state.lip_protrusion.value,
                  false, shape);
    apply_constriction(lip_slot.current(), lip_slot.activation.value, shape);
    const LipProfile lips=lip_profile(state);
    shape.lip_aperture = lips.area(0.0);
    double outward=0;
    for(std::size_t section=oral_sections_;section-->oral_sections_-lip_sections_;){
        const double inward=outward+shape.oral_length[section];
        shape.oral_area[section]=lip_storage_area(lips,outward,inward);
        outward=inward;
    }
    shape.port_area = std::max(state.port_area.value, 0.0) * area_scale_;
    finish(shape);
}

void TractModel::build_vowel(const VowelShape& vowel, TractShape& shape) noexcept {
    geometry_error_=nullptr;
    TongueGesture idle{};
    std::span<const double> residual(vowel.coefficients.data() + tongue_parameter_count, anatomy_parameter_count);
    write_posture(posture_from(vowel.coefficients), idle, idle, residual, vowel.lip_area, vowel.lip_protrusion, true, shape);
    shape.lip_aperture = shape.oral_area.back();
    shape.port_area = 0.0;
    finish(shape);
}

void TractModel::log_area_gradient(const VowelShape& vowel, std::span<double> output) {
    const double step = 0.01;
    for (std::size_t parameter = 0; parameter < shape_parameter_count; ++parameter) {
        VowelShape plus = vowel;
        VowelShape minus = vowel;
        plus.coefficients[parameter] += step;
        minus.coefficients[parameter] -= step;
        build_vowel(plus, scratch_);
        for (std::size_t section = 0; section < oral_sections_; ++section) {
            output[section * shape_parameter_count + parameter] = std::log(scratch_.oral_area[section]);
        }
        build_vowel(minus, scratch_);
        for (std::size_t section = 0; section < oral_sections_; ++section) {
            double high = output[section * shape_parameter_count + parameter];
            double low = std::log(scratch_.oral_area[section]);
            output[section * shape_parameter_count + parameter] = (high - low) / (2.0 * step);
        }
    }
}

Point2 TractModel::required_drag() const noexcept {
    return required_drag_;
}

const TongueShape& TractModel::tongue() const noexcept {
    return tongue_;
}

std::size_t TractModel::oral_sections() const noexcept {
    return oral_sections_;
}

std::size_t TractModel::nasal_sections() const noexcept {
    return nasal_sections_;
}

std::size_t TractModel::port_section() const noexcept {
    return port_section_;
}

std::size_t TractModel::epilarynx_sections() const noexcept {
    return epilarynx_sections_;
}

std::size_t TractModel::lip_sections() const noexcept {
    return lip_sections_;
}

LipProfile TractModel::lip_profile(const ArticulatorState& state) const noexcept {
    const auto& slot=state.slot[static_cast<std::size_t>(ConstrictionSlot::lips)];
    const ConstrictionGesture gesture=slot.current();
    return {.length=double(lip_sections_)*section_length_+std::max(state.lip_protrusion.value,0.)*length_scale_,
            .base_area=std::max(state.lip_area.value,0.)*area_scale_,.activation=slot.activation.value,
            .place=gesture.place_from_lips*length_scale_,.core=constriction_core*length_scale_,
            .width=gesture.width*length_scale_,.target_area=gesture.target_area*area_scale_};
}

double TractModel::length_scale() const noexcept {
    return length_scale_;
}

double TractModel::area_scale() const noexcept {
    return area_scale_;
}

double TractModel::distance_from_lips(std::size_t section) const noexcept {
    return distance_[section];
}

std::size_t TractModel::narrowest_section() const noexcept {
    return narrowest_;
}

double TractModel::narrowest_area() const noexcept {
    return narrowest_area_;
}

std::size_t TractModel::section_at_distance(double reference_distance) const noexcept {
    double wanted = reference_distance * length_scale_;
    std::size_t best = oral_sections_ - 1;
    double best_gap = std::fabs(distance_[best] - wanted);
    for (std::size_t section = 0; section < oral_sections_; ++section) {
        double gap = std::fabs(distance_[section] - wanted);
        if (gap < best_gap) {
            best_gap = gap;
            best = section;
        }
    }
    return best;
}

std::span<const SideCavity> TractModel::side_cavities() const noexcept {
    return std::span<const SideCavity>(side_cavities_);
}

void TractModel::set_floor(bool enabled) noexcept {
    floor_ = enabled;
}

void TractModel::floor_sections(std::size_t& first, std::size_t& count) const noexcept {
    first = oral_sections_;
    count = 0;
    std::size_t lip_start = oral_sections_ - lip_sections_;
    for (std::size_t section = epilarynx_sections_; section < lip_start; ++section) {
        double arc = arc_[section];
        double x = MidsagittalTongue::midline(arc).x;
        bool oral = arc > MidsagittalTongue::pharynx_height + 0.5 * pi * MidsagittalTongue::bend_radius;
        if (oral && x >= MidsagittalTongue::floor_back_x && x < MidsagittalTongue::incisors_x) {
            if (count == 0) {
                first = section;
            }
            count += 1;
        }
    }
}

double TractModel::oral_volume() const noexcept {
    return oral_volume_;
}

}
