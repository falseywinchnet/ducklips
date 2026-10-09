#include "ducklips/air_path.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdlib>

namespace ducklips {
namespace {

double distance(const Position3& a, const Position3& b) noexcept {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    double length = std::sqrt(dx * dx + dy * dy + dz * dz);
    return length;
}

Position3 offset(const Position3& origin, const Position3& relative) noexcept {
    Position3 sum{.x = origin.x + relative.x, .y = origin.y + relative.y, .z = origin.z + relative.z};
    return sum;
}

// Image coordinate along one axis for a source at s in a room [0, size]:
// mirror parity u in {0, 1}, period n. Reflections: |n - u| at the low wall
// plus |n| at the high wall (Allen and Berkley 1979).
struct AxisImage final {
    double coordinate{};
    long reflections{};
};

AxisImage axis_image(double source, double size, long period, long parity) noexcept {
    double mirrored = source;
    if (parity == 1) {
        mirrored = -source;
    }
    AxisImage image{};
    image.coordinate = mirrored + 2.0 * static_cast<double>(period) * size;
    image.reflections = std::labs(period - parity) + std::labs(period);
    return image;
}

}

void AirPath::add_links(const Position3& source, std::size_t radiator, const Position3& receiver,
                        std::size_t microphone, const AirPathOptions& options, const AirProperties& air,
                        SynthesisRate rate) {
    const double spreading = air.density / (4.0 * pi);
    if (!options.room.enabled) {
        double direct = distance(source, receiver);
        links_.push_back(Link{.radiator = radiator, .microphone = microphone,
                              .delay_samples = direct / air.sound_speed * rate.samples_per_second,
                              .gain = spreading / direct, .smoothing = 1.0, .filter_state = 0.0});
        return;
    }
    const RoomGeometry& room = options.room;
    Position3 room_source = offset(room.head, source);
    Position3 room_receiver = offset(room.head, receiver);
    double direct_gain = spreading / distance(room_source, room_receiver);
    long span = static_cast<long>(room.maximum_order);
    for (long px = -span; px <= span; ++px) {
        for (long ux = 0; ux <= 1; ++ux) {
            AxisImage ix = axis_image(room_source.x, room.length, px, ux);
            for (long py = -span; py <= span; ++py) {
                for (long uy = 0; uy <= 1; ++uy) {
                    AxisImage iy = axis_image(room_source.y, room.width, py, uy);
                    for (long pz = -span; pz <= span; ++pz) {
                        for (long uz = 0; uz <= 1; ++uz) {
                            AxisImage iz = axis_image(room_source.z, room.height, pz, uz);
                            long order = ix.reflections + iy.reflections + iz.reflections;
                            if (order > span) {
                                continue;
                            }
                            Position3 image{.x = ix.coordinate, .y = iy.coordinate, .z = iz.coordinate};
                            double path = distance(image, room_receiver);
                            double gain = spreading / path * std::pow(room.wall_reflection, static_cast<double>(order));
                            if (order > 0 && std::fabs(gain) < room.cutoff_ratio * direct_gain) {
                                continue;
                            }
                            double smoothing = 1.0;
                            if (order > 0) {
                                double cutoff = 9000.0 / (1.0 + static_cast<double>(order));
                                smoothing = 1.0 - std::exp(-2.0 * pi * cutoff / rate.samples_per_second);
                            }
                            links_.push_back(Link{.radiator = radiator, .microphone = microphone,
                                                  .delay_samples = path / air.sound_speed * rate.samples_per_second,
                                                  .gain = gain, .smoothing = smoothing, .filter_state = 0.0});
                        }
                    }
                }
            }
        }
    }
}

void AirPath::prepare(const AirPathOptions& options, const AirProperties& air, SynthesisRate rate) {
    links_.clear();
    std::array<Position3, radiator_count> sources{options.lips, options.nostrils, options.throat};
    std::array<Position3, microphone_count> receivers{options.left_microphone, options.right_microphone};
    for (std::size_t radiator = 0; radiator < radiator_count; ++radiator) {
        for (std::size_t microphone = 0; microphone < microphone_count; ++microphone) {
            add_links(sources[radiator], radiator, receivers[microphone], microphone, options, air, rate);
        }
    }
    double longest = 0.0;
    for (std::size_t link = 0; link < links_.size(); ++link) {
        longest = std::max(longest, links_[link].delay_samples);
    }
    std::size_t needed = static_cast<std::size_t>(std::ceil(longest)) + 8;
    ring_size_ = std::bit_ceil(needed);
    ring_mask_ = ring_size_ - 1;
    history_.assign(radiator_count * ring_size_, 0.0);
    write_index_ = 0;
    output_scale_ = 1.0 / options.full_scale_pressure;
    dc_coefficient_ = std::exp(-2.0 * pi * options.dc_block_hz / rate.samples_per_second);
    reset();
}

void AirPath::reset() noexcept {
    for (std::size_t index = 0; index < history_.size(); ++index) {
        history_[index] = 0.0;
    }
    for (std::size_t link = 0; link < links_.size(); ++link) {
        links_[link].filter_state = 0.0;
    }
    for (std::size_t microphone = 0; microphone < microphone_count; ++microphone) {
        dc_input_[microphone] = 0.0;
        dc_output_[microphone] = 0.0;
    }
}

double AirPath::read_delayed(std::size_t radiator, double delay_samples) const noexcept {
    // Cubic Lagrange interpolation between the samples at whole delays i and
    // i + 1, using i - 1 .. i + 2. Delays are at least two samples.
    double whole = std::floor(delay_samples);
    double u = delay_samples - whole;
    std::size_t delay = static_cast<std::size_t>(whole);
    std::size_t base = radiator * ring_size_;
    double newer = history_[base + ((write_index_ - delay + 1) & ring_mask_)];
    double current = history_[base + ((write_index_ - delay) & ring_mask_)];
    double older = history_[base + ((write_index_ - delay - 1) & ring_mask_)];
    double oldest = history_[base + ((write_index_ - delay - 2) & ring_mask_)];
    double c0 = -u * (u - 1.0) * (u - 2.0) / 6.0;
    double c1 = (u + 1.0) * (u - 1.0) * (u - 2.0) / 2.0;
    double c2 = -(u + 1.0) * u * (u - 2.0) / 2.0;
    double c3 = (u + 1.0) * u * (u - 1.0) / 6.0;
    double value = c0 * newer + c1 * current + c2 * older + c3 * oldest;
    return value;
}

void AirPath::process(const std::array<double, radiator_count>& flow_derivative,
                      std::array<double, microphone_count>& microphone) noexcept {
    for (std::size_t radiator = 0; radiator < radiator_count; ++radiator) {
        history_[radiator * ring_size_ + write_index_] = flow_derivative[radiator];
    }
    std::array<double, microphone_count> pressure{0.0, 0.0};
    for (std::size_t index = 0; index < links_.size(); ++index) {
        Link& link = links_[index];
        double arriving = link.gain * read_delayed(link.radiator, link.delay_samples);
        if (link.smoothing < 1.0) {
            link.filter_state += link.smoothing * (arriving - link.filter_state);
            arriving = link.filter_state;
        }
        pressure[link.microphone] += arriving;
    }
    for (std::size_t channel = 0; channel < microphone_count; ++channel) {
        double input = pressure[channel];
        double output = dc_coefficient_ * (dc_output_[channel] + input - dc_input_[channel]);
        dc_input_[channel] = input;
        dc_output_[channel] = output;
        microphone[channel] = output * output_scale_;
    }
    write_index_ = (write_index_ + 1) & ring_mask_;
}

std::size_t AirPath::link_count() const noexcept {
    return links_.size();
}

}
