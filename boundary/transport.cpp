#include "transport.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>

namespace boundary {

namespace {

constexpr double four_pi = 4.0 * ducklips::pi;
constexpr double offset_into_air = 1.0e-5; // m; lifts ray ends off their own surface

Vec3 bilinear_point(const std::array<Vec3, 4>& c, double s, double t) {
    Vec3 bottom = (1.0 - s) * c[0] + s * c[1];
    Vec3 top = (1.0 - s) * c[3] + s * c[2];
    Vec3 result = (1.0 - t) * bottom + t * top;
    return result;
}

// Area element |dP/ds x dP/dt| of the bilinear map at (s, t).
double bilinear_jacobian(const std::array<Vec3, 4>& c, double s, double t) {
    Vec3 along_s = (1.0 - t) * (c[1] - c[0]) + t * (c[2] - c[3]);
    Vec3 along_t = (1.0 - s) * (c[3] - c[0]) + s * (c[2] - c[1]);
    double result = norm(cross(along_s, along_t));
    return result;
}

// Moller-Trumbore. Returns true when the open segment from origin to
// origin + direction crosses the triangle strictly between its ends.
bool segment_hits_triangle(const Vec3& origin, const Vec3& direction, const Vec3& a, const Vec3& b,
                           const Vec3& c) {
    const double tolerance = 1.0e-12;
    Vec3 edge_1 = b - a;
    Vec3 edge_2 = c - a;
    Vec3 p = cross(direction, edge_2);
    double determinant = dot(edge_1, p);
    if (std::fabs(determinant) < tolerance) {
        return false;
    }
    double inverse = 1.0 / determinant;
    Vec3 offset = origin - a;
    double u = dot(offset, p) * inverse;
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    Vec3 q = cross(offset, edge_1);
    double v = dot(direction, q) * inverse;
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    double t = dot(edge_2, q) * inverse;
    bool hit = t > 1.0e-9 && t < 1.0 - 1.0e-9;
    return hit;
}

// Cyclic Jacobi for a dense symmetric matrix (row-major, n x n, destroyed).
// Eigenvector k is column k of vectors.
void jacobi_eigen(std::vector<double>& matrix, std::size_t n, std::vector<double>& values,
                  std::vector<double>& vectors) {
    vectors.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        vectors[i * n + i] = 1.0;
    }
    for (std::size_t sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                off += matrix[i * n + j] * matrix[i * n + j];
            }
        }
        if (off < 1.0e-30) {
            break;
        }
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                double apq = matrix[p * n + q];
                if (std::fabs(apq) < 1.0e-300) {
                    continue;
                }
                double app = matrix[p * n + p];
                double aqq = matrix[q * n + q];
                double theta = 0.5 * (aqq - app) / apq;
                double sign = 1.0;
                if (theta < 0.0) {
                    sign = -1.0;
                }
                double t = sign / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double s = t * c;
                for (std::size_t k = 0; k < n; ++k) {
                    double akp = matrix[k * n + p];
                    double akq = matrix[k * n + q];
                    matrix[k * n + p] = c * akp - s * akq;
                    matrix[k * n + q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    double apk = matrix[p * n + k];
                    double aqk = matrix[q * n + k];
                    matrix[p * n + k] = c * apk - s * aqk;
                    matrix[q * n + k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    double vkp = vectors[k * n + p];
                    double vkq = vectors[k * n + q];
                    vectors[k * n + p] = c * vkp - s * vkq;
                    vectors[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    values.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = matrix[i * n + i];
    }
}

}

bool BoundaryTransport::TapOrder::operator()(const Tap& a, const Tap& b) const noexcept {
    if (a.delay != b.delay) {
        return a.delay < b.delay;
    }
    return a.source < b.source;
}

bool BoundaryTransport::segment_blocked(const Scene& scene, const Vec3& from, const Vec3& to) const {
    Vec3 direction = to - from;
    for (std::size_t s = 0; s < scene.surfaces.size(); ++s) {
        const std::array<Vec3, 4>& c = scene.surfaces[s].corners;
        if (segment_hits_triangle(from, direction, c[0], c[1], c[2])) {
            return true;
        }
        if (segment_hits_triangle(from, direction, c[0], c[2], c[3])) {
            return true;
        }
    }
    return false;
}

void BoundaryTransport::add_interpolated(std::vector<double>& by_delay, double delay, double weight) {
    double clamped = std::max(delay, 1.0);
    double whole = std::floor(clamped);
    double fraction = clamped - whole;
    std::size_t index = static_cast<std::size_t>(whole);
    if (index + 1 >= by_delay.size()) {
        return;
    }
    by_delay[index] += weight * (1.0 - fraction);
    by_delay[index + 1] += weight * fraction;
}

void BoundaryTransport::flush_taps(std::vector<double>& by_delay, std::size_t source, std::vector<Tap>& taps) {
    for (std::size_t d = 0; d < by_delay.size(); ++d) {
        if (by_delay[d] != 0.0) {
            Tap tap{};
            tap.delay = static_cast<std::uint32_t>(d);
            tap.source = static_cast<std::uint32_t>(source);
            tap.weight = by_delay[d];
            taps.push_back(tap);
            by_delay[d] = 0.0;
        }
    }
}

void BoundaryTransport::accumulate_pair(const Scene& scene, const Vec3& receiver, const Vec3& receiver_normal,
                                        bool has_normal, std::size_t source, double factor,
                                        std::vector<double>& pressure_by_delay,
                                        std::vector<double>& acceleration_by_delay) {
    const Patch& patch = scene.patches[source];
    const double c = air_.sound_speed;
    const double rho = air_.density;
    Vec3 to_centroid = patch.centroid - receiver;
    double centroid_distance = norm(to_centroid);
    double edge = std::sqrt(patch.area);
    double wanted = std::ceil(quadrature_resolution_ * edge / std::max(centroid_distance, 1.0e-9));
    std::size_t q = static_cast<std::size_t>(std::clamp(wanted, 2.0, 8.0));
    Vec3 receiver_lifted = receiver;
    if (has_normal) {
        receiver_lifted = receiver - offset_into_air * receiver_normal;
    }
    for (std::size_t a = 0; a < q; ++a) {
        for (std::size_t b = 0; b < q; ++b) {
            double s = (static_cast<double>(a) + 0.5) / static_cast<double>(q);
            double t = (static_cast<double>(b) + 0.5) / static_cast<double>(q);
            Vec3 y = bilinear_point(patch.corners, s, t);
            double weight = bilinear_jacobian(patch.corners, s, t) / static_cast<double>(q * q);
            Vec3 d = y - receiver;
            double r = norm(d);
            if (r < 1.0e-9) {
                continue;
            }
            double cos_b = dot(d, patch.normal_into_wall) / r;
            double facing = cos_b;
            if (visibility_culling_) {
                if (has_normal && dot(d, receiver_normal) > 1.0e-9 * r) {
                    continue; // behind the receiver's wall
                }
                if (cos_b < -1.0e-9) {
                    continue; // the receiver is behind the source's wall
                }
                Vec3 y_lifted = y - offset_into_air * patch.normal_into_wall;
                if (segment_blocked(scene, receiver_lifted, y_lifted)) {
                    statistics_.occluded_subpoints += 1;
                    continue;
                }
                facing = std::max(cos_b, 0.0);
            }
            double delay = r / c * rate_;
            if (delay < 1.5) {
                statistics_.clamped_subpoints += 1;
            }
            double value_gain = factor * weight * facing / (four_pi * r * r);
            double rate_gain = factor * weight * facing / (four_pi * c * r);
            double acceleration_gain = -factor * rho * weight / (four_pi * r);
            add_interpolated(pressure_by_delay, delay, value_gain);
            // Central difference about the delayed instant.
            add_interpolated(pressure_by_delay, delay - 0.5, rate_gain * rate_);
            add_interpolated(pressure_by_delay, delay + 0.5, -rate_gain * rate_);
            add_interpolated(acceleration_by_delay, delay, acceleration_gain);
        }
    }
}

void BoundaryTransport::compile(const Scene& scene, const TransportSettings& settings) {
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
    statistics_ = CompileStatistics{};
    air_ = settings.air;
    quadrature_resolution_ = settings.quadrature_resolution;
    visibility_culling_ = settings.visibility_culling;
    temporal_averaging_ = settings.temporal_averaging;
    rate_ = settings.audio_rate * static_cast<double>(settings.oversampling);
    time_step_ = 1.0 / rate_;
    patch_count_ = scene.patches.size();
    statistics_.patches = patch_count_;

    // Longest path: scene extent and microphone distance.
    double longest = 0.0;
    for (std::size_t i = 0; i < patch_count_; ++i) {
        for (std::size_t k = 0; k < 4; ++k) {
            Vec3 corner = scene.patches[i].corners[k];
            longest = std::max(longest, norm(corner - scene.microphone));
            longest = std::max(longest, norm(corner - scene.patches[0].corners[0]) * 2.0);
        }
    }
    largest_delay_ = static_cast<std::size_t>(std::ceil(longest / air_.sound_speed * rate_)) + 4;
    std::vector<double> pressure_by_delay(largest_delay_ + 2, 0.0);
    std::vector<double> acceleration_by_delay(largest_delay_ + 2, 0.0);
    std::vector<Tap> receiver_pressure{};
    std::vector<Tap> receiver_acceleration{};
    TapOrder order{};

    pressure_begin_.assign(patch_count_ + 1, 0);
    all_acceleration_begin_.assign(patch_count_ + 1, 0);
    pressure_taps_.clear();
    all_acceleration_taps_.clear();
    for (std::size_t j = 0; j < patch_count_; ++j) {
        const Patch& receiver = scene.patches[j];
        TissueProperties properties = tissue_properties(scene.surfaces[receiver.surface].tissue);
        double reflectance = properties.reflectance;
        if (settings.reflectance_override > 0.0) {
            reflectance = settings.reflectance_override;
        }
        double factor = 2.0 * reflectance;
        receiver_pressure.clear();
        receiver_acceleration.clear();
        for (std::size_t i = 0; i < patch_count_; ++i) {
            if (i == j) {
                continue;
            }
            accumulate_pair(scene, receiver.centroid, receiver.normal_into_wall, true, i, factor, pressure_by_delay,
                            acceleration_by_delay);
            std::size_t before = receiver_pressure.size() + receiver_acceleration.size();
            flush_taps(pressure_by_delay, i, receiver_pressure);
            flush_taps(acceleration_by_delay, i, receiver_acceleration);
            if (receiver_pressure.size() + receiver_acceleration.size() > before) {
                statistics_.linked_pairs += 1;
            }
        }
        std::sort(receiver_pressure.begin(), receiver_pressure.end(), order);
        std::sort(receiver_acceleration.begin(), receiver_acceleration.end(), order);
        double static_gain = 0.0;
        for (std::size_t t = 0; t < receiver_pressure.size(); ++t) {
            static_gain += receiver_pressure[t].weight;
        }
        statistics_.largest_static_gain = std::max(statistics_.largest_static_gain, static_gain);
        pressure_taps_.insert(pressure_taps_.end(), receiver_pressure.begin(), receiver_pressure.end());
        all_acceleration_taps_.insert(all_acceleration_taps_.end(), receiver_acceleration.begin(),
                                      receiver_acceleration.end());
        pressure_begin_[j + 1] = pressure_taps_.size();
        all_acceleration_begin_[j + 1] = all_acceleration_taps_.size();
    }
    statistics_.pressure_taps = pressure_taps_.size();

    microphone_pressure_taps_.clear();
    all_microphone_acceleration_taps_.clear();
    Vec3 no_normal{};
    for (std::size_t i = 0; i < patch_count_; ++i) {
        accumulate_pair(scene, scene.microphone, no_normal, false, i, 1.0, pressure_by_delay, acceleration_by_delay);
        flush_taps(pressure_by_delay, i, microphone_pressure_taps_);
        flush_taps(acceleration_by_delay, i, all_microphone_acceleration_taps_);
    }
    std::sort(microphone_pressure_taps_.begin(), microphone_pressure_taps_.end(), order);
    std::sort(all_microphone_acceleration_taps_.begin(), all_microphone_acceleration_taps_.end(), order);

    ring_length_ = std::bit_ceil(largest_delay_ + 4);
    ring_mask_ = ring_length_ - 1;
    statistics_.ring_length = ring_length_;
    pressure_ring_.assign(ring_length_ * patch_count_, 0.0);
    acceleration_ring_.assign(ring_length_ * patch_count_, 0.0);
    forcing_sign_.assign(patch_count_, 0.0);
    injection_sign_.assign(patch_count_, 0.0);
    area_.assign(patch_count_, 0.0);
    for (std::size_t i = 0; i < patch_count_; ++i) {
        area_[i] = scene.patches[i].area;
    }
    if (scene.source_kind == SourceKind::pressure_dipole) {
        for (std::size_t k = 0; k < scene.source_positive.size(); ++k) {
            forcing_sign_[scene.source_positive[k]] = 1.0;
        }
        for (std::size_t k = 0; k < scene.source_negative.size(); ++k) {
            forcing_sign_[scene.source_negative[k]] = -1.0;
        }
    } else {
        for (std::size_t k = 0; k < scene.source_positive.size(); ++k) {
            injection_sign_[scene.source_positive[k]] = 1.0;
        }
        for (std::size_t k = 0; k < scene.source_negative.size(); ++k) {
            injection_sign_[scene.source_negative[k]] = -1.0;
        }
    }

    mode_first_patch_.clear();
    mode_patch_count_.clear();
    mode_shape_offset_.clear();
    mode_shape_.clear();
    mode_omega_.clear();
    mode_damping_.clear();
    if (settings.membranes) {
        compile_membranes(scene);
    }
    mode_position_.assign(mode_omega_.size(), 0.0);
    mode_velocity_.assign(mode_omega_.size(), 0.0);
    membrane_acceleration_.assign(patch_count_, 0.0);
    statistics_.membrane_modes = mode_omega_.size();
    patch_moves_.assign(patch_count_, false);
    if (settings.membranes) {
        for (std::size_t i = 0; i < patch_count_; ++i) {
            patch_moves_[i] = tissue_properties(scene.surfaces[scene.patches[i].surface].tissue).elastic;
        }
    }
    select_active_acceleration_taps();
    step_index_ = 0;
    std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = finished - started;
    statistics_.compile_seconds = elapsed.count();
}

void BoundaryTransport::select_active_acceleration_taps() {
    acceleration_begin_.assign(patch_count_ + 1, 0);
    acceleration_taps_.clear();
    for (std::size_t j = 0; j < patch_count_; ++j) {
        for (std::size_t t = all_acceleration_begin_[j]; t < all_acceleration_begin_[j + 1]; ++t) {
            const Tap& tap = all_acceleration_taps_[t];
            if (patch_moves_[tap.source] || injection_sign_[tap.source] != 0.0) {
                acceleration_taps_.push_back(tap);
            }
        }
        acceleration_begin_[j + 1] = acceleration_taps_.size();
    }
    microphone_acceleration_taps_.clear();
    for (std::size_t t = 0; t < all_microphone_acceleration_taps_.size(); ++t) {
        const Tap& tap = all_microphone_acceleration_taps_[t];
        if (patch_moves_[tap.source] || injection_sign_[tap.source] != 0.0) {
            microphone_acceleration_taps_.push_back(tap);
        }
    }
    acceleration_active_ = !acceleration_taps_.empty() || !microphone_acceleration_taps_.empty();
    statistics_.acceleration_taps = acceleration_taps_.size();
    statistics_.microphone_taps = microphone_pressure_taps_.size() + microphone_acceleration_taps_.size();
}

void BoundaryTransport::compile_membranes(const Scene& scene) {
    std::vector<double> mass{};
    std::vector<double> stiffness{};
    std::vector<double> values{};
    std::vector<double> vectors{};
    for (std::size_t s = 0; s < scene.surfaces.size(); ++s) {
        const Surface& surface = scene.surfaces[s];
        TissueProperties properties = tissue_properties(surface.tissue);
        if (!properties.elastic) {
            continue;
        }
        std::size_t cu = surface.cells_u;
        std::size_t cv = surface.cells_v;
        std::size_t n = cu * cv;
        mass.assign(n, 0.0);
        stiffness.assign(n * n, 0.0);
        for (std::size_t v = 0; v < cv; ++v) {
            for (std::size_t u = 0; u < cu; ++u) {
                std::size_t k = v * cu + u;
                const Patch& patch = scene.patches[surface.first_patch + k];
                double position = (static_cast<double>(u) + 0.5) / static_cast<double>(cu);
                double density = properties.surface_density * (1.0 + properties.density_gradient * (position - 0.5));
                mass[k] = density * patch.area;
                stiffness[k * n + k] += properties.foundation_stiffness * patch.area;
                double edge = std::sqrt(patch.area);
                // Strip coupling along u, edge coupling between strips along v,
                // and fixed outer edges (half-cell spring 2T).
                if (u + 1 < cu) {
                    std::size_t l = k + 1;
                    double spacing = norm(scene.patches[surface.first_patch + l].centroid - patch.centroid);
                    double coupling = properties.tension * edge / spacing;
                    stiffness[k * n + k] += coupling;
                    stiffness[l * n + l] += coupling;
                    stiffness[k * n + l] -= coupling;
                    stiffness[l * n + k] -= coupling;
                } else {
                    stiffness[k * n + k] += 2.0 * properties.tension;
                }
                if (u == 0) {
                    stiffness[k * n + k] += 2.0 * properties.tension;
                }
                if (v + 1 < cv) {
                    std::size_t l = k + cu;
                    double spacing = norm(scene.patches[surface.first_patch + l].centroid - patch.centroid);
                    double coupling = properties.tension * edge / spacing;
                    stiffness[k * n + k] += coupling;
                    stiffness[l * n + l] += coupling;
                    stiffness[k * n + l] -= coupling;
                    stiffness[l * n + k] -= coupling;
                } else {
                    stiffness[k * n + k] += 2.0 * properties.tension;
                }
                if (v == 0) {
                    stiffness[k * n + k] += 2.0 * properties.tension;
                }
            }
        }
        // Symmetric form M^-1/2 K M^-1/2.
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                stiffness[i * n + j] /= std::sqrt(mass[i] * mass[j]);
            }
        }
        jacobi_eigen(stiffness, n, values, vectors);
        for (std::size_t k = 0; k < n; ++k) {
            double omega = std::sqrt(std::max(values[k], 0.0));
            mode_first_patch_.push_back(surface.first_patch);
            mode_patch_count_.push_back(n);
            mode_shape_offset_.push_back(mode_shape_.size());
            for (std::size_t l = 0; l < n; ++l) {
                double shape = vectors[l * n + k] / std::sqrt(mass[l]);
                mode_shape_.push_back(shape);
            }
            mode_omega_.push_back(omega);
            mode_damping_.push_back(2.0 * properties.damping_ratio * omega);
            statistics_.highest_membrane_hz = std::max(statistics_.highest_membrane_hz, omega / (2.0 * ducklips::pi));
        }
    }
}

void BoundaryTransport::reset() noexcept {
    std::fill(pressure_ring_.begin(), pressure_ring_.end(), 0.0);
    std::fill(acceleration_ring_.begin(), acceleration_ring_.end(), 0.0);
    std::fill(mode_position_.begin(), mode_position_.end(), 0.0);
    std::fill(mode_velocity_.begin(), mode_velocity_.end(), 0.0);
    std::fill(membrane_acceleration_.begin(), membrane_acceleration_.end(), 0.0);
    step_index_ = 0;
}

void BoundaryTransport::inject_acceleration_at(std::size_t patch) {
    std::fill(forcing_sign_.begin(), forcing_sign_.end(), 0.0);
    std::fill(injection_sign_.begin(), injection_sign_.end(), 0.0);
    if (patch < patch_count_) {
        injection_sign_[patch] = 1.0;
    }
    select_active_acceleration_taps();
}

double BoundaryTransport::step(double source) noexcept {
    const std::uint64_t n = step_index_;
    const std::size_t count = patch_count_;
    const std::uint64_t mask = ring_mask_;
    const std::size_t write_base = static_cast<std::size_t>(n & mask) * count;
    const std::size_t previous_base = static_cast<std::size_t>((n - 1) & mask) * count;
    double new_weight = 1.0;
    double previous_weight = 0.0;
    if (temporal_averaging_) {
        new_weight = 0.5;
        previous_weight = 0.5;
    }

    // 1. Wall pressures from delayed exchanges (every delay is >= 1 sample,
    //    so the slot being written is never read in this pass).
    for (std::size_t j = 0; j < count; ++j) {
        double sum = 0.0;
        std::size_t end = pressure_begin_[j + 1];
        for (std::size_t t = pressure_begin_[j]; t < end; ++t) {
            const Tap& tap = pressure_taps_[t];
            std::size_t slot = static_cast<std::size_t>((n - tap.delay) & mask);
            sum += tap.weight * pressure_ring_[slot * count + tap.source];
        }
        std::size_t acceleration_end = acceleration_begin_[j];
        if (acceleration_active_) {
            acceleration_end = acceleration_begin_[j + 1];
        }
        for (std::size_t t = acceleration_begin_[j]; t < acceleration_end; ++t) {
            const Tap& tap = acceleration_taps_[t];
            std::size_t slot = static_cast<std::size_t>((n - tap.delay) & mask);
            sum += tap.weight * acceleration_ring_[slot * count + tap.source];
        }
        double fresh = sum + forcing_sign_[j] * source;
        pressure_ring_[write_base + j] = new_weight * fresh + previous_weight * pressure_ring_[previous_base + j];
    }

    // 2. Membranes: the new pressures pluck each surface's modes.
    std::fill(membrane_acceleration_.begin(), membrane_acceleration_.end(), 0.0);
    const double dt = time_step_;
    for (std::size_t k = 0; k < mode_omega_.size(); ++k) {
        std::size_t first = mode_first_patch_[k];
        std::size_t local_count = mode_patch_count_[k];
        std::size_t offset = mode_shape_offset_[k];
        double force = 0.0;
        for (std::size_t l = 0; l < local_count; ++l) {
            force += mode_shape_[offset + l] * area_[first + l] * pressure_ring_[write_base + first + l];
        }
        double omega = mode_omega_[k];
        double acceleration = force - mode_damping_[k] * mode_velocity_[k] - omega * omega * mode_position_[k];
        mode_velocity_[k] += dt * acceleration;
        mode_position_[k] += dt * mode_velocity_[k];
        for (std::size_t l = 0; l < local_count; ++l) {
            membrane_acceleration_[first + l] += mode_shape_[offset + l] * acceleration;
        }
    }
    for (std::size_t j = 0; j < count; ++j) {
        acceleration_ring_[write_base + j] = membrane_acceleration_[j] + injection_sign_[j] * source;
    }

    // 3. Microphone: another receiver, outside, factor 1.
    double microphone = 0.0;
    for (std::size_t t = 0; t < microphone_pressure_taps_.size(); ++t) {
        const Tap& tap = microphone_pressure_taps_[t];
        std::size_t slot = static_cast<std::size_t>((n - tap.delay) & mask);
        microphone += tap.weight * pressure_ring_[slot * count + tap.source];
    }
    std::size_t microphone_acceleration_count = 0;
    if (acceleration_active_) {
        microphone_acceleration_count = microphone_acceleration_taps_.size();
    }
    for (std::size_t t = 0; t < microphone_acceleration_count; ++t) {
        const Tap& tap = microphone_acceleration_taps_[t];
        std::size_t slot = static_cast<std::size_t>((n - tap.delay) & mask);
        microphone += tap.weight * acceleration_ring_[slot * count + tap.source];
    }
    step_index_ = n + 1;
    return microphone;
}

double BoundaryTransport::patch_pressure(std::size_t patch) const noexcept {
    std::size_t slot = static_cast<std::size_t>((step_index_ - 1) & ring_mask_);
    double result = pressure_ring_[slot * patch_count_ + patch];
    return result;
}

double BoundaryTransport::internal_rate() const noexcept {
    return rate_;
}

const CompileStatistics& BoundaryTransport::statistics() const noexcept {
    return statistics_;
}

}
