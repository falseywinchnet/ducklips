#include "tube.hpp"
#include <algorithm>
#include <cmath>

namespace boundary {

void TubeControl::compile(const Scene& scene, double internal_rate, const ducklips::AirProperties& air) {
    air_ = air;
    kind_ = scene.source_kind;
    cell_length_ = scene.tube.cell_length;
    microphone_distance_ = scene.microphone_distance;
    source_face_ = scene.tube.source_face;
    area_ = scene.tube.area;
    const std::size_t cells = area_.size();
    double courant = air_.sound_speed / internal_rate / cell_length_;
    substeps_ = static_cast<std::size_t>(std::ceil(courant / 0.45));
    substeps_ = std::max<std::size_t>(substeps_, 1);
    dt_ = 1.0 / (internal_rate * static_cast<double>(substeps_));
    pressure_.assign(cells, 0.0);
    flow_.assign(cells + 1, 0.0);
    face_inertance_.assign(cells + 1, 0.0);
    const double rho = air_.density;
    for (std::size_t f = 1; f < cells; ++f) {
        face_inertance_[f] = rho * (0.5 * cell_length_ / area_[f - 1] + 0.5 * cell_length_ / area_[f]);
    }
    face_inertance_[cells] = rho * 0.5 * cell_length_ / area_[cells - 1];
    double open_area = area_[cells - 1];
    double radius = std::sqrt(open_area / ducklips::pi);
    radiation_resistance_ = 128.0 * rho * air_.sound_speed / (9.0 * ducklips::pi * ducklips::pi * open_area);
    radiation_inertance_ = 8.0 * rho * radius / (3.0 * ducklips::pi * open_area);
    reset();
}

void TubeControl::reset() noexcept {
    std::fill(pressure_.begin(), pressure_.end(), 0.0);
    std::fill(flow_.begin(), flow_.end(), 0.0);
    radiation_flow_ = 0.0;
    piston_velocity_ = 0.0;
}

double TubeControl::step(double source) noexcept {
    const std::size_t cells = pressure_.size();
    const double rho_c2 = air_.density * air_.sound_speed * air_.sound_speed;
    const double dt = dt_;
    const double decay = 1.0 - dt * damping_;
    double radiated_change = 0.0;
    for (std::size_t sub = 0; sub < substeps_; ++sub) {
        // Faces: momentum.
        if (kind_ == SourceKind::wall_acceleration) {
            piston_velocity_ += dt * source;
            flow_[0] = -area_[0] * piston_velocity_;
        }
        for (std::size_t f = 1; f < cells; ++f) {
            double drive = pressure_[f - 1] - pressure_[f];
            if (kind_ == SourceKind::pressure_dipole && f == source_face_) {
                drive += source;
            }
            flow_[f] = decay * flow_[f] + dt / face_inertance_[f] * drive;
        }
        // Open end: implicit in the radiation resistance, parallel R-L load.
        double k = dt / face_inertance_[cells];
        double old_flow = flow_[cells];
        double new_flow = (old_flow + k * (pressure_[cells - 1] + radiation_resistance_ * radiation_flow_)) /
                          (1.0 + k * radiation_resistance_);
        double load_pressure = radiation_resistance_ * (new_flow - radiation_flow_);
        radiation_flow_ += dt * load_pressure / radiation_inertance_;
        flow_[cells] = new_flow;
        radiated_change += new_flow - old_flow;
        // Cells: continuity.
        for (std::size_t i = 0; i < cells; ++i) {
            double compliance_inverse = rho_c2 / (area_[i] * cell_length_);
            pressure_[i] -= dt * compliance_inverse * (flow_[i + 1] - flow_[i]);
        }
    }
    double flow_derivative = radiated_change / (dt * static_cast<double>(substeps_));
    double microphone = air_.density * flow_derivative / (4.0 * ducklips::pi * microphone_distance_);
    return microphone;
}

}
