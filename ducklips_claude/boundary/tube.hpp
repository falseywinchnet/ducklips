#pragma once
// 1-D control: the same geometry as a plane-wave tube (staggered pressure /
// volume-velocity finite differences), closed at the constriction, radiating
// through a baffled-piston R-L load (Flanagan) at the lips. This is the air-
// column view DuckLips uses, without its other machinery, so the boundary
// model can be compared against the same area function and the same source.
#include "geometry.hpp"
#include "ducklips/physics.hpp"
#include <cstddef>
#include <vector>

namespace boundary {

class TubeControl final {
public:
    // internal_rate is the caller's step rate; the tube substeps internally
    // to keep its Courant number below 0.5. Allocates.
    void compile(const Scene& scene, double internal_rate, const ducklips::AirProperties& air);
    void reset() noexcept;
    // source: series pressure (Pa) at the source face, or the piston
    // acceleration (m/s^2) at the closed end for wall_acceleration scenes.
    // Returns far-field pressure at the scene's microphone distance.
    double step(double source) noexcept;

private:
    ducklips::AirProperties air_{};
    SourceKind kind_{SourceKind::pressure_dipole};
    std::size_t substeps_{1};
    double dt_{0.0};
    double cell_length_{0.0};
    double microphone_distance_{0.15};
    std::size_t source_face_{0};
    std::vector<double> area_{};
    std::vector<double> pressure_{};      // per cell
    std::vector<double> flow_{};          // per face, face 0 = closed end
    std::vector<double> face_inertance_{};
    double radiation_resistance_{0.0};
    double radiation_inertance_{0.0};
    double radiation_flow_{0.0};          // flow through the radiation inertance
    double piston_velocity_{0.0};
    double damping_{250.0};               // 1/s on face flows, a uniform loss stand-in
};

}
