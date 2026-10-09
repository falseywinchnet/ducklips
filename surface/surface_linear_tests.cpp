#include "surface_network.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
using Complex = std::complex<double>;
void require(bool condition, const char* message) {
    if (!condition) { throw std::runtime_error(message); }
}
bool rejects(const ducklips::SurfaceNetwork& network) {
    try { (void)network.linear_snapshot(); }
    catch (const std::logic_error&) { return true; }
    return false;
}
Complex filter(const std::array<double, 3>& b, const std::array<double, 3>& a, Complex q) {
    return (b[0] + q*b[1] + q*q*b[2])/(a[0] + q*a[1] + q*q*a[2]);
}
void check_passive(const ducklips::SurfaceLinearSnapshot& snapshot) {
    require(snapshot.paths.size() % 2 == 0, "paths must be reciprocal pairs");
    const std::array<double, 6> frequencies{48.0, 125.0, 500.0, 2000.0, 6000.0, 12000.0};
    for (double frequency : frequencies) {
        Complex q = std::exp(Complex(0.0, -2.0*3.14159265358979323846*frequency/snapshot.sample_rate));
        for (std::size_t index = 0; index < snapshot.paths.size(); index += 2) {
            const ducklips::SurfaceLinearPath& a = snapshot.paths[index];
            const ducklips::SurfaceLinearPath& b = snapshot.paths[index+1];
            require(a.from == b.to && a.to == b.from && a.admittance == b.admittance,
                    "path reciprocity lost");
            if (a.admittance == 0.0) { continue; }
            Complex h1 = std::pow(q, a.delay)*filter(a.numerator, a.denominator, q);
            Complex h2 = std::pow(q, b.delay)*filter(b.numerator, b.denominator, q);
            Complex denominator = 1.0-h1*h2;
            Complex diagonal = (1.0+h1*h2)/denominator;
            Complex ab = -2.0*h2/denominator;
            Complex ba = -2.0*h1/denominator;
            double minimum = diagonal.real()-std::abs((ab+std::conj(ba))/2.0);
            require(minimum >= -1e-10*std::max(1.0,std::abs(diagonal)), "nonpassive path");
        }
        for (const ducklips::SurfaceLinearPort& port : snapshot.ports) {
            Complex reflection = filter(port.numerator, port.denominator, q);
            Complex impedance = (1.0-reflection)/(1.0+reflection);
            require(impedance.real() >= -1e-10*std::max(1.0,std::abs(impedance)), "nonpassive wall/radiation port");
        }
    }
}
}
int main() {
    try {
        ducklips::SurfaceNetwork network{};
        require(rejects(network), "unprepared snapshot accepted");
        std::vector<double> areas(30, 3e-4);
        ducklips::CavityGeometry geometry = ducklips::single_tube(areas, 0.175, 8);
        geometry.tubes[0].solid = true;
        geometry.tubes[0].bent = true;
        geometry.tubes[0].cross_paths = true;
        geometry.tubes[0].axial_share = 0.6;
        network.build(geometry, 192000.0, ducklips::AirConstants{});
        check_passive(network.linear_snapshot());
        areas[15] *= 0.5;
        network.reshape(areas, 16);
        require(rejects(network), "moving snapshot accepted");
        for (std::size_t sample = 0; sample < 16; ++sample) {
            network.gather(); (void)network.scatter(0.0,0.0,0.0);
        }
        require(rejects(network), "unsettled smoothing snapshot accepted");
        network.reshape(areas, 1);
        network.reset();
        check_passive(network.linear_snapshot());
        for (std::size_t sample = 0; sample < 256; ++sample) {
            network.gather(); (void)network.scatter(sample == 0 ? 1e-6 : 0.0,0.0,0.0);
        }
        require(rejects(network), "nonlinear flow state accepted as zero-flow linearization");
        network.reset();
        check_passive(network.linear_snapshot());
        geometry.tubes[0].edge_tension = 0.1;
        network.build(geometry, 192000.0, ducklips::AirConstants{});
        require(rejects(network), "unsupported membrane coupling silently omitted");
        std::cout << "Frozen snapshot guards, reciprocity, and passive port checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
