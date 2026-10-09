#pragma once
#include "ducklips/acoustic_network.hpp"
#include "ducklips/articulation.hpp"
#include "ducklips/modal.hpp"
#include "ducklips/phonetics.hpp"
#include "ducklips/tree_operator.hpp"
#include "ducklips/voice.hpp"
#include <array>
#include <cstddef>
#include <vector>

namespace ducklips {

inline constexpr std::size_t fitted_formant_count = 3;
using FormantArray = std::array<double, fitted_formant_count>;
using ShapeJacobian = std::array<std::array<double, shape_parameter_count>, fitted_formant_count>;

struct FormantFit final {
    Phone phone{Phone::AX};
    FormantArray target{};   // Hz
    FormantArray achieved{}; // Hz, lossless eigenfrequencies of the fitted posture
    std::size_t iterations{};
    double largest_relative_error{};
};

// Finds a tongue posture (body position, jaw) whose first three oral
// eigenfrequencies match target formants. Levenberg-Marquardt on log
// frequency with Tikhonov regularization toward the neutral posture, inside
// an anatomical box. The Jacobian chains the eigenvalue sensitivity
//   A_i d(lambda_k)/dA_i = sum_branches (A dg/dA) (dphi_k)^2 - lambda_k C_i phi_k(i)^2,
// the kinetic-minus-potential energy density of the mode at section i.
class FormantInverter final {
public:
    void prepare(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections,
                 const AirProperties& air);
    // Oral eigenfrequencies of a posture. Returns false if fewer than three
    // oral modes lie below 6 kHz.
    bool formants(const VowelShape& shape, FormantArray& output);
    // Analytic d(ln F_k)/d(coefficient_m) at the posture.
    bool formant_jacobian(const VowelShape& shape, FormantArray& frequencies, ShapeJacobian& jacobian);
    FormantFit fit(Phone phone, const FormantArray& target, VowelShape& shape);
    [[nodiscard]] const TractModel& model() const noexcept;

private:
    bool solve_posture(const VowelShape& shape, std::array<std::size_t, fitted_formant_count>& oral_modes);

    AirProperties air_{};
    TractModel model_{};
    TractShape tract_{};
    AcousticNetwork network_{};
    TreeOperator matrix_{};
    ModalSolver solver_{};
    ModalBasis basis_{};
    ModalBasis empty_{};
    std::vector<double> area_gradient_{};
};

// Fitted postures for every vowel of one voice.
class VowelLibrary final {
public:
    void build(const VoiceProfile& voice, std::size_t oral_sections, std::size_t nasal_sections,
               const AirProperties& air);
    [[nodiscard]] const VowelShape& shape(Phone phone) const noexcept;
    [[nodiscard]] const std::vector<FormantFit>& fits() const noexcept;

private:
    std::array<VowelShape, phone_count> shapes_{};
    std::vector<FormantFit> fits_{};
};

}
