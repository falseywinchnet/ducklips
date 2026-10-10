#pragma once
#include "ducklips/graph_continuous_glottis.hpp"
namespace ducklips {
struct JointAdvanceInputs final {
    JointGraphGeometry geometry{};
    double reservoir_pressure{},outlet_inertance{-1.};
    std::array<double,4> external_force{};
    std::span<const double> nodal_flow{},series_pressure{};
};
// Return actual old/mid/end geometry and volume-displacement flow for these
// absolute times. Views must remain valid until the next provider call. A
// provider can be called for rejected trials; it must not commit body state.
using JointGeometryProvider=bool (*)(void*,double,double,JointAdvanceInputs&) noexcept;
struct JointAccuracyOptions final {
    // maximum_step bounds a trial; accepted leaves are its two half steps.
    double tolerance{1e-5},maximum_step{1e-6},minimum_step{1e-10};
    std::size_t maximum_trials{4096};
};
struct JointAccuracyStep final {
    double energy_change{},external_work{},physical_loss{},numerical_loss{},energy_residual{};
    double source_port_work{},graph_port_work{},reservoir_volume{},maximum_accepted_error{},minimum_step{};
    std::size_t trials{},rejected{},accepted_steps{};
};
// Step doubling of the entire coupled source/tract system. Accept the two
// half steps, without extrapolation; their physical work receipts remain the
// receipts of the accepted path. Tolerance is local, not a global error bound.
// Topology and loss ownership must remain those prepared by the model. This
// does not add a source-contact transition rule.
enum class JointAccuracyFailure {none,configuration,trial_limit,step_limit};
class JointAccuracy final {
public:
    void prepare(const GraphContinuousGlottis&,std::size_t maximum_trials=4096);
    [[nodiscard]] bool advance(GraphContinuousGlottis&,double first_time,double last_time,
        JointGeometryProvider,void*,const JointAccuracyOptions&) noexcept;
    // Mean pressure and flow over the entire successful request. Mean(P)*
    // mean(Q)*duration is generally NOT its port work; use the work receipt.
    std::vector<double> average_pressure,average_flow;
    [[nodiscard]] std::span<const double> accepted_end_times() const noexcept {
        return {end_times_.data(),last.accepted_steps};
    }
    JointAccuracyStep last{};
    // Attempt diagnostics may change on failure; accepted receipts do not.
    JointAccuracyFailure failure{JointAccuracyFailure::none};
    double attempt_time{},attempt_step{},attempt_error{};
private:
    [[nodiscard]] bool compatible(const GraphContinuousGlottis&) const noexcept;
    [[nodiscard]] double difference(const GraphContinuousGlottis&,const GraphContinuousGlottis&,
        const JointGraphGeometry&) const noexcept;
    void accumulate(const GraphContinuousGlottis&,double,JointAccuracyStep&,
        std::vector<double>&,std::vector<double>&) noexcept;
    GraphContinuousGlottis accepted_,full_,halves_;
    std::vector<double> end_compliance_,end_conductance_,pressure_sum_,flow_sum_,trial_pressure_,trial_flow_,end_times_,trial_times_,pending_ends_;
    std::size_t maximum_trials_{};
};
}
