#include "ducklips/cavity_network.hpp"
#include "ducklips/sources.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <new>
namespace {
std::size_t failures{};
std::atomic<std::size_t> allocations{};
void check(bool okay,const char* label,double value=0) {
    std::printf("%s %s %.12g\n",okay?"PASS":"FAIL",label,value);if (!okay) { ++failures; }
}
}
void* operator new(std::size_t n) { ++allocations;void* p=std::malloc(std::max<std::size_t>(n,1));if (!p) { throw std::bad_alloc(); }return p; }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
int main() {
    const ducklips::AirProperties air{};
    ducklips::CavityNetwork sealed{};sealed.prepare(1,{});
    std::array<double,1> v{1e-4},c{1e-8},p{120},in{0};
    check(sealed.reset(v,c,p),"sealed cavity initialized");v[0]-=1e-5;
    ducklips::CavityStepReport report=sealed.advance(v,c,{},in,.001,air);
    check(report.converged && std::abs(sealed.pressures()[0]-1120)<1e-10,"sealed moving volume follows exact mass balance",sealed.pressures()[0]);
    c[0]*=2;report=sealed.advance(v,c,{},in,.001,air);
    check(report.converged && std::abs(sealed.pressures()[0]-560)<1e-10,"changing compliance preserves stored charge",sealed.pressures()[0]);
    const std::array<ducklips::CavityConnection,3> links{{{0,ducklips::cavity_reservoir},{0,1},{1,ducklips::cavity_reservoir}}};
    std::array<ducklips::CavityPassage,3> ports{{{2.5e-5,1e5,600},{0,0,0},{6e-6,1e5,0}}};
    std::array<double,2> volumes{5e-5,3e-5},caps{5e-9,5e-9},pressure{400,0},inflows{0,0};
    ducklips::CavityNetwork moving{},fixed{};moving.prepare(2,links);fixed.prepare(2,links);
    check(moving.reset(volumes,caps,pressure) && fixed.reset(volumes,caps,pressure),"two sides of stop initialized");
    const std::array<double,2> baseline=volumes;
    bool valid=true;double isolated_error=0,mass_error=0,bridge_flow=0;
    const std::size_t before=allocations.load();
    for (std::size_t i=0;i<100;++i) {
        const std::array<double,2> old_volume=volumes;
        const std::array<double,2> old_pressure{moving.pressures()[0],moving.pressures()[1]};
        volumes[1]-=1e-7;const double dt=.0001;
        report=moving.advance(volumes,caps,ports,inflows,dt,air);
        valid=valid&&report.converged&&fixed.advance(baseline,caps,ports,inflows,dt,air).converged;
        isolated_error=std::max(isolated_error,std::abs(moving.pressures()[0]-fixed.pressures()[0]));
        bridge_flow=std::max(bridge_flow,std::abs(moving.flows()[1]));
        double stored=0;for (std::size_t k=0;k<2;++k) { stored+=caps[k]*(moving.pressures()[k]-old_pressure[k])+volumes[k]-old_volume[k]; }
        mass_error=std::max(mass_error,std::abs(stored+dt*(moving.flows()[0]+moving.flows()[2])));
    }
    check(allocations.load()==before,"pressure stepping performs no allocation");
    check(valid && isolated_error<1e-7 && bridge_flow==0,"front compression cannot cross sealed stop beyond solver tolerance",isolated_error);
    check(moving.pressures()[0]<=600 && moving.pressures()[1]>fixed.pressures()[1],"front compression raises front pressure only",moving.pressures()[1]);
    check(mass_error<1e-16,"all moving-cavity storage and reservoir flow balances",mass_error);
    const std::array<double,2> saved_pressure{moving.pressures()[0],moving.pressures()[1]};
    ports[1].area=std::numeric_limits<double>::quiet_NaN();
    report=moving.advance(volumes,caps,ports,inflows,.001,air);
    check(!report.converged && moving.pressures()[0]==saved_pressure[0] && moving.pressures()[1]==saved_pressure[1],"invalid port leaves pressure state unchanged");
    // At steady flow, independently sum series quadratic and linear losses.
    ports={{{2.5e-5,1e5,600},{8e-6,2e5,0},{6e-6,1e5,0}}};caps={1e-10,1e-10};pressure={0,0};
    check(moving.reset(volumes,caps,pressure),"open series network initialized");valid=true;
    for (std::size_t i=0;i<200;++i) { valid=valid&&moving.advance(volumes,caps,ports,inflows,.001,air).converged; }
    double inverse_area_squared=0,resistance=0;
    for (const ducklips::CavityPassage& port:ports) { inverse_area_squared+=1/(port.area*port.area);resistance+=port.resistance; }
    const double expected=ducklips::orifice_flow(1/std::sqrt(inverse_area_squared),600,resistance,air);
    const double flow_error=std::max({std::abs(moving.flows()[0]+expected),std::abs(moving.flows()[1]-expected),std::abs(moving.flows()[2]-expected)});
    check(valid && flow_error<1e-12,"series network agrees with independent steady pressure-drop solution",flow_error);
    ports[0].reservoir_pressure=-600;check(moving.reset(volumes,caps,pressure),"reverse-pressure network initialized");
    valid=true;for (std::size_t i=0;i<200;++i) { valid=valid&&moving.advance(volumes,caps,ports,inflows,.001,air).converged; }
    check(valid && std::abs(moving.flows()[1]+expected)<1e-12 && moving.pressures()[0]<0,
          "reverse flow obeys same passive pressure law",moving.flows()[1]);
    // Sealed external ports, one open internal connection: dissipate energy.
    ports[0].area=ports[2].area=0;pressure={700,-100};
    check(moving.reset(volumes,caps,pressure),"isolated equalization initialized");
    double previous_energy=.5*(caps[0]*pressure[0]*pressure[0]+caps[1]*pressure[1]*pressure[1]);valid=true;
    for (std::size_t i=0;i<100;++i) {
        report=moving.advance(volumes,caps,ports,inflows,.001,air);valid=valid&&report.converged;
        const double energy=.5*(caps[0]*moving.pressures()[0]*moving.pressures()[0]+caps[1]*moving.pressures()[1]*moving.pressures()[1]);
        valid=valid&&energy<=previous_energy+1e-16;previous_energy=energy;
    }
    check(valid && std::abs(moving.pressures()[0]-300)<1e-7 && std::abs(moving.pressures()[1]-300)<1e-7,
          "passive equalization preserves charge and dissipates energy",previous_energy);
    ducklips::CavityNetwork tree{};
    const std::array<ducklips::CavityConnection,2> fork{{{0,2},{1,2}}};tree.prepare(3,fork);
    const std::array<double,3> tv{1e-5,1e-5,1e-5},tc{1e-10,1e-10,1e-10},tp{400,-200,100},ti{0,0,0};
    const std::array<ducklips::CavityPassage,2> te{{{1e-5,1e6,0},{1e-5,1e6,0}}};
    valid=tree.reset(tv,tc,tp);for (std::size_t i=0;i<200;++i) { valid=valid&&tree.advance(tv,tc,te,ti,.001,air).converged; }
    double tree_error=0;for (double value:tree.pressures()) { tree_error=std::max(tree_error,std::abs(value-100)); }
    check(valid && tree_error<1e-7,"branched cavity tree reaches charge-weighted equilibrium",tree_error);
    return failures?1:0;
}
