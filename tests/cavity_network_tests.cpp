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
    ducklips::CavityNetwork sensitive{};
    const std::array<ducklips::CavityConnection,1> connection{{{0,1}}};sensitive.prepare(2,connection);
    const std::array<double,2> sv{1e-5,1e-5},sc{1e-8,1e-8},sp{273,273},si{1e-8,-1e-8};
    const std::array<ducklips::CavityPassage,1> passage{{{1e-4,.001,0,0,0}}};
    sensitive.reset(sv,sc,sp);report=sensitive.advance(sv,sc,passage,si,1./48000,air);
    check(report.converged && report.maximum_mass_residual<1e-10 && std::abs(sensitive.flows()[0]-1e-8)<1e-10,
          "high conductance converges within pressure representation error",report.maximum_mass_residual);
    check(std::abs(sensitive.pressures()[0]+sensitive.pressures()[1]-546)<1e-12,
          "roundoff-aware solve preserves total isolated charge",sensitive.pressures()[0]+sensitive.pressures()[1]-546);
    // Pressure reversals and reopening mix linear and directional quadratic
    // passages. With fixed volume and no nodal source, pressures must remain
    // within the initial/reservoir bounds, regardless of the Newton path.
    constexpr std::size_t count=12;ducklips::CavityNetwork switching{};
    std::vector<ducklips::CavityConnection> slinks{{0,ducklips::cavity_reservoir}};
    std::vector<ducklips::CavityPassage> sports{{2.5e-5,1e5,600}};
    for(std::size_t i=0;i+1<count;++i){slinks.push_back({i,i+1});sports.push_back({1e-4,i%2?.001:1e5,0,i%2?0.:1.,i%2?1.:0.});}
    slinks.push_back({count-1,ducklips::cavity_reservoir});sports.push_back({1e-4,1e5,0});
    std::vector<double> vv(count,1e-5),cc(count),pp(count),ii(count,0);
    for(std::size_t i=0;i<count;++i){cc[i]=double(i+1)*1e-10;pp[i]=i%2?600:0;}
    switching.prepare(count,slinks);valid=switching.reset(vv,cc,pp);mass_error=0;double excess=0;
    for(std::size_t step=0;step<100;++step){
        sports[0].reservoir_pressure=(step/10)%2?-600:600;sports[5].area=step%20>=10 && step%20<15?0:1e-4;
        double old_charge=0;for(std::size_t i=0;i<count;++i){old_charge+=cc[i]*switching.pressures()[i];}
        report=switching.advance(vv,cc,sports,ii,.001,air);valid=report.converged&&valid;
        double charge=0;for(std::size_t i=0;i<count;++i){charge+=cc[i]*switching.pressures()[i];excess=std::max(excess,std::abs(switching.pressures()[i])-600);}
        mass_error=std::max(mass_error,std::abs(charge-old_charge+.001*(switching.flows().front()+switching.flows().back())));
    }
    check(valid,"reversing and reopening passive network converges");
    check(excess<1e-7,"reversing and reopening passive network obeys pressure bounds",excess);
    check(mass_error<1e-13,"reversing and reopening passive network conserves mass",mass_error);
    // A four-cell ring has a genuine cycle: its antisymmetric mode has
    // Laplacian eigenvalue 2/R. Deleting either route changes this solution.
    ducklips::CavityNetwork ring{};
    const std::array<ducklips::CavityConnection,4> ring_links{{{0,1},{1,2},{2,3},{3,0}}};
    std::array<ducklips::CavityPassage,4> ring_ports{};
    for(ducklips::CavityPassage& port:ring_ports){port={1e-4,1e6,0,0,0};}
    const std::array<double,4> rv{1e-5,1e-5,1e-5,1e-5},rc{1e-10,1e-10,1e-10,1e-10},rp{723,123,-477,123},ri{};
    ring.prepare(4,ring_links);valid=ring.reset(rv,rc,rp);
    const std::size_t ring_allocations=allocations.load();
    report=ring.advance(rv,rc,ring_ports,ri,.001,air);
    const double amplitude=600/(1+.001*2/(1e6*1e-10));double ring_error=0;
    const std::array<double,4> ring_expected{123+amplitude,123,123-amplitude,123};
    for(std::size_t i=0;i<4;++i){ring_error=std::max(ring_error,std::abs(ring.pressures()[i]-ring_expected[i]));}
    check(valid && report.converged && ring_error<1e-8,"cyclic pressure retains both paths and matches analytical mode",ring_error);
    check(allocations.load()==ring_allocations,"cyclic pressure stepping performs no allocation");
    // Close both links at node zero: compressing that chamber cannot affect
    // the still-connected remaining three. No topology rebuilding is needed.
    ring_ports[0].area=ring_ports[3].area=0;
    const std::array<double,4> uniform{100,100,100,100};std::array<double,4> compressed=rv;
    valid=ring.reset(rv,rc,uniform);compressed[0]-=1e-8;
    report=ring.advance(compressed,rc,ring_ports,ri,.001,air);double other_error=0;
    for(std::size_t i=1;i<4;++i){other_error=std::max(other_error,std::abs(ring.pressures()[i]-100));}
    check(valid && report.converged && std::abs(ring.pressures()[0]-200)<1e-8 && other_error<1e-8
        && ring.flows()[0]==0 && ring.flows()[3]==0,"closed cyclic contacts isolate compression exactly",other_error);
    // Reopening dissipates energy without creating or destroying stored mass.
    ring_ports[0].area=ring_ports[3].area=1e-4;valid=true;mass_error=0;previous_energy=INFINITY;
    for(std::size_t step=0;step<40;++step){
        report=ring.advance(compressed,rc,ring_ports,ri,.001,air);valid=valid&&report.converged;
        double charge=0,energy=0;for(double value:ring.pressures()){charge+=1e-10*value;energy+=.5e-10*value*value;}
        mass_error=std::max(mass_error,std::abs(charge-5e-8));valid=valid&&energy<=previous_energy+1e-16;previous_energy=energy;
    }
    check(valid && mass_error<1e-16,"reopened cycle dissipates energy and preserves isolated mass",mass_error);
    ducklips::CavityNetwork disjoint{};disjoint.prepare(3,{});
    std::array<double,3> disjoint_volume=tv;disjoint_volume[1]-=1e-8;
    valid=disjoint.reset(tv,tc,tp);report=disjoint.advance(disjoint_volume,tc,{},ti,.001,air);
    check(valid && report.converged && std::abs(disjoint.pressures()[0]-400)<1e-8 && std::abs(disjoint.pressures()[1]+100)<1e-8,
        "disconnected graph keeps independent pressure states");
    // Multiple physical passages between two cavities are parallel conductances.
    ducklips::CavityNetwork parallel{};
    const std::array<ducklips::CavityConnection,2> parallel_links{{{0,1},{0,1}}};
    const std::array<ducklips::CavityPassage,2> parallel_ports{{{1e-4,1e6,0,0,0},{1e-4,1e6,0,0,0}}};
    const std::array<double,2> pv{1e-5,1e-5},pc{1e-10,1e-10},ppar{600,0},pi{};
    parallel.prepare(2,parallel_links);valid=parallel.reset(pv,pc,ppar);report=parallel.advance(pv,pc,parallel_ports,pi,.001,air);
    const double parallel_expected=300+300/(1+4*.001/(1e6*1e-10));
    check(valid && report.converged && std::abs(parallel.pressures()[0]-parallel_expected)<1e-8,
        "parallel passages both participate in pressure equalization",parallel.pressures()[0]-parallel_expected);
    // A vanished finite-volume cell has zero storage, not an epsilon cavity.
    // Its old air must leave through an interval-open passage in this step.
    ducklips::CavityNetwork transport{};transport.prepare(2,connection);
    const double bulk=air.density*air.sound_speed*air.sound_speed;
    const std::array<double,2> oldv{1e-6,2e-6},oldc{oldv[0]/bulk,oldv[1]/bulk},oldp{500,100};
    std::array<double,2> newv{0,3e-6},newc{0,newv[1]/bulk};
    std::array<ducklips::CavityPassage,1> exit_port{{{1e-5,1e6,0,0,0}}};
    valid=transport.reset(oldv,oldc,oldp);const double oldmass=oldv[0]+oldv[1]+oldc[0]*oldp[0]+oldc[1]*oldp[1];
    const std::size_t transport_allocations=allocations.load();
    report=transport.advance(newv,newc,exit_port,pi,.001,air);
    const double expected_outflow=(oldv[0]+oldc[0]*oldp[0])/.001;
    check(valid && report.converged && std::abs(transport.flows()[0]-expected_outflow)<1e-12
        && std::abs(newv[1]+newc[1]*transport.pressures()[1]-oldmass)<1e-16,
        "disappearing cell expels all its air without artificial storage",transport.flows()[0]-expected_outflow);
    check(allocations.load()==transport_allocations,"vanishing-cell step allocates no storage");
    exit_port[0].area=0;const double retained=transport.pressures()[1];
    report=transport.advance(newv,newc,exit_port,pi,.001,air);
    check(report.converged && transport.pressures()[0]==0 && transport.pressures()[1]==retained && transport.flows()[0]==0,
        "absent sealed cell has no residual pressure or stored air");
    // Reopening a new cell fills it by flow, not by copying another pressure.
    exit_port[0].area=1e-5;newv=oldv;newc=oldc;
    const double charge=oldmass-newv[0]-newv[1],gdt=.001/exit_port[0].resistance;
    const double born_pressure=(gdt*charge/newc[1]-newv[0])/(newc[0]+gdt*(1+newc[0]/newc[1]));
    report=transport.advance(newv,newc,exit_port,pi,.001,air);
    check(report.converged && std::abs(transport.pressures()[0]-born_pressure)<1e-7
        && std::abs(newv[0]+newv[1]+newc[0]*transport.pressures()[0]+newc[1]*transport.pressures()[1]-oldmass)<1e-16,
        "appearing cell fills conservatively with the independently solved pressure",transport.pressures()[0]-born_pressure);
    ducklips::CavityNetwork impossible{};impossible.prepare(1,{});
    std::array<double,1> iv{1e-6},ic{iv[0]/bulk},ip{500},iz{};
    valid=impossible.reset(iv,ic,ip);report=impossible.advance(iz,iz,{},iz,.001,air);
    check(valid && !report.converged && impossible.pressures()[0]==500,
        "sealed disappearing air is rejected transactionally rather than deleted");
    report=impossible.advance(iv,ic,{},iz,.001,air);
    check(report.converged && impossible.pressures()[0]==500,"rejected disappearance leaves original stored mass intact");
    const std::array<ducklips::CavityConnection,3> tri_links{{{0,1},{1,2},{2,0}}};
    std::array<ducklips::CavityPassage,3> tri_ports{{{0,0,0,0,0},{1e-4,1e6,0,0,0},{0,0,0,0,0}}};
    ducklips::CavityNetwork chord{};chord.prepare(3,tri_links);
    std::array<double,3> cv{1e-5,1e-6,1e-5},ccap{1e-10,1e-11,1e-10},cp{};
    valid=chord.reset(cv,ccap,cp);cv[1]=ccap[1]=0;report=chord.advance(cv,ccap,tri_ports,ti,.001,air);
    check(valid && report.converged && std::abs(chord.flows()[1]-.001)<1e-12 && std::abs(chord.pressures()[2]-10000)<1e-6,
        "zero-storage cell can drain through an edge outside the preconditioner tree",chord.pressures()[2]);
    // A moving wedge has V=a*t^2 and interval conductance b*t. Its
    // prescribed geometry and opening admit a finite constant pressure.
    // This distinguishes consistent moving geometry from independently
    // interpolating cell volume and aperture across a birth/death event.
    double wedge_error=0,wedge_mass_error=0;valid=true;
    for(std::size_t subdivisions:{1u,2u,4u,8u,16u}){
        for(bool closing:{false,true}){
            ducklips::CavityNetwork wedge{};
            const std::array<ducklips::CavityConnection,1> vent{{{0,ducklips::cavity_reservoir}}};wedge.prepare(1,vent);
            const double acoef=1e-4,bcoef=1e-6,duration=.01;
            const double expected_pressure=closing?2*acoef/(bcoef-2*acoef/bulk):-2*acoef/(bcoef+2*acoef/bulk);
            std::array<double,1> wv{closing?acoef*duration*duration:0},wc{wv[0]/bulk},wp{closing?expected_pressure:0},wi{};
            const double initial_mass=wv[0]+wc[0]*wp[0];double integrated_outflow=0;
            valid=wedge.reset(wv,wc,wp)&&valid;
            for(std::size_t step=1;step<=subdivisions;++step){
                const double t0=duration*double(step-1)/double(subdivisions),t1=duration*double(step)/double(subdivisions);
                const double old_height=closing?duration-t0:t0,new_height=closing?duration-t1:t1;
                const double average_height=.5*(old_height+new_height);wv[0]=acoef*new_height*new_height;wc[0]=wv[0]/bulk;
                const std::array<ducklips::CavityPassage,1> opening{{{1e-5*average_height/duration,1/(bcoef*average_height),0,0,0}}};
                report=wedge.advance(wv,wc,opening,wi,t1-t0,air);valid=valid&&report.converged;
                wedge_error=std::max(wedge_error,std::abs(wedge.pressures()[0]-expected_pressure));
                integrated_outflow+=(t1-t0)*wedge.flows()[0];
                wedge_mass_error=std::max(wedge_mass_error,std::abs(wv[0]+wc[0]*wedge.pressures()[0]+integrated_outflow-initial_mass));
            }
        }
    }
    check(valid && wedge_error<1e-7,"coherent moving-wedge birth and disappearance retain finite exact pressure under timestep refinement",wedge_error);
    check(wedge_mass_error<1e-18,"moving wedge exchanges exactly its stored air with the reservoir",wedge_mass_error);
    // The same disappearing-cell conservation must hold through nonlinear
    // passage loss, where the drained amount fixes Q independently of Newton.
    valid=transport.reset(oldv,oldc,oldp);newv={0,3e-6};newc={0,newv[1]/bulk};exit_port[0]={1e-5,1e6,0,1,2};
    report=transport.advance(newv,newc,exit_port,pi,.001,air);
    const double receiving_pressure=(oldmass-newv[1])/newc[1];
    const double transport_pressure=receiving_pressure+1e6*expected_outflow+air.density*expected_outflow*expected_outflow/(2e-10);
    check(valid && report.converged && std::abs(transport.flows()[0]-expected_outflow)<1e-12
        && std::abs(transport.pressures()[0]-transport_pressure)<1e-5,
        "nonlinear disappearing-cell outflow matches independent mass and pressure-drop equations",transport.pressures()[0]-transport_pressure);
    return failures?1:0;
}
