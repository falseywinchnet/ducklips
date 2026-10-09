#include "ducklips/moving_tree.hpp"
#include "ducklips/moving_mode.hpp"
#include "ducklips/canonical_speech.hpp"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr double pi=3.14159265358979323846;
void require(bool okay,const char* message) { if (!okay) { throw std::runtime_error(message); } }
void read(std::ifstream& in,void* target,std::size_t bytes) {
    in.read(static_cast<char*>(target),static_cast<std::streamsize>(bytes));require(bool(in),"truncated fixture");
}
void fixture(const char* input,const char* output) {
    std::ifstream in(input,std::ios::binary);std::uint64_t header[3]{};double dt{};
    read(in,header,sizeof(header));read(in,&dt,sizeof(dt));
    require(header[0]>0 && header[0]<=4096 && header[1]<=16384 && header[2]<=10000000,"invalid fixture size");
    require(std::isfinite(dt) && dt>0,"invalid fixture period");
    const std::size_t n=header[0],m=header[1],frames=header[2];
    std::vector<ducklips::TreeEdge> edges(m);
    for (ducklips::TreeEdge& edge:edges) {
        std::uint64_t pair[2];read(in,pair,sizeof(pair));edge={pair[0],pair[1]};
    }
    ducklips::MovingTreeReference solver;solver.prepare(n,edges);
    read(in,solver.charge.data(),n*sizeof(double));read(in,solver.momentum.data(),m*sizeof(double));
    std::vector<double> c1(n),c2(n),g1(m),g2(m);
    std::ofstream out(output,std::ios::binary);require(bool(out),"cannot open fixture output");
    for (std::size_t frame=0;frame<frames;++frame) {
        read(in,c1.data(),n*sizeof(double));read(in,c2.data(),n*sizeof(double));
        read(in,g1.data(),m*sizeof(double));read(in,g2.data(),m*sizeof(double));
        for (double c:c1) { require(c>0 && std::isfinite(c),"invalid C1"); }
        for (double c:c2) { require(c>0 && std::isfinite(c),"invalid C2"); }
        for (double g:g1) { require(g>=0 && std::isfinite(g),"invalid G1"); }
        for (double g:g2) { require(g>=0 && std::isfinite(g),"invalid G2"); }
        solver.advance(c1,c2,g1,g2,dt);
        out.write(reinterpret_cast<const char*>(solver.charge.data()),static_cast<std::streamsize>(n*sizeof(double)));
        out.write(reinterpret_cast<const char*>(solver.momentum.data()),static_cast<std::streamsize>(m*sizeof(double)));
    }
    require(bool(out),"failed fixture output");
}
void benchmark(int substeps) {
    constexpr std::size_t n=96;
    const std::size_t frames=48000*5*static_cast<std::size_t>(substeps);
    std::vector<ducklips::TreeEdge> edges;
    for (std::size_t i=1;i<n;++i) { edges.push_back({(i-1)/2,i}); }
    for (std::size_t i=0;i<n;++i) { edges.push_back({i,ducklips::ground_node}); }
    ducklips::MovingTreeReference solver;solver.prepare(n,edges);
    std::vector<double> c1(n,5e-11),c2(n,5e-11),g1(edges.size(),.004),g2(edges.size(),.004);
    solver.charge[0]=1e-8;
    const std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    for (std::size_t frame=0;frame<frames;++frame) {
        // Bounded moving metrics. Include stage-metric filling and every
        // factorization; excludes articulation, sources, losses and microphone.
        const double position=1+.2*std::sin(2*pi*double(frame)/(4800*substeps));
        for (std::size_t i=0;i<n;++i) {
            c1[i]=5e-11*(1+double(i%7)*.1)*position;
            c2[i]=c1[i]*(1+1e-5*std::cos(double(frame)*.01/substeps));
        }
        for (std::size_t i=0;i<edges.size();++i) {
            g1[i]=.004*(1+double(i%5)*.1)/position;g2[i]=g1[i]*1.00001;
        }
        solver.advance(c1,c2,g1,g2,1.0/(48000*substeps));
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    double checksum=0;for (double q:solver.charge) { checksum+=q; }
    require(std::isfinite(checksum),"nonfinite benchmark");
    std::cout<<std::setprecision(17)<<"{\"nodes\":"<<n<<",\"branches\":"<<edges.size()
             <<",\"substeps\":"<<substeps<<",\"states\":"<<n+edges.size()<<",\"simulated_seconds\":5,\"wall_seconds\":"<<seconds
             <<",\"realtime_multiple\":"<<5/seconds<<",\"checksum\":"<<checksum<<"}\n";
}
void tests() {
    ducklips::MovingTreeReference solver;
    const ducklips::TreeEdge edges[]={{0,1},{0,1},{1,ducklips::ground_node},{0,ducklips::ground_node}};
    solver.prepare(2,edges);
    solver.charge={.3,.1};solver.momentum={.2,.1,.4,-.1};
    const std::vector<double> c{.7,1.1},g{1.2,.5,.8,.3};
    double initial=0;
    for (std::size_t i=0;i<c.size();++i) { initial+=solver.charge[i]*solver.charge[i]/c[i]/2; }
    for (std::size_t i=0;i<g.size();++i) { initial+=g[i]*solver.momentum[i]*solver.momentum[i]/2; }
    for (int i=0;i<10000;++i) { solver.advance(c,c,g,g,.01); }
    double energy=0;
    for (std::size_t i=0;i<c.size();++i) { energy+=solver.charge[i]*solver.charge[i]/c[i]/2; }
    for (std::size_t i=0;i<g.size();++i) { energy+=g[i]*solver.momentum[i]*solver.momentum[i]/2; }
    require(std::abs(energy/initial-1)<1e-10,"stationary energy drift");
    const ducklips::TreeEdge ground[]={{0,ducklips::ground_node}};
    const double omega=2*pi*1234,rate=std::log(4.0)/.1,duration=.002;
    ducklips::UniformAreaModeStep exact;exact.prepare(omega,0,rate,duration);
    const std::complex<double> expected=exact.apply({.3,.7});
    double previous_error=0;
    for (int steps:{32,64,128}) {
        solver.prepare(1,ground);solver.charge[0]=.3;solver.momentum[0]=.7/omega;
        const double dt=duration/steps;
        for (int i=0;i<steps;++i) {
            const double s1=std::exp(rate*(i+.5-std::sqrt(3.0)/6)*dt);
            const double s2=std::exp(rate*(i+.5+std::sqrt(3.0)/6)*dt);
            const double c1[]={s1},c2[]={s2},g1[]={omega*omega*s1},g2[]={omega*omega*s2};
            solver.advance(c1,c2,g1,g2,dt);
        }
        const double scale=std::exp(rate*duration/2);
        const std::complex<double> actual(solver.charge[0]/scale,solver.momentum[0]*omega*scale);
        const double error=std::abs(actual-expected);
        if (previous_error>0) { require(previous_error/error>15,"moving uniform limit is not fourth order"); }
        previous_error=error;
    }
    require(previous_error<4e-6,"uniform reference mismatch");
    // Affine pressure and volume-flow forcing against a closed-form oscillator.
    solver.prepare(1,ground);
    const double compliance[]={1},conductance[]={omega*omega},nodal[]={.2},series[]={.3/omega};
    for (int i=0;i<1024;++i) { solver.advance(compliance,compliance,conductance,conductance,duration/1024,nodal,series); }
    const std::complex<double> pole(0,omega),drive(.2,.3);
    const std::complex<double> forced=drive*(std::exp(pole*duration)-std::complex<double>(1,0))/pole;
    require(std::abs(std::complex<double>(solver.charge[0],omega*solver.momentum[0])-forced)<1e-12,"forced oscillator mismatch");
    // A held source must exchange work with the Gauss average, not the
    // endpoint trapezoid. Include acoustics far above the output Nyquist.
    for (double w:{1000.,25000.,200000.,400000.}) {
        solver.prepare(1,ground);const double metric[]={w*w},input[]={.4};
        solver.advance(compliance,compliance,metric,metric,1.0/48000,input);
        const double energy=(solver.charge[0]*solver.charge[0]+w*w*solver.momentum[0]*solver.momentum[0])/2;
        const double work=input[0]*solver.average_pressure[0]/48000;
        require(std::abs(energy-work)<1e-22 && work>=0,"source port does not conserve discrete work");
    }
    // Source work and passive decay including viscothermal/radiation memory.
    // The other power port is a series pressure driving branch volume flow.
    // Its average must account for the same transported volume as charge.
    solver.prepare(2,std::vector<ducklips::TreeEdge>{{0,1}});
    const double pair_c[]={2e-10,3e-10},pair_g[]={.02},pair_source[]={25};
    const double pair_dt=1./48000;
    solver.advance(pair_c,pair_c,pair_g,pair_g,pair_dt,{},pair_source);
    const double transported=pair_dt*solver.average_flow[0];
    require(std::abs(solver.charge[0]+transported)<1e-24 &&
            std::abs(solver.charge[1]-transported)<1e-24,"branch flow is not the transported volume");
    const double pair_energy=.5*(solver.charge[0]*solver.charge[0]/pair_c[0]+
        solver.charge[1]*solver.charge[1]/pair_c[1]+pair_g[0]*solver.momentum[0]*solver.momentum[0]);
    require(std::abs(pair_energy-pair_source[0]*transported)<1e-20,
            "series source port does not conserve discrete work");
    ducklips::AcousticNetwork network;network.prepare(4,2,1,{});
    ducklips::TractShape shape;shape.oral_area.assign(4,2e-4);shape.oral_length.assign(4,.04);
    shape.nasal_area.assign(2,1e-4);shape.nasal_length.assign(2,.04);shape.port_area=3e-5;
    const ducklips::AirProperties air{};network.assemble(shape,air,{});
    ducklips::CanonicalSpeech acoustic;acoustic.prepare(network,air);
    acoustic.free_step(network,network,0,0,1.0/48000,8,0,1,0,0);
    const double input=1e-5,source_work=input*input*acoustic.load_impedance()/48000;
    acoustic.inject(input);double previous_energy=acoustic.energy();
    require(previous_energy>0 && source_work>=previous_energy*(1-1e-12),"causal loss source port generates energy");
    const double starting_energy=previous_energy;
    for (int i=0;i<500;++i) {
        acoustic.free_step(network,network,0,0,1.0/48000,8,0,1,0,0);
        const double now=acoustic.energy();
        require(now<=previous_energy*(1+1e-11),"unforced causal loss memory is not passive");previous_energy=now;
    }
    require(previous_energy<starting_energy*.8,"causal memory did not dissipate energy");
    bool rejected=false;
    try { const ducklips::TreeEdge cyclic[]={{0,1},{1,2},{2,0}};solver.prepare(3,cyclic); }
    catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"cyclic topology accepted");
    std::cout<<"Moving tree stationary energy, parallel branches, fourth-order uniform reference and topology checks passed\n";
}
}
int main(int argc,char** argv) {
    try {
        if (argc==4 && std::string(argv[1])=="--fixture") { fixture(argv[2],argv[3]); }
        else if (argc==2 && std::string(argv[1])=="--benchmark") { benchmark(1); }
        else if (argc==3 && std::string(argv[1])=="--benchmark") {
            const int substeps=std::stoi(argv[2]);require(substeps>=1 && substeps<=16,"invalid substeps");benchmark(substeps);
        }
        else if (argc==1) { tests(); }
        else { throw std::invalid_argument("moving_tree_probe [--fixture input output | --benchmark]"); }
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
