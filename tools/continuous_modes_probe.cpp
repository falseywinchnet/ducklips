#include "ducklips/continuous_modes.hpp"
#include "ducklips/moving_mode.hpp"
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <new>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace { std::atomic<std::size_t> allocations{0}; }
void* operator new(std::size_t count) {
    ++allocations;
    void* result=std::malloc(std::max<std::size_t>(count,1));
    if (!result) { throw std::bad_alloc(); }
    return result;
}
void* operator new[](std::size_t count) { return ::operator new(count); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }

namespace {
bool use_gauss=false;
void advance(ducklips::ContinuousModesReference& solver,const ducklips::ModalMotionSample& a,const ducklips::ModalMotionSample& b,const ducklips::ModalMotionSample& c,double dt) {
    if (use_gauss) { solver.advance_gauss(a,b,c,dt); } else { solver.advance(a,b,c,dt); }
}
constexpr double pi=3.14159265358979323846;
void require(bool okay,const char* message) { if (!okay) { throw std::runtime_error(message); } }
struct Motion final {
    std::vector<double> omega,pressure,flow,diagonal;
    explicit Motion(std::size_t n):omega(n),pressure(n*n),flow(n*n),diagonal(n) {}
    ducklips::ModalMotionSample view() const { return {omega,pressure,flow}; }
    ducklips::ModalMotionSample canonical_view() const { return {omega,pressure,{},diagonal}; }
};
void read(std::ifstream& in,void* target,std::size_t bytes) {
    in.read(static_cast<char*>(target),static_cast<std::streamsize>(bytes));require(bool(in),"truncated fixture");
}
void fixture(const char* input,const char* output,bool canonical=false) {
    std::ifstream in(input,std::ios::binary);std::uint64_t header[2];double dt;
    read(in,header,sizeof(header));read(in,&dt,sizeof(dt));
    require(header[0]>0 && header[0]<=256 && header[1]<=10000000 && std::isfinite(dt) && dt>0,"invalid fixture header");
    const std::size_t n=header[0],frames=header[1];
    ducklips::ContinuousModesReference solver;solver.prepare(n);
    read(in,solver.state.data(),2*n*sizeof(double));
    Motion first(n),middle(n),last(n);
    std::ofstream out(output,std::ios::binary);require(bool(out),"cannot open output");
    for (std::size_t frame=0;frame<frames;++frame) {
        for (Motion* motion:{&first,&middle,&last}) {
            read(in,motion->omega.data(),n*sizeof(double));
            read(in,motion->pressure.data(),n*n*sizeof(double));read(in,motion->flow.data(),n*n*sizeof(double));
            for (double w:motion->omega) { require(std::isfinite(w) && (canonical ? w>0 : w>=0),"invalid frequency"); }
            for (double v:motion->pressure) { require(std::isfinite(v),"invalid pressure transport"); }
            for (double v:motion->flow) { require(std::isfinite(v),"invalid flow transport"); }
        }
        if (canonical) {
            for (Motion* motion:{&first,&middle,&last}) {
                for (std::size_t i=0;i<n;++i) { motion->diagonal[i]=motion->flow[i*n+i]+motion->pressure[i*n+i]; }
            }
            advance(solver,first.canonical_view(),middle.canonical_view(),last.canonical_view(),dt);
        } else { advance(solver,first.view(),middle.view(),last.view(),dt); }
        out.write(reinterpret_cast<const char*>(solver.state.data()),static_cast<std::streamsize>(2*n*sizeof(double)));
    }
    require(bool(out),"failed output");
}
void dense_rotation_test(std::size_t n) {
    // Dense orthogonal change of coordinates applied to independently known
    // pair rotations. All acoustic frequencies coincide, so acoustic phase
    // and spatial rotation commute. Exercises SIMD widths and odd tails.
    Motion motion(n);std::vector<double> h(n*n),q(n*n),v(n),expected(2*n),paired(2*n);
    double length=0;for (std::size_t i=0;i<n;++i) { v[i]=double(i+1);length+=v[i]*v[i]; }
    for (double& value:v) { value/=std::sqrt(length); }
    for (std::size_t i=0;i<n;++i) {
        motion.omega[i]=2*pi*1234;
        for (std::size_t j=0;j<n;++j) { h[i*n+j]=(i==j ? 1.0 : 0.0)-2*v[i]*v[j]; }
    }
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<n;++j) {
            for (std::size_t k=0;k+1<n;k+=2) {
                const double rate=double(k+1);
                q[i*n+j]+=rate*(-h[i*n+k]*h[j*n+k+1]+h[i*n+k+1]*h[j*n+k]);
            }
        }
    }
    motion.pressure=q;motion.flow=q;
    ducklips::ContinuousModesReference solver;solver.prepare(n);
    for (std::size_t i=0;i<2*n;++i) { solver.state[i]=.02*std::sin(double(i+1)); }
    const std::vector<double> initial=solver.state;
    for (std::size_t block=0;block<2;++block) {
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t j=0;j<n;++j) { paired[block*n+i]+=h[i*n+j]*solver.state[block*n+j]; }
        }
        for (std::size_t i=0;i+1<n;i+=2) {
            const double c=std::cos(double(i+1)*.01),s=std::sin(double(i+1)*.01);
            const double a=paired[block*n+i],b=paired[block*n+i+1];
            paired[block*n+i]=c*a-s*b;paired[block*n+i+1]=s*a+c*b;
        }
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t j=0;j<n;++j) { expected[block*n+i]+=h[i*n+j]*paired[block*n+j]; }
        }
    }
    for (std::size_t i=0;i<n;++i) {
        const double c=std::cos(motion.omega[i]*.01),s=std::sin(motion.omega[i]*.01);
        const double a=expected[i],b=expected[n+i];expected[i]=c*a-s*b;expected[n+i]=s*a+c*b;
    }
    const std::size_t before=allocations.load();
    for (int i=0;i<480;++i) { advance(solver,motion.view(),motion.view(),motion.view(),1.0/48000); }
    require(allocations.load()==before,"continuous acoustic step allocates");
    for (std::size_t i=0;i<2*n;++i) { require(std::abs(solver.state[i]-expected[i])<1e-10,"dense rotation mismatch"); }
    const std::vector<double> full=solver.state;
    solver.state=initial;
    const std::size_t canonical_before=allocations.load();
    for (int i=0;i<480;++i) { advance(solver,motion.canonical_view(),motion.canonical_view(),motion.canonical_view(),1.0/48000); }
    require(allocations.load()==canonical_before,"canonical acoustic step allocates");
    for (std::size_t i=0;i<2*n;++i) {
        require(std::abs(solver.state[i]-expected[i])<1e-10,"dense canonical rotation mismatch");
        require(std::abs(solver.state[i]-full[i])<1e-12,"canonical/full transport diverged");
    }
}
void tests() {
    for (std::size_t n:{8,9,20,32}) { dense_rotation_test(n); }
    ducklips::ContinuousModesReference solver;solver.prepare(1);
    Motion motion(1);motion.omega[0]=2*pi*12000;
    solver.state={.3,.7};
    for (int i=0;i<48000;++i) { advance(solver,motion.view(),motion.view(),motion.view(),1.0/48000); }
    require(std::abs(std::complex<double>(solver.state[0],solver.state[1])-std::complex<double>(.3,.7))<1e-10,
            "stationary acoustic rotation accumulated phase bias");
    // Independently known phase integral of a quadratic chirp; zero coupling.
    solver.state={.3,.7};Motion first(1),middle(1),last(1);
    const double duration=.1,dt=1.0/48000;
    for (int i=0;i<4800;++i) {
        int stage=0;
        for (Motion* sample:{&first,&middle,&last}) {
            const double t=(i+.5*stage)*dt;
            sample->omega[0]=2*pi*(300+4000*t+7000*t*t);++stage;
        }
        advance(solver,first.view(),middle.view(),last.view(),dt);
    }
    const double phase=2*pi*(300*duration+2000*duration*duration+7000*duration*duration*duration/3);
    const std::complex<double> expected=std::exp(std::complex<double>(0,phase))*std::complex<double>(.3,.7);
    require(std::abs(std::complex<double>(solver.state[0],solver.state[1])-expected)<1e-11,"quadratic frequency integration failed");
    // Verify both sides of the bounded phase-reuse threshold and its direct
    // trigonometric fallback against the independently integrated phase.
    for (double correction:{-0.001001,-0.001,0.0,0.001,0.001001,0.03}) {
        first.omega[0]=12000;middle.omega[0]=13000;last.omega[0]=12000+4*correction/dt;
        solver.state={.3,.7};advance(solver,first.view(),middle.view(),last.view(),dt);
        const double whole=dt*(first.omega[0]+4*middle.omega[0]+last.omega[0])/6;
        const std::complex<double> exact_phase=std::exp(std::complex<double>(0,whole))*std::complex<double>(.3,.7);
        require(std::abs(std::complex<double>(solver.state[0],solver.state[1])-exact_phase)<1e-14,"phase reuse/fallback mismatch");
    }
    // Uniform moving metric has an independently derived exact propagator.
    motion.pressure[0]=-7;motion.flow[0]=7;motion.omega[0]=2*pi*1234;
    solver.state={.3,.7};
    for (int i=0;i<4800;++i) { advance(solver,motion.view(),motion.view(),motion.view(),dt); }
    ducklips::UniformAreaModeStep exact;exact.prepare(motion.omega[0],0,14,duration);
    require(std::abs(std::complex<double>(solver.state[0],solver.state[1])-exact.apply({.3,.7}))<1e-7,
            "uniform moving transport mismatch");
    solver.state={.3,.7};motion.diagonal[0]=motion.flow[0]+motion.pressure[0];
    for (int i=0;i<4800;++i) { advance(solver,motion.canonical_view(),motion.canonical_view(),motion.canonical_view(),dt); }
    require(std::abs(std::complex<double>(solver.state[0],solver.state[1])-exact.apply({.3,.7}))<1e-7,
            "canonical uniform moving transport mismatch");
    std::cout<<"Continuous modes dense rotations, allocation check, stationary phase, quadratic chirp and uniform metric reference passed\n";
}
void benchmark(std::size_t n,bool canonical=false) {
    require(n>0 && n<=256,"invalid benchmark mode count");
    ducklips::ContinuousModesReference solver;solver.prepare(n);solver.state[0]=.3;solver.state[n]=.7;
    Motion first(n),middle(n),last(n),base(n);
    for (std::size_t i=0;i<n;++i) { base.omega[i]=2*pi*(200+double(i)*10000/double(n)); }
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<n;++j) {
            const double difference=double(i)-double(j);
            base.pressure[i*n+j]=-2*std::sin(difference)/(1+std::abs(difference));
            base.flow[i*n+j]=base.pressure[i*n+j]*base.omega[i]/base.omega[j];
        }
    }
    const std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    for (int frame=0;frame<48000*5;++frame) {
        int stage=0;
        for (Motion* sample:{&first,&middle,&last}) {
            const double t=(frame+.5*stage)/48000.0;
            const double modulation=1+.1*std::sin(2*pi*3*t);
            for (std::size_t i=0;i<n;++i) { sample->omega[i]=base.omega[i]*modulation; }
            for (std::size_t i=0;i<n*n;++i) {
                sample->pressure[i]=base.pressure[i]*modulation;
                if (!canonical) { sample->flow[i]=base.flow[i]*modulation; }
            }
            ++stage;
        }
        if (canonical) { advance(solver,first.canonical_view(),middle.canonical_view(),last.canonical_view(),1.0/48000); }
        else { advance(solver,first.view(),middle.view(),last.view(),1.0/48000); }
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    double checksum=0;for (double v:solver.state) { checksum+=v; }
    require(std::isfinite(checksum),"nonfinite benchmark");
    std::cout<<std::setprecision(17)<<"{\"canonical\":"<<(canonical ? "true" : "false")<<",\"modes\":"<<n<<",\"states\":"<<2*n<<",\"simulated_seconds\":5,\"wall_seconds\":"<<seconds
             <<",\"realtime_multiple\":"<<5/seconds<<",\"checksum\":"<<checksum<<"}\n";
}
}
int main(int argc,char** argv) {
    try {
        if (argc==1) { tests();use_gauss=true;tests(); }
        else if (argc==4 && std::string(argv[1])=="--fixture") { fixture(argv[2],argv[3]); }
        else if (argc==4 && std::string(argv[1])=="--canonical-fixture") { fixture(argv[2],argv[3],true); }
        else if (argc==3 && std::string(argv[1])=="--canonical-benchmark") { benchmark(std::stoul(argv[2]),true); }
        else if (argc==3 && std::string(argv[1])=="--benchmark") { benchmark(std::stoul(argv[2])); }
        else { throw std::invalid_argument("continuous_modes_probe [--fixture input output | --benchmark modes]"); }
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
