#include "ducklips/motion_segment.hpp"
#include <chrono>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <new>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace { std::atomic<std::size_t> allocations{0}; }
void* operator new(std::size_t count) {
    ++allocations;void* result=std::malloc(std::max<std::size_t>(count,1));
    if (!result) { throw std::bad_alloc(); }return result;
}
void* operator new[](std::size_t count) { return ::operator new(count); }
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
namespace {
constexpr double pi=3.14159265358979323846;
void require(bool okay,const char* message) { if (!okay) { throw std::runtime_error(message); } }
struct Knot final {
    std::vector<double> omega,rate,pressure;
    explicit Knot(std::size_t n):omega(n),rate(n),pressure(n*n) {}
    ducklips::CanonicalMotionKnot view() const { return {omega,rate,pressure}; }
};
void read(std::ifstream& in,void* target,std::size_t bytes) {
    in.read(static_cast<char*>(target),static_cast<std::streamsize>(bytes));require(bool(in),"truncated motion segment fixture");
}
void fixture(const char* input,const char* output) {
    std::ifstream in(input,std::ios::binary);std::uint64_t header[3];double interval;
    read(in,header,sizeof(header));read(in,&interval,sizeof(interval));
    const std::size_t n=header[0],cells=header[1],samples=header[2];
    require(n>0 && n<=128 && cells>0 && cells<=10000 && samples>0 && samples<=4096 && std::isfinite(interval) && interval>0,"invalid segment fixture");
    ducklips::ContinuousModesReference solver;solver.prepare(n);
    read(in,solver.state.data(),2*n*sizeof(double));
    std::vector<Knot> knots;knots.reserve(cells+5);
    for (std::size_t i=0;i<cells+5;++i) {
        knots.emplace_back(n);Knot& knot=knots.back();
        read(in,knot.omega.data(),n*sizeof(double));read(in,knot.rate.data(),n*sizeof(double));read(in,knot.pressure.data(),n*n*sizeof(double));
    }
    ducklips::CanonicalMotionSegment segment;segment.prepare(n);
    ducklips::CanonicalMotionFrame first,middle,last;first.prepare(n);middle.prepare(n);last.prepare(n);
    std::ofstream out(output,std::ios::binary);require(bool(out),"cannot open segment output");
    for (std::size_t cell=0;cell<cells;++cell) {
        require(segment.configure({knots[cell].view(),knots[cell+1].view(),knots[cell+2].view(),knots[cell+3].view(),knots[cell+4].view(),knots[cell+5].view()},interval),"invalid control knots");
        require(segment.evaluate(0,first),"invalid initial frame");
        for (std::size_t sample=0;sample<samples;++sample) {
            require(segment.evaluate((double(sample)+.5)/samples,middle) && segment.evaluate(double(sample+1)/samples,last),"invalid interpolated frame");
            solver.advance_gauss(first.view(),middle.view(),last.view(),interval/samples);
            std::swap(first,last);
            out.write(reinterpret_cast<const char*>(solver.state.data()),static_cast<std::streamsize>(2*n*sizeof(double)));
        }
    }
    require(bool(out),"failed segment output");
}
void tests() {
    std::array<Knot,6> knots{Knot(1),Knot(1),Knot(1),Knot(1),Knot(1),Knot(1)};
    for (std::size_t i=0;i<6;++i) {
        const double x=double(i)-2;
        knots[i].omega[0]=1000+x*(2+x*(3+x*(4+x*(5+6*x))));
        knots[i].rate[0]=(2+x*(6+x*(12+x*(20+30*x))))/.01;
        knots[i].pressure[0]=.3+x*(.2+x*(-.1+.03*x));
    }
    ducklips::CanonicalMotionSegment segment;segment.prepare(1);
    ducklips::CanonicalMotionFrame frame;frame.prepare(1);
    require(!segment.evaluate(0,frame),"unconfigured segment accepted");
    require(segment.configure({knots[0].view(),knots[1].view(),knots[2].view(),knots[3].view(),knots[4].view(),knots[5].view()},.01),"polynomial setup failed");
    for (int i=0;i<=20;++i) {
        const double x=i/20.0,w=1000+x*(2+x*(3+x*(4+x*(5+6*x))));
        const double rate=(2+x*(6+x*(12+x*(20+30*x))))/.01;
        require(segment.evaluate(x,frame),"polynomial evaluation rejected");
        require(std::abs(frame.angular_frequency[0]-w)<1e-11 && std::abs(frame.flow_diagonal[0]-rate/w)<1e-10 &&
                std::abs(frame.pressure_transport[0]-(.3+x*(.2+x*(-.1+.03*x))))<1e-14,"polynomial reproduction failed");
    }
    const double old=frame.angular_frequency[0];
    require(!segment.evaluate(-.1,frame) && frame.angular_frequency[0]==old,"invalid time changed output");
    require(!segment.configure({knots[0].view(),knots[1].view(),knots[2].view(),knots[3].view(),knots[4].view(),knots[5].view()},std::numeric_limits<double>::infinity()),"invalid interval accepted");
    require(segment.evaluate(1,frame) && frame.angular_frequency[0]==old,"bad configuration lost previous segment");
    std::array<Knot,7> curved{Knot(1),Knot(1),Knot(1),Knot(1),Knot(1),Knot(1),Knot(1)};
    for (std::size_t i=0;i<7;++i) {
        const double x=double(i)-2;curved[i].omega[0]=1000+.1*std::sin(x);curved[i].rate[0]=10*std::cos(x);
        curved[i].pressure[0]=std::sin(.7*x)+.03*x*x*x*x;
    }
    require(segment.configure({curved[0].view(),curved[1].view(),curved[2].view(),curved[3].view(),curved[4].view(),curved[5].view()},.01),"left segment rejected");
    require(segment.evaluate(1-1e-5,frame),"left slope rejected");const double near_left=frame.pressure_transport[0];
    require(segment.evaluate(1,frame),"left endpoint rejected");const double left=frame.pressure_transport[0];
    require(segment.configure({curved[1].view(),curved[2].view(),curved[3].view(),curved[4].view(),curved[5].view(),curved[6].view()},.01),"right segment rejected");
    require(segment.evaluate(0,frame),"right endpoint rejected");const double right=frame.pressure_transport[0];
    require(segment.evaluate(1e-5,frame),"right slope rejected");const double near_right=frame.pressure_transport[0];
    require(std::abs(left-right)<1e-13 && std::abs((left-near_left-near_right+right)/1e-5)<1e-4,"pressure value/slope jump at control knot");
    for (Knot& knot:knots) { knot.omega[0]=1;knot.rate[0]=0; }
    knots[2].rate[0]=-10000;
    require(segment.configure({knots[0].view(),knots[1].view(),knots[2].view(),knots[3].view(),knots[4].view(),knots[5].view()},.01),"overshoot fixture setup rejected");
    const double saved_frequency=frame.angular_frequency[0],saved_pressure=frame.pressure_transport[0];
    require(!segment.evaluate(.5,frame) && frame.angular_frequency[0]==saved_frequency && frame.pressure_transport[0]==saved_pressure,"nonpositive interpolated frequency changed frame");
    std::cout<<"Motion segment polynomial reproduction, knot slope continuity, physical phase derivative and transactional rejection passed\n";
}
void fill(Knot& knot,std::size_t n,double t,const std::vector<double>& base) {
    const double modulation=1+.1*std::sin(2*pi*3*t),rate=.1*2*pi*3*std::cos(2*pi*3*t);
    for (std::size_t i=0;i<n;++i) {
        const double omega=2*pi*(200+double(i)*10000/double(n));knot.omega[i]=omega*modulation;knot.rate[i]=omega*rate;
    }
    for (std::size_t i=0;i<n*n;++i) { knot.pressure[i]=base[i]*modulation; }
}
void benchmark(std::size_t n) {
    require(n>0 && n<=128,"invalid benchmark size");
    std::vector<double> base(n*n);
    for (std::size_t i=0;i<n;++i) {
        for (std::size_t j=0;j<n;++j) { const double d=double(i)-double(j);base[i*n+j]=-2*std::sin(d)/(1+std::abs(d)); }
    }
    std::array<Knot,6> knots{Knot(n),Knot(n),Knot(n),Knot(n),Knot(n),Knot(n)};
    for (std::size_t i=0;i<6;++i) { fill(knots[i],n,(double(i)-2)*.001,base); }
    ducklips::ContinuousModesReference solver;solver.prepare(n);solver.state[0]=.3;solver.state[n]=.7;
    ducklips::CanonicalMotionSegment segment;segment.prepare(n);
    ducklips::CanonicalMotionFrame first,middle,last;first.prepare(n);middle.prepare(n);last.prepare(n);
    const std::size_t before=allocations.load();
    const std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    for (std::size_t cell=0;cell<5000;++cell) {
        const std::size_t head=cell%6;
        require(segment.configure({knots[head].view(),knots[(head+1)%6].view(),knots[(head+2)%6].view(),knots[(head+3)%6].view(),knots[(head+4)%6].view(),knots[(head+5)%6].view()},.001),"benchmark configuration rejected");
        require(segment.evaluate(0,first),"benchmark frame rejected");
        for (int sample=0;sample<48;++sample) {
            require(segment.evaluate((sample+.5)/48,middle) && segment.evaluate((sample+1.0)/48,last),"benchmark interpolation rejected");
            solver.advance_gauss(first.view(),middle.view(),last.view(),1.0/48000);
            std::swap(first,last);
        }
        fill(knots[head],n,(double(cell)+4)*.001,base);
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    const std::size_t heap_allocations=allocations.load()-before;
    require(heap_allocations==0,"segment interpolation/stepping allocates");
    double checksum=0;for (double v:solver.state) { checksum+=v; }require(std::isfinite(checksum),"nonfinite benchmark");
    std::cout<<std::setprecision(17)<<"{\"modes\":"<<n<<",\"simulated_seconds\":5,\"control_rate_hz\":1000,\"wall_seconds\":"<<seconds
             <<",\"realtime_multiple\":"<<5/seconds<<",\"heap_allocations\":"<<heap_allocations<<",\"checksum\":"<<checksum<<"}\n";
}
}
int main(int argc,char** argv) {
    try {
        if (argc==1) { tests(); }
        else if (argc==4 && std::string(argv[1])=="--fixture") { fixture(argv[2],argv[3]); }
        else if (argc==3 && std::string(argv[1])=="--benchmark") { benchmark(std::stoul(argv[2])); }
        else { throw std::invalid_argument("motion_segment_probe [--fixture input output | --benchmark modes]"); }
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
