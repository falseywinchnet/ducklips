#include "ducklips/lip_profile.hpp"
#include "ducklips/moving_graph.hpp"
#include "ducklips/physics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
namespace {
using ducklips::LipProfile;
bool read(LipProfile& p){return std::scanf("%lf %lf %lf %lf %lf %lf %lf",&p.length,&p.base_area,&p.activation,&p.place,&p.core,&p.width,&p.target_area)==7 && p.valid();}
bool connected(const LipProfile& p){
    return p.minimum_area()>0 || p.place-p.core<=0 || p.place+p.core>=p.length;
}
LipProfile between(const LipProfile& a,const LipProfile& b,double t){
    return {std::lerp(a.length,b.length,t),std::lerp(a.base_area,b.base_area,t),std::lerp(a.activation,b.activation,t),
            std::lerp(a.place,b.place,t),std::lerp(a.core,b.core,t),std::lerp(a.width,b.width,t),std::lerp(a.target_area,b.target_area,t)};
}
void values(const char* name,const std::vector<double>& v){std::printf("\"%s\":[",name);for(std::size_t i=0;i<v.size();++i){std::printf("%s%.17g",i?",":"",v[i]);}std::putchar(']');}
struct Metrics {
    std::vector<double> volume,centre,compliance,conductance,midpoint;
    explicit Metrics(std::size_t n):volume(n),centre(n),compliance(n),conductance(n+1),midpoint(n){}
    bool build(const LipProfile& p){
        const std::size_t n=volume.size();const ducklips::AirProperties air{};
        for(std::size_t i=0;i<n;++i){
            const double a=p.length*double(n-i-1)/double(n),b=p.length*double(n-i)/double(n);
            if(!p.measure(a,b,volume[i],centre[i])){return false;}
            compliance[i]=volume[i]/(air.density*air.sound_speed*air.sound_speed);midpoint[i]=p.area(.5*(a+b));
        }
        if(!p.conductance(centre[0],p.length,air.density,conductance[0])){return false;}
        for(std::size_t i=1;i<n;++i){if(!p.conductance(centre[i],centre[i-1],air.density,conductance[i])){return false;}}
        return p.conductance(0,centre.back(),air.density,conductance.back());
    }
    void write()const{std::putchar('{');values("volume",volume);std::putchar(',');values("centroid",centre);std::putchar(',');values("compliance",compliance);std::putchar(',');values("conductance",conductance);std::putchar(',');values("midpoint_area",midpoint);std::putchar('}');}
};
void state(const ducklips::MovingGraphReference& acoustic,const Metrics& m){
    std::printf(",\"charge\":");std::putchar('[');for(std::size_t i=0;i<acoustic.charge.size();++i){std::printf("%s%.17g",i?",":"",acoustic.charge[i]);}std::putchar(']');
    std::putchar(',');values("momentum",acoustic.momentum);std::putchar(',');values("average_flow",acoustic.average_flow);
    std::vector<double> pressure(m.volume.size());for(std::size_t i=0;i<pressure.size();++i){pressure[i]=acoustic.charge[i]/m.compliance[i];}
    std::putchar(',');values("pressure",pressure);
}
bool tests(){
    LipProfile p{.length=.0109375,.base_area=.0003,.activation=0,.place=0,.core=.0015,.width=.0045,.target_area=0};
    Metrics m(4);if(!m.build(p)){return false;}
    for(std::size_t i=0;i<4;++i){if(std::abs(m.volume[i]/(p.base_area*p.length/4)-1)>1e-14){return false;}}
    p.activation=1;if(!m.build(p) || m.midpoint.back()!=0 || m.volume.back()<=0 || m.conductance.back()!=0 || m.conductance[3]<=0){return false;}
    const double delta=p.length/4-p.core,u=delta/p.width;
    const double exact=p.base_area*p.width*(u-std::sqrt(ducklips::pi)/2*std::erf(u));
    if(std::abs(m.volume.back()/exact-1)>1e-12){return false;}
    // A tiny positive shoulder remains air; no subtraction-to-zero shortcut.
    const double x=p.core+1e-10;if(!(p.area(x)>0)){return false;}
    double v=0,c=0,g=0;if(!p.measure(p.core,x,v,c) || v<=0){return false;}
    if(!p.conductance(p.core,p.length,1.14,g) || g!=0){return false;}
    p.activation=0;if(!p.conductance(0,p.length,1.14,g) || std::abs(g/(p.base_area/(1.14*p.length))-1)>1e-14){return false;}
    p.activation=1;p.target_area=1e-18;
    if(!p.conductance(0,p.length,1.14,g) || g<=0 || g>=p.target_area/(1.14*p.core)){return false;}
    p.activation=1-1e-9;p.target_area=0;
    if(p.minimum_area()<=0 || !p.conductance(0,p.length,1.14,g) || g<=0){return false;}
    p.width=0;if(p.measure(0,p.length,v,c)){return false;}
    return true;
}
}
int main(int argc,char** argv){
    if(argc==1){if(!tests()){std::fprintf(stderr,"lip profile test failed\n");return 1;}std::puts("PASS uniform tube, partial contact volume, exact sealed port, tiny positive shoulder and aperture, invalid input");return 0;}
    if(argc!=2){return 2;}const std::string mode=argv[1];
    if(mode=="--metrics"){
        std::size_t n=0;while(std::scanf("%zu",&n)==1){LipProfile p;if(!n || n>256 || !read(p)){return 2;}Metrics m(n);if(!m.build(p)){return 1;}m.write();std::puts("");}return 0;
    }
    if(mode!="--sequence"){return 2;}
    // n intervals subdivisions initial_pressure inlet_pressure; first profile;
    // each interval: duration then end profile. Exterior pressure is zero.
    std::size_t n=0,intervals=0,subdivisions=0;double initial=0,inlet=0;
    if(std::scanf("%zu %zu %zu %lf %lf",&n,&intervals,&subdivisions,&initial,&inlet)!=5 || !n || n>64 || !intervals || intervals>100000 || !subdivisions || subdivisions>8192 || !std::isfinite(initial) || !std::isfinite(inlet)){return 2;}
    LipProfile p;if(!read(p)){return 2;}
    if(!connected(p)){std::fprintf(stderr,"interior contact requires separate air-region storage\n");return 1;}
    Metrics end(n),first(n),second(n);if(!end.build(p)){return 1;}
    ducklips::MovingGraphReference acoustic;std::vector<ducklips::TreeEdge> edges{{0,ducklips::ground_node}};
    for(std::size_t i=1;i<n;++i){edges.push_back({i-1,i});}edges.push_back({n-1,ducklips::ground_node});acoustic.prepare(n,edges);
    std::vector<double> drive(n),series(n+1,0.);series[0]=-inlet;
    for(std::size_t i=0;i<n;++i){if(end.volume[i]<=0){std::fprintf(stderr,"empty section requires state transport\n");return 1;}acoustic.charge[i]=end.compliance[i]*initial;}
    std::printf("{\"kind\":\"initial\",\"geometry\":");end.write();state(acoustic,end);std::puts("}");
    for(std::size_t interval=0;interval<intervals;++interval){
        double duration=0;LipProfile q;if(std::scanf("%lf",&duration)!=1 || !std::isfinite(duration) || duration<=0 || !read(q)){return 2;}
        if(!connected(q)){std::fprintf(stderr,"interior contact requires separate air-region storage\n");return 1;}
        std::vector<double> cuts;cuts.reserve(subdivisions+2);
        for(std::size_t k=0;k<=subdivisions;++k){cuts.push_back(double(k)/double(subdivisions));}
        if(p.target_area==0 && q.target_area==0 && (p.activation<1)!=(q.activation<1)){
            const double contact=(1-p.activation)/(q.activation-p.activation);
            if(contact>0 && contact<1){cuts.push_back(contact);}
        }
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(std::size_t k=0;k+1<cuts.size();++k){
            const double a=cuts[k],b=cuts[k+1],dt=duration*(b-a);
            const auto old=end.volume;
            const LipProfile p1=between(p,q,a+.21132486540518711775*(b-a)),p2=between(p,q,a+.78867513459481288225*(b-a)),pe=between(p,q,b);
            if(!connected(p1) || !connected(p2) || !connected(pe)){std::fprintf(stderr,"interior contact requires separate air-region storage\n");return 1;}
            if(!first.build(p1) || !second.build(p2) || !end.build(pe)){std::fprintf(stderr,"lip geometry integration failed\n");return 1;}
            for(std::size_t i=0;i<n;++i){if(end.volume[i]<=0){std::fprintf(stderr,"empty section requires state transport\n");return 1;}drive[i]=-(end.volume[i]-old[i])/dt;}
            for(std::size_t i=0;i<n+1;++i){if(first.conductance[i]==0 && second.conductance[i]==0){acoustic.momentum[i]=0;}}
            if(!acoustic.advance(first.compliance,second.compliance,first.conductance,second.conductance,dt,drive,series)){return 1;}
            for(std::size_t i=0;i<n+1;++i){if(end.conductance[i]==0){acoustic.momentum[i]=0;}}
            std::printf("{\"interval\":%zu,\"start_fraction\":%.17g,\"fraction\":%.17g,\"dt\":%.17g,\"stage1\":",interval,a,b,dt);first.write();std::printf(",\"stage2\":");second.write();std::printf(",\"geometry\":");end.write();state(acoustic,end);std::puts("}");
        }
        p=q;
    }
}
