#include "ducklips/contour_pressure_loss.hpp"
#include <cstdio>
#include <stdexcept>
#include <string>
namespace {template<std::size_t N> void array(const std::array<double,N>& a){std::printf("[");for(std::size_t i=0;i<N;++i){std::printf("%s%.17g",i?",":"",a[i]);}std::printf("]");}}
int main(int argc,char** argv){try{
 const auto profile=argc>1 && std::string(argv[1])=="linear"?ducklips::ContourPressureProfile::linear_span:ducklips::ContourPressureProfile::flat_caps;ducklips::ContourPressureKinetics k;ducklips::ContourPressureDissipation d;ducklips::ContourPressureStorage s;if(!k.configure({},{},profile) || !d.configure({},{},profile) || !s.configure({},{},profile)){throw std::runtime_error("configure");}
 double a,b,u0,u1;while(std::scanf("%lf%lf%lf%lf",&a,&b,&u0,&u1)==4){ducklips::ContourPressureFlow f;ducklips::ContourPressureLoss r;ducklips::ContourPressureState g;bool storage_ok=s.evaluate({a,b},{u0,u1},g),flow_ok=storage_ok && k.prepare({a,b},{u0,u1},f),loss_ok=flow_ok && d.prepare(f,r);bool ok=storage_ok && flow_ok && loss_ok;if(!ok){std::printf("{\"accepted\":false,\"storage\":%s,\"flow\":%s,\"loss\":%s}\n",storage_ok?"true":"false",flow_ok?"true":"false",loss_ok?"true":"false");continue;}std::printf("{\"accepted\":true,\"mass\":");array(g.mass);std::printf(",\"volume\":%.17g,\"energy\":%.17g,\"coenergy\":%.17g,\"area_effort\":",g.volume,g.energy,g.coenergy);array(g.area_effort);std::printf(",\"compliance\":");array(g.compliance);std::printf(",\"mass_area\":");array(g.mass_area);std::printf(",\"inertia\":");array(f.inertia);std::printf(",\"resistance\":");array(r.resistance);std::printf("}\n");}return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
