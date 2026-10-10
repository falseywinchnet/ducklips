#include "ducklips/contour_pressure_flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace {
using Matrix=std::array<double,25>;using Vector=std::array<double,5>;
constexpr double velocity_scale[]={1e-4,1e-4,1e-4,.001,.001};
bool solve(Matrix a,Vector b,Vector& x){
    for(std::size_t j=0;j<5;++j){std::size_t pivot=j;for(std::size_t i=j+1;i<5;++i){if(std::abs(a[5*i+j])>std::abs(a[5*pivot+j])){pivot=i;}}
        if(!std::isfinite(a[5*pivot+j]) || a[5*pivot+j]==0){return false;}
        for(std::size_t k=j;k<5;++k){std::swap(a[5*j+k],a[5*pivot+k]);}std::swap(b[j],b[pivot]);
        for(std::size_t i=j+1;i<5;++i){const double ratio=a[5*i+j]/a[5*j+j];for(std::size_t k=j+1;k<5;++k){a[5*i+k]-=ratio*a[5*j+k];}b[i]-=ratio*b[j];}}
    for(std::size_t j=5;j-->0;){double value=b[j];for(std::size_t k=j+1;k<5;++k){value-=a[5*j+k]*x[k];}x[j]=value/a[5*j+j];}return true;
}
Matrix metric(const ducklips::ContourPressureFlow& flow,const ducklips::FoldParameters& f){
    Matrix result=flow.inertia;result[18]+=f.lower_mass/(2*f.length*f.length);result[24]+=f.upper_mass/(2*f.length*f.length);
    for(std::size_t i=0;i<5;++i){for(std::size_t j=0;j<5;++j){result[5*i+j]*=velocity_scale[i]*velocity_scale[j]/1e-6;}}return result;
}
double quadratic(const Vector& p,const Matrix& inverse){double result=0;for(std::size_t i=0;i<5;++i){for(std::size_t j=0;j<5;++j){result+=.5*p[i]*inverse[5*i+j]*p[j];}}return result;}
}
int main(int argc,char** argv){try{
    const bool relative_mass=argc==2 && std::string(argv[1])=="--relative-mass";
    if(argc!=1 && !relative_mass){throw std::runtime_error("unknown option");}
    ducklips::FoldParameters f;ducklips::ContourPressureStorage storage;ducklips::ContourPressureKinetics kinetics;
    if(!storage.configure(f,{},ducklips::ContourPressureProfile::linear_span) || !kinetics.configure(f,{},ducklips::ContourPressureProfile::linear_span)){throw std::runtime_error("configure");}
    double factor;
    while(std::scanf("%lf",&factor)==1){std::array<double,9> y;for(double& v:y){if(std::scanf("%lf",&v)!=1){throw std::runtime_error("input");}}
        if(factor<=0 || !std::isfinite(factor)){throw std::runtime_error("factor");}
        ducklips::ContourPressureState gas;ducklips::ContourPressureFlow flow;
        if(!storage.evaluate({y[2],y[3]},{y[0],y[1]},gas) || !kinetics.prepare({y[2],y[3]},{y[0],y[1]},flow)){throw std::runtime_error("initial configuration");}
        const Matrix original=metric(flow,f);Vector momentum{};
        for(std::size_t i=0;i<5;++i){for(std::size_t j=0;j<5;++j){momentum[i]+=original[5*i+j]*y[4+j]/velocity_scale[j];}}
        const std::array<double,4> z={gas.mass[0]/1e-8,gas.mass[1]/1e-8,y[2]/1e-6,y[3]/1e-6};
        const auto inverse=[&](const std::array<double,4>& state){
            ducklips::ContourPressureState recovered;ducklips::ContourPressureFlow current;
            const std::array<double,2> area={state[2]*1e-6,state[3]*1e-6},mass={state[0]*1e-8,state[1]*1e-8};
            if(!storage.recover(area,mass,{y[0],y[1]},recovered) || !kinetics.prepare(area,recovered.effort,current)){throw std::runtime_error("perturbed configuration");}
            Matrix m=metric(current,f);Matrix result{};
            const std::array<bool,5> active={area[0]>0,true,area[1]>0,true,true};
            for(std::size_t j=0;j<5;++j){if(!active[j]){for(std::size_t i=0;i<5;++i){m[5*i+j]=m[5*j+i]=0.;}m[5*j+j]=1.;}}
            for(std::size_t j=0;j<5;++j){if(!active[j]){continue;}Vector right{},column{};right[j]=1.;if(!solve(m,right,column)){throw std::runtime_error("inverse");}for(std::size_t i=0;i<5;++i){result[5*i+j]=active[i]?column[i]:0.;}}
            return result;
        };
        std::array<double,4> gradient{},steps{};
        for(std::size_t k=0;k<4;++k){double h=factor*2e-5*(relative_mass && k<2?std::abs(z[k]):std::max(1.,std::abs(z[k])));h=std::min(h,.1*std::abs(z[k]));if(h<=0){throw std::runtime_error("zero coordinate");}
            auto plus=z,minus=z;plus[k]+=h;minus[k]-=h;Matrix a=inverse(plus),b=inverse(minus);
            for(std::size_t j=0;j<25;++j){a[j]=(a[j]-b[j])/(2*h);}gradient[k]=quadratic(momentum,a)*(k<2?100.:1.);steps[k]=h*(k<2?1e-8:1e-6);}
        std::printf("{\"gradient\":[%.17g,%.17g,%.17g,%.17g],\"steps\":[%.17g,%.17g,%.17g,%.17g]}\n",gradient[0],gradient[1],gradient[2],gradient[3],steps[0],steps[1],steps[2],steps[3]);
    }return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
