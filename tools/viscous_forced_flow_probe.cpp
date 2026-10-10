#include "ducklips/slit_viscous.hpp"
#include <cmath>
#include <cstdio>
int main(){
    double inertia,gap,pressure,dt;unsigned steps;
    while(std::scanf("%lf %lf %lf %lf %u",&inertia,&gap,&pressure,&dt,&steps)==5){
        ducklips::SlitViscous model;const ducklips::AirProperties air{};const double g=1/inertia;
        if(!steps || steps>1000000 || !std::isfinite(pressure) || !model.prepare(air,inertia,g,gap)){return 2;}
        const double resistance=3*inertia*(air.viscosity/air.density)/(gap*gap),exact=pressure/resistance;
        if(!model.equilibrate(exact)){return 2;}double p=exact/g,flow=0,residual=0;
        for(unsigned i=0;i<steps;++i){
            const double before=.5*g*p*p+model.energy();
            if(!model.advance(p,dt/2)){return 1;}double heat=model.last.dissipated_energy;
            flow=g*(p+.5*dt*pressure);p+=dt*pressure;
            if(!model.advance(p,dt/2)){return 1;}heat+=model.last.dissipated_energy;
            const double after=.5*g*p*p+model.energy();residual=after-before+heat-dt*pressure*flow;
        }
        std::printf("%.17g %.17g %.17g %.17g\n",flow,g*p,exact,residual);
    }
}
