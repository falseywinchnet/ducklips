#include "ducklips/moving_mode.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
constexpr double pi=3.14159265358979323846;
void require(bool okay,const char* message) { if (!okay) { throw std::runtime_error(message); } }
}
int main() {
    try {
        ducklips::UniformAreaModeStep step;
        const std::complex<double> initial(.3,.7);
        for (double rate : {0.0,13.862943611198906,-13.862943611198906,2*2*pi*1234,20000.0}) {
            const double omega=2*pi*1234,dt=1.0/48000;
            step.prepare(omega,0,rate,dt);
            // Determinant is one for conservative area motion: phase-space
            // volume stays fixed even when instantaneous acoustic energy changes.
            require(std::abs(step.aa*step.bb-step.ab*step.ba-1)<1e-12,"moving step lost canonical determinant");
            std::complex<double> state=initial;
            ducklips::UniformAreaModeStep whole;
            const int count=rate>10000 ? 12 : 4800;
            whole.prepare(omega,0,rate,count*dt);
            for (int i=0;i<count;++i) { state=step.apply(state); }
            require(std::abs(state-whole.apply(initial))<1e-9,"sample subdivision changes exact trajectory");
            // Independent fine midpoint integration of wall power over one
            // step, using a second-order Taylor derivative reference at tiny h.
            const int divisions=1000;
            const double h=dt/divisions;
            std::complex<double> x=initial;
            double work=0;
            for (int i=0;i<divisions;++i) {
                const double a=x.real(),b=x.imag();
                const std::complex<double> derivative(-rate*a/2-omega*b,omega*a+rate*b/2);
                const double da=derivative.real(),db=derivative.imag();
                const std::complex<double> acceleration(-rate*da/2-omega*db,omega*da+rate*db/2);
                const std::complex<double> middle=x+derivative*(h/2)+acceleration*(h*h/8);
                work+=h*rate*(middle.imag()*middle.imag()-middle.real()*middle.real())/2;
                x+=h*derivative+(h*h/2)*acceleration;
            }
            const double energy_change=(std::norm(step.apply(initial))-std::norm(initial))/2;
            require(std::abs(energy_change-work)<1e-8,"energy change does not equal moving-wall work");
        }
        step.prepare(2*pi*1000,37,0,1.0/48000);
        const std::complex<double> expected=std::exp(std::complex<double>(-37,2*pi*1000)/48000.0)*initial;
        require(std::abs(step.apply(initial)-expected)<1e-15,"stationary step differs from existing modal dynamics");
        step.prepare(0,0,0,0);require(step.apply(initial)==initial,"zero step is not identity");
        bool rejected=false;
        try { step.prepare(100,0,std::numeric_limits<double>::infinity(),.1); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected && step.apply(initial)==initial,"invalid preparation is not transactional");
        std::cout<<"Uniform moving-area semigroup, canonical determinant, wall work and static limit passed\n";
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
