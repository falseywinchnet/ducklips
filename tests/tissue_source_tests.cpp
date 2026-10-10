#include "ducklips/tissue_source.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool condition,const char* message){if(!condition){throw std::runtime_error(message);}}
int main(){
    ducklips::TissueSource source;source.reset();
    constexpr double dt=1./96000.;
    for(int i=0;i<9600;++i){
        require(source.advance(0,0,dt)==0,"zero breath must not create flow");
        require(source.kinetic_energy()==0 && source.elastic_energy()==0,"resting tissue must stay at rest without pressure");
    }
    for(int i=0;i<28800;++i){source.advance(600*(1-std::exp(-i*dt/.012)),0,dt);}
    const double initial=source.kinetic_energy()+source.elastic_energy();
    require(initial>1e-10,"breath must transfer energy to tissue");
    double maximum=initial;
    for(int i=0;i<96000;++i){
        require(source.advance(0,0,dt)==0,"zero pressure drop must not prescribe an acoustic pulse");
        const double energy=source.kinetic_energy()+source.elastic_energy();
        require(std::isfinite(energy),"mechanical energy must remain finite");
        maximum=std::max(maximum,energy);
    }
    const double final=source.kinetic_energy()+source.elastic_energy();
    require(maximum<initial*1.03,"unforced tissue must not generate mechanical energy");
    require(final<initial*1e-4,"tissue must recoil and dissipate after pressure removal");
    const auto old_gap=source.gap();
    const auto old_velocity=source.velocity();
    source.set_adduction(-.0006);
    for(std::size_t i=0;i<ducklips::TissueSource::count;++i){
        require(std::abs(source.gap()[i]-old_gap[i])<1e-15,"adduction must not teleport tissue");
        require(source.velocity()[i]==old_velocity[i],"adduction must preserve momentum");
    }
    // A stationary open leak is not phonation. Require the configured,
    // initially apposed body to develop recurring open/closed flow under
    // constant pressure, without a periodic driver.
    ducklips::TissueMaterial apposed;apposed.rest_gap=0;
    source.reset(apposed);
    double low=1,high=0,previous_flow=0;int openings=0;
    for(int i=0;i<48000;++i){
        const double flow=source.advance(600*(1-std::exp(-i*dt/.012)),0,dt);
        if(i>=33600){
            low=std::min(low,flow);high=std::max(high,flow);
            if(previous_flow<1e-6 && flow>=1e-6){++openings;}
        }
        previous_flow=flow;
    }
    require(low<1e-8 && high>1e-4 && openings>=10,"apposed body must sustain pressure-driven phonation with closure");
    std::cout<<"zero-input rest, pressure-driven motion, unforced recoil and dissipation passed\n"
        <<"energy initial="<<initial<<" maximum="<<maximum<<" final="<<final<<'\n'
        <<"constant-pressure opening cycles="<<openings<<" flow range="<<low<<" to "<<high<<'\n';
}
