#include "ducklips/tissue_source.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "ducklips/motor_hold.hpp"

void require(bool condition,const char* message){if(!condition){throw std::runtime_error(message);}}
int main(){
    // A held deformation must actually release stress. Changing random
    // rate changes the trajectory but cannot supply mechanical energy.
    ducklips::TissueRelaxation memory;memory.reset(.1,.1,.6,42);
    ducklips::TissueRelaxation::Coordinates held{};held[0]=1e-6;
    double correction=0;bool varying_rate=false;const auto initial_rates=memory.rates();
    for(int step=0;step<1000;++step){
        memory.advance(held,.001);
        const double next=memory.energy_correction(held);
        require(next<=correction+1e-18,"stochastic relaxation cannot add energy at fixed strain");correction=next;
        varying_rate=varying_rate || memory.rates()!=initial_rates;
    }
    require(correction<0 && varying_rate,"held strain must relax with changing positive rates");
    memory.reset(.1,.1,0,42);memory.advance(held,.1);
    const auto strain=memory.strain(held);double expected=0;
    for(double value:strain){expected+=.1*value*value*(std::exp(-2)-1);}
    require(std::abs(memory.energy_correction(held)-expected)<1e-16,"constant-rate relaxation must match analytic held-strain energy");
    ducklips::MotorPool pool;pool.reset(42);
    require(pool.drive_for_activation(0)==0,"zero requested activation must not recruit a motor unit");
    for(int step=0;step<500;++step){pool.advance(.6,.001);}
    const double held_force=pool.activation();const auto spike_count=pool.spikes();
    for(int step=0;step<1000;++step){pool.advance(0,.001);}
    require(held_force>0 && pool.activation()<held_force*1e-8,"muscle force must relax after firing stops");
    require(pool.spikes()==spike_count,"disabled muscle must not fire");
    ducklips::MotorHold reference_hold,perturbed_hold;
    reference_hold.reset(.0003,.00025,44.86,19);perturbed_hold.reset(.0003,.00025,44.86,19);
    for(int step=0;step<10;++step){
        reference_hold.advance(.0003);perturbed_hold.advance(.0005);
        require(reference_hold.drive()==perturbed_hold.drive(),"feedback must not see a perturbation before its delay");
    }
    for(int step=0;step<90;++step){reference_hold.advance(.0003);perturbed_hold.advance(.0005);}
    require(perturbed_hold.drive()>reference_hold.drive(),"delayed opening error must recruit greater closing effort");
    for(int step=0;step<400;++step){perturbed_hold.advance(.0008,false);}
    require(perturbed_hold.drive()<0 && perturbed_hold.force()<0,"opening gesture must recruit an opposing muscle force");
    std::cout<<"stochastic stress relaxation, analytic held-strain limit, motor twitch release and delayed posture correction passed\n";
    ducklips::TissueSource source;source.reset();
    constexpr double dt=1./96000.;
    require(source.advance(0,600,dt)<0,"reversed pressure must produce reversed finite quasi-steady flow");
    source.reset();
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
