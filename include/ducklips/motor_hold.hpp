#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>

namespace ducklips {
// Experimental motor pool. A discharge initiates a continuous alpha twitch;
// successive twitches overlap. Firing intervals are irregular, not periodic.
// 15/35 ms twitch times and 0.2 interval CV follow Manriquez et al. (2017).
// Unit count, recruitment thresholds and maximum force are model assumptions.
class MotorPool final {
public:
    static constexpr std::size_t units=64;
    void reset(std::uint64_t seed){
        random_.seed(seed);normal_.reset();
        for(std::size_t i=0;i<units;++i){
            phase_[i]=std::generate_canonical<double,53>(random_);
            interval_[i]=next_interval();
        }
        fast_x_=fast_y_=slow_x_=slow_y_=0;spikes_=0;
    }
    [[nodiscard]] double expected_activation(double drive) const noexcept{
        double sum=0;
        for(std::size_t i=0;i<units;++i){sum+=firing_rate(drive,i)/40;}
        return sum/units;
    }
    [[nodiscard]] double drive_for_activation(double activation) const noexcept{
        if(activation<=0){return 0;}
        if(activation>=1){return 1;}
        double low=0,high=1;
        for(int i=0;i<12;++i){const double mid=.5*(low+high);if(expected_activation(mid)<activation){low=mid;}else{high=mid;}}
        return .5*(low+high);
    }
    double advance(double drive,double dt){
        const double fast_decay=std::exp(-dt/.015),slow_decay=std::exp(-dt/.035);
        fast_y_=(fast_y_+dt/.015*fast_x_)*fast_decay;fast_x_*=fast_decay;
        slow_y_=(slow_y_+dt/.035*slow_x_)*slow_decay;slow_x_*=slow_decay;
        for(std::size_t i=0;i<units;++i){
            const double rate=firing_rate(drive,i);
            if(rate==0){continue;}
            const double old=phase_[i];phase_[i]+=rate*dt;
            if(phase_[i]>=interval_[i]){
                const double remaining=dt-(interval_[i]-old)/rate;
                const double tau=i%10==0?.035:.015;
                // Per-unit force scale gives a full-drive mean of 1. It is
                // fixed analytically, never normalized to a rendered clip.
                const double jump=1./(units*40*tau);
                const double tail=jump*std::exp(-remaining/tau);
                if(i%10==0){slow_x_+=tail;slow_y_+=remaining/tau*tail;}
                else{fast_x_+=tail;fast_y_+=remaining/tau*tail;}
                phase_[i]-=interval_[i];interval_[i]=next_interval();++spikes_;
            }
        }
        return activation();
    }
    [[nodiscard]] double activation() const noexcept{return fast_y_+slow_y_;}
    [[nodiscard]] std::uint64_t spikes() const noexcept{return spikes_;}
private:
    static double firing_rate(double drive,std::size_t unit) noexcept{
        const double threshold=.5*double(unit)/units;
        return drive<=threshold?0:8+32*std::clamp((drive-threshold)/(1-threshold),0.,1.);
    }
    double next_interval(){return std::max(.2,1+.2*normal_(random_));}
    std::array<double,units> phase_{},interval_{};
    double fast_x_{},fast_y_{},slow_x_{},slow_y_{};
    std::uint64_t spikes_{};std::mt19937_64 random_{1};std::normal_distribution<double> normal_{};
};

// A deliberately small posture-control hypothesis, not an anatomical reflex
// model. It sees a 20 ms filtered shape after 30 ms transport delay. An error
// changes motor drive; it never clamps the tissue or specifies vibration F0.
class MotorHold final {
public:
    void set_target_gap(double target) noexcept{target_=target;}
    void reset(double passive_gap,double target_gap,double stiffness,std::uint64_t seed){
        pool_.reset(seed);opener_.reset(seed+101);rest_=passive_gap;target_=target_gap;stiffness_=stiffness;
        maximum_force_=stiffness*.00025;sensed_=passive_gap;delayed_.fill(passive_gap);
        index_=0;integral_=0;drive_=force_=0;
    }
    double advance(double gap,bool enabled=true){
        constexpr double dt=.001;
        sensed_+=-std::expm1(-dt/.020)*(gap-sensed_);
        const double delayed=delayed_[index_];delayed_[index_]=sensed_;index_=(index_+1)%delayed_.size();
        const double requested=enabled?target_:.0012;
        const double error=delayed-requested;
        const double feedforward=.5*stiffness_*(rest_-requested);
        const double proportional=.25*stiffness_*error;
        const double candidate=integral_+.5*stiffness_*error*dt/.2;
        const double effort=feedforward+proportional+candidate;
        if((effort>=-maximum_force_ && effort<=maximum_force_)||(effort< -maximum_force_ && error>0)||(effort>maximum_force_ && error<0)){integral_=candidate;}
        const double activation=std::clamp((feedforward+proportional+integral_)/maximum_force_,-1.,1.);
        drive_=activation>=0?pool_.drive_for_activation(activation):-opener_.drive_for_activation(-activation);
        // Opposing opening and closing muscles carry independent twitch
        // histories. A new gesture changes firing, never existing force.
        force_=maximum_force_*(pool_.advance(std::max(0.,drive_),dt)-opener_.advance(std::max(0.,-drive_),dt));
        return force_;
    }
    [[nodiscard]] double force() const noexcept{return force_;}
    [[nodiscard]] double drive() const noexcept{return drive_;}
    [[nodiscard]] double sensed_gap() const noexcept{return sensed_;}
    [[nodiscard]] std::uint64_t spikes() const noexcept{return pool_.spikes()+opener_.spikes();}
private:
    MotorPool pool_{},opener_{};std::array<double,30> delayed_{};std::size_t index_{};
    double rest_{},target_{},stiffness_{},maximum_force_{},sensed_{},integral_{},drive_{},force_{};
};
}
