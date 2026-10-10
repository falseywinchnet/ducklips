#pragma once
#include "ducklips/tissue_body_data.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>

namespace ducklips {
// Strain-history relaxation in material coordinates. Kmaterial = L L^T.
// Each positive stochastic rate dissipates stored branch energy at fixed
// deformation. No random force, waveform, phase or output gain is supplied.
// The noise statistics are an explicit hypothesis, not measured physiology.
class TissueRelaxation final {
public:
    static constexpr std::size_t size=compiled_tissue::modes;
    using Coordinates=std::array<double,size>;
    void reset(double fraction=.1,double seconds=.1,double log_rate_sd=.6,std::uint64_t seed=1){
        if(!std::isfinite(fraction)||fraction<0||fraction>=1||!std::isfinite(seconds)||seconds<=0
            ||!std::isfinite(log_rate_sd)||log_rate_sd<0){throw std::invalid_argument("invalid tissue relaxation");}
        fraction_=fraction;seconds_=seconds;sd_=log_rate_sd;history_.fill(0);clock_=0;
        random_.seed(seed);normal_.reset();
        for(std::size_t i=0;i<size;++i){log_rate_[i]=sd_*normal_(random_);update_rate(i);}
    }
    [[nodiscard]] Coordinates strain(const Coordinates& deformation) const noexcept{
        Coordinates value{};
        for(std::size_t j=0;j<size;++j){
            for(std::size_t i=j;i<size;++i){value[j]+=compiled_tissue::material_factor[i*size+j]*deformation[i];}
        }
        return value;
    }
    [[nodiscard]] Coordinates relieved_force() const noexcept{
        Coordinates force{};
        for(std::size_t i=0;i<size;++i){
            for(std::size_t j=0;j<=i;++j){force[i]+=fraction_*compiled_tissue::material_factor[i*size+j]*history_[j];}
        }
        return force;
    }
    void advance(const Coordinates& deformation,double dt){
        const Coordinates current=strain(deformation);
        // The slow stochastic process runs at 1 kHz; stress history remains
        // continuous at the tissue time step. This avoids sample-rate-sized
        // random kicks and repeated audio-rate random-number generation.
        double remaining=dt;
        while(remaining>0){
            const double step=std::min(remaining,.001-clock_);
            for(std::size_t i=0;i<size;++i){
                history_[i]+=-std::expm1(-step*rate_[i])*(current[i]-history_[i]);
            }
            clock_+=step;remaining-=step;
            if(clock_>=.001-1e-14){
                const double retain=std::exp(-.001/seconds_);
                const double innovation=sd_*std::sqrt(1-retain*retain);
                for(std::size_t i=0;i<size;++i){log_rate_[i]=retain*log_rate_[i]+innovation*normal_(random_);update_rate(i);}
                clock_=0;
            }
        }
    }
    // Correction to e^T Kinstant e, including both folds. The remaining
    // material has (1-fraction) stiffness; prestress stiffness is untouched.
    [[nodiscard]] double energy_correction(const Coordinates& deformation) const noexcept{
        const Coordinates current=strain(deformation);double result=0;
        for(std::size_t i=0;i<size;++i){
            const double offset=current[i]-history_[i];
            result+=fraction_*(offset*offset-current[i]*current[i]);
        }
        return result;
    }
    [[nodiscard]] const Coordinates& rates() const noexcept{return rate_;}
private:
    void update_rate(std::size_t i){rate_[i]=std::exp(log_rate_[i]-.5*sd_*sd_)/seconds_;}
    Coordinates history_{},log_rate_{},rate_{};
    double fraction_{},seconds_{.1},sd_{},clock_{};
    std::mt19937_64 random_{1};std::normal_distribution<double> normal_{};
};
}
