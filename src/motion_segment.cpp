#include "ducklips/motion_segment.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace ducklips {
void CanonicalMotionFrame::prepare(std::size_t modes) {
    if (modes==0 || modes>4096) { throw std::invalid_argument("invalid motion frame size"); }
    angular_frequency.resize(modes);flow_diagonal.resize(modes);pressure_transport.resize(modes*modes);
}
ModalMotionSample CanonicalMotionFrame::view() const noexcept { return {angular_frequency,pressure_transport,{},flow_diagonal}; }
void CanonicalMotionSegment::prepare(std::size_t modes) {
    if (modes==0 || modes>4096) { throw std::invalid_argument("invalid motion segment size"); }
    modes_=modes;configured_=false;
    pressure_.resize(4*modes*modes);next_pressure_.resize(pressure_.size());
    frequency_.resize(6*modes);next_frequency_.resize(frequency_.size());
    evaluated_frequency_.resize(modes);evaluated_diagonal_.resize(modes);
}
bool CanonicalMotionSegment::configure(const std::array<CanonicalMotionKnot,6>& knots,double interval) noexcept {
    if (modes_==0 || !std::isfinite(interval) || interval<=0) { return false; }
    for (const CanonicalMotionKnot& knot:knots) {
        if (knot.angular_frequency.size()!=modes_ || knot.frequency_rate.size()!=modes_ || knot.pressure_transport.size()!=modes_*modes_) { return false; }
        for (double value:knot.angular_frequency) { if (!std::isfinite(value) || value<=0) { return false; } }
        for (double value:knot.frequency_rate) { if (!std::isfinite(value)) { return false; } }
        for (double value:knot.pressure_transport) { if (!std::isfinite(value)) { return false; } }
    }
    for (std::size_t i=0;i<modes_;++i) {
        const double a=knots[1].angular_frequency[i],b=interval*knots[1].frequency_rate[i];
        const double c=knots[2].angular_frequency[i],d=interval*knots[2].frequency_rate[i];
        const double e=knots[3].angular_frequency[i],f=interval*knots[3].frequency_rate[i];
        next_frequency_[i]=c;next_frequency_[modes_+i]=d;
        next_frequency_[2*modes_+i]=a+.25*b-2*c+e-.25*f;
        next_frequency_[3*modes_+i]=-1.25*a-.25*b-2*d+1.25*e-.25*f;
        next_frequency_[4*modes_+i]=-.5*a-.25*b+c-.5*e+.25*f;
        next_frequency_[5*modes_+i]=.75*a+.25*b+d-.75*e+.25*f;
    }
    const std::size_t count=modes_*modes_;
    for (std::size_t i=0;i<count;++i) {
        const double previous2=knots[0].pressure_transport[i],previous=knots[1].pressure_transport[i];
        const double first=knots[2].pressure_transport[i],last=knots[3].pressure_transport[i];
        const double next=knots[4].pressure_transport[i],next2=knots[5].pressure_transport[i];
        const double first_rate=(previous2-8*previous+8*last-next)/12;
        const double last_rate=(previous-8*first+8*next-next2)/12;
        next_pressure_[i]=first;next_pressure_[count+i]=first_rate;
        next_pressure_[2*count+i]=3*(last-first)-2*first_rate-last_rate;
        next_pressure_[3*count+i]=2*(first-last)+first_rate+last_rate;
        // Certify the entire [0,1] polynomial once. Every Horner partial sum
        // is bounded by this absolute coefficient sum. Keep rounding margin.
        const double bound=std::abs(next_pressure_[i])+std::abs(next_pressure_[count+i])+
                           std::abs(next_pressure_[2*count+i])+std::abs(next_pressure_[3*count+i]);
        if (!std::isfinite(bound) || bound>std::numeric_limits<double>::max()/2) { return false; }
    }
    for (double v:next_frequency_) { if (!std::isfinite(v)) { return false; } }
    pressure_.swap(next_pressure_);frequency_.swap(next_frequency_);interval_=interval;configured_=true;
    return true;
}
bool CanonicalMotionSegment::evaluate(double fraction,CanonicalMotionFrame& frame) noexcept {
    const std::size_t count=modes_*modes_;
    if (!configured_ || !std::isfinite(fraction) || fraction<0 || fraction>1 || frame.angular_frequency.size()!=modes_ ||
        frame.flow_diagonal.size()!=modes_ || frame.pressure_transport.size()!=count) { return false; }
    for (std::size_t i=0;i<modes_;++i) {
        double omega=frequency_[5*modes_+i],derivative=0;
        for (int degree=4;degree>=0;--degree) {
            derivative=derivative*fraction+omega;
            omega=omega*fraction+frequency_[static_cast<std::size_t>(degree)*modes_+i];
        }
        const double diagonal=derivative/interval_/omega;
        if (!std::isfinite(omega) || omega<=0 || !std::isfinite(diagonal)) { return false; }
        evaluated_frequency_[i]=omega;evaluated_diagonal_[i]=diagonal;
    }
    for (std::size_t i=0;i<count;++i) {
        const double value=((pressure_[3*count+i]*fraction+pressure_[2*count+i])*fraction+pressure_[count+i])*fraction+pressure_[i];
        frame.pressure_transport[i]=value;
    }
    std::copy(evaluated_frequency_.begin(),evaluated_frequency_.end(),frame.angular_frequency.begin());
    std::copy(evaluated_diagonal_.begin(),evaluated_diagonal_.end(),frame.flow_diagonal.begin());
    return true;
}
}
