#include "ducklips/somatic_control.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
constexpr std::array<const char*,SomaticControl::count> names{
    "tongue_x","tongue_y","jaw","blade","root","lip_area","lip_protrusion",
    "lips_activation","lips_place","lips_area","lips_width","lips_tip",
    "tip_activation","tip_place","tip_area","tip_width","tip_down",
    "body_activation","body_place","body_area","body_width","body_tip",
    "tongue_drag_x","tongue_drag_y","velum","glottal_area","breath_pressure","laryngeal_pitch_intent","effort"};
struct Range {double scale,low,high;std::size_t group;};
constexpr std::array<Range,SomaticControl::count> ranges{{
    {2,0,8,0},{2,2,9,0},{1,0,1,1},{1,-1,1,0},{1,-1,1,0},
    {3e-4,0,.002,1},{.02,-.005,.04,1},
    {1,0,1.5,1},{.03,0,.175,1},{1e-4,0,.002,1},{.01,.0001,.06,1},{1,0,1,1},
    {1,0,1.5,0},{.03,0,.175,0},{1e-4,0,.002,0},{.01,.0001,.06,0},{1,0,1,0},
    {1,0,1.5,0},{.03,0,.175,0},{1e-4,0,.002,0},{.01,.0001,.06,0},{1,0,1,0},
    {.5,-2,2,0},{.5,-2,2,0},{1e-4,0,.001,2},
    {1e-5,-1e-5,1e-4,3},{600,0,2500,4},{120,30,500,3},{1,.3,2,4}}};
}
const char* SomaticControl::name(std::size_t i) noexcept{return names[i];}
std::array<CriticallyDampedTracker*,SomaticControl::count> SomaticControl::coordinates(ArticulatorState& b) noexcept{
    std::array<CriticallyDampedTracker*,count> p{};std::size_t n=0;
    for(std::size_t i=0;i<5;++i){p[n++]=&b.shape[i];}
    p[n++]=&b.lip_area;p[n++]=&b.lip_protrusion;
    for(auto& slot:b.slot){p[n++]=&slot.activation;p[n++]=&slot.place;p[n++]=&slot.target_area;p[n++]=&slot.width;p[n++]=&slot.tip_down;}
    p[n++]=&b.body_offset_x;p[n++]=&b.body_offset_y;p[n++]=&b.port_area;
    p[n++]=&b.glottal_rest_area;p[n++]=&b.lung_pressure;p[n++]=&b.pitch;p[n++]=&b.effort;
    return p;
}
double SomaticControl::noise(std::uint64_t& state) noexcept{
    state^=state>>12;state^=state<<25;state^=state>>27;
    const double uniform=double((state*2685821657736338717ULL)>>11)*0x1.0p-53;
    return std::sqrt(12.)*(uniform-.5);
}
double SomaticControl::drive(double activation) noexcept{
    // Shared inverse recruitment table is constructed during reset, not in
    // the audio callback. Avoid solving 29 identical recruitment maps per ms.
    static const std::array<double,257> table=[]{
        std::array<double,257> t{};MotorPool pool;
        for(std::size_t i=0;i<t.size();++i){t[i]=pool.drive_for_activation(double(i)/256);}
        return t;
    }();
    const double x=std::clamp(activation,0.,1.)*256;
    const auto i=std::min(std::size_t(x),std::size_t(255));
    return std::lerp(table[i],table[i+1],x-double(i));
}
void SomaticControl::validate(const Settings& settings){
    if(!std::isfinite(settings.noise)||settings.noise<0||settings.noise>1
        ||!std::isfinite(settings.tone)||settings.tone<0||settings.tone>.25
        ||!std::isfinite(settings.feedback)||settings.feedback<0||settings.feedback>2
        ||!std::isfinite(settings.relaxation)||settings.relaxation<0||settings.relaxation>=1
        ||!std::isfinite(settings.respiratory_compliance)||settings.respiratory_compliance<=0
        ||!std::isfinite(settings.breath_volume)||settings.breath_volume<=0){throw std::invalid_argument("invalid somatic control settings");}
}
void SomaticControl::reset(ArticulatorState& body,std::uint64_t seed,const Settings& settings){
    validate(settings);
    settings_=settings;random_=seed?seed:1;common_.fill(0);loads_.fill(0);breath_=settings.breath_volume;
    const auto p=coordinates(body);(void)drive(.5);
    for(std::size_t i=0;i<count;++i){
        auto& s=state_[i];s=Coordinate{};s.origin=p[i]->value;s.scale=ranges[i].scale;
        s.low=ranges[i].low;s.high=ranges[i].high;s.group=ranges[i].group;
        s.previous=s.current=s.origin;s.random=seed+1009*(i+1);if(s.random==0){s.random=1;}
        s.positive.reset(s.random);s.negative.reset(s.random+503);s.delayed.fill(0);
        observations_[i]={.target=p[i]->target,.position=s.origin};
    }
}
void SomaticControl::tick(ArticulatorState& body,double pressure,double flow){
    constexpr double dt=.001;
    const auto p=coordinates(body);
    const double common_decay=std::exp(-dt/.15);
    for(double& value:common_){value=common_decay*value+std::sqrt(1-common_decay*common_decay)*noise(random_);}
    breath_=std::clamp(breath_-flow*dt,0.,settings_.breath_volume);
    for(std::size_t i=0;i<count;++i){
        auto& s=state_[i];auto& tracker=*p[i];
        double target=(std::clamp(tracker.target,s.low,s.high)-s.origin)/s.scale;
        double rate=tracker.angular_rate>0?tracker.angular_rate:40;
        if(tracker.speed_limit>0){rate=std::min(rate,tracker.speed_limit);}
        s.rate=std::clamp(rate,8.,180.);
        const double local_decay=std::exp(-dt/.12);
        s.rate_noise=local_decay*s.rate_noise+std::sqrt(1-local_decay*local_decay)*noise(s.random);
        const double relaxation_rate=std::exp(.12*settings_.noise*(s.rate_noise+common_[s.group]))/.2;
        s.memory+=-std::expm1(-dt*relaxation_rate)*(s.q-s.memory);
        s.sense+=-std::expm1(-dt/.02)*(s.q-s.sense);
        const double delayed=s.delayed[s.cursor];s.delayed[s.cursor]=s.sense;s.cursor=(s.cursor+1)%s.delayed.size();
        const double error=target-delayed;
        const double integral=s.integral+settings_.feedback*.15*error*dt/.2;
        double requested=(1-settings_.relaxation)*target+settings_.feedback*.15*error+integral;
        constexpr double capacity=8;
        if(std::abs(requested)<capacity || requested*error<0){s.integral=std::clamp(integral,-2.,2.);}
        requested=std::clamp(requested,-capacity*(1-settings_.tone),capacity*(1-settings_.tone));
        // Shared neural fluctuation slightly changes effort; each coordinate
        // also has independent motor-unit firing and relaxing strain memory.
        requested*=1+.002*settings_.noise*common_[s.group];
        const double activation=std::abs(requested)/capacity;
        // Antagonists maintain a small opposing tone even at neutral intent.
        // Their means cancel; independent overlapping twitches do not.
        // Slot activation is a gesture recruitment signal, not resting muscle
        // tone. No selected contact means no tonic recruitment of that slot.
        // The underlying tongue, jaw and lip posture still has opposing tone.
        const bool gesture_activation=i==7||i==12||i==17;
        const double tone=settings_.tone*(gesture_activation?std::clamp(tracker.target,0.,1.):1.);
        const double pd=drive(tone+(requested>0?activation:0));
        const double nd=drive(tone+(requested<0?activation:0));
        const double pool=capacity*(s.positive.advance(pd,dt)-s.negative.advance(nd,dt));
        // Mean-field recruitment carries the unresolved population; the
        // finite pool contributes a configurable fraction of fluctuation.
        // Both evolve twitch force. No random position/output offset exists.
        const double mean=capacity*(s.positive.expected_activation(pd)-s.negative.expected_activation(nd));
        const double decay=std::exp(-dt/.017);
        s.twitch=requested+decay*((s.twitch-requested)+dt/.017*(s.twitch_input-requested));
        s.twitch_input=requested+(s.twitch_input-requested)*decay;
        s.mean_twitch=mean+decay*((s.mean_twitch-mean)+dt/.017*(s.mean_input-mean));
        s.mean_input=mean+(s.mean_input-mean)*decay;
        const double muscle=s.twitch+settings_.noise*(pool-s.mean_twitch);
        double load=loads_[i];
        // Small effective pressure-load projections. These adapters need
        // anatomical identification; units here are normalized coordinates.
        const double normalized_pressure=pressure/600;
        if(i==5){load+=.012*normalized_pressure;}
        if(i==2){load+=.004*normalized_pressure;}
        if(i==0||i==1||i==24){load+=.002*normalized_pressure;}
        // Passive material + overlapping active effort + explicit load.
        // Exact held-force damped update keeps high cue rates bounded.
        const double equilibrium=muscle+settings_.relaxation*s.memory+load;
        const double displacement=s.q-equilibrium;
        const double slope=s.v+s.rate*displacement;
        const double body_decay=std::exp(-s.rate*dt);
        s.q=equilibrium+body_decay*(displacement+dt*slope);
        s.v=body_decay*(s.v-s.rate*dt*slope);
        if(i==26){s.q-=flow*dt/(settings_.respiratory_compliance*s.scale);}
        double position=s.origin+s.scale*s.q;
        const double upper=i==26?std::min(s.high,breath_/settings_.respiratory_compliance):s.high;
        if(position<s.low||position>upper){position=std::clamp(position,s.low,upper);s.q=(position-s.origin)/s.scale;s.v=0;}
        s.current=position;tracker.value=position;tracker.velocity=s.scale*s.v;tracker.command=position;
        observations_[i]={tracker.target,position,tracker.velocity,s.memory,muscle,s.origin+s.scale*s.sense,s.integral,load};
    }
}
void SomaticControl::relocate_idle_gesture(ArticulatorState& body,std::size_t slot) noexcept{
    const auto p=coordinates(body);
    for(std::size_t i=8+5*slot;i<12+5*slot;++i){
        auto& s=state_[i];
        // begin() placed a new idle reference in the tracker. Its previous
        // commanded reference is the origin plus the normalized target kept
        // by the last observation. Translate, retaining q/v/stress/twitches.
        const double delta=p[i]->target-observations_[i].target;
        s.origin+=delta;s.previous+=delta;s.current+=delta;
        p[i]->value=s.current;p[i]->velocity=s.scale*s.v;p[i]->command=s.current;
        observations_[i].target=p[i]->target;observations_[i].position+=delta;observations_[i].sensed+=delta;
    }
}
void SomaticControl::advance(ArticulatorState& body,double seconds,double pressure,double flow){
    const auto steps=std::llround(seconds/.001);
    if(steps<1||std::abs(seconds-double(steps)*.001)>1e-9){throw std::invalid_argument("somatic control requires whole millisecond blocks");}
    for(auto& s:state_){s.previous=s.current;}
    for(long long i=0;i<steps;++i){tick(body,pressure,flow);}
    // Anatomical residual coefficients and microphone gain are not muscles.
    for(std::size_t i=5;i<shape_parameter_count;++i){body.shape[i].advance(seconds);}
}
void SomaticControl::physiology(ArticulatorState& body,double fraction) const noexcept{
    const auto p=coordinates(body);
    for(std::size_t i=25;i<count;++i){p[i]->value=std::lerp(state_[i].previous,state_[i].current,std::clamp(fraction,0.,1.));}
}
}
