#include "ducklips/somatic_control.hpp"
#include "ducklips/synthesizer.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>

void require(bool condition,const char* message){if(!condition){throw std::runtime_error(message);}}
int main(){
    ducklips::ArticulatorState body;
    body.place_vowel(ducklips::VowelShape{});body.pitch.place(120);body.effort.place(1);
    ducklips::SomaticControl regulated,uncorrected;
    auto reference=body;
    ducklips::SomaticControl::Settings settings;settings.noise=0;
    regulated.reset(body,17,settings);settings.feedback=0;uncorrected.reset(reference,17,settings);
    regulated.set_load(0,.08);uncorrected.set_load(0,.08);
    const double desired=body.shape[0].value;
    for(int i=0;i<4000;++i){regulated.advance(body,.001,0,0);uncorrected.advance(reference,.001,0,0);}
    const double corrected_error=std::abs(body.shape[0].value-desired),open_error=std::abs(reference.shape[0].value-desired);
    require(corrected_error<.3*open_error,"delayed adaptation must reduce a sustained load error");
    // Replaying the same body/seed gives identical trajectories; unrelated
    // seeds must change the holding effort without changing the intention.
    ducklips::SomaticControl a,b,c;
    auto ba=body,bb=body,bc=body;settings.feedback=1;settings.noise=.08;
    a.reset(ba,3,settings);b.reset(bb,3,settings);c.reset(bc,4,settings);
    ba.shape[0].target+=.2;bb.shape[0].target+=.2;bc.shape[0].target+=.2;
    double differing=0;
    for(int i=0;i<1000;++i){
        a.advance(ba,.001,0,1e-4);b.advance(bb,.001,0,1e-4);c.advance(bc,.001,0,1e-4);
        require(ba.shape[0].value==bb.shape[0].value,"seeded control must replay exactly");
        differing+=std::abs(ba.shape[0].value-bc.shape[0].value);
        for(const auto& observation:a.observations()){require(std::isfinite(observation.position)&&std::isfinite(observation.muscle),"all controls must remain finite");}
    }
    require(differing>1e-6,"independent motor histories must affect motion");
    require(std::abs(a.remaining_breath()-(.003-.0001))<1e-12,"exhaled flow must deplete the finite breath supply");
    require(ba.shape[5].value==body.shape[5].value,"fixed anatomy must not wander");
    ducklips::ArticulatorState neutral,noiseless;
    neutral.place_vowel(ducklips::VowelShape{});neutral.pitch.place(120);neutral.effort.place(1);
    noiseless=neutral;
    ducklips::SomaticControl held,quiet;
    held.reset(neutral,31,settings);settings.noise=0;quiet.reset(noiseless,31,settings);
    double minimum_pitch=120,maximum_pitch=120;
    for(int i=0;i<2000;++i){
        held.advance(neutral,.001,0,0);quiet.advance(noiseless,.001,0,0);
        minimum_pitch=std::min(minimum_pitch,neutral.pitch.value);maximum_pitch=std::max(maximum_pitch,neutral.pitch.value);
        require(noiseless.pitch.value==120,"balanced tone without fluctuation must maintain neutral posture");
    }
    require(maximum_pitch-minimum_pitch>1e-4 && maximum_pitch-minimum_pitch<2,"neutral holding must have small nonzero motor variation");
    // A client's audio chunk boundaries cannot reset or alter neural history.
    ducklips::SynthesisOptions options;options.local_flow=true;options.living_control=true;
    options.conservative_tongue=true;options.oral_sections=24;options.nasal_sections=9;options.local_substeps=1;
    auto voice=ducklips::baritone_voice();voice.aspiration_gain=voice.frication_gain=0;
    ducklips::Synthesizer whole,chunked;whole.prepare(voice,options);chunked.prepare(voice,options);
    ducklips::Cue breath;breath.kind=ducklips::CueKind::lungs;breath.value=600;breath.angular_rate=60;
    require(whole.enqueue(breath)&&chunked.enqueue(breath),"breath cue must fit");
    ducklips::Cue contact;contact.kind=ducklips::CueKind::constriction;contact.slot=ducklips::ConstrictionSlot::tongue_body;
    contact.gesture.place_from_lips=.06;contact.gesture.width=.009;contact.gesture.target_area=1e-4;
    contact.value=.15;contact.angular_rate=60;
    require(whole.enqueue(contact)&&chunked.enqueue(contact),"new contact cue must fit");
    std::vector<float> contiguous(19200),pieces(19200);whole.render(contiguous);
    for(std::size_t offset=0;offset<pieces.size();){
        const std::size_t size=std::min(std::size_t(146),pieces.size()-offset);
        chunked.render(std::span<float>(pieces).subspan(offset,size));offset+=size;
    }
    if(whole.statistics().airflow_failures||whole.statistics().geometry_failures){
        std::cerr<<"integration failure: air="<<whole.statistics().airflow_failures<<" at "<<whole.statistics().airflow_failure_time
            <<" geometry="<<whole.statistics().geometry_failures<<" at "<<whole.statistics().geometry_failure_time
            <<" reason="<<(whole.statistics().geometry_error?whole.statistics().geometry_error:"unknown")<<'\n';
    }
    require(whole.statistics().airflow_failures==0 && whole.statistics().geometry_failures==0,"integrated control must render without a physical failure");
    require(contiguous==pieces,"stream chunk boundaries must preserve identical audio");
    require(*std::max_element(contiguous.begin(),contiguous.end())>0,"integrated control must produce nonzero audio");
    ducklips::Synthesizer initial;initial.prepare(voice,options);
    ducklips::Cue seal;seal.kind=ducklips::CueKind::constriction;seal.slot=ducklips::ConstrictionSlot::lips;
    seal.gesture.width=.005;seal.value=1.25;seal.angular_rate=0;
    require(initial.enqueue(seal),"initial lip seal must fit");
    std::array<float,96> first_block{};initial.render(first_block);
    require(initial.somatic_control().observations()[7].position>1.2,"time-zero preparatory lip closure must be established before sound generation");
    options.somatic.tone=-1;bool rejected=false;
    try{whole.prepare(voice,options);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"invalid control settings must be rejected before noexcept reset");
    std::cout<<"29-channel control: load error "<<open_error<<" -> "<<corrected_error
        <<"; neutral pitch range "<<minimum_pitch<<".."<<maximum_pitch
        <<"; seed replay, independent motor histories, finite breath, fixed anatomy and streaming passed\n";
}
