// Geometry-only replay of planner cues through the native articulators and tract.
// Non-geometric cues are intentionally ignored; validate against a real speech probe.
#include "ducklips/body_config.hpp"
#include "ducklips/planner.hpp"
#include "ducklips/synthesizer.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <string>
int main(int argc,char** argv) {
 if(argc!=4){std::fprintf(stderr,"usage: ducklips_tongue_reach_audit BODY DICTIONARY TEXT > trajectory.csv\n");return 2;}
 ducklips::VoiceProfile voice=ducklips::baritone_voice();std::string error;
 if(!ducklips::load_body_config(argv[1],voice,error)){std::fprintf(stderr,"%s\n",error.c_str());return 2;}
 ducklips::SynthesisOptions options;ducklips::Synthesizer engine;engine.prepare(voice,options);
 ducklips::UtterancePlanner planner;planner.prepare(voice,engine.vowels(),options.rate);if(!planner.load_dictionary(argv[2])){std::fprintf(stderr,"dictionary has no definitions\n");return 2;}
 const std::string text=argv[3];
 const ducklips::UtterancePlan plan=planner.plan_text(text);
 if(plan.phones.empty()){std::fprintf(stderr,"text has no spoken phones\n");return 2;}
 const ducklips::SampleCount block=static_cast<ducklips::SampleCount>(options.control_block);
 const double block_seconds=ducklips::seconds_from_samples(options.rate,block);
 ducklips::ArticulatorState state;state.apply_speeds(voice.speeds);state.place_vowel(engine.vowels().shape(ducklips::Phone::AX));
 ducklips::TractModel model;model.prepare(voice,options.oral_sections,options.nasal_sections);model.set_conservative_body(true);
 ducklips::TractShape airway;model.allocate(airway);std::size_t cue=0,phone=0;
 std::puts("time,phone,body_activation,body_place,body_gap,body_width,tip_activation,tip_place,tip_gap,tip_down,posture_x,posture_y,jaw,blade,offset_x,offset_y,root_x,root_y,tip_x,tip_y,miss,overlap,drag_x,drag_y,attachment_motion");
 for(ducklips::SampleCount sample=0;sample<ducklips::samples_from_seconds(options.rate,plan.duration);sample+=block){
  while(cue<plan.cues.size() && plan.cues[cue].time<sample+block){
   const ducklips::Cue& c=plan.cues[cue++];
   if(c.kind==ducklips::CueKind::vowel){if(c.angular_rate<=0)state.place_vowel(c.vowel);else state.target_vowel(c.vowel,c.angular_rate);}
   if(c.kind==ducklips::CueKind::constriction)state.slot[static_cast<std::size_t>(c.slot)].begin(c.gesture,c.value,c.angular_rate);
   if(c.kind==ducklips::CueKind::release){ducklips::CriticallyDampedTracker& a=state.slot[static_cast<std::size_t>(c.slot)].activation;if(c.angular_rate<=0)a.place(c.value);else{a.target=c.value;a.angular_rate=c.angular_rate;}}
  }
  state.advance_geometry(block_seconds);model.build(state,airway);
  if(model.geometry_error()){fprintf(stderr,"%.6f %s\n",ducklips::seconds_from_samples(options.rate,sample),model.geometry_error());return 1;}
  const ducklips::TongueShape& s=model.tongue();
  const ducklips::Point2 pull=model.required_drag();state.body_offset_x.target=pull.x;state.body_offset_y.target=pull.y;
  const double time=ducklips::seconds_from_samples(options.rate,sample);while(phone+1<plan.phones.size() && plan.phones[phone+1].start<=time)++phone;
  const ducklips::ConstrictionSlotState& b=state.slot[2];const ducklips::ConstrictionSlotState& t=state.slot[1];
  const double bg=ducklips::MidsagittalTongue::gap_from_area(b.target_area.value*1e4,17.5-b.place.value*100);
  const double tg=ducklips::MidsagittalTongue::gap_from_area(t.target_area.value*1e4,17.5-t.place.value*100);
  printf("%.6f,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",time,ducklips::phone_spec(plan.phones[phone].phone).name,b.activation.value,b.place.value*100,bg,b.width.value*100,t.activation.value,t.place.value*100,tg,t.tip_down.value,state.shape[0].value,state.shape[1].value,state.shape[2].value,state.shape[3].value,state.body_offset_x.value,state.body_offset_y.value,s.blade_root.x,s.blade_root.y,s.tip.x,s.tip.y,s.tip_target_error,s.body_wall_penetration,pull.x,pull.y,s.attachment_motion);
 }
}
