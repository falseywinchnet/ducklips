// Geometry-only replay of planner cues through the native articulators and tract.
// Source targets may be exported separately; this tool does not render audio.
#include "ducklips/body_config.hpp"
#include "ducklips/planner.hpp"
#include "ducklips/synthesizer.hpp"
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <string>
int main(int argc,char** argv) {
 if(argc<4 || argc%2!=0){std::fprintf(stderr,"usage: ducklips_tongue_reach_audit BODY DICTIONARY TEXT [--partition SECTION] [--airway CSV] [--lip-profile CSV] [--source-cues JSON] > trajectory.csv\n");return 2;}
 std::size_t partition=0;
 std::ofstream airway_file;
 std::ofstream lip_file;
 std::ofstream source_file;
 for(int arg=4;arg<argc;arg+=2){
  if(std::string(argv[arg])=="--partition"){
   const std::string value=argv[arg+1];const auto parsed=std::from_chars(value.data(),value.data()+value.size(),partition);
   if(parsed.ec!=std::errc{} || parsed.ptr!=value.data()+value.size() || partition==0){std::fprintf(stderr,"invalid partition section\n");return 2;}
  }else if(std::string(argv[arg])=="--airway" && !airway_file.is_open()){
   airway_file.open(argv[arg+1]);if(!airway_file){std::fprintf(stderr,"cannot write airway CSV\n");return 2;}
  }else if(std::string(argv[arg])=="--lip-profile" && !lip_file.is_open()){
   lip_file.open(argv[arg+1]);if(!lip_file){return 2;}
  }else if(std::string(argv[arg])=="--source-cues" && !source_file.is_open()){
   source_file.open(argv[arg+1]);if(!source_file){return 2;}
  }else{std::fprintf(stderr,"unknown or repeated option\n");return 2;}
 }
 ducklips::VoiceProfile voice=ducklips::baritone_voice();std::string error;
 if(!ducklips::load_body_config(argv[1],voice,error)){std::fprintf(stderr,"%s\n",error.c_str());return 2;}
 ducklips::SynthesisOptions options;ducklips::Synthesizer engine;engine.prepare(voice,options);
 if(partition>=options.oral_sections-1){std::fprintf(stderr,"partition needs oral sections on both sides\n");return 2;}
 ducklips::UtterancePlanner planner;planner.prepare(voice,engine.vowels(),options.rate);if(!planner.load_dictionary(argv[2])){std::fprintf(stderr,"dictionary has no definitions\n");return 2;}
 const std::string text=argv[3];
 const ducklips::UtterancePlan plan=planner.plan_text(text);
 if(plan.phones.empty()){std::fprintf(stderr,"text has no spoken phones\n");return 2;}
 if(source_file.is_open()){
  // Export the planner's physical targets, not the older renderer's source
  // waveform, measured pitch, gain-shaped audio or smoothed state samples.
  source_file<<std::setprecision(17)<<"{\"duration_s\":"<<plan.duration
   <<",\"initial\":{\"glottal_rest_area_m2\":"<<voice.abducted_area
   <<",\"lung_pressure_Pa\":0,\"pitch_Hz\":"<<voice.base_pitch
   <<",\"effort\":1,\"loudness\":1,\"velum_area_m2\":0},\"cues\":[";
  bool first=true;
  for(const auto& c:plan.cues){
   const char* kind=nullptr;
   switch(c.kind){
    case ducklips::CueKind::glottis:kind="glottal_rest_area_m2";break;
    case ducklips::CueKind::lungs:kind="lung_pressure_Pa";break;
    case ducklips::CueKind::pitch:kind="pitch_Hz";break;
    case ducklips::CueKind::effort:kind="effort";break;
    case ducklips::CueKind::loudness:kind="loudness";break;
    case ducklips::CueKind::velum:kind="velum_area_m2";break;
    default:break;
   }
   if(!kind){continue;}
   source_file<<(first?"":",")<<"{\"time_s\":"<<ducklips::seconds_from_samples(options.rate,c.time)
    <<",\"kind\":\""<<kind<<"\",\"value\":"<<c.value<<",\"angular_rate\":"<<c.angular_rate<<'}';first=false;
  }
  source_file<<"],\"phones\":[";first=true;
  for(const auto& p:plan.phones){
   source_file<<(first?"":",")<<"{\"name\":\""<<ducklips::phone_spec(p.phone).name
    <<"\",\"stress\":"<<p.stress<<",\"start_s\":"<<p.start<<",\"end_s\":"<<p.end<<'}';first=false;
  }
  source_file<<"]}\n";
 }
 const ducklips::SampleCount block=static_cast<ducklips::SampleCount>(options.control_block);
 const double block_seconds=ducklips::seconds_from_samples(options.rate,block);
 ducklips::ArticulatorState state;state.apply_speeds(voice.speeds);state.place_vowel(engine.vowels().shape(ducklips::Phone::AX));
 ducklips::TractModel model;model.prepare(voice,options.oral_sections,options.nasal_sections);model.set_conservative_body(true);
 ducklips::TractShape airway;model.allocate(airway);std::size_t cue=0,phone=0;
 if(lip_file.is_open()){lip_file<<std::setprecision(17)<<"time,length,base_area,activation,place,core,width,target_area,sections\n";}
 if(airway_file.is_open()){
  airway_file<<std::setprecision(17)<<"time,narrowest_section,oral_volume";
  for(std::size_t i=0;i<airway.oral_area.size();++i){airway_file<<",area_"<<i<<",length_"<<i;}
  airway_file<<'\n';
 }
 std::printf("time,phone,body_activation,body_place,body_gap,body_width,tip_activation,tip_place,tip_gap,tip_down,posture_x,posture_y,jaw,blade,offset_x,offset_y,root_x,root_y,tip_x,tip_y,miss,overlap,drag_x,drag_y,attachment_motion,posture_root");
 if(partition){std::printf(",partition_section,upstream_volume,downstream_volume,neck_area,front_area,neck_length");}
 std::puts("");
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
  if(lip_file.is_open()){
   const auto p=model.lip_profile(state);
   lip_file<<time<<','<<p.length<<','<<p.base_area<<','<<p.activation<<','<<p.place<<','<<p.core<<','<<p.width<<','<<p.target_area<<','<<model.lip_sections()<<'\n';
  }
  if(airway_file.is_open()){
   airway_file<<time<<','<<model.narrowest_section()<<','<<model.oral_volume();
   for(std::size_t i=0;i<airway.oral_area.size();++i){airway_file<<','<<airway.oral_area[i]<<','<<airway.oral_length[i];}
   airway_file<<'\n';
  }
  const ducklips::ConstrictionSlotState& b=state.slot[2];const ducklips::ConstrictionSlotState& t=state.slot[1];
  const double bg=ducklips::MidsagittalTongue::gap_from_area(b.target_area.value*1e4,17.5-b.place.value*100);
  const double tg=ducklips::MidsagittalTongue::gap_from_area(t.target_area.value*1e4,17.5-t.place.value*100);
  printf("%.6f,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g",time,ducklips::phone_spec(plan.phones[phone].phone).name,b.activation.value,b.place.value*100,bg,b.width.value*100,t.activation.value,t.place.value*100,tg,t.tip_down.value,state.shape[0].value,state.shape[1].value,state.shape[2].value,state.shape[3].value,state.body_offset_x.value,state.body_offset_y.value,s.blade_root.x,s.blade_root.y,s.tip.x,s.tip.y,s.tip_target_error,s.body_wall_penetration,pull.x,pull.y,s.attachment_motion);
  printf(",%.12g",state.shape[ducklips::root_parameter].value);
  if(partition){
   // Fixed spatial cut at the centre of the requested section. This is a
   // diagnostic partition, not a production rule for finding air chambers.
   double upstream=0,downstream=0;
   for(std::size_t i=0;i<airway.oral_area.size();++i){
    const double volume=airway.oral_area[i]*airway.oral_length[i];
    if(i<partition)upstream+=volume;else if(i>partition)downstream+=volume;else{upstream+=volume/2;downstream+=volume/2;}
   }
   double front=airway.oral_area[partition+1];
   for(std::size_t i=partition+2;i<airway.oral_area.size();++i)front=std::min(front,airway.oral_area[i]);
   const double neck=std::min(airway.oral_area[partition-1],airway.oral_area[partition]);
   printf(",%zu,%.17g,%.17g,%.17g,%.17g,%.17g",partition,upstream,downstream,neck,front,.5*(airway.oral_length[partition-1]+airway.oral_length[partition]));
  }
  std::puts("");
 }
 if(airway_file.is_open() && !airway_file){std::fprintf(stderr,"failed writing airway CSV\n");return 1;}
 if(lip_file.is_open() && !lip_file){return 1;}
}
