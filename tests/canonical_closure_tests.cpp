#include "ducklips/canonical_speech.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

int main(){
    const ducklips::AirProperties air{};
    const ducklips::NetworkLosses losses{};
    const double dt=1./48000;
    ducklips::TractShape shape{};
    shape.oral_area.assign(4,1e-4);shape.oral_length.assign(4,.04);
    shape.nasal_area.assign(2,1e-4);shape.nasal_length.assign(2,.04);
    ducklips::AcousticNetwork closed{},open{};
    closed.prepare(4,2,1,{});open.prepare(4,2,1,{});
    open.assemble(shape,air,losses);
    shape.oral_area[1]=0;closed.assemble(shape,air,losses);
    ducklips::CanonicalSpeech engine{};engine.prepare(closed,air);
    // A pressure source across absent air moves no air and stores no energy.
    // Turn the source off before reopening. No pressure or flow remains to
    // drive a release, and the reopening geometry alone must stay silent.
    for(std::size_t i=0;i<2400;++i){
        engine.free_step(closed,closed,0,0,dt,8,0,0,0,500);engine.inject(0);
    }
    const double closed_energy=engine.energy();double peak_energy=0,peak_radiation=0;
    for(std::size_t i=0;i<480;++i){
        if(i<48){engine.free_step(closed,open,double(i)/48,1./48,dt,8,0,0,0,0);}
        else{engine.free_step(open,open,0,0,dt,8,0,0,0,0);}
        engine.inject(0);
        peak_energy=std::max(peak_energy,engine.energy());
        peak_radiation=std::max(peak_radiation,std::abs(engine.flow_derivative(open.lip_branch(),dt)));
    }
    const bool isolated=closed_energy==0 && peak_energy==0 && peak_radiation==0;
    std::printf("%s no latent source history through sealed passage: closed_energy %.12g reopened_energy %.12g lip_acceleration %.12g\n",
                isolated?"PASS":"FAIL",closed_energy,peak_energy,peak_radiation);
    // The same source must drive an existing air path; closure handling must
    // not erase all acoustic history or suppress pressure sources generally.
    engine.reset();
    for(std::size_t i=0;i<48;++i){engine.free_step(open,open,0,0,dt,8,0,0,0,500);engine.inject(0);}
    const double driven_energy=engine.energy();
    engine.free_step(open,open,0,0,dt,8,0,0,0,0);engine.inject(0);
    const bool retained=std::isfinite(driven_energy) && driven_energy>0 && engine.energy()>0;
    std::printf("%s open passage retains driven acoustic state: energy %.12g\n",retained?"PASS":"FAIL",driven_energy);
    return isolated && retained?0:1;
}
