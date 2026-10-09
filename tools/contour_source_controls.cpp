#include "ducklips/contour_glottal_storage.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
namespace {void require(bool v,const char* message){if(!v){throw std::runtime_error(message);}}}
int main(){try{
    ducklips::FoldParameters f;f.asymmetry=1.17;ducklips::AirProperties air;ducklips::FoldState initial;initial.rest_area=3e-6;
    initial.lower_position={1e-6,-2e-6};initial.upper_position={-3e-6,4e-6};initial.lower_velocity={.001,-.002};initial.upper_velocity={.003,-.001};
    ducklips::ContourGlottalStorage source;require(source.prepare(f,air,.9,initial,{500.,100.}),"asymmetric prepare");source.momentum={.001,-.0005,.0004};
    const std::array<double,4> force={1e-4,-2e-4,3e-5,-4e-5};double e0=0;require(source.energy(e0),"initial energy");
    double work=0,loss=0,max_balance=0,max_force=0;const double dt=2e-6;
    for(int i=0;i<1000;++i){const double rest=3e-6+2e-7*std::sin((i+1)*dt*1000);
        require(source.advance(dt,300.,20.,1e6,rest,force),"asymmetric advance");
        work+=source.last.reservoir_work-source.last.load_work+source.last.actuator_work+source.last.external_work;loss+=source.last.damping_loss+source.last.passage_loss;
        double e=0;require(source.energy(e),"energy");max_balance=std::max(max_balance,std::abs(e-e0+loss-work));
        for(double r:source.last.force_residual){max_force=std::max(max_force,std::abs(r));}}
    require(max_balance<1e-13 && max_force<1e-9,"asymmetric external work");
    const auto saved=source;auto unchanged=[&](){return source.state.mass==saved.state.mass && source.momentum==saved.momentum && source.state.folds.lower_position==saved.state.folds.lower_position
        && source.state.folds.upper_position==saved.state.folds.upper_position && source.state.folds.lower_velocity==saved.state.folds.lower_velocity && source.state.folds.upper_velocity==saved.state.folds.upper_velocity
        && source.state.folds.rest_area==saved.state.folds.rest_area && source.last.energy_residual==saved.last.energy_residual;};
    const double nan=std::numeric_limits<double>::quiet_NaN();
    require(!source.advance(0,300,0,0,3e-6) && unchanged(),"zero time rejection");
    require(!source.advance(dt,300,0,-1,3e-6) && unchanged(),"negative load rejection");
    require(!source.advance(dt,300,0,0,3e-6,{nan,0,0,0}) && unchanged(),"force rejection");
    ducklips::FoldState sealed;sealed.rest_area=0;require(!source.prepare(f,air,1.,sealed,{}) && unchanged(),"contact preparation rejection");
    // No-drive, no-motion equilibrium remains an exact accepted state.
    initial={};initial.rest_area=3e-6;require(source.prepare(f,air,1.,initial,{}),"equilibrium prepare");const auto equilibrium=source.state;
    require(source.advance(dt,0,0,0,3e-6),"equilibrium advance");require(source.state.mass==equilibrium.mass && source.state.folds.lower_position==equilibrium.folds.lower_position && source.momentum==std::array<double,3>{},"equilibrium moved");
    std::printf("PASS asymmetric force/rest drive: work %.17g J force %.17g N; rejection and equilibrium\n",max_balance,max_force);return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
