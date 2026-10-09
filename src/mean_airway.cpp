#include "ducklips/mean_airway.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
constexpr double closed_area=1e-12;
// Bernoulli acceleration from first to second area, with no recovery on
// expansion. Expressed relative to a reference area to avoid huge 1/A^2.
double contraction(double first,double second,double reference) noexcept {
    if(second>=first){return 0;}
    const double a=reference/second,b=reference/first;return a*a-b*b;
}
}
void MeanAirway::prepare(const AcousticNetwork& layout,const TractShape& shape,
                        const AirProperties& air,double oral_wall_compliance){
    if(!std::isfinite(oral_wall_compliance) || oral_wall_compliance<=0){throw std::invalid_argument("mean airway needs positive wall compliance");}
    air_=air;oral_=shape.oral_area.size();nasal_=shape.nasal_area.size();
    const std::size_t n=layout.node_count();
    old_area_.resize(n);new_area_.resize(n);area_.resize(n);
    old_length_.resize(n);new_length_.resize(n);length_.resize(n);
    volume_.resize(n);compliance_.resize(n);wall_compliance_.assign(n,0);
    inflow_.assign(n,0);initial_pressure_.assign(n,0);oral_ports_.resize(oral_);
    links_.clear();links_.push_back({0,cavity_reservoir,false,false});
    std::vector<CavityConnection> connections{{0,cavity_reservoir}};
    for(const AcousticBranch& branch:layout.branches){
        if(branch.kind==BranchKind::wall){continue;}
        const std::size_t index=links_.size();
        links_.push_back({branch.from,branch.to,branch.kind==BranchKind::velar_port,
                         branch.kind==BranchKind::side && branch.from<oral_+nasal_});
        connections.push_back({branch.from,branch.to});
        if(branch.kind==BranchKind::tract){oral_ports_[branch.from]=index;}
        if(branch.kind==BranchKind::lip_radiation){lip_port_=index;oral_ports_.back()=index;}
        if(branch.kind==BranchKind::nostril_radiation){nasal_port_=index;}
    }
    passages_.resize(links_.size());network_.prepare(n,connections);
    double oral_length=0;for(double value:shape.oral_length){oral_length+=value;}
    for(std::size_t i=0;i<oral_;++i){wall_compliance_[i]=oral_wall_compliance*shape.oral_length[i]/oral_length;}
    reset();target(shape);
}
void MeanAirway::reset() noexcept {ready_=false;}
void MeanAirway::copy_shape(const TractShape& shape,std::vector<double>& area,std::vector<double>& length) noexcept {
    for(std::size_t i=0;i<oral_;++i){area[i]=shape.oral_area[i];length[i]=shape.oral_length[i];}
    for(std::size_t i=0;i<nasal_;++i){area[oral_+i]=shape.nasal_area[i];length[oral_+i]=shape.nasal_length[i];}
    for(std::size_t i=0;i<shape.side_area.size();++i){area[oral_+nasal_+i]=shape.side_area[i];length[oral_+nasal_+i]=shape.side_length[i];}
}
void MeanAirway::target(const TractShape& shape) noexcept {
    old_area_=new_area_;old_length_=new_length_;old_port_area_=new_port_area_;old_port_length_=new_port_length_;
    copy_shape(shape,new_area_,new_length_);new_port_area_=shape.port_area;new_port_length_=shape.port_length;
    if(!ready_){
        old_area_=new_area_;old_length_=new_length_;old_port_area_=new_port_area_;old_port_length_=new_port_length_;
        build(0);ready_=network_.reset(volume_,compliance_,initial_pressure_);
    }
}
void MeanAirway::build(double fraction) noexcept {
    const double bulk=air_.density*air_.sound_speed*air_.sound_speed;
    for(std::size_t i=0;i<area_.size();++i){
        area_[i]=old_area_[i]+fraction*(new_area_[i]-old_area_[i]);
        length_[i]=old_length_[i]+fraction*(new_length_[i]-old_length_[i]);
        const double first=old_area_[i]*old_length_[i],last=new_area_[i]*new_length_[i];
        volume_[i]=first+fraction*(last-first);
        compliance_[i]=volume_[i]/bulk+wall_compliance_[i];
    }
    const double port_area=old_port_area_+fraction*(new_port_area_-old_port_area_);
    const double port_length=old_port_length_+fraction*(new_port_length_-old_port_length_);
    for(std::size_t e=1;e<links_.size();++e){
        const Link& link=links_[e];CavityPassage& p=passages_[e];p={};
        const double first=area_[link.from];
        const double second=link.to==cavity_reservoir?first:area_[link.to];
        const double neck=link.port?std::min({first,second,port_area}):std::min(first,second);
        if(neck<=closed_area){continue;}
        p.area=neck;
        double resistance=link.side_entrance?0:length_[link.from]/(2*first*first);
        if(link.to!=cavity_reservoir){resistance+=length_[link.to]/(2*second*second);}
        if(link.port){resistance+=port_length/(port_area*port_area);}
        p.resistance=8*pi*air_.viscosity*resistance;
        p.forward_loss=contraction(first,second,neck);p.reverse_loss=contraction(second,first,neck);
        if(link.port){
            p.forward_loss=contraction(first,port_area,neck)+contraction(port_area,second,neck);
            p.reverse_loss=contraction(second,port_area,neck)+contraction(port_area,first,neck);
        }
        // At an open mouth/nose, static pressure tends to ambient. The
        // acceleration through the lip/nostril contraction was already
        // counted inside the tube; do not add a second outlet-orifice loss.
    }
}
CavityStepReport MeanAirway::advance(double fraction,double glottal_area,double inlet_pressure,
                                    double inlet_resistance,double dt) noexcept {
    if(!ready_ || !std::isfinite(fraction) || fraction<0 || fraction>1){return {};}
    build(fraction);
    inlet_viscosity_=area_[0]>closed_area?4*pi*air_.viscosity*length_[0]/(area_[0]*area_[0]):0;
    passages_[0]={area_[0]>closed_area?glottal_area:0,inlet_resistance+inlet_viscosity_,inlet_pressure};
    return network_.advance(volume_,compliance_,passages_,inflow_,dt,air_);
}
double MeanAirway::oral_flow(std::size_t section) const noexcept {return network_.flows()[oral_ports_[section]];}
}
