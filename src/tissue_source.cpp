#include "ducklips/tissue_source.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ducklips {
namespace {
// Fraction of a rectangular surface patch whose linearly reconstructed gap
// g + ax*x + az*z is positive. This is a geometric area integral, with no
// chosen smoothing length. u and v are the half-ranges across the patch.
double exposed_fraction(double gap,double u,double v) noexcept {
    if(gap>=u+v){return 1;}
    if(gap<=-u-v){return 0;}
    if(std::min(u,v)<1e-15){return std::clamp(.5+gap/(2*std::max(u,v)),0.,1.);}
    const double t=gap+u+v;
    const double a=std::max(0.,t-2*u),b=std::max(0.,t-2*v);
    return std::clamp((t*t-a*a-b*b)/(8*u*v),0.,1.);
}
}
void TissueSource::reset(const TissueMaterial& material){
    material_=material;
    relaxation_.reset(material.relaxation_fraction,material.relaxation_seconds,material.relaxation_log_rate_sd,material.relaxation_seed);
    // The mode basis belongs to one exact configured body. Never silently
    // interpret a different material using these reference-body modes.
    if(material.length!=compiled_tissue::length || material.height!=compiled_tissue::height
        || material.cover_depth!=compiled_tissue::cover_depth || material.body_depth!=compiled_tissue::body_depth
        || material.density!=compiled_tissue::density || material.cover_shear!=compiled_tissue::cover_shear
        || material.body_shear!=compiled_tissue::body_shear || material.viscosity!=compiled_tissue::cover_viscosity
        || material.body_viscosity!=compiled_tissue::body_viscosity
        || material.longitudinal_stress!=compiled_tissue::longitudinal_stress
        || material.body_longitudinal_stress!=compiled_tissue::body_longitudinal_stress){
        throw std::invalid_argument("elastic body changed: recompile data/tissue_body.json");
    }
    if(!std::isfinite(material.rest_gap)||!std::isfinite(material.convergence)
        || !std::isfinite(material.contact_depth)||material.contact_depth<=0){
        throw std::invalid_argument("invalid tissue posture/contact");
    }
    dx_=material.length/(nx-1);dz_=material.height/(nz-1);reference_gap_=material.rest_gap;
    displacement_.fill(0);velocity_.fill(0);equilibrium_.fill(0);force_.fill(0);
    pressure_.fill(0);contact_.fill(0);flow_=0;holding_force_=holding_work_=0;
    for(std::size_t z=0;z<nz;++z){
        for(std::size_t x=0;x<nx;++x){
            const double along=double(x)/(nx-1);
            rest_[z*nx+x]=material.rest_gap*4*along*(1-along)
                +2*material.convergence*(material.height-z*dz_);
        }
    }
    tension_scale_=1;
    for(std::size_t i=0;i<modes;++i){for(std::size_t j=0;j<modes;++j){
        double material=0;
        for(std::size_t k=0;k<modes;++k){material+=compiled_tissue::material_factor[i*modes+k]*compiled_tissue::material_factor[j*modes+k];}
        geometric_stiffness_[i*modes+j]=(i==j?compiled_tissue::stiffness[i]:0)-material;
    }}
    geometry();
}

double TissueSource::posture_stiffness() noexcept{
    double result=0;
    for(std::size_t i=0;i<modes;++i){result+=compiled_tissue::posture[i]*compiled_tissue::stiffness[i]*compiled_tissue::posture[i];}
    return result;
}

void TissueSource::set_adduction(double central_rest_gap){
    if(!std::isfinite(central_rest_gap)){throw std::invalid_argument("nonfinite adduction");}
    const double shift=.5*(central_rest_gap-reference_gap_);
    for(std::size_t mode=0;mode<modes;++mode){equilibrium_[mode]=shift*compiled_tissue::posture[mode];}
    material_.rest_gap=central_rest_gap;
    // The actual displacement and velocity are untouched. Elastic stress
    // performs the posturing work during subsequent physical time steps.
}

void TissueSource::geometry(){
    minimum_area_=1;
    for(std::size_t z=0;z<nz;++z){
        area_[z]=0;
        for(std::size_t x=0;x<nx;++x){
            const std::size_t i=z*nx+x;
            double surface=0,inside=0,velocity=0;
            for(std::size_t mode=0;mode<modes;++mode){
                const double normal=compiled_tissue::surface_basis[i*modes+mode];
                surface+=normal*displacement_[mode];velocity+=normal*velocity_[mode];
                inside+=compiled_tissue::inside_basis[i*modes+mode]*displacement_[mode];
            }
            cover_[i]=surface;body_[i]=inside;cover_v_[i]=velocity;
            gap_[i]=rest_[i]+2*surface;
            if(!std::isfinite(surface)||std::abs(surface)>.01){throw std::runtime_error("tissue displacement outside experiment range");}
            const double weight=(x==0||x+1==nx)?.5:1;
            area_[z]+=weight*dx_*std::max(0.,gap_[i]);
        }
        minimum_area_=std::min(minimum_area_,area_[z]);
    }
}

void TissueSource::aerodynamics(double lung,double outlet){
    // Reverse the axial coordinate when the quasi-steady pressure drop reverses.
    const bool reverse=lung<outlet;
    const double upstream=reverse?outlet:lung,downstream=reverse?lung:outlet;
    const auto at=[reverse](std::size_t z){return reverse?nz-1-z:z;};
    std::array<double,nz> ordered_area{};
    for(std::size_t z=0;z<nz;++z){ordered_area[z]=area_[at(z)];}
    // Laminar slit resistance: channels across x are in parallel, slices
    // along z are in series. Contacted patches carry no throughflow.
    std::size_t throat=0;
    for(std::size_t z=0;z<nz;++z){
        if(ordered_area[z]<ordered_area[throat]){throat=z;}
        double cube_integral=0;
        for(std::size_t x=0;x<nx;++x){
            const double gap=std::max(0.,gap_[at(z)*nx+x]);
            const double weight=(x==0||x+1==nx)?.5:1;
            cube_integral+=weight*dx_*gap*gap*gap;
        }
        resistance_[z]=12*material_.air_viscosity/std::max(cube_integral,1e-24);
    }
    // Separation is a geometric closure approximation, not fitted pulse
    // timing: first divergent section at 1.2 times the throat area.
    std::size_t separation=nz-1;
    double fraction=1;
    double jet_area=ordered_area[nz-1];
    for(std::size_t z=throat+1;z<nz;++z){
        if(ordered_area[z]>=1.2*ordered_area[throat]){
            separation=z;
            jet_area=1.2*ordered_area[throat];
            fraction=(jet_area-ordered_area[z-1])/std::max(ordered_area[z]-ordered_area[z-1],1e-24);
            break;
        }
    }
    double total_resistance=0;
    for(std::size_t z=1;z<=separation;++z){
        const double part=z==separation?fraction:1;
        total_resistance+=dz_*part*(resistance_[z-1]+.5*part*(resistance_[z]-resistance_[z-1]));
    }
    const double delta=std::abs(lung-outlet);
    jet_area=std::max(jet_area,1e-12);
    const double beta=.5*material_.air_density/(jet_area*jet_area);
    flow_=minimum_area_<=1e-12?0:2*delta/(total_resistance+std::sqrt(total_resistance*total_resistance+4*beta*delta));
    const double magnitude=flow_;
    if(reverse){flow_=-flow_;}
    double upstream_resistance=0;
    bool upstream_open=true;
    for(std::size_t z=0;z<nz;++z){
        if(z){upstream_resistance+=.5*dz_*(resistance_[z-1]+resistance_[z]);}
        if(ordered_area[z]<=1e-12){upstream_open=false;}
        double pressure=downstream;
        if((z<separation || (z==separation && fraction==1)) && upstream_open){
            const double velocity=magnitude/std::max(ordered_area[z],1e-12);
            pressure=upstream-upstream_resistance*magnitude-.5*material_.air_density*velocity*velocity;
        }
        for(std::size_t x=0;x<nx;++x){
            const std::size_t physical_z=at(z);
            const std::size_t i=physical_z*nx+x;
            const std::size_t left=x==0?x:x-1,right=x+1==nx?x:x+1;
            const std::size_t below=physical_z==0?physical_z:physical_z-1,above=physical_z+1==nz?physical_z:physical_z+1;
            const double u=.25*std::abs(gap_[physical_z*nx+right]-gap_[physical_z*nx+left]);
            const double v=.25*std::abs(gap_[above*nx+x]-gap_[below*nx+x]);
            pressure_[i]=exposed_fraction(gap_[i],u,v)*pressure;
        }
    }
}

double TissueSource::advance(double lung,double outlet,double dt){
    aerodynamics(lung,outlet);
    force_.fill(0);
    for(std::size_t i=0;i<count;++i){
        const double penetration=std::max(0.,-.5*gap_[i]);
        const double strain=penetration/material_.contact_depth;
        // Compliant contact is separate from bulk elasticity. The normal
        // modulus 3G is the incompressible uniaxial small-strain limit.
        contact_[i]=3*material_.cover_shear*strain*(1+strain*strain)
            +3*material_.viscosity/material_.contact_depth*strain*std::max(0.,-cover_v_[i]);
        const double force=(pressure_[i]+contact_[i])*compiled_tissue::surface_area[i];
        for(std::size_t mode=0;mode<modes;++mode){force_[mode]+=compiled_tissue::surface_basis[i*modes+mode]*force;}
    }
    std::array<double,modes> acceleration{};
    const auto relieved=relaxation_.relieved_force();
    const double holding_shift=-holding_force_/posture_stiffness();
    for(std::size_t mode=0;mode<modes;++mode){
        double damping=0;
        for(std::size_t other=0;other<modes;++other){damping+=compiled_tissue::damping[mode*modes+other]*velocity_[other];}
        force_[mode]+=holding_shift*compiled_tissue::stiffness[mode]*compiled_tissue::posture[mode];
        double tension_force=0;
        if(tension_scale_!=1){for(std::size_t other=0;other<modes;++other){tension_force+=(tension_scale_-1)*geometric_stiffness_[mode*modes+other]*(displacement_[other]-equilibrium_[other]);}}
        acceleration[mode]=force_[mode]-tension_force-compiled_tissue::stiffness[mode]*(displacement_[mode]-equilibrium_[mode])-damping+relieved[mode];
    }
    for(std::size_t mode=0;mode<modes;++mode){
        velocity_[mode]+=dt*acceleration[mode];displacement_[mode]+=dt*velocity_[mode];
        holding_work_+=2*holding_shift*compiled_tissue::stiffness[mode]*compiled_tissue::posture[mode]*velocity_[mode]*dt;
    }
    TissueRelaxation::Coordinates offset{};
    for(std::size_t mode=0;mode<modes;++mode){offset[mode]=displacement_[mode]-equilibrium_[mode];}
    relaxation_.advance(offset,dt);
    geometry();
    return flow_;
}

double TissueSource::kinetic_energy() const noexcept{
    double energy=0;
    for(double velocity:velocity_){energy+=velocity*velocity;}
    // Unit modal masses, .5 per fold, two symmetric folds.
    return energy;
}

double TissueSource::elastic_energy() const noexcept{
    TissueRelaxation::Coordinates deformation{};
    for(std::size_t mode=0;mode<modes;++mode){deformation[mode]=displacement_[mode]-equilibrium_[mode];}
    double energy=relaxation_.energy_correction(deformation);
    for(std::size_t mode=0;mode<modes;++mode){
        const double offset=displacement_[mode]-equilibrium_[mode];
        energy+=compiled_tissue::stiffness[mode]*offset*offset;
    }
    for(std::size_t i=0;i<modes;++i){for(std::size_t j=0;j<modes;++j){energy+=(tension_scale_-1)*deformation[i]*geometric_stiffness_[i*modes+j]*deformation[j];}}
    for(std::size_t i=0;i<count;++i){
        const double strain=std::max(0.,-.5*gap_[i])/material_.contact_depth;
        energy+=6*compiled_tissue::surface_area[i]*material_.cover_shear*material_.contact_depth
            *(.5*strain*strain+.25*strain*strain*strain*strain);
    }
    return energy;
}
}
