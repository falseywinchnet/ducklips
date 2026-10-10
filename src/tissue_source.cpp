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
    for(double value:{material.length,material.height,material.cover_depth,material.body_depth,
            material.density,material.cover_shear,material.body_shear,material.contact_depth}){
        if(!std::isfinite(value)||value<=0){throw std::invalid_argument("invalid tissue material");}
    }
    dx_=material.length/(nx-1);dz_=material.height/(nz-1);
    cover_.fill(0);body_.fill(0);cover_v_.fill(0);body_v_.fill(0);
    rotation_.fill(0);rotation_v_.fill(0);rotation_a_.fill(0);shear_.fill(0);
    pressure_.fill(0);contact_.fill(0);flow_=0;
    for(std::size_t z=0;z<nz;++z){
        for(std::size_t x=0;x<nx;++x){
            const double along=double(x)/(nx-1);
            // Adducted attachment points and a slightly convergent medial
            // surface. This is a rest geometry, not a vibration mode shape.
            rest_[z*nx+x]=material.rest_gap*4*along*(1-along)
                +2*material.convergence*(material.height-z*dz_);
        }
    }
    geometry();
}

void TissueSource::set_adduction(double central_rest_gap){
    if(!std::isfinite(central_rest_gap)){throw std::invalid_argument("nonfinite adduction");}
    const double difference=material_.rest_gap-central_rest_gap;
    for(std::size_t z=0;z<nz;++z){
        for(std::size_t x=0;x<nx;++x){
            const std::size_t i=z*nx+x;
            const double along=double(x)/(nx-1);
            const double shift=.5*difference*4*along*(1-along);
            rest_[i]-=2*shift;cover_[i]+=shift;body_[i]+=shift;
        }
    }
    material_.rest_gap=central_rest_gap;
    geometry();
}

void TissueSource::geometry(){
    minimum_area_=1;
    for(std::size_t z=0;z<nz;++z){
        area_[z]=0;
        for(std::size_t x=0;x<nx;++x){
            const std::size_t i=z*nx+x;
            gap_[i]=rest_[i]+2*cover_[i];
            const double weight=(x==0||x+1==nx)?.5:1;
            area_[z]+=weight*dx_*std::max(0.,gap_[i]);
        }
        minimum_area_=std::min(minimum_area_,area_[z]);
    }
}

void TissueSource::aerodynamics(double lung,double outlet){
    // Laminar slit resistance: channels across x are in parallel, slices
    // along z are in series. Contacted patches carry no throughflow.
    std::size_t throat=0;
    for(std::size_t z=0;z<nz;++z){
        if(area_[z]<area_[throat]){throat=z;}
        double cube_integral=0;
        for(std::size_t x=0;x<nx;++x){
            const double gap=std::max(0.,gap_[z*nx+x]);
            const double weight=(x==0||x+1==nx)?.5:1;
            cube_integral+=weight*dx_*gap*gap*gap;
        }
        resistance_[z]=12*material_.air_viscosity/std::max(cube_integral,1e-24);
    }
    // Separation is a geometric closure approximation, not fitted pulse
    // timing: first divergent section at 1.2 times the throat area.
    std::size_t separation=nz-1;
    double fraction=1;
    double jet_area=area_[nz-1];
    for(std::size_t z=throat+1;z<nz;++z){
        if(area_[z]>=1.2*area_[throat]){
            separation=z;
            jet_area=1.2*area_[throat];
            fraction=(jet_area-area_[z-1])/std::max(area_[z]-area_[z-1],1e-24);
            break;
        }
    }
    double total_resistance=0;
    for(std::size_t z=1;z<=separation;++z){
        const double part=z==separation?fraction:1;
        total_resistance+=dz_*part*(resistance_[z-1]+.5*part*(resistance_[z]-resistance_[z-1]));
    }
    const double delta=lung-outlet;
    if(delta<0){throw std::runtime_error("reverse flow not implemented in source experiment");}
    jet_area=std::max(jet_area,1e-12);
    const double beta=.5*material_.air_density/(jet_area*jet_area);
    flow_=minimum_area_<=1e-12?0:2*delta/(total_resistance+std::sqrt(total_resistance*total_resistance+4*beta*delta));
    double upstream_resistance=0;
    bool upstream_open=true;
    for(std::size_t z=0;z<nz;++z){
        if(z){upstream_resistance+=.5*dz_*(resistance_[z-1]+resistance_[z]);}
        if(area_[z]<=1e-12){upstream_open=false;}
        double pressure=outlet;
        if((z<separation || (z==separation && fraction==1)) && upstream_open){
            const double velocity=flow_/std::max(area_[z],1e-12);
            pressure=lung-upstream_resistance*flow_-.5*material_.air_density*velocity*velocity;
        }
        for(std::size_t x=0;x<nx;++x){
            const std::size_t i=z*nx+x;
            const std::size_t left=x==0?x:x-1,right=x+1==nx?x:x+1;
            const std::size_t below=z==0?z:z-1,above=z+1==nz?z:z+1;
            const double u=.25*std::abs(gap_[z*nx+right]-gap_[z*nx+left]);
            const double v=.25*std::abs(gap_[above*nx+x]-gap_[below*nx+x]);
            pressure_[i]=exposed_fraction(gap_[i],u,v)*pressure;
        }
    }
}

double TissueSource::advance(double lung,double outlet,double dt){
    aerodynamics(lung,outlet);
    const double cover_mass=material_.density*material_.cover_depth;
    // Linear displacement across backing depth gives integral rho*(r/d)^2
    // dr = rho*d/3. The attached lateral boundary has zero displacement.
    const double body_mass=material_.density*material_.body_depth/3;
    const double coupling=material_.cover_shear/material_.cover_depth;
    const double support=material_.body_shear/material_.body_depth;
    const double relative_drag=material_.viscosity/material_.cover_depth;
    const double support_drag=material_.body_viscosity/material_.body_depth;
    const double rotation_inertia=material_.density*material_.cover_depth*material_.cover_depth*material_.cover_depth/12;
    for(std::size_t z=0;z+1<nz;++z){
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            const double strain=(cover_[i+nx]-cover_[i])/dz_-rotation_[i];
            const double rate=(cover_v_[i+nx]-cover_v_[i])/dz_-rotation_v_[i];
            shear_[i]=material_.cover_depth*(material_.cover_shear*strain+material_.viscosity*rate);
            rotation_a_[i]=shear_[i]/rotation_inertia;
        }
    }
    for(std::size_t z=0;z<nz;++z){
        const std::size_t below=z==0?1:z-1,above=z+1==nz?nz-2:z+1;
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            const double cxx=(cover_[i-1]-2*cover_[i]+cover_[i+1])/(dx_*dx_);
            const double bxx=(body_[i-1]-2*body_[i]+body_[i+1])/(dx_*dx_);
            const double bzz=(body_[below*nx+x]-2*body_[i]+body_[above*nx+x])/(dz_*dz_);
            const double vxx=(cover_v_[i-1]-2*cover_v_[i]+cover_v_[i+1])/(dx_*dx_);
            const double bvxx=(body_v_[i-1]-2*body_v_[i]+body_v_[i+1])/(dx_*dx_);
            const double bvzz=(body_v_[below*nx+x]-2*body_v_[i]+body_v_[above*nx+x])/(dz_*dz_);
            const double relative=cover_[i]-body_[i];
            const double strain=relative/material_.cover_depth;
            // Positive stored energy .5*k*r² + .25*k*r⁴/d²: strain stiffening
            // of the compliant cover, not a waveform-shaping nonlinearity.
            const double coupling_pressure=coupling*relative*(1+strain*strain)
                +relative_drag*(cover_v_[i]-body_v_[i]);
            const double penetration=std::max(0.,-.5*gap_[i]);
            const double contact_strain=penetration/material_.contact_depth;
            // Distributed compliant contact: pressure and damping approach
            // zero with penetration; no binary damping jump at first touch.
            contact_[i]=material_.cover_shear*contact_strain*(1+contact_strain*contact_strain)
                +relative_drag*contact_strain*std::max(0.,-cover_v_[i]);
            const double lower_shear=z==0?0:shear_[(z-1)*nx+x];
            const double upper_shear=z+1==nz?0:shear_[i];
            const double vertical_force=(upper_shear-lower_shear)/(dz_*((z==0||z+1==nz)?.5:1));
            cover_a_[i]=(material_.cover_depth*material_.longitudinal_stress*cxx+vertical_force
                +material_.viscosity*material_.cover_depth*vxx
                -coupling_pressure+pressure_[i]+contact_[i])/cover_mass;
            body_a_[i]=(material_.body_depth/3*(material_.body_longitudinal_stress*bxx+material_.body_shear*bzz)
                +material_.body_viscosity*material_.body_depth/3*(bvxx+bvzz)
                +coupling_pressure-support*body_[i]-support_drag*body_v_[i])/body_mass;
        }
    }
    // Local explicit mechanical update. Work scales with surface samples;
    // there is no global matrix factorization or nonlinear root iteration.
    for(std::size_t z=0;z<nz;++z){
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            cover_v_[i]+=dt*cover_a_[i];body_v_[i]+=dt*body_a_[i];
            cover_[i]+=dt*cover_v_[i];body_[i]+=dt*body_v_[i];
            if(z+1<nz){rotation_v_[i]+=dt*rotation_a_[i];rotation_[i]+=dt*rotation_v_[i];}
            if(!std::isfinite(cover_[i])||std::abs(cover_[i])>.01){throw std::runtime_error("tissue experiment diverged");}
        }
    }
    // Keep observable geometry at the new tissue state. It also supplies
    // the next step's aerodynamic geometry without a redundant evaluation.
    geometry();
    return flow_;
}

double TissueSource::kinetic_energy() const noexcept{
    double energy=0;
    for(std::size_t z=0;z<nz;++z){
        const double weight=(z==0||z+1==nz)?.5:1;
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            energy+=weight*dx_*dz_*material_.density*(material_.cover_depth*cover_v_[i]*cover_v_[i]
                +material_.body_depth/3*body_v_[i]*body_v_[i]);
            if(z+1<nz){energy+=dx_*dz_*material_.density*material_.cover_depth*material_.cover_depth*material_.cover_depth/12*rotation_v_[i]*rotation_v_[i];}
        }
    }
    // .5 per fold, multiplied by the two symmetric folds.
    return energy;
}

double TissueSource::elastic_energy() const noexcept{
    double energy=0;
    const double coupling=material_.cover_shear/material_.cover_depth;
    const double support=material_.body_shear/material_.body_depth;
    for(std::size_t z=0;z<nz;++z){
        const double wz=(z==0||z+1==nz)?.5:1;
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            const double r=cover_[i]-body_[i],strain=r/material_.cover_depth;
            const double penetration=std::max(0.,-.5*(rest_[i]+2*cover_[i]));
            const double contact_strain=penetration/material_.contact_depth;
            energy+=dx_*dz_*wz*(.5*coupling*r*r+.25*coupling*r*r*strain*strain
                +.5*support*body_[i]*body_[i]
                +material_.cover_shear*material_.contact_depth*(.5*contact_strain*contact_strain
                    +.25*contact_strain*contact_strain*contact_strain*contact_strain));
        }
        for(std::size_t x=0;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            const double c=(cover_[i+1]-cover_[i])/dx_,b=(body_[i+1]-body_[i])/dx_;
            energy+=.5*dx_*dz_*wz*(material_.cover_depth*material_.longitudinal_stress*c*c
                +material_.body_depth/3*material_.body_longitudinal_stress*b*b);
        }
    }
    for(std::size_t z=0;z+1<nz;++z){
        for(std::size_t x=1;x+1<nx;++x){
            const std::size_t i=z*nx+x;
            const double c=(cover_[i+nx]-cover_[i])/dz_-rotation_[i],b=(body_[i+nx]-body_[i])/dz_;
            energy+=.5*dx_*dz_*(material_.cover_depth*material_.cover_shear*c*c
                +material_.body_depth/3*material_.body_shear*b*b);
        }
    }
    return 2*energy;
}
}
