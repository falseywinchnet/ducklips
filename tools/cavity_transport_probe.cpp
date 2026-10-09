// Persistent Eulerian cells with disappearing/reappearing air volume.
// This is a mean-pressure transport reference, not acoustic speech.
// stdin: nodes edges states; edges: from to (-1 reservoir);
// initial nodes: V C p; then each step: dt, nodes (V C inflow),
// edges (A R reservoir_pressure forward_loss reverse_loss).
#include "ducklips/cavity_network.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>
int main(){
    std::size_t n=0,m=0,states=0;
    if(std::scanf("%zu %zu %zu",&n,&m,&states)!=3 || n==0 || n>65536 || m>262144 || states<2 || states>100000){return 2;}
    std::vector<ducklips::CavityConnection> links(m);
    for(ducklips::CavityConnection& link:links){
        long long to=0;if(std::scanf("%zu %lld",&link.from,&to)!=2 || to < -1){return 2;}
        link.to=to==-1?ducklips::cavity_reservoir:static_cast<std::size_t>(to);
    }
    std::vector<double> v(n),c(n),p(n),in(n);std::vector<ducklips::CavityPassage> ports(m);
    for(std::size_t i=0;i<n;++i){if(std::scanf("%lf %lf %lf",&v[i],&c[i],&p[i])!=3){return 2;}}
    ducklips::CavityNetwork network{};try{network.prepare(n,links);}catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 2;}
    if(!network.reset(v,c,p)){return 2;}
    std::puts("step,node,pressure");
    for(std::size_t i=0;i<n;++i){std::printf("0,%zu,%.17g\n",i,network.pressures()[i]);}
    const ducklips::AirProperties air{};double largest_residual=0;std::size_t iterations=0;
    for(std::size_t step=1;step<states;++step){
        double dt=0;if(std::scanf("%lf",&dt)!=1){return 2;}
        for(std::size_t i=0;i<n;++i){if(std::scanf("%lf %lf %lf",&v[i],&c[i],&in[i])!=3){return 2;}}
        for(ducklips::CavityPassage& port:ports){
            if(std::scanf("%lf %lf %lf %lf %lf",&port.area,&port.resistance,&port.reservoir_pressure,&port.forward_loss,&port.reverse_loss)!=5){return 2;}
        }
        const ducklips::CavityStepReport r=network.advance(v,c,ports,in,dt,air);
        if(!r.converged){std::fprintf(stderr,"transport failed step %zu at Newton %zu residual %.17g\n",step,r.iterations,r.maximum_mass_residual);return 1;}
        largest_residual=std::max(largest_residual,r.maximum_mass_residual);iterations=std::max(iterations,r.iterations);
        for(std::size_t i=0;i<n;++i){std::printf("%zu,%zu,%.17g\n",step,i,network.pressures()[i]);}
    }
    std::fprintf(stderr,"nodes=%zu edges=%zu states=%zu max_residual=%.17g max_Newton=%zu\n",n,m,states,largest_residual,iterations);
}
