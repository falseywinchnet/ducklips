// Offline pressure integration for an arbitrary, fixed air-space graph.
// stdin: node_count edge_count step_count dt
//        V C initial_pressure inflow (one row per node)
//        from to area resistance reservoir_pressure forward_loss reverse_loss
//        (one row per edge; to=-1 denotes reservoir).
// stdout: kind,index,value (pressure per node, flow per edge).
#include "ducklips/cavity_network.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
int main() {
    std::size_t n=0,m=0,steps=0;double dt=0;
    if(std::scanf("%zu %zu %zu %lf",&n,&m,&steps,&dt)!=4 || n==0 || n>65536 || m>262144 || steps==0 || steps>100000
        || !std::isfinite(dt) || dt<=0){return 2;}
    std::vector<double> v(n),c(n),p(n),in(n);
    for(std::size_t i=0;i<n;++i){if(std::scanf("%lf %lf %lf %lf",&v[i],&c[i],&p[i],&in[i])!=4){return 2;}}
    std::vector<ducklips::CavityConnection> links(m);std::vector<ducklips::CavityPassage> passages(m);
    for(std::size_t e=0;e<m;++e){
        long long to=0;ducklips::CavityPassage& a=passages[e];
        if(std::scanf("%zu %lld %lf %lf %lf %lf %lf",&links[e].from,&to,&a.area,&a.resistance,&a.reservoir_pressure,&a.forward_loss,&a.reverse_loss)!=7 || to < -1){return 2;}
        links[e].to=to==-1?ducklips::cavity_reservoir:static_cast<std::size_t>(to);
    }
    ducklips::CavityNetwork network{};
    try{network.prepare(n,links);}catch(const std::exception& error){std::fprintf(stderr,"%s\n",error.what());return 2;}
    if(!network.reset(v,c,p)){return 2;}const ducklips::AirProperties air{};
    double largest_residual=0;std::size_t most_iterations=0;
    for(std::size_t step=0;step<steps;++step){
        const ducklips::CavityStepReport r=network.advance(v,c,passages,in,dt,air);
        if(!r.converged){std::fprintf(stderr,"failed step %zu, Newton iteration %zu, residual %.17g\n",step,r.iterations,r.maximum_mass_residual);return 1;}
        largest_residual=std::max(largest_residual,r.maximum_mass_residual);most_iterations=std::max(most_iterations,r.iterations);
    }
    std::puts("kind,index,value");
    for(std::size_t i=0;i<n;++i){std::printf("pressure,%zu,%.17g\n",i,network.pressures()[i]);}
    for(std::size_t e=0;e<m;++e){std::printf("flow,%zu,%.17g\n",e,network.flows()[e]);}
    std::fprintf(stderr,"nodes=%zu edges=%zu steps=%zu max_residual=%.17g max_Newton=%zu\n",n,m,steps,largest_residual,most_iterations);
}
