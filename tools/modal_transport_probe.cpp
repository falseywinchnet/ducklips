#include "ducklips/modal_transport.hpp"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
constexpr double pi=3.14159265358979323846;
void require(bool okay,const char* message) { if (!okay) { throw std::runtime_error(message); } }
void read(std::ifstream& in,void* target,std::size_t bytes) {
    in.read(static_cast<char*>(target),static_cast<std::streamsize>(bytes));require(bool(in),"truncated transport fixture");
}
void fixture(const char* input,const char* output) {
    std::ifstream in(input,std::ios::binary);std::uint64_t header[4];read(in,header,sizeof(header));
    const std::size_t m=header[0],n=header[1],b=header[2],frames=header[3];
    require(m>0 && m<=n && n<=4096 && b<=16384 && frames<=100000,"invalid fixture size");
    ducklips::AcousticNetwork network;network.nodes.resize(n);network.branches.resize(b);
    for (ducklips::AcousticBranch& branch:network.branches) {
        std::uint64_t edge[2];read(in,edge,sizeof(edge));branch.from=edge[0];branch.to=edge[1];
    }
    ducklips::ModalBasis basis;basis.prepare(n,m);basis.mode_count=m;
    ducklips::ModalTransportReference transport;transport.prepare(m,n,b);
    std::vector<double> dc(n),dg(b);
    std::ofstream out(output,std::ios::binary);require(bool(out),"cannot open transport output");
    for (std::size_t frame=0;frame<frames;++frame) {
        read(in,basis.angular_frequency.data(),m*sizeof(double));read(in,basis.pressure_shape.data(),m*n*sizeof(double));
        read(in,dc.data(),n*sizeof(double));read(in,dg.data(),b*sizeof(double));
        require(transport.derive(network,basis,dc,dg),"rejected transport fixture");
        const ducklips::ModalMotionSample sample=transport.view();
        for (std::span<const double> values:{transport.frequency_rate(),sample.pressure_transport,sample.flow_transport}) {
            out.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(values.size()*sizeof(double)));
        }
    }
    require(bool(out),"failed transport output");
}
void tests() {
    ducklips::AcousticNetwork network;network.nodes.resize(1);network.branches.resize(1);
    network.branches[0].from=0;network.branches[0].to=ducklips::ground_node;
    ducklips::ModalBasis basis;basis.prepare(1,1);basis.mode_count=1;basis.pressure_shape[0]=1;
    basis.angular_frequency[0]=100;
    ducklips::ModalTransportReference transport;transport.prepare(1,1,1);
    const double dc[]={14},dg[]={140000};
    require(transport.derive(network,basis,dc,dg),"uniform transport rejected");
    require(std::abs(transport.view().pressure_transport[0]+7)<1e-12 &&
            std::abs(transport.view().flow_transport[0]-7)<1e-12 && transport.frequency_rate()[0]==0,"uniform transport mismatch");
    const double dc2[]={-2},dg2[]={30000};
    require(transport.derive(network,basis,dc2,dg2),"nonuniform metric rejected");
    require(transport.view().pressure_transport[0]==1 && transport.view().flow_transport[0]==1.5 &&
            transport.frequency_rate()[0]==250,"independent compliance/inertance derivatives wrong");
    require(transport.canonical_view().flow_transport.empty() && transport.canonical_view().flow_diagonal[0]==2.5,
            "canonical coefficient view is inconsistent");
    basis.angular_frequency[0]=0;
    require(!transport.derive(network,basis,dc,dg) && transport.view().flow_transport[0]==1.5 &&
            transport.view().angular_frequency[0]==100,"failed transport replaced last valid coefficients");
    network.nodes.resize(2);basis.prepare(2,2);basis.mode_count=2;basis.angular_frequency={100,100};
    transport.prepare(2,2,1);const double dc3[]={0,0};
    require(!transport.derive(network,basis,dc3,dg),"degenerate modes silently accepted");
    std::cout<<"Physical modal transport uniform limit, independent metrics, transaction and degeneracy checks passed\n";
}
void benchmark(std::size_t modes) {
    constexpr std::size_t nodes=96,frames=5000;
    require(modes>0 && modes<=nodes,"invalid mode count");
    ducklips::AcousticNetwork network;network.nodes.resize(nodes);network.branches.resize(2*nodes-1);
    for (std::size_t i=1;i<nodes;++i) { network.branches[i-1].from=i-1;network.branches[i-1].to=i; }
    for (std::size_t i=0;i<nodes;++i) { network.branches[nodes-1+i].from=i;network.branches[nodes-1+i].to=ducklips::ground_node; }
    ducklips::ModalBasis basis;basis.prepare(nodes,modes);basis.mode_count=modes;
    constexpr double compliance=5e-11;
    for (std::size_t i=0;i<modes;++i) {
        basis.angular_frequency[i]=2*pi*(200+double(i)*10000/double(modes));
        for (std::size_t node=0;node<nodes;++node) {
            basis.pressure_shape[i*nodes+node]=std::sqrt(2.0/nodes/compliance)*std::cos(pi*(node+.5)*(i+1)/nodes);
        }
    }
    std::vector<double> dc(nodes),dg(network.branches.size());
    ducklips::ModalTransportReference transport;transport.prepare(modes,nodes,dg.size());
    const std::chrono::steady_clock::time_point begin=std::chrono::steady_clock::now();
    for (std::size_t frame=0;frame<frames;++frame) {
        const double rate=std::sin(2*pi*double(frame)/1000);
        for (std::size_t i=0;i<nodes;++i) { dc[i]=compliance*(double(i%7)-3)*rate; }
        for (std::size_t i=0;i<dg.size();++i) { dg[i]=.004*(double(i%5)-2)*rate; }
        require(transport.derive(network,basis,dc,dg),"transport benchmark rejected");
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    std::cout<<std::setprecision(17)<<"{\"modes\":"<<modes<<",\"nodes\":"<<nodes<<",\"updates\":"<<frames
             <<",\"updates_per_simulated_second\":1000,\"wall_seconds\":"<<seconds<<",\"equivalent_realtime_multiple\":"<<5/seconds<<"}\n";
}
}
int main(int argc,char** argv) {
    try {
        if (argc==1) { tests(); }
        else if (argc==4 && std::string(argv[1])=="--fixture") { fixture(argv[2],argv[3]); }
        else if (argc==3 && std::string(argv[1])=="--benchmark") { benchmark(std::stoul(argv[2])); }
        else { throw std::invalid_argument("modal_transport_probe [--fixture input output | --benchmark modes]"); }
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
