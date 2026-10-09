#include "frozen_modes.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition,const char* message) { if (!condition) { throw std::runtime_error(message); } }
void write(const std::filesystem::path& path,const std::string& text) {
    std::ofstream out(path);out<<text;require(static_cast<bool>(out),"cannot write test model");
}
const std::string valid="DUCKLIPS_FROZEN_MODES_V1 192000 2 1 2\n2 0.3\n0.8 0.1 1 0.2 1.4 -0.8 0.2 0.1\n0.5 0 0.25 0 0.7 0 0.9 0\n";
}
int main() {
    const std::filesystem::path path=std::filesystem::temp_directory_path()/
        ("ducklips-modes-test-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        write(path,valid);
        ducklips::FrozenModes general,voiced;
        general.load(path.string());voiced.load(path.string());
        std::complex<double> pair=0;double real=0;
        for (int iteration=0;iteration<4096;++iteration) {
            const double input=std::sin(iteration*0.37)+0.3*std::cos(iteration*0.111);
            const double pressure=(std::complex<double>(1.4,-0.8)*pair).real()+0.7*real+2*input;
            const double radiation=(std::complex<double>(0.2,0.1)*pair).real()+0.9*real+0.3*input;
            double output[2]{};
            general.tick(std::span<const double>(&input,1),output);
            require(std::abs(output[0]-pressure)<1e-12 && std::abs(output[1]-radiation)<1e-12,"generic sample timing or complex pair incorrect");
            require(std::abs(voiced.free_pressure(0)+voiced.direct(0,0)*input-pressure)<1e-12,"coupling pressure incorrect");
            require(std::abs(voiced.tick_glottis(input)-radiation)<1e-12,"specialized tick disagrees");
            pair=std::complex<double>(0.8,0.1)*pair+std::complex<double>(1,0.2)*input;
            real=0.5*real+0.25*input;
        }
        general.reset();require(general.free_pressure(0)==0 && general.free_pressure(1)==0,"reset retained energy");
        bool rejected=false;
        try { double output[1];const double input=0;general.tick(std::span<const double>(&input,1),output); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"wrong port count accepted");
        for (const std::string& bad : {std::string(""),valid+"garbage",std::string("DUCKLIPS_FROZEN_MODES_V1 192000 1 1 2\n0 0\n1.1 0 1 0 1 0 1 0\n"),std::string("DUCKLIPS_FROZEN_MODES_V1 192000 1 1 2\n0 0\n0.8 -0.1 1 0 1 0 1 0\n"),std::string("DUCKLIPS_FROZEN_MODES_V1 nan 1 1 2\n")}) {
            write(path,bad);rejected=false;
            try { general.load(path.string()); } catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"invalid modal data accepted");
            require(general.modes()==2 && general.direct(0,0)==2,"failed load changed prepared model");
        }
        std::filesystem::remove(path);
        std::cout<<"Frozen real/complex recurrence, coupling, reset and transactional parser passed\n";
    } catch (const std::exception& error) {
        std::filesystem::remove(path);std::cerr<<error.what()<<'\n';return 1;
    }
}
