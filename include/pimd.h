/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

// standard C++ headers
#include <fstream>

// uovie headers
#include "simu_para.h"
#include "mol_geom.h"

#ifndef PHY_CONST_SHORTHAND
#define PHY_CONST_SHORTHAND
constexpr double h_bar = uovie::phy_const::red_Planck_const;
constexpr double k = uovie::phy_const::Boltzmann_const;
#endif // !PHY_CONST_SHORTHAND

namespace uovie {
namespace pimd {
    
    class pimd_via_nhc {
    public:
        pimd_via_nhc() = default;
        pimd_via_nhc(std::ofstream& _out, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int _nbead, const int _nchain) :
            out(_out), bsp(_bsp), sys(_sys), nbead(_nbead), nchain(_nchain) { }

        void pimd_implement();

    private:
        std::ofstream& out;
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const int nbead;
        const int nchain;

        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const double& T = sys.temperature;
        const int& M = nchain;
    };

} // !pimd
} // !uovie
#endif // !PIMD_H_