// standard C++ headers
#include <vector>
#include <fstream>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "phy_const.h"
#include "pimd.h"
#include "thermostat/nhc.h"

#ifndef PHY_CONST_SHORTHAND
#define PHY_CONST_SHORTHAND
constexpr double h_bar = uovie::phy_const::red_Planck_const;
constexpr double k = uovie::phy_const::Boltzmann_const;
#endif // !PHY_CONST_SHORTHAND

namespace uovie {
namespace pimd {

    void implement(std::ofstream& out, const Global::basic_simu_para& bsp,
        const Global::system& sys, const int nbead, const int nchain)
    {
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const double& T = sys.temperature;
        const int& M = nchain;

        Eigen::ArrayXXd zero_array = Eigen::ArrayXXf::Zero(d * N, nbead);
        
        /* NHC Preparation */
        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        thermostat::nhc::nhc_procedure nhc_proce(bsp, sys, zero_array, tfs, 4);
        //std::cout << "NHC is ready, and PIMD starts." << std::endl;
        /* NHC Preparation */

        nhc_proce.implement(out);
    }

} // !pimd
} // !uovie