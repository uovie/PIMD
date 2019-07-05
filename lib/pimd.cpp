// standard C++ headers
#include <vector>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "phy_const.h"
#include "pimd.h"
#include "thermostat/nhc.h"

namespace uovie {
namespace pimd {

    void pimd_via_nhc::pimd_implement()
    {
        Eigen::ArrayXXd ZeroXXd(d * N, nbead); // ZeroXXd.setZero();
        Eigen::ArrayXXd ZeroXXXd(d * N * M, nbead); // ZeroXXXd.setZero();
        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        thermostat::nhc::nhc_procedure_for_pimd nhc_proce(bsp, sys, tfs, ZeroXXd, ZeroXXXd);
        nhc_proce.implement(out);
    }

} // !pimd
} // !uovie