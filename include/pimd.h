/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

namespace uovie {
namespace pimd {
    
    void implement(std::ofstream& out, const Global::basic_simu_para& bsp,
        const Global::system& sys, const int nbead, const int nchain);

} // !pimd
} // !uovie
#endif // !PIMD_H_