/* Nose-Hoover Chain */
#ifndef NOSE_HOOVER_CHAIN_H_
#define NOSE_HOOVER_CHAIN_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <cmath>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "simu_para.h"
#include "phy_const.h"
#include "mol_geom.h"

#ifndef PHY_CONST_SHORTHAND
#define PHY_CONST_SHORTHAND
constexpr double h_bar = uovie::phy_const::red_Planck_const;
constexpr double k = uovie::phy_const::Boltzmann_const;
#endif // !PHY_CONST_SHORTHAND

namespace uovie {
namespace thermostat {
namespace nhc {

    /*** ================================================== ***/
    /*** Thermostat Factorization Scheme                    ***/
    /*** ================================================== ***/

    // Suzuki¨CYoshida scheme
    class thermo_factor_scheme {
    public:
        thermo_factor_scheme() = default;
        thermo_factor_scheme(const int nsy, const int nff): n_sy(nsy), n_ff(nff) {
            if (nsy == 3)
                weight = { 1.351207191959658, -1.702414383919316, 1.351207191959658 };
            else if (nsy == 7)
                weight = { 0.784513610477560, 0.235573213359357, -1.17767998417887,
                    1.315186320683906, -1.17767998417887, 0.235573213359357, 0.784513610477560 };
            else
                throw "unsupported Suzuki-Yoshida scheme";
        }

        int nsy() const { return n_sy; }
        int nff() const { return n_ff; }
        double w(const int& i) const { return weight[i]; }

    private:
        int n_sy;                       // the number of Suzuki-Yoshida weights
        int n_ff;                       // the number of terms in the further factorization
        std::vector<double> weight;     // Suzuki-Yoshida weights
    };

    /*** ================================================== ***/
    /*** NHC Procedure for pimd                             ***/
    /*** ================================================== ***/

    class nhc_procedure_for_pimd {
    public:
        nhc_procedure_for_pimd() = default;
        nhc_procedure_for_pimd(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const thermo_factor_scheme& _tfs, const Eigen::ArrayXXd& _ZeroXXd, const Eigen::ArrayXXd& _ZeroXXXd) :
            bsp(_bsp), sys(_sys), tfs(_tfs), ZeroXXd(_ZeroXXd), ZeroXXXd(_ZeroXXXd) { }
        
        void initialize();
        void implement();
        void implement(std::ofstream& out);
        //double sys_ene() const { return syst_energy; }

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const thermo_factor_scheme& tfs;
        const Eigen::ArrayXXd& ZeroXXd;
        const Eigen::ArrayXXd& ZeroXXXd;
        const int nbead = ZeroXXd.cols();
        const int nchain = ZeroXXXd.rows() / ZeroXXd.rows();    // extented dimension

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double k = uovie::phy_const::Boltzmann_const;
        const double& T = sys.temperature;
        const int& M = nchain;
        const double fic_omega = k * T * sqrt(nbead) / h_bar;

        Eigen::ArrayXXd m = ZeroXXd;            // real mass
        Eigen::ArrayXXd q = ZeroXXd;            // real position
        Eigen::ArrayXXd m_tilde = ZeroXXd;      // transformed mass
        Eigen::ArrayXXd r = ZeroXXd;            // transformed position
        Eigen::ArrayXXd s = ZeroXXd;            // fictitious momentum
        Eigen::ArrayXXd F = ZeroXXd;            // fictitious force
        
        Eigen::ArrayXXd mu = ZeroXXXd;          // extented mass
        Eigen::ArrayXXd eta = ZeroXXXd;         // extented position
        Eigen::ArrayXXd theta = ZeroXXXd;       // extented momentum
        Eigen::ArrayXXd Gamma = ZeroXXXd;       // thermostat force

        Eigen::ArrayXXd kine_energy = ZeroXXd;
        Eigen::ArrayXXd pote_energy = ZeroXXd;
        Eigen::ArrayXXd ther_energy = ZeroXXd;
        Eigen::ArrayXXd cons_energy = ZeroXXd;

        double sys_tot_energy = 0;
        double cano_prob_dens = 0;

        void stag_trans();
        void inve_stag_trans();

        void calc_physic_force();
        void calc_thermo_force(const int& j);
        void physic_propagate();
        void thermo_propagate();
        void calc_syco_energy();

        void print_nhc_procedure_title(std::ofstream& out);
        void print_nhc_procedure_data(std::ofstream& out, double& t);
    };

} // !nhc
} // !thermostat
} // !uovie
#endif // !NOSE_HOOVER_CHAIN_H_