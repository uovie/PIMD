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
    /*** Thermostat Variables                               ***/
    /*** ================================================== ***/

    class thermo_vari {
    public:
        thermo_vari() = default;
        thermo_vari(double m, double e, double t):
            mu(m), eta(e), theta(t) { }

        double mu;                      // extented mass
        double eta;                     // extented position
        double theta;                   // extented momentum
        double Gamma;                   // thermostat force
    };

    
    void thermo_vari_generator(std::vector<thermo_vari>& tmvs,
        const int M, const double& T, const double tau);
    
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
    /*** NHC Procedure                                      ***/
    /*** ================================================== ***/

    Eigen::ArrayXXd stag_trans(Eigen::ArrayXXd& cart_pos_coll, Eigen::ArrayXXd& tran_pos_coll);
    Eigen::ArrayXXd inve_stag_trans(Eigen::ArrayXXd& cart_pos_coll, Eigen::ArrayXXd& tran_pos_coll);

    std::mt19937 nhc_tmv_mte(27);
    std::mt19937 car_mom_mte(36);

    class nhc_procedure {
    public:
        nhc_procedure() = default;
        nhc_procedure(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const Eigen::ArrayXXd& _zero, const thermo_factor_scheme& _tfs, const int _M) :
            bsp(_bsp), sys(_sys), zero(_zero), tfs(_tfs), M(_M) { }
        
        void initialize();
        void implement();
        void implement(std::ofstream& out);
        //double sys_ene() const { return syst_energy; }

    protected:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const Eigen::ArrayXXd& zero;
        const thermo_factor_scheme& tfs;
        const int M;                        // extented dimension
        const int nbead = zero.cols();

        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const double k = uovie::phy_const::Boltzmann_const;
        const double& T = sys.temperature;
        const double fic_omega = k * T * sqrt(nbead) / h_bar;

        
        Eigen::ArrayXXd cpc = zero;         // cart_pos_coll
        Eigen::ArrayXXd tpc = zero;         // tran_pos_coll
        Eigen::ArrayXXd rmc = zero;         // real_mas_coll
        Eigen::ArrayXXd tmc = zero;         // tran_mas_coll
        Eigen::ArrayXXd fmc = zero;         // fict_mom_coll
        Eigen::ArrayXXd F = zero;
        
        std::vector<Eigen::ArrayXXd> mu;          // extented mass
        std::vector<Eigen::ArrayXXd> eta;         // extented position
        std::vector<Eigen::ArrayXXd> theta;       // extented momentum
        std::vector<Eigen::ArrayXXd> Gamma;       // thermostat force

        Eigen::ArrayXXd kine_energy = zero;
        Eigen::ArrayXXd pote_energy = zero;
        Eigen::ArrayXXd ther_energy = zero;
        Eigen::ArrayXXd cons_energy = zero;

        double sys_tot_energy = 0;
        double cano_prob_dens = 0;

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