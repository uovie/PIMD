// Nose-Hoover Chain
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <random>
#include <chrono>

#include "simu_para.h"
#include "thermostat/nhc.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    /*** ================================================== ***/
    /*** Thermostat Variables Generator                     ***/
    /*** ================================================== ***/

    void thermo_vari_generator(std::vector<thermo_vari>& tmvs,
        const int M, const double& T, const double tau)
    {
        double mu = k * T * pow(tau, 2);
        std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu) };
        
        for (int j = 0; j < M; j++) {
            tmvs.push_back(thermostat::nhc::thermo_vari(mu, 0, ndrm(nhc_tmv_mte)));
        }
    }

    /*** ================================================== ***/
    /*** Staging Transformation                             ***/
    /*** ================================================== ***/

    Eigen::ArrayXXd stag_trans(Eigen::ArrayXXd& cpc, Eigen::ArrayXXd& tpc)
    {
        tpc.col(0) = cpc.col(0);
        for (int bi = 1; bi < tpc.cols() - 1; bi++)
            tpc.col(bi) = -(1 / (bi + 1)) * cpc.col(bi) + cpc.col(bi) - (bi / (bi + 1)) * cpc.col(bi + 1);
        tpc.col(tpc.cols() - 1) = cpc.col(tpc.cols() - 1) - cpc.col(0);
    }

    /*** ================================================== ***/
    /*** Inverse Staging Transformation                     ***/
    /*** ================================================== ***/

    Eigen::ArrayXXd inve_stag_trans(Eigen::ArrayXXd& cpc, Eigen::ArrayXXd& tpc)
    {
        cpc.col(0) = tpc.col(0);
        cpc.col(cpc.cols() - 1) = tpc.col(tpc.cols() - 1) + tpc.col(0);
        for (int bi = cpc.cols() - 2; bi > 0; bi--)
            cpc.col(bi) = (1 / (bi + 1)) * tpc.col(0) + tpc.col(bi) + (bi / (bi + 1)) * cpc.col(bi + 1);
    }

    /*** ================================================== ***/
    /*** NHC Procedure Base Class Member Functions          ***/
    /*** ================================================== ***/

    void nhc_procedure::initialize()
    {
        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    cpc(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    rmc(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < cpc.cols(); ci++) {
            cpc.col(ci) = cpc.col(0);
            rmc.col(ci) = rmc.col(0);
        }

        // initialize staging transformed positions
        tpc = cpc;
        stag_trans(cpc, tpc);

        // initialize staging transformed masses
        tmc.col(0) = rmc.col(0);
        for (auto ci = 1; ci < tmc.cols(); ci++)
            tmc.col(ci) = ((ci + 1) / ci) * rmc.col(ci);

        // initialize fictition momenta
        for (auto ri = 0; ri < fmc.rows(); ri++) {
            for (auto ci = 0; ci < fmc.cols(); ci++){
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * rmc(ri, ci)) };
                fmc(ri, ci) = ndrm(car_mom_mte);
            }
        }

        // initialize extented masses and extented momenta
        double tau = 1;
        for (auto j = 0; j < M; j++) {
            mu[j] = k * T * pow(tau, 2) * (zero + 1);
            for (auto ri = 0; ri < fmc.rows(); ri++) {
                for (auto ci = 0; ci < fmc.cols(); ci++) {
                    std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu[j](ri, ci)) };
                    theta[j](ri, ci) = ndrm(nhc_tmv_mte);
                }
            }
        }

    }

    // calculate physical forces (simple harmonic oscillator)
    void nhc_procedure::calc_physic_force()
    {
        double omega = 1;

        F.col(0) *= 0;
        for (auto cj = 0; cj < F.cols(); cj++)
            F.col(0) -= rmc.col(0) * pow(omega, 2) * cpc.col(cj);
        F.col(0) /= nbead;

        for (auto ci = 1; ci < F.cols(); ci++)
            F.col(ci) = -tmc.col(ci) * pow(fic_omega, 2) * tpc.col(ci)
            - (rmc.col(ci) * pow(omega, 2) * cpc.col(ci)
                + ((ci - 1) / ci) * rmc.col(ci) * pow(omega, 2) * cpc.col(ci - 1)) / F.cols();
    }

    void nhc_procedure::calc_thermo_force(const int& j)
    {
        if (j == 0) {
            kine_energy = fmc.pow(2) / (2 * rmc);
            Gamma[0] = 2 * kine_energy - k * T;
        }
        else
            Gamma[j] = theta[j - 1] * theta[j - 1] / mu[j - 1] - k * T;
    }

    void nhc_procedure::physic_propagate()
    {
        calc_physic_force();
        fmc += bsp.time_step_size * F / 2;
        tpc += bsp.time_step_size * fmc / tmc;
        calc_physic_force();
        fmc += bsp.time_step_size * F / 2;
    }

    void nhc_procedure::thermo_propagate()
    {
        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * bsp.time_step_size / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta[M - 1] += tmp_delta * Gamma[M - 1] / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta[j] *= exp(-1 * tmp_delta * theta[j + 1]
                        / (8 * mu[j + 1]));
                    calc_thermo_force(j);
                    theta[j] += tmp_delta * Gamma[j] / 4;
                    theta[j] *= exp(-1 * tmp_delta * theta[j + 1]
                        / (8 * mu[j + 1]));
                }
                for (int j = 0; j < M; j++)
                    eta[j] += tmp_delta * theta[j] / (2 * mu[j]);

                Eigen::ArrayXXd momen_scale = (-1 * tmp_delta * theta[0] / (2 * mu[0])).exp();
                fmc *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta[j] *= exp(-1 * tmp_delta * theta[j + 1]
                        / (8 * mu[j + 1]));
                    calc_thermo_force(j);
                    theta[j] += tmp_delta * Gamma[j] / 4;
                    theta[j] *= exp(-1 * tmp_delta * theta[j + 1]
                        / (8 * mu[j + 1]));
                }
                calc_thermo_force(M - 1);
                theta[M - 1] += tmp_delta * Gamma[M - 1] / 4;
            }
        }
    }

    void nhc_procedure::calc_syco_energy()
    {
        double omega = 1;

        kine_energy = fmc.pow(2) / (2 * tmc);

        pote_energy.col(0) *= 0;
        for (auto cj = 0; cj < pote_energy.cols(); cj++)
            pote_energy.col(0) += 0.5 * rmc.col(0) * pow(omega, 2) * cpc.col(cj).pow(2);
        pote_energy.col(0) /= nbead;

        for (auto ci = 1; ci < pote_energy.cols(); ci++)
            pote_energy.col(ci) += 0.5 * tmc.col(ci) * pow(fic_omega, 2) * tpc.col(ci).pow(2)
            + (0.5 * rmc.col(ci) * pow(omega, 2) * cpc.col(ci).pow(2)
                + 0.5 * ((ci - 1) / ci) * rmc.col(ci) * pow(omega, 2) * cpc.col(ci - 1).pow(2)) / nbead;

        ther_energy *= 0;
        for (int j = 0; j < M; j++)
            ther_energy += theta[j].pow(2) / (2 * mu[j]) + k * T * eta[j];

        cons_energy = kine_energy + pote_energy + ther_energy;

        double sys_kine_energy = (0.5 * fmc.pow(2) / rmc).sum();
        double sys_pote_energy = (0.5 * rmc * pow(omega, 2) * cpc.pow(2)).sum();

        sys_tot_energy = sys_kine_energy + sys_kine_energy;
        cano_prob_dens = exp(-sys_tot_energy / (k * T));
    }

    void nhc_procedure::print_nhc_procedure_title(std::ofstream& out) {
        std::cout << "\nNHC Procedure:\n   Time" << "  " << "position" << "  " << "momentum" << "  "
            << "cons_energy" << "  " << "cano_prob_dens";
        out << "\nNHC Procedure:\n   Time" << "  " << "position" << "  " << "momentum" << "  "
            << "cons_energy" << "  " << "cano_prob_dens";
    }

    void nhc_procedure::print_nhc_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::setw(15) << tpc(0, 0) << std::setw(15) << fmc(0, 0)
            << std::setw(15) << cons_energy << std::setw(15) << cano_prob_dens;
        out << std::setw(15) << tpc(0, 0) << std::setw(15) << fmc(0, 0)
            << std::setw(15) << cons_energy << std::setw(15) << cano_prob_dens;
    }

    void nhc_procedure::implement() {
        for (double t = 0; t <= bsp.run_time; t += bsp.time_step_size) {
            thermo_propagate();
            physic_propagate();
            thermo_propagate();
            t += bsp.time_step_size;
        }
    }
    
    void nhc_procedure::implement(std::ofstream& out) {
        double t = 0;
        int ctr = 0;
        
        print_nhc_procedure_title(out);
        calc_syco_energy();
        print_nhc_procedure_data(out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            thermo_propagate();
            physic_propagate();
            thermo_propagate();
            inve_stag_trans(cpc, tpc);
            
            if (ctr == bsp.data_coll_peri) {
                calc_syco_energy();
                print_nhc_procedure_data(out, t);
                ctr = 0;
            }
            ctr++;
            t += bsp.time_step_size;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of NHC procedure (s): " << time_elapsed.count() << std::endl;
    }

} // !nhc
} // !thermostat
} // !uovie