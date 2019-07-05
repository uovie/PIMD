// Nose-Hoover Chain
#include <iostream>
#include <iomanip>
#include <cassert>
#include <chrono>
#include <random>

// uovie headers
#include "thermostat/nhc.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    std::mt19937 nhc_tmv_mte(27);
    std::mt19937 car_mom_mte(36);

    /*** ================================================== ***/
    /*** NHC Procedure for PIMD Class Member Functions      ***/
    /*** ================================================== ***/

    // initialization
    void nhc_procedure_for_pimd::initialize()
    {
        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    q(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    m(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < q.cols(); ci++) {
            q.col(ci) = q.col(0);
            m.col(ci) = m.col(0);
        }

        // initialize staging transformed positions
        r = q;
        stag_trans();

        // initialize staging transformed masses
        m_tilde.col(0) = m.col(0);
        for (auto ci = 1; ci < m_tilde.cols(); ci++)
            m_tilde.col(ci) = ((ci + 1.0) / ci) * m.col(ci);

        // initialize fictition momenta
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++){
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(car_mom_mte);
            }
        }

        // initialize extented masses and extented momenta
        double tau = 1;
        mu = k * T * pow(tau, 2) * (ZeroXXXd + 1);
        for (auto ri = 0; ri < theta.rows(); ri++) {
            for (auto ci = 0; ci < theta.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu(ri, ci)) };
                theta(ri, ci) = ndrm(nhc_tmv_mte);
            }
        }

    }

    // staging transformation
    void nhc_procedure_for_pimd::stag_trans()
    {
        r.col(0) = q.col(0);
        for (int ci = 1; ci < r.cols() - 1; ci++)
            r.col(ci) = -(1 / (ci + 1.0)) * q.col(ci) + q.col(ci) - (ci / (ci + 1.0)) * q.col(ci + 1);
        r.col(r.cols() - 1) = q.col(r.cols() - 1) - q.col(0);
    }

    // inverse staging transformation
    void nhc_procedure_for_pimd::inve_stag_trans()
    {
        q.col(0) = r.col(0);
        q.col(q.cols() - 1) = r.col(r.cols() - 1) + r.col(0);
        for (int ci = q.cols() - 2; ci > 0; ci--)
            q.col(ci) = (1 / (ci + 1.0)) * r.col(0) + r.col(ci) + (ci / (ci + 1.0)) * q.col(ci + 1);
    }

    // calculate physical forces (simple harmonic oscillator)
    void nhc_procedure_for_pimd::calc_physic_force()
    {
        double omega = 1;

        F.col(0).setZero();
        for (auto cj = 0; cj < F.cols(); cj++)
            F.col(0) -= m.col(cj) * pow(omega, 2) * q.col(cj);
        F.col(0) /= nbead;

        for (auto ci = 1; ci < F.cols(); ci++)
            F.col(ci) = -1 * m_tilde.col(ci) * pow(fic_omega, 2) * r.col(ci)
            - (m.col(ci) * pow(omega, 2) * q.col(ci)
                + ((ci - 1.0) / ci) * m.col(ci - 1) * pow(omega, 2) * q.col(ci - 1)) / nbead;
    }

    void nhc_procedure_for_pimd::calc_thermo_force(const int& j)
    {
        if (j == 0)
            Gamma.block(0, 0, dof, nbead) = s.pow(2) / m_tilde - k * T;
        else
            Gamma.block(j * dof, 0, dof, nbead) = theta.block((j - 1) * dof, 0, dof, nbead).pow(2)
            / mu.block((j - 1) * dof, 0, dof, nbead) - k * T;
    }

    void nhc_procedure_for_pimd::physic_propagate()
    {
        calc_physic_force();
        s += Dt * F / 2;
        r += Dt * s / m_tilde;
        inve_stag_trans();
        calc_physic_force();
        s += Dt * F / 2;
    }

    void nhc_procedure_for_pimd::thermo_propagate()
    {
        for (auto wi = 0; wi < tfs.nsy(); wi++) {
            double tmp_delta = tfs.w(wi) * Dt / tfs.nff();
            for (auto ni = 0; ni < tfs.nff(); ni++) {

                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;

                for (int j = M - 2; j >= 0; j--) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                eta += tmp_delta * theta / (2 * mu);

                Eigen::ArrayXXd momen_scale = (-1 * tmp_delta * theta.block(0, 0, dof, nbead)
                    / (2 * mu.block(0, 0, dof, nbead))).exp();
                s *= momen_scale;

                for (int j = 0; j < M - 1; j++) {
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                    calc_thermo_force(j);
                    theta.block(j * dof, 0, dof, nbead) += tmp_delta * Gamma.block(j * dof, 0, dof, nbead) / 4;
                    theta.block(j * dof, 0, dof, nbead) *= (-1 * tmp_delta * theta.block((j + 1) * dof, 0, dof, nbead)
                        / (8 * mu.block((j + 1) * dof, 0, dof, nbead))).exp();
                }
                calc_thermo_force(M - 1);
                theta.block((M - 1) * dof, 0, dof, nbead) += tmp_delta * Gamma.block((M - 1) * dof, 0, dof, nbead) / 4;
            }
        }
    }

    void nhc_procedure_for_pimd::calc_syco_energy()
    {
        double omega = 1;

        kine_energy = s.pow(2) / (2 * m_tilde);

        pote_energy.col(0).setZero();
        for (auto cj = 0; cj < pote_energy.cols(); cj++)
            pote_energy.col(0) += 0.5 * m.col(0) * pow(omega, 2) * q.col(cj).pow(2);
        pote_energy.col(0) /= nbead;

        for (auto ci = 1; ci < pote_energy.cols(); ci++)
            pote_energy.col(ci) += 0.5 * m_tilde.col(ci) * pow(fic_omega, 2) * r.col(ci).pow(2)
            + (0.5 * m.col(ci) * pow(omega, 2) * q.col(ci).pow(2)
                + 0.5 * ((ci - 1.0) / ci) * m.col(ci) * pow(omega, 2) * q.col(ci - 1).pow(2)) / nbead;

        ther_energy.setZero();
        for (int j = 0; j < M; j++)
            ther_energy += theta.block(j * dof, 0, dof, nbead).pow(2) / (2 * mu.block(j * dof, 0, dof, nbead))
            + k * T * eta.block(j * dof, 0, dof, nbead);

        cons_energy = kine_energy + pote_energy + ther_energy;

        double sys_kine_energy = (0.5 * s.pow(2) / m).sum();
        double sys_pote_energy = (0.5 * m * pow(omega, 2) * q.pow(2)).sum();

        sys_tot_energy = sys_kine_energy + sys_kine_energy;
        cano_prob_dens = exp(-sys_tot_energy / (k * T));
    }

    void nhc_procedure_for_pimd::print_nhc_procedure_title(std::ofstream& out) {
        std::cout << "\nNHC Procedure:\n   Time" << "  " << "position" << "  " << "momentum" << "  "
            << "cons_energy" << "  " << "cano_prob_dens";
        out << "\nNHC Procedure:\n   Time" << "  " << "position" << "  " << "momentum" << "  "
            << "cons_energy" << "  " << "cano_prob_dens";
    }

    void nhc_procedure_for_pimd::print_nhc_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::setw(15) << r(0, 0) << std::setw(15) << s(0, 0)
            << std::setw(15) << cons_energy(0, 0) << std::setw(15) << cano_prob_dens;
        out << std::setw(15) << r(0, 0) << std::setw(15) << s(0, 0)
            << std::setw(15) << cons_energy(0, 0) << std::setw(15) << cano_prob_dens;
    }

    void nhc_procedure_for_pimd::implement() {
        for (double t = 0; t <= bsp.run_time; t += Dt) {
            thermo_propagate();
            physic_propagate();
            thermo_propagate();
            t += Dt;
        }
    }
    
    void nhc_procedure_for_pimd::implement(std::ofstream& out)
    {
        initialize();

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
            
            if (ctr == bsp.data_coll_peri) {
                calc_syco_energy();
                print_nhc_procedure_data(out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of NHC procedure (s): " << time_elapsed.count() << std::endl;
    }

} // !nhc
} // !thermostat
} // !uovie