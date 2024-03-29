// Nose-Hoover Chain
#include <iostream>
#include <iomanip>
#include <cassert>
#include <chrono>
#include <random>

// uovie headers
#include "thermostat/nhc.h"
#include "model.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    std::mt19937 s_mte(36);
    std::mt19937 theta_mte(27);

    /*** ================================================== ***/
    /*** NHC Procedure for PIMD Class Member Functions      ***/
    /*** ================================================== ***/

    // initialization
    void nhc_procedure_for_pimd::initialize()
    {
        // resize arrays
        m.resize(d * N, nbead);
        q.resize(d * N, nbead);
        m_tilde.resize(d * N, nbead);
        r.resize(d * N, nbead);
        s.resize(d * N, nbead);
        F.resize(d * N, nbead);

        mu.resize(d * N * M, nbead);
        eta.resize(d * N * M, nbead);
        theta.resize(d * N * M, nbead);
        Gamma.resize(d * N * M, nbead);

        kine_energy.resize(d * N, nbead);
        pote_energy.resize(d * N, nbead);
        ther_energy.resize(d * N, nbead);
        cons_energy.resize(d * N, nbead);

        // initialize cartisian positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi, 0) = sys.molecules[mi].atoms[ai].m;
                    q(vi, 0) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }
        for (auto ci = 1; ci < q.cols(); ci++) {
            q.col(ci) = q.col(0);
            m.col(ci) = m.col(0);
        }

        // initialize staging transformed masses
        m_tilde.col(0) = m.col(0);
        for (auto ci = 1; ci < m_tilde.cols(); ci++)
            m_tilde.col(ci) = ((ci + 1.0) / ci) * m.col(ci);

        m_bar = m_tilde;
        m_bar.col(0).setZero();

        // initialize staging transformed positions
        stag_trans();

        // initialize fictition momenta
        for (auto ri = 0; ri < s.rows(); ri++) {
            for (auto ci = 0; ci < s.cols(); ci++){
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri, ci)) };
                s(ri, ci) = ndrm(s_mte);
            }
        }

        // initialize extented masses and extented momenta
        double tau = 1;
        mu = k * T * pow(tau, 2) * Eigen::ArrayXXd::Constant(d * N * M, nbead, 1);
        for (auto ri = 0; ri < theta.rows(); ri++) {
            for (auto ci = 0; ci < theta.cols(); ci++) {
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * mu(ri, ci)) };
                theta(ri, ci) = ndrm(theta_mte);
            }
        }

        // initialize extented positions
        eta.setZero();

    }

    // staging transformation
    void nhc_procedure_for_pimd::stag_trans()
    {
        r.col(0) = q.col(0);
        for (int ci = 1; ci < r.cols() - 1; ci++)
            r.col(ci) = -(1 / (ci + 1.0)) * q.col(0) + q.col(ci) - (ci / (ci + 1.0)) * q.col(ci + 1);
        r.col(r.cols() - 1) = q.col(r.cols() - 1) - q.col(0);
    }

    // inverse staging transformation
    void nhc_procedure_for_pimd::inve_stag_trans()
    {
        q.col(q.cols() - 1) = r.col(q.cols() - 1) + r.col(0);
        for (int ci = q.cols() - 2; ci > 0; ci--)
            q.col(ci) = (1 / (ci + 1.0)) * r.col(0) + r.col(ci) + (ci / (ci + 1.0)) * q.col(ci + 1);
        q.col(0) = r.col(0);
    }

    // calculate physical forces
    void nhc_procedure_for_pimd::calc_physic_force()
    {
        Eigen::ArrayXXd pV_pq = Eigen::ArrayXXd::Zero(d * N, nbead); // \frac{\partial V}{\partial q}
        Eigen::ArrayXXd pV_pr = Eigen::ArrayXXd::Zero(d * N, nbead); // \frac{\partial V}{\partial r}

        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pV_pq = m * pow(HO.ome(), 2) * q;
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ci = 0; ci < pV_pq.cols(); ci++) {
                for (auto ni = 0; ni < N; ni++) {
                    for (auto nj = 0; nj < N; nj++) {
                        if (ni == nj) continue;
                        pV_pq.block(d * ni, ci, d, 1) -= LJ.F(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
                    }
                }
            }
        }

        for (auto cj = 0; cj < pV_pr.cols(); cj++)
            pV_pr.col(0) += pV_pq.col(cj);
        for (auto ci = 1; ci < pV_pr.cols(); ci++)
            pV_pr.col(ci) = pV_pq.col(ci) + ((ci - 1.0) / ci) * pV_pr.col(ci - 1);
        F = -1 * m_bar * pow(fic_omega, 2) * r - pV_pr / nbead;

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

    void nhc_procedure_for_pimd::calc_cons_quant()
    {
        kine_energy = s.pow(2) / (2 * m_tilde);

        
        pote_energy.setZero();
        if (sys.model_type == "HO") { // harmonic oscillator
            model::harmonic_oscilator HO(sys.model_para[0]);
            pote_energy = 0.5 * m * pow(HO.ome(), 2) * q.pow(2);
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ci = 0; ci < pote_energy.cols(); ci++) {
                for (auto ni = 0; ni < N; ni++) {
                    for (auto nj = 0; nj < N; nj++) {
                        if (ni == nj) continue;
                        pote_energy.block(d * ni, ci, d, 1) += LJ.V(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
                    }
                }
            }
            pote_energy /= 6;
        }
        pote_energy = 0.5 * m_bar * pow(fic_omega, 2) * r.pow(2) + pote_energy / nbead;

        ther_energy.setZero();
        for (int j = 0; j < M; j++)
            ther_energy += theta.block(j * dof, 0, dof, nbead).pow(2) / (2 * mu.block(j * dof, 0, dof, nbead))
            + k * T * eta.block(j * dof, 0, dof, nbead);

        cons_energy = kine_energy + pote_energy + ther_energy;
    }

    void nhc_procedure_for_pimd::calc_prim_estor()
    {
        prim_kine_estor = 0;
        for (int bi = 0; bi < nbead - 1; bi++)
            prim_kine_estor += (m.col(bi) * (q.col(bi + 1) - q.col(bi)).pow(2)).sum();
        prim_kine_estor += (m.col(nbead - 1) * (q.col(0) - q.col(nbead - 1)).pow(2)).sum();
        prim_kine_estor *= -1 * nbead / (2 * pow(h_bar * beta, 2));
        prim_kine_estor += dof * nbead / (2 * beta);

        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pote_energy = 0.5 * m * pow(HO.ome(), 2) * q.pow(2);
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ci = 0; ci < pote_energy.cols(); ci++) {
                for (auto ni = 0; ni < N; ni++) {
                    for (auto nj = 0; nj < N; nj++) {
                        if (ni == nj) continue;
                        pote_energy.block(d * ni, ci, d, 1) += LJ.V(q.block(d * ni, ci, d, 1), q.block(d * nj, ci, d, 1));
                    }
                }
            }
            pote_energy /= 2;
        }
        prim_pote_estor = pote_energy.sum() / nbead;

        //prim_pres_estor = N*nbead/(beta * V)
    }

    void nhc_procedure_for_pimd::print_nhc_procedure_title(std::ofstream& out) {
        std::cout << "\nNHC Procedure for PIMD:\n   Time" << "            " << "position" << "            " << "momentum"
            << "          " << "cons_energy" << "        " << "prim_kine_estor" << "     " << "prim_pote_estor";
        out << "\nNHC Procedure for PIMD:\n   Time" << "            " << "position" << "            " << "momentum"
            << "          " << "cons_energy" << "        " << "prim_kine_estor" << "     " << "prim_pote_estor";
    }

    void nhc_procedure_for_pimd::print_nhc_procedure_data(std::ofstream& out, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        std::cout << std::scientific << std::setw(20) << r(0, 0) << std::setw(20) << s(0, 0) << std::setw(20)
            << cons_energy.sum() << std::setw(20) << prim_kine_estor << std::setw(20) << prim_pote_estor;
        out << std::scientific << std::setw(20) << r(0, 0) << std::setw(20) << s(0, 0) << std::setw(20)
            << cons_energy.sum() << std::setw(20) << prim_kine_estor << std::setw(20) << prim_pote_estor;
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
        calc_cons_quant();
        calc_prim_estor();
        print_nhc_procedure_data(out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            // NHC numerical evolution
            thermo_propagate();
            physic_propagate();
            thermo_propagate();
            
            if (ctr == bsp.data_coll_peri) {
                calc_cons_quant();
                calc_prim_estor();
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