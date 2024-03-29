// standard C++ headers
#include <iostream>
#include <sstream>
#include <random>
#include <cstdlib>
#include <cassert>

// Piano headers
#include "process.h"
#include "phy_const.h"
#include "atom_data.h"

#ifndef PHY_CONST_SHORTHAND
#define PHY_CONST_SHORTHAND
constexpr double h_bar = uovie::phy_const::red_Planck_const;
constexpr double k = uovie::phy_const::Boltzmann_const;
#endif // !PHY_CONST_SHORTHAND

using namespace uovie::Global;

void process::open(const std::string& filename) {
    // check the extension
    if (filename.rfind(".vie") == std::string::npos)
        throw std::invalid_argument("Only .vie file is a valid input file.");

    in.open(filename);

    if (in.fail()) {
        std::cout << "Can not open the file " << filename << '.' << std::endl;
        exit(EXIT_FAILURE);
    }

    for (auto it = filename.begin(); it != filename.end() - 4; it++)
        fn_no_ex += *it;

    // open output file
    out.open(fn_no_ex + ".out");
}

void process::read() {

    constexpr double angstrom_to_bohr = 1e-10 / phy_const::a_u_length;

    std::cout << "Read Infomation from " << fn_no_ex + ".vie" << std::endl;
    in >> job >> bsp.run_time >> bsp.step_size >> bsp.data_coll_peri;
    in >> sys.dimension >> sys.volume >> sys.temperature >> sys.pressure;
    in >> sys.model_type;

    if (sys.model_type == "HO") {
        double tmp_data;
        in >> tmp_data;
        sys.model_para.push_back(tmp_data);
    }
    else if (sys.model_type == "LJ") {
        double tmp_data;
        for (int i = 0; i < 2; i++) {
            in >> tmp_data;
            sys.model_para.push_back(tmp_data);
        }
    }else
        throw "unsupported model";

    const double& T = sys.temperature;
    
    int mol_ctr = 0; int& mi = mol_ctr;
    std::mt19937 mte(36);
    std::string line;
    molecule tmp_mole;
    
    while (getline(in, line)) {
        std::string str;
        if (line.size() >= 5) {
            for (int a = 0; a < 5; a++)
                str.push_back(line[a]);
        }
        if (str == "*****") {
            int natom;
            in >> natom;
            sys.molecules.push_back(tmp_mole);
            for (int ai = 0; ai < natom; ai++) {
                atom tmp_atom;
                sys.molecules[mi].atoms.push_back(tmp_atom);
                in >> sys.molecules[mi].atoms[ai].symbol;
                for (int di = 0; di < sys.dimension; di++) {
                    double tmp_q;
                    sys.molecules[mi].atoms[ai].q.push_back(tmp_q);
                    in >> sys.molecules[mi].atoms[ai].q[di];
                    sys.molecules[mi].atoms[ai].q[di] /= phy_const::a_u_length * 1e10;
                }
            }
            // identify atomic numbers and assign atomic masses, momentums, forces
            for (int ai = 0; ai < natom; ai++) {
                for (int ei = 0; ei < atom_data::element.size(); ei++) {
                    if (sys.molecules[mi].atoms[ai].symbol == atom_data::element[ei]) {
                        sys.molecules[mi].atoms[ai].atomic_number = ei + 1;
                        sys.molecules[mi].atoms[ai].m = atom_data::atomic_mass_kg[ei] / phy_const::a_u_mass;
                        break;
                    }
                }
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * sys.molecules[mi].atoms[ai].m) };
                for (int di = 0; di < sys.dimension; di++) {
                    sys.molecules[mi].atoms[ai].p.push_back(ndrm(mte));
                    sys.molecules[mi].atoms[ai].F.push_back(0);
                }
            }
            mol_ctr++;
        }
    }

    sys.num_part = 0;
    for (auto mi = 0; mi < sys.molecules.size(); mi++)
        sys.num_part += sys.molecules[mi].atoms.size();

}

void process::print() {
    std::cout << "\nNormal termination. Congratulations!" << std::endl;
}

void process::close() {
    // close files
    in.close();
    out.close();
}