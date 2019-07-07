// Atomic data
#ifndef ATOM_DATA_H_
#define ATOM_DATA_H_

#include <vector>
#include <string>

namespace uovie {
namespace atom_data {

#define ELEMENT_NUM_MAX 118
    
    static const std::vector<std::string> element{"H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
        "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca", "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co",
        "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr", "Nb", "Mo", "Tc", "Ru",
        "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt",
        "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" , "Np", "Pu", "Am",
        "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
        "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

    static const std::vector<double> atomic_mass_kg{ // -1 represents an uncertain data
        9.10938356e-31, -1, -1, -1, -1, -1, -1, -1, -1, 3.35e-26
    };
    
    
    static const std::vector<double> atomic_weight{ // -1 represents an uncertain data
        1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305,
        26.982, 28.085, 30.974, 32.06, 35.45, 39.88, 39.098, 40.078, 44.956, 47.867, 50.942,
        51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.63, 74.922, 78.971,
        79.904, 83.798, 85.468, 87.62, 88.906, 91.224, 92.906, 95.95, -1, 101.07, 102.91, 106.42,
        107.87, 112.41, 114.82, 118.71, 121.76, 127.6, 126.9, 131.29, 132.91, 137.33, 138.91,
        140.12, 140.91, 144.24, -1, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93,
        173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
        204.38, 207.2, 208.98, -1, -1, -1, -1, -1, -1, 232.04, 231.04, 238.03, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

    static const std::vector<double> atomic_mass{ // unit: Da or u // -1 represents an uncertain data
        1/*1.00782503223*/, 4.00260325413, 7.0160034366, 9.012183065, 11.00930536, 12.0000000,
        14.00307400443, 15.99491461957, 18.99840316273, 19.9924401762, 22.9897692820, 23.985041697,
        26.98153853, 27.97692653465, 30.97376199842, 31.9720711744, 34.968852682, 39.9623831237,
        38.9637064864, 39.962590863, 44.95590828, 47.94794198, 50.94395704, 51.94050623, 54.93804391,
        55.93493633, 58.93319429, 57.93534241, 62.92959772, 63.92914201, 68.9255735, 73.921177761,
        74.92159457, 79.9165218, 78.9183376, 83.9114977282, 84.9117897379, 87.9056125, 88.9058403,
        89.9046977, 92.9063730, 97.90540482, 97.9072124, 101.9043441, 102.9054980, 105.9034804,
        106.9050916, 113.90336509, 114.903878776, 119.90220163, 120.9038120, 129.906222748, 126.9044719,
        131.9041550856, 132.9054519610, 137.90524700, 138.9063563, 139.9054431, 140.9076576,
        141.9077290, 144.9127559, 151.9197397, 152.9212380, 157.9241123, 158.9253547, 163.9291819,
        164.9303288, 165.9302995, 168.9342179, 173.9388664, 174.9407752, 179.9465570, 180.9479958,
        183.95093092, 186.9557501, 191.9614770, 192.9629216, 194.9647917, 196.96656879, 201.97064340,
        204.9744278, 207.9766525, 208.9803991, 208.9824308, 209.9871479, 222.0175782, 223.0197360,
        226.0254103, 227.0277523, 232.0380558, 231.0358842, 238.0507884, 237.0481736, 244.0642053,
        243.0613813, 247.0703541, 247.0703073, 251.0795886, 252.082980, 257.0951061, 258.0984315,
        259.10103, 266.11983, 267.12179, 268.12567, 271.13393, 270.13336, 269.13375, 278.15631,
        281.16451, 282.16912, 285.17712, 286.18221, 289.19042, 289.19363, 293.20449, 294.21046, -1 };

}   // Atom_Data
}   // uovie

#endif