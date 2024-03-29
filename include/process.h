// Process
#ifndef PROCESS_H_
#define PROCESS_H_

// standard C++ headers
#include <fstream>
#include <string>

// uovie headers
#include "mol_geom.h"
#include "simu_para.h"

namespace uovie {
namespace Global {

    class process {
    public:
        std::string fn_no_ex;   // filename no extension
        std::ifstream in;
        std::ofstream out;
        std::string job;
        basic_simu_para bsp;
        system sys;

        void open(const std::string& filename);
        void read();
        void print();
        void close();
    };
    
} // !Global
} // !uovie
#endif // !PROCESS_H_