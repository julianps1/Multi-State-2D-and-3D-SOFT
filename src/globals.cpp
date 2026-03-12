#include "globals.hpp"
#include "grid.hpp"
#include "wave.hpp"
#include <cmath>
#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>

static std::unordered_map<std::string,double*> double_params = {
    {"xmin", &grid::xmin},
    {"ymin", &grid::ymin},
    {"xmax", &grid::xmax},
    {"ymax", &grid::ymax},
    {"dx",   &grid::dx},
    {"dy",   &grid::dy},
    {"an0",  &grid::an0},
    {"alfa", &wave::alfa},
    {"beta", &wave::beta},
    {"q0",   &wave::q0},
    {"y0",   &wave::y0},
    {"px0",   &wave::px0},
    {"py0",   &wave::py0},
    {"dt",   &wave::dt},
};

static std::unordered_map<std::string,int*> int_params = {
    {"tsub",   &wave::tsub},
    {"tmax",   &wave::tmax},
    {"nwpackets",   &wave::nwpackets},
}; 

namespace dyn {

real ak2[ni][nj];
real xg[ni];
real yg[nj];

real v[ni][nj][ns][ns];
real vcos[ni][nj];
real vsin[ni][nj];

cplx c[ns];
cplx vd[ni][nj][ns];

cplx im;
cplx psi[ni][nj][ns];
cplx ksi[ni][nj][ns];
cplx psi0[ni][nj][ns]; //For now grid size is fixed

void init()
 {
    im = {0.0, 1.0};
    constexpr real pi = 2.0 * std::asin(1.0);
    // …other setup, maybe read from file…
 };

void read_input(const std::string &fname)
 {
    std::ifstream ifs(fname);
    if (!ifs) throw std::runtime_error("cannot open " + fname);

    std::string name;
    double val;
    while (ifs >> name >> val) {
     auto itd = double_params.find(name);
      if (itd != double_params.end()) {
          *(itd->second) = val;
           }
      else {
          auto iti = int_params.find(name);
          if (iti != int_params.end()) {
              *(iti->second) = static_cast<int>(val);
           }
      else {
            std::cerr << "Unknown parameter\n";
           } 
      }
    }
  }
}
