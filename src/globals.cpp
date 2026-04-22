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
    {"alfa", &wave::alfa},
    {"beta", &wave::beta},
    {"q0",   &wave::q0},
    {"y0",   &wave::y0},
    {"px0",   &wave::px0},
    {"py0",   &wave::py0},
    {"dt",   &wave::dt},
    {"h0",   &wave::h0},
    {"h1",   &wave::h1},
    {"h2",   &wave::h2},
    {"rCAP", &grid::rCAP},
    {"kCAP", &grid::kCAP},
};

static std::unordered_map<std::string,int*> int_params = {
    {"tsub",   &wave::tsub},
    {"tmax",   &wave::tmax},
    {"nwpackets",   &wave::nwpackets},
    {"istate", &grid::istate},
    {"pot_name", &grid::pot_name},
    {"nmom", &grid::nmom},
    {"DN", &grid::DN},
}; 

void init(const std::string &fname)
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
