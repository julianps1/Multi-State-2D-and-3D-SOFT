#include <iostream>
#include <complex>
#include <cmath>
#include "wave.hpp"

namespace wave {
    double alfa,beta,q0,y0,px0,py0,h0,h1,h2,amx,amc,amy,dt,dtsub,dt2;
    int tmax,tsub,nwpackets;
}

//Initialize the masses
namespace wave {

    void init_mass_2D(double h0, double h1)
    //Initialize a generic x and y mass for a 2D system
    {

        amx = h0;
        amy = h1;

        std::cout << "MASS OF PARTICLE X: " << amx << std::endl;
        std::cout << "MASS OF PARTICLE Y: " << amy << std::endl;
    }
    
    void init_mass(double h0, double h1, double h2)
    //For atomic system like NaHF 
    {
        double ame = 1836.15;

        double amf  = ame * h1;
        double amh  = ame * h2;
        double amna = ame * h0;

        amy = amna * (amh + amf) / (amna + amh + amf); //Center of mass
        //amx = amna * amh / (amna + amh); // Reduced mass of NaH
        amx = amh * amf / (amh + amf); // Reduced mass of HF

        std::cout << "REDUCED MASS OF HF " << amx << std::endl;
        std::cout << "CENTER OF MASS " << amy << std::endl;
    }

}
