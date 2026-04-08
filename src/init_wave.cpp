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

        double amh  = ame * 1.008 * h1;
        double amd  = ame * 1.008 * h2;
        double amet = ame * h0;

        amc = amet * (amh + amd) / (amet + amh + amd);
        amx = amh * amd / (amh + amd);
        amy = amx;

        std::cout << "REDUCED MASS OF H2 " << amx << std::endl;
        std::cout << "MASS OF THE COM " << amc << std::endl;
    }

}
