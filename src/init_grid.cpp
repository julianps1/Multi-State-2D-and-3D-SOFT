#include <iostream>
#include "wave.hpp"


//Initialize the masses
namespace wave {
    
    void init_mass(double h0, double h1, double h2)
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
