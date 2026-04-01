#pragma once
#include <string>

namespace wave {
    extern double alfa, beta, q0, y0, px0, py0, dt, dt2, dtsub, h0, h1, h2, amx, amc, amy;
    extern int tmax,tsub,nwpackets;

    void init_mass(double h0, double h1, double h2);
    void init_mass_electron();
}
