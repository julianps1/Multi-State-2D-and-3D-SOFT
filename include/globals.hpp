#pragma once
#include <complex>

constexpr int ni = 256, nj = 256, ns = 2;
constexpr char defaultInputFile[] = "INq";

namespace dyn {
    using real = double;
    using cplx = std::complex<double>;

    extern real ak2[ni][nj];
    extern real xg[ni];
    extern real yg[nj];

    extern real v[ni][nj][ns][ns];
    extern real vcos[ni][nj];
    extern real vsin[ni][nj];

    extern cplx c[ns];
    extern cplx vd[ni][nj][ns];

    extern cplx im;
    extern cplx psi[ni][nj][ns];
    extern cplx ksi[ni][nj][ns];
    extern cplx psi0[ni][nj][ns];

    void init();            // set im, compute pi, N2, …
    void read_input(const std::string &fname);
}

