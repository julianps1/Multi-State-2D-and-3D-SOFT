#include "globals.hpp"
#include <cmath>

namespace wave {

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
 }
}

