#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include "FFT2D.hpp"
#include <fstream>
#include <string>
#include <filesystem>
#include <cmath>
#include <fftw3.h>

namespace grid {


void kinetic_step(fft::FFT2D &f2d, GridData &gridd, int n)
{
    double N = 1.0/(ni*nj);

    for (int i = 0; i<ni; ++i)
    {
    for (int j = 0; j<nj; ++j)
    {
	gridd.tmp[i][j] = gridd.psi[i][j][n];
    }
    }

    fft::run(f2d, gridd.tmp, 1); //Forward FFT, output goes to tspi

    for (int i = 0; i<ni; ++i)
    {
    for (int j = 0; j<nj; ++j)
    {
	gridd.tpsi[i][j] = gridd.tpsi[i][j]*gridd.ak2[i][j];
    }
    }

    fft::run(f2d, gridd.tpsi, -1); //Backward FFT, output goes to tspi

    for (int i = 0; i<ni; ++i)
    {
    for (int j = 0; j<nj; ++j)
    {
        gridd.tpsi[i][j] = gridd.tpsi[i][j] * N;
    }
    }
}

cplx overlapgg(const GridData &gridd, int n, const cplx wf2[][nj])
{
    cplx z = {0.0, 0.0};
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
      z += std::conj(gridd.psi[i][j][n]) * wf2[i][j];
    }
    }

    return z;
}

cplx overlapggs(const GridData &gridd, int n, const cplx wf2[][nj][ns])
{
    cplx z = {0.0, 0.0};
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
      z += std::conj(gridd.psi[i][j][n]) * wf2[i][j][n];
    }
    }

    return z;
}

void Ham(fft::FFT2D &f2d, GridData &gridd, int pot_name)
{
  pick_pot(pot_name, gridd);
 
  cplx zkin = {0.0, 0.0};
  cplx zpot = {0.0, 0.0}; 
	
  for (int n = 0; n<ns; ++n) {
   kinetic_step(f2d, gridd, n);
   
   zkin = overlapgg(gridd, n, gridd.tpsi);
   zpot = overlapggs(gridd, n, gridd.vpsi);


   std::cout << "Printing Initial Energy on State," << n << '\n'
	     << "<K> = " << zkin
             << "<V> = " << zpot
	     << "<H> = " << zkin+zpot << '\n'; 
    }

}

}

