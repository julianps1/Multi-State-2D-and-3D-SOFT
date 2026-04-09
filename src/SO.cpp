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

void apply_kinetic_phase(cplx wf[][nj], const GridData &gridd, double ts)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        cplx phase = std::exp(im * ts * gridd.ak2[i][j]);
         wf[i][j] *= phase;

    }
    }
}
void apply_potential_phase(GridData &gridd, double ts)
{
    //For a rotation matrix U:=[cos,-sin,sin,cos]
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        
        cplx z1 = (gridd.vcos[i][j] * gridd.psi[i][j][0] - gridd.vsin[i][j] * gridd.psi[i][j][1]);
        z1 = z1 * std::exp(im * gridd.vd[i][j][0] * ts);
        cplx z2 = (gridd.vsin[i][j] * gridd.psi[i][j][0] + gridd.vcos[i][j] * gridd.psi[i][j][1]);
        z2 = z2 * std::exp(im * gridd.vd[i][j][1] * ts);
        
        gridd.psi[i][j][0] = gridd.vcos[i][j] * z1 + gridd.vsin[i][j] * z2;
        gridd.psi[i][j][1] = -gridd.vsin[i][j] * z1 + gridd.vcos[i][j] * z2;

    }
    }
}
void apply_1_state_potential_phase(GridData &gridd, double ts)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        
        gridd.psi[i][j][0] = gridd.psi[i][j][0] * std::exp(im * gridd.v[i][j][0][0] * ts);

    }
    }
}
void split(double ts, double ts2, GridData &gridd, fft::FFT2D &f2d) //Propagate the wavefunction (gridd.psi) for one time step using the split-operator method
{
    double N = 1.0/(ni*nj);
    for (int n = 0; n < ns; ++n)
    {
        for (int i = 0; i<ni; ++i)
        {
        for (int j = 0; j<nj; ++j)
         {
	    gridd.tmp[i][j] = gridd.psi[i][j][n];
         }
         }
    
        fft::run(f2d, gridd.tmp, 1); //Forward FFT, output goes to tspi
        apply_kinetic_phase(gridd.tpsi, gridd, ts2);
        fft::run(f2d, gridd.tpsi, -1); //Backward FFT, output goes to tspi

        for (int i = 0; i<ni; ++i)
        {
        for (int j = 0; j<nj; ++j)
        {        
            gridd.psi[i][j][n] = gridd.tpsi[i][j] * N;
        }
        }
    }

    if (ns == 1)
    {
        apply_1_state_potential_phase(gridd, ts); //Only using first state as this is a single state model
    }
    else
    {
        apply_potential_phase(gridd, ts); //Using both states for potential phase
    }

    for (int n = 0; n < ns; ++n)
    {
        for (int i = 0; i<ni; ++i)
        {
        for (int j = 0; j<nj; ++j)
         {
	    gridd.tmp[i][j] = gridd.psi[i][j][n];
         }
         }
    
        fft::run(f2d, gridd.tmp, 1); //Forward FFT, output goes to tspi
        apply_kinetic_phase(gridd.tpsi, gridd, ts2);
        fft::run(f2d, gridd.tpsi, -1); //Backward FFT, output goes to tspi

        for (int i = 0; i<ni; ++i)
        {
        for (int j = 0; j<nj; ++j)
        {        
            gridd.psi[i][j][n] = gridd.tpsi[i][j] * N;
        }
        }
    }

}
}