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

void apply_kinetic_phase(cplx wf[][nj], const GridData &gridd)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {

        cplx phase = gridd.expk[i][j];
         wf[i][j] *= phase;

    }
    }
}
void apply_potential_phase(GridData &gridd)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        cplx adiabatic[ns] = {};
        cplx diabatic[ns] = {};

        for (int adi = 0; adi < ns; ++adi)
        {
            for (int dia = 0; dia < ns; ++dia)
            {
                adiabatic[adi] += gridd.vvec[i][j][dia][adi] * gridd.psi[i][j][dia];
            }
            adiabatic[adi] *= gridd.expv[i][j][adi];
        }

        for (int dia = 0; dia < ns; ++dia)
        {
            for (int adi = 0; adi < ns; ++adi)
            {
                diabatic[dia] += gridd.vvec[i][j][dia][adi] * adiabatic[adi];
            }
            gridd.psi[i][j][dia] = diabatic[dia];
        }

    }
    }
}
void apply_1_state_potential_phase(GridData &gridd)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        cplx phase = gridd.expv[i][j][0];
        gridd.psi[i][j][0] *= phase;

    }
    }
}
void init_exp(double ts, double ts2, GridData &gridd)
{
    for (int i = 0; i < ni; ++i)
    {
    for (int j = 0; j < nj; ++j)
    {
        if (ns == 1)
        {
        gridd.expv[i][j][0] = std::exp(im * gridd.v[i][j][0][0] * ts);
        }
        else
        {
        for (int n = 0; n < ns; ++n)
        {
        gridd.expv[i][j][n] = std::exp(im * gridd.vd[i][j][n] * ts);
        }
        }
        gridd.expk[i][j] = std::exp(im * ts2 * gridd.ak2[i][j]);
    }
    }
}
void split(GridData &gridd, fft::FFT2D &f2d) //Propagate the wavefunction (gridd.psi) for one time step using the split-operator method
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
        apply_kinetic_phase(gridd.tpsi, gridd);
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
        apply_1_state_potential_phase(gridd); //Only using first state as this is a single state model
    }
    else
    {
        apply_potential_phase(gridd); //Using both states for potential phase
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
        apply_kinetic_phase(gridd.tpsi, gridd);
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
