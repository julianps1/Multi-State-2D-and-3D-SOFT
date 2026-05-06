#pragma once
#include <string>
#include <array>
#include "globals.hpp"
#include "grid.hpp"

namespace grid 
{
    void print_psi(int iwa, GridData &gridd);
    void print_init_psi(GridData &gridd);
    void compute_vpsi(GridData &gridd);
    void write_potential(GridData &gridd);
    void corr(GridData &gridd, double t, int istate);
    void crosscorr(GridData &gridd, double t, int istate);     
    void RxnWeight(GridData &gridd, double rCAPx, double rCAPy);
    void RxnProb(GridData &gridd, double t, double ts);
    void compute_populations(GridData &gridd, double t);
    void compute_adipopulations(GridData &gridd, double t);
    void compute_moments(GridData &gridd, double t, int nmom);
    cplx overlapgg(const GridData &gridd, int n, const cplx wf2[][nj]);
    cplx overlapggs(const GridData &gridd, int n, const cplx wf2[][nj][ns]);
    void Ham(fft::FFT2D &f2d, GridData &gridd);
    void compute_flux(GridData &gridd, const double rCAPx, const double rCAPy, double t);
    double cent_diff(GridData &gridd, int ns, int i, int j, int xy);
}