#pragma once
#include <string>
#include <array>
#include "globals.hpp"

namespace fft {
    struct FFT2D;
}

namespace grid {
    extern double xmin, ymin, xmax, ymax, dx, dy, an0; 
    extern int istate,iwa,pot_name;

    struct GridData
	{
	double ak2[ni][nj];

	double xg[ni];
	double yg[nj];

	double v[ni][nj][ns][ns];
	double vcos[ni][nj];
	double vsin[ni][nj];
    double wt[ni];

	cplx c[ns];
    cplx cc[ns];
	cplx vd[ni][nj][ns];

	cplx psi[ni][nj][ns];
	cplx ksi[ni][nj][ns];
	cplx vpsi[ni][nj][ns];
	cplx psi0[ni][nj][ns];
    cplx psiref[ni][nj][ns];
    cplx tpsi[ni][nj];
	cplx tmp[ni][nj];
	};
    
     void init_grid(double xmin, double ymin, double xmax, double ymax, GridData &gridd, double& an0);
     void init_psi(int istate, GridData &gridd);
     void init_psiref(int istate, GridData &gridd);
     void print_psi(int iwa, GridData &gridd);
     void print_init_psi(GridData &gridd);
     void pot_psi(GridData &gridd);
     void pot_Ferretti(GridData &gridd);
     void pick_pot(int pot_name, GridData &gridd);
     void pot_diag(GridData &gridd);
     void write_potential(GridData &gridd);
     void apply_kinetic_phase(cplx wf[][nj], const GridData &gridd, double ts);
     void apply_potential_phase(GridData &gridd, double ts);
     void split(double ts, double ts2, GridData &gridd, fft::FFT2D &f2d);
     void potential_step(GridData &gridd, int n);
     void corr(GridData &gridd, double t);
     void crosscorr(GridData &gridd, double t);     
     void RxnWeight(GridData &gridd);
     void RxnProb(GridData &gridd, double t);
     void kinetic_step(fft::FFT2D &f2d, GridData &gridd, int n);
     cplx overlapgg(const GridData &gridd, int n, const cplx wf2[][nj]);
     cplx overlapggs(const GridData &gridd, int n, const cplx wf2[][nj][ns]);
     void Ham(fft::FFT2D &f2d, GridData &gridd, int pot_name);
}

extern grid::GridData gridd;
