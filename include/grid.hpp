#pragma once
#include <string>
#include <array>
#include "globals.hpp"

namespace fft {
    struct FFT2D;
}

namespace grid {
    extern double xmin, ymin, xmax, ymax, dx, dy,an0; 
    extern int istate,iwa;

    struct GridData
	{
	double ak2[ni][nj];

	double xg[ni];
	double yg[nj];

	double v[ni][nj][ns][ns];
	double vcos[ni][nj];
	double vsin[ni][nj];

	cplx c[ns];
	cplx vd[ni][nj][ns];

	cplx psi[ni][nj][ns];
	cplx ksi[ni][nj][ns];
	cplx vpsi[ni][nj][ns];
	cplx psi0[ni][nj][ns];
        cplx tpsi[ni][nj];
	cplx tmp[ni][nj];
	};
    
     void init_grid(double xmin, double ymin, double xmax, double ymax, GridData &gridd, double& an0);
     void init_psi(int istate, GridData &gridd);
     void print_psi(int iwa, GridData &gridd);
     void pot_psi(GridData &gridd);
     void pot_diag(GridData &gridd);
     void write_potential(GridData &gridd);
     void kinetic_step(fft::FFT2D &f2d, GridData &gridd, int n);
     cplx overlapgg(const GridData &gridd, int n, const cplx wf2[][nj]);
     cplx overlapggs(const GridData &gridd, int n, const cplx wf2[][nj][ns]);
     void Ham(fft::FFT2D &f2d, GridData &gridd);
}

extern grid::GridData gridd;
