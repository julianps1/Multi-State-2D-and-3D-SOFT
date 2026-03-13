#pragma once
#include <string>
#include <array>
#include "globals.hpp"

namespace grid {
    extern double xmin, ymin, xmax, ymax, dx, dy; 

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
	cplx psi0[ni][nj][ns];
	};
    
     void init_grid(GridData &gridd);
}

extern grid::GridData gridd;
