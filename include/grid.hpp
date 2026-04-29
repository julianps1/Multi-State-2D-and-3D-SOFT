#pragma once
#include <string>
#include <array>
#include "globals.hpp"

namespace fft {
    struct FFT2D;
}

namespace grid {
    extern double xmin, ymin, xmax, ymax, dx, dy, an0, rCAPx, rCAPy, kCAP; 
    extern int istate,iwa,pot_name,nmom;

    struct GridData
	{
    double ak2[ni][nj];

	double xg[ni];
	double yg[nj];

	cplx v[ni][nj][ns][ns];
	double vvec[ni][nj][ns][ns]; // diabatic row, adiabatic eigenvector column
    double wt[ni][nj];

	cplx vd[ni][nj][ns];

	cplx psi[ni][nj][ns];
	cplx ksi[ni][nj][ns];
	cplx vpsi[ni][nj][ns];
	cplx psi0[ni][nj][ns];
    cplx psiref[ni][nj][ns];
    cplx tpsi[ni][nj];
	cplx tmp[ni][nj];
    cplx expv[ni][nj][ns];
    cplx expk[ni][nj];

	};
 
    //Potential and wf on grid
     void init_grid(double xmin, double ymin, double xmax, double ymax, GridData &gridd, double& an0);
     void init_psi(int istate, GridData &gridd);
     void init_psiref(int istate, GridData &gridd);
	 void init_CAP(GridData &gridd, double rCAPx, double rCAPy, double kCAP);
     void pot_diag(GridData &gridd);
     
    //Potentials
     void pot_LandauSomb(GridData &gridd);
     void pot_Ferretti(GridData &gridd);
     void pot_doublewell(GridData &gridd);
     void load_surface_file(const std::string &path, GridData &gridd, int n1, int n2);
     void write_nafh_input_file(const std::string &path);
     void run_nafh_driver(const std::string &surface_dir);
     void pot_NaFH(GridData &gridd);
     void pot_h2p(GridData &gridd);
     void pick_pot(int pot_name, GridData &gridd);


    //SO
     void apply_kinetic_phase(cplx wf[][nj], const GridData &gridd);
     void apply_potential_phase(GridData &gridd);
     void apply_1_state_potential_phase(GridData &gridd);
     void init_exp(double ts, double ts2, GridData &gridd);
     void split(GridData &gridd, fft::FFT2D &f2d);
     void potential_step(GridData &gridd, int n);
     void kinetic_step(fft::FFT2D &f2d, GridData &gridd, int n);
     
}

extern grid::GridData gridd;
