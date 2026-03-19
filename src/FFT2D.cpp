#include "FFT2D.hpp"
#include "globals.hpp"
#include "grid.hpp"
#include "wave.hpp"
#include <fftw3.h> 
#include <iostream>
#include <cmath>

fft::FFT2D f2d(ni, nj);

namespace fft {
FFT2D::FFT2D(int ni_, int nj_) : ni(ni_), nj(nj_)
{
    data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ni * nj);

    forward = fftw_plan_dft_2d(ni, nj, data, data,
                              FFTW_FORWARD, FFTW_MEASURE);

    backward = fftw_plan_dft_2d(ni, nj, data, data,
                               FFTW_BACKWARD, FFTW_MEASURE);
}

FFT2D::~FFT2D()
{
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(data);
}


void run(FFT2D& f2d, cplx wf[][nj], int call)
{
    // copy in
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
        {
            int idx = i*nj + j;
            f2d.data[idx][0] = std::real(wf[i][j]);
            f2d.data[idx][1] = std::imag(wf[i][j]);
        }
    if (call == 1) {fftw_execute(f2d.forward);}
    else if (call == -1) {fftw_execute(f2d.backward);}
    else {std::cout <<"ERROR: Unkown Call in FFT" << '\n';}

    // copy back
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
        {
            int idx = i*nj + j;
            gridd.tpsi[i][j] = {f2d.data[idx][0], f2d.data[idx][1]};
        }
}
}
