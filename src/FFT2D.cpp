#include "FFT2D.hpp"
#include <fftw3.h> 
#include <iostream>
#include <cmath>


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

namespace fft {
void run(FFT2D& fft, GridData& gridd, int ns, int call)
{
    // copy in
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
        {
            int idx = i*nj + j;
            fft.data[idx][0] = std::real(gridd.psi[i][j][ns]);
            fft.data[idx][1] = std::imag(gridd.psi[i][j][ns]);
        }
    if (call == 1) {fftw_execute(fft.forward);}
    else if (call == -1) {fftw_execute(fft.backward);}
    else {std::cout <<"ERROR: Unkown Call in FFT" << '\n';}

    // copy back
    for(int i=0;i<ni;i++)
        for(int j=0;j<nj;j++)
        {
            int idx = i*nj + j;
            gridd.psi[i][j][ns] = {fft.data[idx][0], fft.data[idx][1]};
        }
}
}
