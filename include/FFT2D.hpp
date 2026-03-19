#pragma once
#include <fftw3.h>
#include "grid.hpp"
#include "globals.hpp"

namespace fft {
struct FFT2D
{

    int ni,nj;
    fftw_complex* data;
    fftw_plan forward;
    fftw_plan backward;

    FFT2D(int ni_, int nj_);
    ~FFT2D();

   
};
    void run(FFT2D& f2d, cplx wf[][nj], int call);
}

 extern fft::FFT2D f2d;
