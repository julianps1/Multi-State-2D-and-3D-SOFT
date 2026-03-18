#pragma once
#include <fftw3.h>

struct FFT2D
{
    int ni,nj; 
    ftw_complex* data;
    fftw_plan forward;
    fftw_plan backward;

    FFT2D(int ni_, int nj_);
    ~FFT2D();
};

namespace fft {
		void run(FFT2D& fft, GridData& gridd, int call);
}
