#pragma once
#include <complex>

constexpr int ni = 256, nj = 256, ns = 2; //for now fixed grid and states
constexpr char defaultInputFile[] = "INq";

using cplx = std::complex<double>;
constexpr cplx im = {0.0, 1.0}; //define imaginary unit 

void init(const std::string &fname); //read input, compute pi

