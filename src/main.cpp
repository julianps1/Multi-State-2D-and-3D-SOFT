#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include "FFT2D.hpp"
#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <fftw3.h>

int main()
{
    // read whatever input you have;
    init(defaultInputFile);
    wave::dtsub = wave::dt/wave::tsub;
    std::filesystem::create_directory("output");
    wave::dt2 = wave::dtsub * 0.5;
    
    std::cout << "wave parameters:\n"
              << "  alfa = " << wave::alfa
              << "  beta = " << wave::beta
              << "  q0   = " << wave::q0 
              << "  y0   = " << wave::y0 << '\n'
              << "  px0  = " << wave::px0
              << "  py0  = " << wave::py0
              << "  istate = " << grid::istate 
              << "  nwpackets  = " << wave::nwpackets << '\n';
    std::cout << "Propogation parameters:\n"
              << "  dt   = " << wave::dt 
              << "  tsub  = " << wave::tsub
              << "  tmax  = " << wave::tmax 
              << "  dtsub = " << wave::dtsub
              << "  dt2   = " << wave::dt2 << '\n';
    std::cout << "grid extents:\n"
              << "  xmin = " << grid::xmin
              << "  xmax = " << grid::xmax 
              << "  ymin = " << grid::ymin
              << "  ymax = " << grid::ymax << '\n';
    std::cout << "Initial masses:\n"
              << "  Na   = " << wave::h0 
              << "  H1   = " << wave::h1
              << "  H2   = " << wave::h2  << '\n';   
    //Initialize
    fft::FFT2D(ni,nj); //Initialize FFT
    
    wave::init_mass(wave::h0,wave::h1,wave::h2);
    grid::init_grid(grid::xmin,grid::ymin,grid::xmax,grid::ymax, gridd, grid::an0);
    grid::init_psi(grid::istate, gridd);
    grid::iwa = 0;
    grid::print_psi(grid::iwa, gridd);
    grid::Ham(f2d, gridd); //Initialize potential and Initial E
    grid::pot_diag(gridd);
    grid::write_potential(gridd);


    //Time loop
    double t = 0.0;
    for (int it = 0; it < wave::tmax; ++it)
    {
        for (int isub = 0; isub < wave::tsub; ++isub)
        {
            grid::split(wave::dtsub, wave::dt2, gridd, f2d);
            t += wave::dtsub;
        }
        grid::corr(gridd, t);

        if (it % wave::nwpackets == 0)
        {
            grid::iwa += 1;
            grid::print_psi(grid::iwa, gridd);
        }        
    }    
    return 0;
}
