#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include "FFT2D.hpp"
#include <iostream>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fftw3.h>

int main()
{
 
    double cpu_start = std::clock();
    auto wall_start = std::chrono::steady_clock::now();
 
    // read whatever input you have;
    init(defaultInputFile);
    wave::dtsub = wave::dt/wave::tsub;
    std::filesystem::remove_all("output");
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
    std::cout << "Potential:\n"
              << "  pot_name = " << grid::pot_name << '\n';
    std::cout << "Pi = " << pi << '\n';
    //Initialize
    fft::FFT2D(ni,nj); //Initialize FFT
    
    //wave::init_mass(wave::h0,wave::h1,wave::h2); //For tri atomic system like NaHF
    wave::init_mass_electron(); //For electronic dynamics
    grid::init_grid(grid::xmin,grid::ymin,grid::xmax,grid::ymax, gridd, grid::an0);
    grid::init_psi(grid::istate, gridd);
    grid::init_psiref(grid::istate, gridd);
    grid::iwa = 0;
    grid::print_psi(grid::iwa, gridd);
    grid::print_init_psi(gridd); //Print the reference functions for correlations
    std::cout << "Initial energy...\n";
    grid::Ham(f2d, gridd, grid::pot_name); //Initialize potential and Initial E
    //grid::pot_diag(gridd);
    grid::write_potential(gridd);

    //DOUBE WELL INITIALIZE WEIGHTS FOR RXN PROB
    grid::RxnWeight(gridd);

    //Time loop
    double t = 0.0;
    // Progress bar setup
    int total_steps = wave::tmax;
    int update_step = std::max(1, total_steps / 100); // update ~every 1%
    const int barWidth = 50;
    std::cout << "Simulation progress:" << std::endl;
    for (int it = 0; it < total_steps; ++it)
    {
        if ((it % update_step) == 0 || it == total_steps - 1)
        {
            int pct = (100 * it) / std::max(1, total_steps);
            int pos = (pct * barWidth) / 100;
            std::cout << "\r[";
            for (int i = 0; i < pos; ++i) std::cout << '=';
            for (int i = pos; i < barWidth; ++i) std::cout << ' ';
            std::cout << "] " << std::setw(3) << pct << "% (" << it << "/" << total_steps << ")";
            std::cout << std::flush;
        }
        for (int isub = 0; isub < wave::tsub; ++isub)
        {
            grid::split(wave::dtsub, wave::dt2, gridd, f2d);
            t += wave::dtsub;
        }
        grid::corr(gridd, t);
        grid::crosscorr(gridd, t);
        grid::RxnProb(gridd, t);

        if (it % wave::nwpackets == 0)
        {
            grid::iwa += 1;
            grid::print_psi(grid::iwa, gridd);
        }        
    }   
    
    //Recompute final energy
    std::cout << "Final energy...\n";
    grid::Ham(f2d, gridd, grid::pot_name);

    //Check Norm
    double norm = 0.0;
    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
            for (int n = 0; n < ns; ++n)
            {
                norm += std::norm(gridd.psi[i][j][n]);
            }
        }
    }
    std::cout << std::endl;
    std::cout << "Final Norm = " << norm << '\n';

    // Timing end and report
    auto wall_end = std::chrono::steady_clock::now();
    std::clock_t cpu_end = std::clock();
    double wall_secs = std::chrono::duration<double>(wall_end - wall_start).count();
    double cpu_secs = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Wall time: " << wall_secs << " s\n";
    std::cout << "CPU time:  " << cpu_secs << " s\n";
    return 0;
}
