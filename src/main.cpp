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
    std::filesystem::create_directory("output/snapshots");
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
              << "  h0   = " << wave::h0 
              << "  h1   = " << wave::h1
              << "  h2   = " << wave::h2  << '\n';   
    //Initialize
     fft::FFT2D(ni,nj); //Initialize FFT
    //Initialize masses, if 3 masses given will use triatomicx model, otherwise will use generic 2D mass
    if (wave::h2 > 0.0)
        {
        wave::init_mass(wave::h0, wave::h1, wave::h2);
        }
    else
        {
        wave::init_mass_2D(wave::h0, wave::h1);
        }    


    //Initialize grid and wavefunction
     grid::init_grid(grid::xmin,grid::ymin,grid::xmax,grid::ymax, gridd, grid::an0);
     grid::init_psi(grid::istate, gridd);
     grid::init_psiref(grid::istate, gridd);
     grid::iwa = 0;
     grid::print_psi(grid::iwa, gridd);
     grid::print_init_psi(gridd); //Print the reference functions for correlations
     std::cout << "Initial energy...\n";
     grid::pick_pot(grid::pot_name, gridd);
     grid::Ham(f2d, gridd); //Initialize potential and Initial E
     if (ns == 2)
        {
        std::cout << "2-state system: Diagonalizing Potential... \n";  
        grid::pot_diag(gridd);
        }
     grid::write_potential(gridd);

    //DOUBE WELL INITIALIZE WEIGHTS FOR RXN PROB
    grid::RxnWeight(gridd);

    //Time loop
    double t = 0.0;

    //Print initial values
    grid::corr(gridd, t); //Autocorrelation
    grid::crosscorr(gridd, t); //Cross correlation with psi_ref
    grid::compute_populations(gridd, t); //Compute populations in each state
    grid::RxnProb(gridd, t); //Reaction probability for double well system
    if(ns == 2) { grid::compute_adipopulations(gridd, t); } //Compute adiabatic populations for 2 state system, may throw error if ns!= 2
    
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

        grid::corr(gridd, t); //Autocorrelation
        grid::crosscorr(gridd, t); //Cross correlation with psi_ref
        grid::compute_populations(gridd, t); //Compute populations in each state
        grid::RxnProb(gridd, t); //Reaction probability for double well system
        if(ns == 2) { grid::compute_adipopulations(gridd, t); } //Compute adiabatic populations for 2 state system, may throw error if ns!= 2

        if (it % wave::nwpackets == 0)
        {
            grid::iwa += 1;
            grid::print_psi(grid::iwa, gridd);
        }        
    }   
    
    //Recompute final energy
    std::cout << "\n";
    std::cout << "Final energy...\n";
    grid::Ham(f2d, gridd);

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
    std::cout << "Calculation completed successfully :) \n";
    return 0;
}
