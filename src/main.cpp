#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    // read whatever input you have;
    init(defaultInputFile);
    wave::dtsub = wave::dt/wave::tsub;
    wave::dt2 = wave::dtsub * 0.5;

    std::cout << "wave parameters:\n"
              << "  alfa = " << wave::alfa
              << "  beta = " << wave::beta
              << "  q0   = " << wave::q0 
              << "  y0   = " << wave::y0 << '\n'
              << "  px0  = " << wave::px0
              << "  py0  = " << wave::py0 
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
    wave::init_mass(wave::h0,wave::h1,wave::h2);

    return 0;
}
