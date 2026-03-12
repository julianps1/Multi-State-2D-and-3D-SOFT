#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    // read whatever input you have;
    dyn::read_input(defaultInputFile);

    std::cout << "wave parameters:\n"
              << "  alfa = " << wave::alfa
              << "  beta = " << wave::beta
              << "  q0   = " << wave::q0 
              << "  y0   = " << wave::y0 << '\n'
              << "  px0  = " << wave::px0
              << "  py0  = " << wave::py0 
              << "  dt   = " << wave::dt << '\n'
              << "  tsub  = " << wave::tsub
              << "  tmax  = " << wave::tmax 
              << "  nwpackets  = " << wave::nwpackets << '\n';
    std::cout << "grid extents:\n"
              << "  xmin = " << grid::xmin
              << "  xmax = " << grid::xmax 
              << "  ymin = " << grid::ymin
              << "  ymax = " << grid::ymax << '\n'
              << "  dx   = " << grid::dx
              << "  dy   = " << grid::dy  
              << "  an0  = " << grid::an0 << '\n';
    return 0;
}
