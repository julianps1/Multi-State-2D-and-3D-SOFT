#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <iostream>

int main(int argc, char* argv[])
{
    // read whatever input you have; these functions
    // are defined in grid.cpp / wave.cpp
    wave::read(defaultInputFile);
    grid::read(defaultInputFile);

    std::cout << "wave parameters:\n"
              << "  alfa = " << wave::alfa
              << "  beta = " << wave::beta
              << "  q0   = " << wave::q0 << '\n';
    // …print the rest…

    std::cout << "grid extents:\n"
              << "  xmin = " << grid::xmin
              << "  xmax = " << grid::xmax << '\n';
    // …etc.

    return 0;
}
