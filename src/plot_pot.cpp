#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>
#include "grid.hpp"

namespace grid {

void write_potential(GridData &gridd)
{
    std::filesystem::create_directory("output");

    std::string filename = "output/v_" + std::to_string(iwa) + ".dat";
    std::ofstream outfile(filename);

    if (!outfile)
    {
        std::cerr << "Error opening " << filename << "\n";
        return;
    }

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
		outfile << gridd.xg[i] << " "
        		<< gridd.yg[j] << " "
		        << gridd.v[i][j][0][0] << " "
		        << gridd.v[i][j][1][1] << " "
		        << gridd.v[i][j][0][1] << " "
		        << std::real(gridd.vd[i][j][0]) << " "
		        << std::imag(gridd.vd[i][j][0]) << " "
		        << std::real(gridd.vd[i][j][1]) << " "
		        << std::imag(gridd.vd[i][j][1]) << " "
		        << gridd.vcos[i][j] << " "
		        << gridd.vsin[i][j] << "\n";
        }

        // blank line for gnuplot grid formatting
        outfile << "\n";
    }
}

}
