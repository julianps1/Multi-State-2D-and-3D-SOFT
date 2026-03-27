#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <cmath>
#include <fftw3.h>

#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include "FFT2D.hpp"


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

// print wf denisty to file
	void print_psi(int iwa, GridData &gridd)
          {
	   std::string filename = "output/wf_" + std::to_string(iwa) + ".dat";
	   std::ofstream outfile(filename);

		for (int i = 0; i < ni; ++i)
		{

		    for (int j = 0; j < nj; ++j)
    		    {
			double rho = 0.0;

			//for (int n = 0; n < ns; ++n)
		        //{
			// rho = rho + std::norm(gridd.psi[i][j][n]);
			//}
			//if (rho > 1e-6)
			// {
				 outfile << gridd.xg[i] << " "
        			   << gridd.yg[j] << " "
			           << std::norm(gridd.psi[i][j][0]) << " "
			           << std::norm(gridd.psi[i][j][1])
			           << "\n";
			 //} 
    		     }
                    // blank line after each i row
                       outfile << "\n";
                  }
	    }

    void corr(GridData &gridd, double t)
    {

    for (int n = 0; n < ns; ++n)
    {
        cplx sum = 0.0;
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                sum += std::conj(gridd.psi0[i][j][n]) * gridd.psi[i][j][n];
            }
        }
        gridd.c[n] = sum;
    }

        std::string filename = "output/correl.dat";
        std::ofstream outfile(filename, std::ios::app);
        outfile << t << " " << std::real(gridd.c[0]) 
                     << " " << std::imag(gridd.c[0]) 
                     << " " << std::real(gridd.c[1])
                     << " " << std::imag(gridd.c[1])
                     << '\n';

    }

}