#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>


namespace grid {
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
}
