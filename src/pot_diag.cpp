#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <string>
#include <filesystem>
#include <cmath>

namespace grid {

	void pot_diag(GridData &gridd)
          {

		std::filesystem::create_directory("output");

		    std::string filename = "output/theta.dat";
		    std::ofstream outfile(filename);

		    if (!outfile)
		    {
		        std::cerr << "Error opening " << filename << "\n";
		        return;
		    }

		constexpr double eps = 1e-12;

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
    			double v11 = gridd.v[i][j][0][0];
			double v22 = gridd.v[i][j][1][1];
			double v12 = gridd.v[i][j][0][1];
			
			//diagonalize
			
			double dv = v11-v22;
			if (std::abs(dv) < eps) { dv = eps; }
			//double tgh = 2 * v12/dv; 
                        //double theta = 0.5 * std::atan(tgh);
                        double theta = 0.5 * std::atan2(2.0 * v12, dv);

			//check theta is smooth
			outfile << gridd.xg[i] << " "
                                << gridd.yg[j] << " "
                                << theta << "\n";

			gridd.vcos[i][j] = std::cos(theta);
			gridd.vsin[i][j] = std::sin(theta);

			double sqrtv12 = std::sqrt(dv * dv + 4.0 * v12 * v12);
			gridd.vd[i][j][0] = 0.5 * (v11 + v22 - sqrtv12);
			gridd.vd[i][j][1] = 0.5 * (v11 + v22 + sqrtv12);
		}
		// blank line for gnuplot grid formatting
	        outfile << "\n";
		}


		//Add imaginary potential
		
		constexpr double rim = 0.0;
		constexpr double dim = 0.0;

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
			double t1 = std::abs(gridd.xg[i]) - rim;
			if (t1 > 0) 
			 {
			  gridd.vd[i][j][0] = gridd.vd[i][j][0] - im * dim * t1 * t1;
			  gridd.vd[i][j][1] = gridd.vd[i][j][1] - im * dim * t1 * t1;
			 }
			double t2 = std::abs(gridd.yg[j]) - rim;
			if (t2 > 0) double theta = 0.5 * std::atan2(2.0 * v12, dv);
			 {
			  gridd.vd[i][j][0] = gridd.vd[i][j][0] - im * dim * t2 * t2;
			  gridd.vd[i][j][1] = gridd.vd[i][j][1] - im * dim * t2 * t2;
			 }
			
	    	}
		}
		
	   }
}
