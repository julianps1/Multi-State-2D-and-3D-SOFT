#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <cmath>
#include <filesystem>
#include <vector>

grid::GridData gridd;

namespace grid {
	double xmin, ymin, xmax, ymax, dx, dy, an0;
    int istate, iwa, pot_name;

	double stitch_theta_branch(double theta, double reference)
		{
			
			double tm = theta - reference;
			int j = std::round(- tm / (0.5 * pi));
			return theta + j * 0.5 * pi;
		}

	void init_psi(int istate, GridData &gridd)
       {	
		double cc = 0.0;

		for (int i = 0; i < ni; ++i)
		{
		    double x = xmin + dx * i;
		    double xx = x - wave::q0;
		    gridd.xg[i] = x;

		    for (int j = 0; j < nj; ++j)
    		{
		        double y = ymin + dy * j;
		        double yy = y - wave::y0;
			gridd.yg[j] = y;
		        
			gridd.psi[i][j][istate] =
		            std::exp(-wave::alfa*xx*xx - wave::beta*yy*yy
                	     + im*wave::px0*xx
	                     + im*wave::py0*yy)
        		     * std::sqrt(an0);

		        cc += std::norm(gridd.psi[i][j][istate]);

		        gridd.psi0[i][j][istate] = gridd.psi[i][j][istate];
    		}
		}
	    std::cout << "Norm =" << cc << '\n';
	   }
	void init_psiref(int istate, GridData &gridd)
	//Modify depending on model of interest
	   	{
			//Mirror image of psi for double well problem
			for (int i = 0; i < ni; ++i)
			{
				int idx = ni - 1 - i;
				for (int j = 0; j < nj; ++j)
				{
					int indxj = nj - 1 - j;
					gridd.psiref[i][j][istate] = gridd.psi0[idx][indxj][istate];
				}
			}
		}

	void init_grid(double xmin, double ymin, double xmax, double ymax, GridData &gridd, double& an0)
        {
        
		double alfa = wave::alfa;
		double beta = wave::beta;
		double amx = wave::amx;
		double amy = wave::amy;
		double akj,aki;
		double pi = ::pi;


        	dx = (xmax-xmin)/(ni-1);
	        dy = (ymax-ymin)/(nj-1);

		an0 = 2*std::sqrt(alfa*beta)/pi*dx*dy;

		double consI = 2.0*pi/dx/ni;
		double consJ = 2.0*pi/dy/nj;
		int ni2 = ni/2;
		int nj2 = nj/2;

		for(int j = 0; j < nj; j++)
		 {
		    akj = j * consJ;

		    if (j >= nj2)
		        akj = -(nj - j) * consJ;

		  for(int i = 0; i < ni; i++)
    		   {
		        aki = i * consI;

		        if (i >= ni2)
		            aki = -(ni - i) * consI;

		        gridd.ak2[i][j] = (aki*aki/amx + akj*akj/amy) * 0.5;
		    }
		  }

		std::cout <<"Grid Spacing \n"
			  <<" dx = " <<dx
			  <<" dy = " <<dy << '\n';
	   }

	   	void pot_diag(GridData &gridd)
        {
		
    	std::string filename = "output/theta.dat";
    	std::ofstream outfile(filename, std::ios::app);

		constexpr double eps = 1e-12;

		std::vector<std::vector<double>> theta_grid(ni, std::vector<double>(nj, 0.0));

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
    			double v11 = gridd.v[i][j][0][0];
		    	double v22 = gridd.v[i][j][1][1];
			    double v12 = gridd.v[i][j][0][1];
			
			//diagonalize
			//Matrix is [cos,-sin,sin,cos]

			//NOTE: When v12 is 0, atan2 returns 0 BUT Pi is another solution and the needed one to mirror the adibatic states, so need to switch to another branch hence subtracting pi/2 (Note theta from atan2 already halved)
			//NOTE: atan2 has  discontinutiy related to Pi to -Pi flip when sign of first argument changes, hence multiplying by -1 when v12 is positive 
			double dv = v11-v22;
			if (std::abs(dv) <  eps) {dv = eps;}

			double theta = 0.5 * std::atan(-2 * v12 / dv);
			theta_grid[i][j] = theta;
			
			if (1==2)
			{
			if (dv > 0 && v12 < 0)
			{			
			double theta = 0.5 * std::atan2(2 * v12 , dv) ;
			theta_grid[i][j] = theta - 0.5 * pi; 
			}
			else if (dv > 0 && v12 > 0)
			{
			double theta = 0.5 * std::atan2(2 * v12 , dv) ;
			theta_grid[i][j] = theta + 0.5 * pi;
			}
			else if (dv < 0 && v12 < 0)
			{
			double theta = 0.5 * std::atan2(-2 * v12 , dv);
			theta_grid[i][j] = theta + 0.5 * pi; 
			}
			else if (dv < 0 && v12 > 0)
			{
			double theta = 0.5 * std::atan2(-2 * v12 , dv);
			theta_grid[i][j] = theta - 0.5 * pi;
			}

			//double sqrtv12 = std::sqrt(dv * dv + 4.0 * v12 * v12);
			//gridd.vd[i][j][0] = 0.5 * (v11 + v22 - sqrtv12);
			//gridd.vd[i][j][1] = 0.5 * (v11 + v22 + sqrtv12);
			}
			}
		}
		


		//Compute vcos, vsin, and vd matrices
		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
			{
				
			double v11 = gridd.v[i][j][0][0];
			double v22 = gridd.v[i][j][1][1];

			//Because atan=0 even when states need to be mirrored (i.e. v00>v11) we can set that here with sin and cos
			if (v11 > v22)
			{
				gridd.vcos[i][j] = std::cos(theta_grid[i][j] + 0.5 * pi);
				gridd.vsin[i][j] = std::sin(theta_grid[i][j] + 0.5 * pi);
			}
			else
			{
				gridd.vcos[i][j] = std::cos(theta_grid[i][j]);
				gridd.vsin[i][j] = std::sin(theta_grid[i][j]);
			}

			double vd12 = 2 * gridd.vcos[i][j] * gridd.vcos[i][j] * gridd.v[i][j][0][1] + gridd.vsin[i][j] * gridd.vcos[i][j] * (gridd.v[i][j][0][0] - gridd.v[i][j][1][1]) - gridd.v[i][j][0][1];

			gridd.vd[i][j][0] = gridd.vcos[i][j]*gridd.vcos[i][j]*gridd.v[i][j][0][0] + gridd.vsin[i][j]*gridd.vsin[i][j]*gridd.v[i][j][1][1] - 2*gridd.vsin[i][j]*gridd.vcos[i][j]*gridd.v[i][j][0][1];
			gridd.vd[i][j][1] = gridd.vcos[i][j]*gridd.vcos[i][j]*gridd.v[i][j][1][1] + gridd.vsin[i][j]*gridd.vsin[i][j]*gridd.v[i][j][0][0] + 2*gridd.vsin[i][j]*gridd.vcos[i][j]*gridd.v[i][j][0][1];

			outfile << gridd.xg[i] << " " << gridd.yg[j] << " " << vd12 << "\n";
			}
			outfile << "\n"; //Blank line after each i row for gnuplot grid formatting
		}
	   
		//For debugging
		//int i = ni/2; int j = nj/2;
		//std::cout << "Printing transformed matrix:\n"
		//		  << " x = " << gridd.xg[i] << " y = " << gridd.yg[j] << '\n'
		//		  << " vd1' = " << gridd.vd[i][j][0] <<'\n'
		//		  << " vd2' = "<< gridd.vd[i][j][1] << '\n'
		//		  << " v12' = " << 2 * gridd.vcos[i][j] * gridd.vcos[i][j] * gridd.v[i][j][0][1] + gridd.vsin[i][j] * gridd.vcos[i][j] * (gridd.v[i][j][0][0] - gridd.v[i][j][1][1]) - gridd.v[i][j][0][1] << '\n'
		//		  << " va1' = " << gridd.vcos[i][j] * gridd.vcos[i][j] * gridd.v[i][j][0][0] + gridd.vsin[i][j] * gridd.vsin[i][j] * gridd.v[i][j][1][1] - 2 * gridd.vsin[i][j] * gridd.vcos[i][j] * gridd.v[i][j][0][1] << '\n'
		//		  << " va2' = " << gridd.vcos[i][j] * gridd.vcos[i][j] * gridd.v[i][j][1][1] + gridd.vsin[i][j] * gridd.vsin[i][j] * gridd.v[i][j][0][0] + 2 * gridd.vsin[i][j] * gridd.vcos[i][j] * gridd.v[i][j][0][1] << '\n';
		}
	void compute_vpsi(GridData &gridd)
	{
		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
			{
			for (int n1 = 0; n1 < ns; ++n1)
			{ cplx z(0.0, 0.0);
			for (int n2 = 0; n2 < ns; ++n2)
			 {
			 z = z+gridd.v[i][j][n1][n2]*gridd.psi[i][j][n2];
			 }
			 gridd.vpsi[i][j][n1] = z;
			}
			}
		}
	}
}
