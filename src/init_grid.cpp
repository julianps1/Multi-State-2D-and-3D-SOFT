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
	double xmin, ymin, xmax, ymax, dx, dy, an0, rCAPx, rCAPy, kCAP;
    int istate, iwa, pot_name, nmom;

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

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
    			double v11 = std::real(gridd.v[i][j][0][0]);
		    	double v22 = std::real(gridd.v[i][j][1][1]);
			    double v12 = std::real(gridd.v[i][j][0][1]);
			
			//diagonalize
			//Matrix is [cos,-sin,sin,cos]

			double dv = v11-v22;
			if (std::abs(dv) <  eps) {dv = eps;}

			double theta = 0.5 * std::atan(-2 * v12 / dv);
			

			if (v11 > v22) //Need to perform a 180 degree rotation when the diabatic states are not properly ordered
 			 {
				gridd.vcos[i][j] = std::cos(theta + 0.5 * pi);
				gridd.vsin[i][j] = std::sin(theta + 0.5 * pi);
			 }
			else
			 {
				gridd.vcos[i][j] = std::cos(theta);
				gridd.vsin[i][j] = std::sin(theta);
			 }
			
			double vd12 = 2 * gridd.vcos[i][j] * gridd.vcos[i][j] * v12 + gridd.vsin[i][j] * gridd.vcos[i][j] * (v11 - v22) - v12;
			gridd.vd[i][j][0] = gridd.vcos[i][j]*gridd.vcos[i][j]*v11 + gridd.vsin[i][j]*gridd.vsin[i][j]*v22 - 2*gridd.vsin[i][j]*gridd.vcos[i][j]*v12;
			gridd.vd[i][j][1] = gridd.vcos[i][j]*gridd.vcos[i][j]*v22 + gridd.vsin[i][j]*gridd.vsin[i][j]*v11 + 2*gridd.vsin[i][j]*gridd.vcos[i][j]*v12;

			outfile << gridd.xg[i] << " " << gridd.yg[j] << " " << vd12 << "\n";
			
			}
			outfile << "\n"; //Blank line after each i row for gnuplot grid formatting
		}
	}

	void init_CAP(GridData &gridd, double rCAPx, double rCAPy, double kCAP)
	{
		if (ns == 1)
			{
				for (int i = 0; i < ni; ++i)
				{
				for (int j = 0; j < nj; ++j)
				{
					double t1 = std::abs(gridd.xg[i]) - rCAPx;
					if (t1 > 0)
					{
						gridd.v[i][j][0][0] += - im * kCAP * t1 * t1;
					}
					double t2 = std::abs(gridd.yg[j]) - rCAPy;
					if (t2 > 0)
					{
						gridd.v[i][j][0][0] += - im * kCAP * t2 * t2;
					}
				}
				}
			}
		else
			{
				for (int i = 0; i < ni; ++i)
				{
				for (int j = 0; j < nj; ++j)
				{
					double t1 = std::abs(gridd.xg[i]) - rCAPx;
					if (t1 > 0)
					{
						gridd.vd[i][j][0] += - im * kCAP * t1 *t1;
						gridd.vd[i][j][1] += - im * kCAP * t1 *t1;
					}
					double t2 = std::abs(gridd.yg[j]) - rCAPy;
					if (t2 > 0)
					{
						gridd.vd[i][j][0] += - im * kCAP * t2 *t2;
						gridd.vd[i][j][1] += - im * kCAP * t2 *t2;
					}
				}
				}
			}
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
