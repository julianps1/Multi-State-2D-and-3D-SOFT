#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include "linalg.hpp"
#include <fstream>
#include <cmath>
#include <filesystem>
#include <vector>

grid::GridData gridd;

namespace grid {
	double xmin, ymin, xmax, ymax, dx, dy, an0, rCAPx, rCAPy, kCAP;
	double coupling_scale = 1.0;
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
		outfile << "#x,y,max_offdiag_after_diag\n";

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
			std::vector<double> potential(ns * ns);
			for (int n1 = 0; n1 < ns; ++n1)
			{
				for (int n2 = 0; n2 < ns; ++n2)
				{
					if (n1 == n2)
					{
						potential[n1 * ns + n2] = std::real(gridd.v[i][j][n1][n2]);
					}
					else
					{
						//ASSUMING SYMMETRIC COUPLING
						const double v12 = std::real(gridd.v[i][j][n1][n2]);
						potential[n1 * ns + n2] = (v12);
					}
				}
			}

			const auto eig = linalg::diagonalize_symmetric(potential, ns);

			for (int n = 0; n < ns; ++n)
			{
				gridd.vd[i][j][n] = eig.eigenvalues[n];
			}

			for (int row = 0; row < ns; ++row)
			{
				for (int col = 0; col < ns; ++col)
				{
					gridd.vvec[i][j][row][col] = eig.eigenvector(row, col);
				}
			}

			double max_offdiag = 0.0;
			for (int a = 0; a < ns; ++a)
			{
				for (int b = 0; b < ns; ++b)
				{
					if (a == b)
					{
						continue;
					}

					double value = 0.0;
					for (int row = 0; row < ns; ++row)
					{
						for (int col = 0; col < ns; ++col)
						{
							value += gridd.vvec[i][j][row][a] *
								potential[row * ns + col] *
								gridd.vvec[i][j][col][b];
						}
					}

					const double offdiag = std::abs(value);
					if (offdiag > max_offdiag)
					{
						max_offdiag = offdiag;
					}
				}
			}

			outfile << gridd.xg[i] << " " << gridd.yg[j] << " " << max_offdiag << "\n";
			
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
						for (int n = 0; n < ns; ++n)
						{
							gridd.vd[i][j][n] += - im * kCAP * t1 * t1;
						}
					}
					double t2 = std::abs(gridd.yg[j]) - rCAPy;
					if (t2 > 0)
					{
						for (int n = 0; n < ns; ++n)
						{
							gridd.vd[i][j][n] += - im * kCAP * t2 * t2;
						}
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
