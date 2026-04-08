#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <cmath>

namespace grid {
	void pot_doublewell(GridData &gridd)
	{
		constexpr double eta = 1/1.3544;
		constexpr double b =  1.0;
		constexpr double c =  1.0;
		constexpr double d = 1.0; //Linear Coupling

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
			{
			
			gridd.v[i][j][0][0] = 0.0625 * eta * (gridd.xg[i]*gridd.xg[i]) * (gridd.xg[i]*gridd.xg[i]) - 0.5 * b * (gridd.xg[i]*gridd.xg[i]) + 0.5 * c * gridd.yg[j]*gridd.yg[j] + 0.5 * d * gridd.xg[i]*gridd.yg[j]; //Ground PES Linear Coupling
		//	gridd.v[i][j][1][1] = 0.0625 * eta * (gridd.xg[i]*gridd.xg[i]) * (gridd.xg[i]*gridd.xg[i]) - 0.5 * b * (gridd.xg[i]*gridd.xg[i]) + 0.5 * c * gridd.yg[j]*gridd.yg[j] + 0.5 * d * gridd.xg[i]*gridd.yg[j]; //Excited PES Linear Coupling
		//	gridd.v[i][j][1][0] = 0; //Coupling
		//	gridd.v[i][j][0][1] = gridd.v[i][j][1][0]; //Symmetric Coupling

			
			
			}
		}
	   }

	void pot_Ferretti(GridData &gridd)
	{
		//Model parameters
		constexpr double Kx = 0.02;
		constexpr double Ky = 0.1;
		constexpr double X1 = 4.0;
		constexpr double X2 = 3.0;
		constexpr double X3 = 3.0;
		constexpr double alpha = 3.0;
		constexpr double beta = 1.5;
		constexpr double gamma = 0.04;
		constexpr double Delta = 0.01; 

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
			{
			
			gridd.v[i][j][0][0] = 0.5 * Kx * (gridd.xg[i]-X1)*(gridd.xg[i]-X1) + 0.5 * Ky * gridd.yg[j]*gridd.yg[j]; //Ground PES
			gridd.v[i][j][1][1] = 0.5 * Kx * (gridd.xg[i]-X2)*(gridd.xg[i]-X2) + 0.5 * Ky * gridd.yg[j]*gridd.yg[j] + Delta; //Excited PES
			gridd.v[i][j][1][0] = gamma * gridd.yg[j] * exp(-alpha * (gridd.xg[i]-X3)*(gridd.xg[i]-X3)  - beta * gridd.yg[j]*gridd.yg[j]); //Coupling
			gridd.v[i][j][0][1] = gridd.v[i][j][1][0]; //Symmetric Coupling



			}

		}
	   }

	void pot_psi(GridData &gridd)
          {
		//constexpr double va = 0.17;
		constexpr double re = 1.4;

		for (int i = 0; i < ni; ++i)
		{
		for (int j = 0; j < nj; ++j)
    		{
			double r = std::sqrt(gridd.xg[i]*gridd.xg[i]+gridd.yg[j]*gridd.yg[j]);
			
			for (int n1 = 0; n1 < ns; ++n1)
			{
			for (int n2 = 0; n2 < ns; ++n2)
			{
				if (n1 == n2)
				 { gridd.v[i][j][n1][n2] = 0.5 * (r-re) * (r-re); //single PES
				 }
				else
				 { gridd.v[i][j][n1][n2] = 0.0; // single uncoupled PES
				 }
			}
			}

    		}
		}
	   }
	   void pick_pot(int pot_name, GridData &gridd)
	   {
		switch (pot_name)
		{
		case 0:
			pot_psi(gridd);
			break;
		case 1:
			pot_Ferretti(gridd);
			break;
		case 2:
			pot_doublewell(gridd);
			break;
		default:
			std::cerr << "Invalid potential name: " << pot_name << "\n";
			break;
		}
    	}
}