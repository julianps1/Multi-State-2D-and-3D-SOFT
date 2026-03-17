#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <cmath>

namespace grid {

	void pot_psi(GridData &gridd)
          {
		constexpr double va = 0.17;
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
