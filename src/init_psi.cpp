#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <cmath>

namespace grid {

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
}
