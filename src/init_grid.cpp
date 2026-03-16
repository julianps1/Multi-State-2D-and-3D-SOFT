#include <iostream>
#include "grid.hpp"
#include "wave.hpp"
#include "globals.hpp"
#include <fstream>
#include <cmath>

grid::GridData gridd;

namespace grid {
	double xmin, ymin, xmax, ymax, dx, dy, an0;
        int istate, iwa;
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
}
