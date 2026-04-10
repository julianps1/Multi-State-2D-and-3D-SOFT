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

void compute_populations(GridData &gridd, double t)
{
    std::string filename = "output/populations.dat";
    std::ofstream outfile(filename, std::ios::app);
    //outfile << "#t,pop0,pop1\n"; //File header

    outfile << t << " ";
    for (int n = 0; n < ns; ++n)
    {
        double pop = 0.0;
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                pop += std::norm(gridd.psi[i][j][n]);
            }
        }
        outfile << pop << " ";
    }

    outfile << "\n"; //New line after all states appended
}

void compute_adipopulations(GridData &gridd, double t)
{
    std::string filename = "output/adipopulations.dat";
    std::ofstream outfile(filename, std::ios::app);
    //outfile << "#t,pop0,pop1\n"; //File header

    //FOR TWO STATES ONLY


    outfile << t << " ";
        double pop0 = 0.0;
        double pop1 = 0.0;
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
               cplx z1 = (gridd.vcos[i][j] * gridd.psi[i][j][0] - gridd.vsin[i][j] * gridd.psi[i][j][1]);
               cplx z2 = (gridd.vsin[i][j] * gridd.psi[i][j][0] + gridd.vcos[i][j] * gridd.psi[i][j][1]);
            
               pop0 += std::norm(z1);
               pop1 += std::norm(z2);
            }
        }
        outfile << pop0 << " " << pop1 << " ";

    outfile << "\n"; //New line after all states appended
}

void write_potential(GridData &gridd)
{
    
    std::string filename = "output/v_" + std::to_string(iwa) + ".dat";
    std::ofstream outfile(filename);
    
    if (ns == 1) 
    {
        outfile << "#x,y,v00,Re(vd0),Im(vd0)\n"; //File header

    for (int i = 0; i < ni; ++i)
    {
        for (int j = 0; j < nj; ++j)
        {
		outfile << gridd.xg[i] << " "
        		<< gridd.yg[j] << " "
		        << gridd.v[i][j][0][0] << "\n";
        }

        // blank line for gnuplot grid formatting
        outfile << "\n";
    } 
    } 
    else 
    {
    outfile << "#x,y,v00,v11,v01,Re(vd0),Im(vd0),Re(vd1),Im(vd1),vcos,vsin\n"; //File header
    
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
}

//Print reference Psi's
	void print_init_psi(GridData &gridd)
          {
	   std::string filename = "output/psiref.dat";
	   std::ofstream outfile(filename, std::ios::app);
        outfile << "#x,y,|psi_0|^2,|psi_ref|^2\n"; //File header

        double gnorm = (1.0/dx/dy); //Normalization factor for grid density (since dx*dy included in normalization constants for psi, need to divide by it here to get correct density values)
		
        for (int i = 0; i < ni; ++i)
		{

		    for (int j = 0; j < nj; ++j)
    		    {

				 outfile << gridd.xg[i] << " "
        			   << gridd.yg[j] << " "
			           << std::norm(gridd.psi0[i][j][istate])*gnorm << " "
			           << std::norm(gridd.psiref[i][j][istate])*gnorm
			           << "\n";
    		     }
                    // blank line after each i row
                       outfile << "\n";
                  }
	    }

// print wf denisty to file
	void print_psi(int iwa, GridData &gridd)
        {
	   std::string filename = "output/snapshots/wf_" + std::to_string(iwa) + ".dat";
	   std::ofstream outfile(filename, std::ios::app);
        outfile << "#x,y,|psi0|^2,|psi1|^2\n"; //File header

        double gnorm = (1.0/dx/dy); //Normalization factor for grid density (since dx*dy included in normalization constants for psi, need to divide by it here to get correct density values)
		
		for (int i = 0; i < ni; ++i)
 		 {

		    for (int j = 0; j < nj; ++j)
    		    {
				 outfile << gridd.xg[i] << " "
           			     << gridd.yg[j] << " ";
                for (int n = 0; n < ns; ++n)      
                 {              
				 outfile << std::norm(gridd.psi[i][j][n])*gnorm << " ";
			     }
                 outfile << "\n"; //blank line after all states appended
                }
                    // blank line after each i row
                       outfile << "\n";
         }
	    }

    void corr(GridData &gridd, double t)
    {

    std::string filename = "output/correl.dat";
    std::ofstream outfile(filename, std::ios::app);
    //outfile << "#t,Re(cn),Im(cn),|cn|^2\n"; //File header

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

        outfile << t << " " << std::real(gridd.c[n]) 
                     << " " << std::imag(gridd.c[n]) 
                     << " " << std::norm(gridd.c[n]); 

    }

    outfile << '\n'; //Blank line after all states appended

    }

    void crosscorr(GridData &gridd, double t)
    {

    std::string filename = "output/crosscorrel.dat";
    std::ofstream outfile(filename, std::ios::app);
    //outfile << "#t,Re(cn),Im(cn),|cn|^2\n"; //File header
    
    for (int n = 0; n < ns; ++n)
    {
        cplx sum = 0.0;
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                sum += std::conj(gridd.psiref[i][j][n]) * gridd.psi[i][j][n];
            }
        }
        gridd.c[n] = sum;

        outfile << t << " " << std::real(gridd.c[n]) 
                     << " " << std::imag(gridd.c[n]) 
                     << " " << std::norm(gridd.c[n]); 
    }
    }
    void RxnWeight(GridData &gridd)
    //Define the reaction region with a weight function (for double well system)
    {
        for (int i=0; i < ni; ++i)
        {
            if (gridd.xg[i]>0)
            {
                gridd.wt[i] = 1;
            }
            else
            {
                gridd.wt[i] = 0; 
            }
        }
    }

    void RxnProb(GridData &gridd, double t)
    //For the double well system
    {

    for (int n = 0; n < ns; ++n)
    {
        cplx sum = 0.0;
        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                sum += std::conj(gridd.psi[i][j][n]) * gridd.psi[i][j][n] * gridd.wt[i]; //Sum with weighting funciton (1 in reaction region 0 otherwise)
            }
        }
        gridd.c[n] = sum;
    }

        std::string filename = "output/RxnProb.dat";
        std::ofstream outfile(filename, std::ios::app);
        outfile << t << " " << std::real(gridd.c[0]) 
                     << '\n';

    }
}