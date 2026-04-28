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

void compute_moments(GridData &gridd, double t, int nmom)
{
    std::string filename = "output/moments.dat";
    std::ofstream outfile(filename, std::ios::app);
    if (t == 0) outfile << "#t,<x>,<y>,<x^2>,<xy>,<y^2>..." << '\n'; //File header

    if (ns != 1) {std::cout << "WARNING: MOMENTS ONLY IMPLEMENTED FOR GROUND STATE \n"; return;}

    double cmom = 0.0; 

    outfile << t << " ";
    for (int nm = 1; nm < nmom + 1; ++nm)
    {
        for (int nx = 0; nx < nm + 1; ++nx)
        {
            cmom = 0.0; 
            for (int i = 0; i < ni; ++i)
            {
            for (int j = 0; j < nj; ++j)
            {
                cmom += std::norm(gridd.psi[i][j][0]) * std::pow(gridd.xg[i], nm - nx) * std::pow(gridd.yg[j], nx); //Just ground state for now
            }
            }
            outfile << cmom << " ";
        }
    }
    outfile << "\n"; //New line after all moments appended

}

void compute_populations(GridData &gridd, double t)
{
    std::string filename = "output/populations.dat";
    std::ofstream outfile(filename, std::ios::app);
    if (t == 0) outfile << "#t,pop0,pop1\n"; //File header

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
    if (t == 0) outfile << "#t,pop0,pop1\n"; //File header

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
		        << std::real(gridd.v[i][j][0][0]) << " " 
                << std::imag(gridd.v[i][j][0][0]) << "\n";
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
		        << std::real(gridd.v[i][j][0][0]) << " "
		        << std::real(gridd.v[i][j][1][1]) << " "
		        << std::real(gridd.v[i][j][0][1]) << " "
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

    outfile << '\n'; //Blank line after all states appended

    }
    void RxnWeight(GridData &gridd, double rCUTx, double rCUTy)
    //Define the reaction/survival region with a weight function
    //rCUTx,y is the cutoff in x,y (can be different, if you want to include all values of a given dimension set the cutoff to the max)
    {
        
        std::string filename = "output/RxnWeight.dat";
        std::ofstream outfile(filename, std::ios::app);
       
        
        for (int i=0; i < ni; ++i)
        {
            for (int j=0; j < nj; ++j)
            {
				if (gridd.xg[i] > rCUTx or gridd.yg[j] > rCUTy)
				{
					gridd.wt[i][j] = 0;
				}
                else
                {
                    gridd.wt[i][j] = 1;
                }
            }
        }

        for (int i=0; i < ni; ++i)
        {
        for (int j=0; j < nj; ++j)
        {
            outfile << gridd.xg[i] << " " << gridd.yg[j] << " " << gridd.wt[i][j] <<'\n';
        }
        outfile << '\n'; //Blank line for gnuplot
        }
    }

    void RxnProb(GridData &gridd, double t, double ts)
    //Compute flux and reaction/survival probabilities (based on how weights are defined)
    {
    (void)ts;

    std::string filename1 = "output/RxnProb.dat";
    std::ofstream outfile1(filename1, std::ios::app);
    if (ts == 0) outfile1 << "#t,surv_prob,rxn_prob\n"; //File header

    //cplx flux = 0.0;
    outfile1 << t << " ";
    //outfile2 << t << " ";
    
    for (int n = 0; n < ns; ++n)
    {
        cplx sum = 0.0;


        for (int i = 0; i < ni; ++i)
        {
            for (int j = 0; j < nj; ++j)
            {
                sum += std::conj(gridd.psi[i][j][n]) * gridd.psi[i][j][n] * gridd.wt[i][j]; //Sum with weighting funciton (1 in reaction region 0 otherwise)
            }
        }

        //if(t > 0) flux += (sum - gridd.rc[n])/ts; //Compute flux as difference from last point over all states
        gridd.rc[n] = sum;

    }

        outfile1 << std::real(gridd.rc[0]) + std::real(gridd.rc[1]) << " " << 1 - std::real(gridd.rc[0]) - std::real(gridd.rc[1]) << " "; //Survival probability is 1 - reaction probability (for two states, just sum over both)
        //outfile2 << std::real(flux) << " ";

        outfile1 << '\n';
        //outfile2 << '\n';

    }

    double cent_diff(GridData &gridd, int ns, int i, int j, int xy)
    {
        
        double p1 = 0.0;
        double p2 = 0.0;
        double h = 1/dy/dx; //psi is stored as psi*sqrt(dx*dy), need to remove dy,dx from density
        
        switch (xy)
        {
        case 1: //Diff in x
            p1 = h * std::norm(gridd.psi[i+1][j][ns])/dx;
            p2 = h * std::norm(gridd.psi[i-1][j][ns])/dx;
            break;
       
        case 2: //Diff in y
            p1 = h * std::norm(gridd.psi[i][j+1][ns])/dy;
            p2 = h * std::norm(gridd.psi[i][j-1][ns])/dy;
            break;

        default:
            break;
        }

        double cdiff = 0.5 * (p1 - p2);

        return cdiff;

    }


    void compute_flux(GridData &gridd, const double rCAPx, const double rCAPy, double t)
    {

        int rxi = round((rCAPx - xmin)/dx);
        int ryi = round((rCAPy - ymin)/dy);

        //Flux is over a total area (in this case a line along a given direction) for now total length for rectangular borders this may need to be modified
        double Ay = ymax - ymin;
        double Ax = xmax - xmin;


        int fluxflagx = 1;
        int fluxflagy = 1;
        double fluxxt = 0.0;
        double fluxyt = 0.0;
        double fluxx = 0.0;
        double fluxy = 0.0;

        if (rCAPx == xmax) fluxflagx = 0;
        if (rCAPy == ymax) fluxflagy = 0;

        //std::cout << "PRINTING rxi, " << gridd.xg[rxi] << '\n';


        std::string filename = "output/Flux.dat";
        std::ofstream outfile(filename, std::ios::app);
        if (t == 0) outfile << "#flux in x for state0..ns, flux in y for state 0..ns, flux in x on all states, flux in y on all states, total flux \n";
        outfile << t << " ";

        for (int n = 0; n < ns; ++n)
        {

            fluxx = 0.0;

            for (int j = 0; j < nj; ++j)
            {
                fluxx += fluxflagx * cent_diff(gridd, n, rxi, j, 1); //Flux in x, summing all cent diffs in x over y 
            }

            outfile << fluxx/Ay << " ";
            fluxxt += fluxx; 

        }

        for (int n = 0; n < ns; ++n)
        {
            fluxy = 0.0;        
            
            for (int i = 0; i < ni; ++i)
            {
                fluxy += fluxflagy * cent_diff(gridd, n, i, ryi, 2); //Flux in y, summing all cent diffs in y over x 
            }

            outfile << fluxy/Ax << " ";
            fluxyt += fluxy;
        
        }

        outfile << fluxxt/Ax << " " << fluxyt/Ay << " " << fluxxt/Ax + fluxyt/Ay << '\n';

    }
}
