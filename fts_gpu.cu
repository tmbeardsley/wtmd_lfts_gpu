
//------------------------------------------------------------
// GPU version of the FTS code for a diblock copolymer melt
// Note that lengths are expressed in units of R0=a*N^0.5
//------------------------------------------------------------

//#include<math.h>        // math subroutines
#include<stdlib.h>      // standard library
//#include<time.h>        // required to seed random number generator
//#include <cuda.h>       // required for GPUs
//#include"random.h"      // generates Gaussian random numbers
//#include "diblock.h"
#include <iostream>
// #include "anderson.h"
// #include "metadynamics.h"
#include "wtmd_simulation.h"

using namespace std;

//------------------------------------------------------------
// Global CPU variables:
//
// m[i] = number of mesh points in i'th dimension
// M = m[0]*m[1]*m[2] = total number of mesh points
// Mk = m[0]*m[1]*(m[2]/2+1) = total mesh points in k-space
// NA and NB = number of monomers in A and B blocks
// g[k] = FT of Boltzmann weight for bond potential divided by M 
// noise[r] = Gaussian noise for Langevin step
//------------------------------------------------------------
// Global GPU variables:
//
// qr[r] = q_i(r) interspliced with q^\dagger_i(r) 
// q1[i] = pointers to q_i(r) 
// q2[i] = pointers to q^\dagger_{N+1-i}(r) 
// qk[k] = Fourier transform (FT) of qr[r] times M/V
// w_gpu[4*M] = GPU copy of w[4*M] in main
// q1_N[r] = q_N(r)
// h[2*M] = stores Boltzmann weights hA[r] and hB[r] sequentially
// g_gpu[k] = GPU copy of g[k]
// qr_to_qk = plan for forward transform from qr[r] to qk[k] 
// qk_to_qr = plan for forward transform from qk[r] to qr[k] 
// wr_to_wk = plan for forward transform from w-[r] to wk[k] 
//------------------------------------------------------------




//============================================================
// Key variables in main:
//
// L[i] = size of the simulation box in the i'th dimension
// V = L[0]*L[1]*L[2] = volume of the simulation box
// w[4*M] = array containing N*w-, N*w+, phi-, phi+
// r = (x*m[1]+y)*m[2]+z = array position for (x,y,z)
// N = total number of monomers 
// chi_b = bare chi times N
// C = sqrt(Nbar) = dimensionless concentration
// sigma = standard deviation of the random noise times N
// dt = size of the Langevin step times N
// Hf = Hamiltonian in units of nkT
// lnQ = log of the single-chain partition function
// S[k] = structure function
// K[k] = modulus of the wavevector
// wt[k] = weighting of the wavevectors (equals 1 or 2)
// equil_its = number of equilibration steps
// sim_its = number of simulation steps
// sample_freq = frequency that observables are sampled
// wk[k] = fourier transform of w
// wk_gpu[k] = GPU copy of wk[k]

//------------------------------------------------------------
int main ()
{
    time_t t;

    wtmd_simulation *wtmd_sim = new wtmd_simulation("input", "bias_in");
    
    cout << "Starting Equilibrating..." << endl;
    t = time(NULL);
    wtmd_sim -> equilibrate();
    cout << "Equlibration time = " << time(NULL)-t << "secs" << endl << endl;

    cout << "Starting Statistics..." << endl;
    t = time(NULL);
    wtmd_sim -> statistics();
    cout << "Statistics time = " << time(NULL)-t << "secs" << endl << endl;

    cout.precision(6);
    cout << wtmd_sim->getH() << endl;

    delete wtmd_sim;
    return 0;
}
