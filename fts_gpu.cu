// ########################################################################
// Creates a new wtmd_simulation(...) object, passing input files to 
// its contructor and subsequently accessing its public methods to 
// equilibrate the system, gather statistics from the equilibrated system,
// and finally output the system's energy
// ########################################################################

#include <stdlib.h>
#include <iostream>
#include "wtmd_simulation.h"
using namespace std;


int main ()
{
    time_t t;

    // New wtmd_simulation object with input file name specified
    wtmd_simulation *wtmd_sim = new wtmd_simulation("input", "bias_in");

    // Time the equilibration period
    cout << "Starting Equilibrating..." << endl;
    t = time(NULL);
    wtmd_sim -> equilibrate();
    cout << "Equlibration time = " << time(NULL)-t << "secs" << endl << endl;

    // Time the statistics gathering period
    cout << "Starting Statistics..." << endl;
    t = time(NULL);
    wtmd_sim -> statistics();
    cout << "Statistics time = " << time(NULL)-t << "secs" << endl << endl;

    // Output the final energy of the system
    cout.precision(6);
    cout << wtmd_sim->getH() << endl;

    delete wtmd_sim;
    return 0;
}
