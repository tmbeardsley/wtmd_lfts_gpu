#include <stdlib.h>
#include <iostream>
#include "wtmd_simulation.h"
using namespace std;


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
