// #######################################################################
// Creates a new wtmd_simulation(...) object, passing in input file to 
// its contructor and subsequently accessing its public methods to 
// equilibrate the system, gather statistics from the equilibrated system
// and finally output the system's energy
// #######################################################################

#include <string>
#include <iostream>
#include "wtmd_simulation.h"
#include <chrono>

using namespace std;

int main (int argc, char *argv[])
{
    // Get input file name from command-line argument
    if (argc != 3) {
        cout << "Please supply input file names as command-line arguments: <lfts_input_file> <bias_input_file>" << endl << endl;
        exit(1);
    }
    string inputFile(argv[1]);
    string biasFile(argv[2]);

    // New wtmd_simulation object with input file name specified
    // and 512 threads per block on the gpu
    wtmd_simulation *wtmd_sim = new wtmd_simulation(inputFile, biasFile, 512);
    
    // Time the equilibration period
    cout << "Starting Equilibration..." << endl;
    auto start = std::chrono::steady_clock::now();
    wtmd_sim -> equilibrate();
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Equlibration time = " << duration.count() << "secs" << endl << endl;

    // Time the statistics gathering period
    cout << "Starting Statistics..." << endl;
    start = std::chrono::steady_clock::now();
    wtmd_sim -> statistics();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Statistics time = " << duration.count() << "secs" << endl << endl;

    // Output the final energy of the system
    cout.precision(6);
    cout << wtmd_sim->getH() << endl;

    delete wtmd_sim;
    return 0;
}
