
#pragma once
#include <string>
#include <fstream>
#include <iostream>

class wtmd_params {

    double ell_;            // Constant, l, used in definition of order parameter
    double kc_;             // Wave vector cutoff (constant) in order parameter
    double sigma_Psi_;      // Width of Gaussians added to bias potential
    double DT_;             // Controls the rate at which the amplitude of Guassians is reduced
    double Psi_min_;        // Lowest value of Psi for which the bias potential is collected
    int mPsi_;              // Total number of mesh points covering a range of psi
    double dPsi_;           // Distance between adjacent mesh points in psi
    int update_freq_;       // Number of Langevin steps after which the bias field is updated
    int read_bias_;         // Flag indicating whether an initial bias potential should be read from file

    public:
        wtmd_params(std::string biasFile) {
            read_Bias_Params(biasFile);
        }

        ~wtmd_params() {

        }

        // getters
        double ell() { return ell_; }
        double kc() { return kc_; }
        double sigma_Psi() { return sigma_Psi_; }
        double DT() { return DT_; }
        double Psi_min() { return Psi_min_; }
        int mPsi() { return mPsi_; }
        double dPsi() { return dPsi_; }
        int update_freq() { return update_freq_; }
        int read_bias() { return read_bias_; }

        // Print the bias parameters to std output
        void printInputParams() {
            std::cout << "ell = "           << ell_         << std::endl;
            std::cout << "kc = "            << kc_          << std::endl;
            std::cout << "sigma_Psi = "     << sigma_Psi_   << std::endl;
            std::cout << "DT = "            << DT_          << std::endl;
            std::cout << "Psi_min = "       << Psi_min_     << std::endl;
            std::cout << "mPsi = "          << mPsi_        << std::endl;
            std::cout << "dPsi = "          << dPsi_        << std::endl;
            std::cout << "update_freq = "   << update_freq_ << std::endl;
            std::cout << "read_bias = "     << read_bias_   << std::endl;
        }

        // Save the bias parameters and fields to file in the same format as the bias input file
        void saveBiasParams(std::string fileName, bool append=false) 
        {
            std::ofstream outstream;
            if (append) outstream.open(fileName,std::ios_base::app);
            else outstream.open(fileName);
            outstream << std::fixed;

            // Output the bias parameters
            outstream   << ell_ << " " << kc_ << " " << sigma_Psi_ << " " << DT_ << " "  
                        << Psi_min_ << " " << mPsi_ << " " << dPsi_ << " " << update_freq_ << " " << 1 
                        << std::endl;
            outstream.close();
        }

    
    private:
        // Read the parameters from the bias input file (first line)
        void read_Bias_Params(std::string fileName) {
            std::ifstream inFile;
            inFile.open(fileName);
            inFile >> ell_ >> kc_ >> sigma_Psi_ >> DT_ >> Psi_min_ >> mPsi_ >> dPsi_ >> update_freq_ >> read_bias_;
            inFile.close();
        }

};