#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include <array>
#include <cuba.h>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iomanip>



class PhysicsCalculator
{
private:
    // Constants
    static constexpr double alphaem = 1.0 / 137.036;
    static constexpr double thetac = 3.0e-2; // rad [value for FCC-ee]; search for adequate values for other colliders
    static constexpr double me = 5.11e-4;    // in GeV
    static constexpr double SqrtS = 50.0;    // GeV [value for FCC-ee: 92.]
    static constexpr double S = SqrtS * SqrtS;
    static constexpr double PI = 3.14159265358979323846;

    // Variable extremes
    static constexpr double Q20 = 3.0; // GeV
    static constexpr double Q2M = S - 4 * Q20;
    static constexpr double csimax = 1.0;
    static constexpr double ymax = 1.0;

    // Cost factor
    static constexpr double Cost = 3 * alphaem * alphaem * alphaem / (4 * PI);

    // Member variables for current calculation
    double xB, y, csi, Q2;
    double zeta, zetamin, zetamax;
    double shat;
    double Q2min, Q2max;
    double fgl, DLfgl;


public:
    // Constructor
    PhysicsCalculator() : xB(0), y(0), csi(0), Q2(0), zeta(0), zetamin(0), zetamax(0), shat(0), Q2min(0), Q2max(0), fgl(0), DLfgl(0) {}

    // Set variables 
    void setVariables(double xB_val, double y_val, double csi_val)
    {
        xB = xB_val;
        y = y_val;
        csi = csi_val;
        Q2 = xB * y * S;

        shat = (csi - xB) * y * S;

        // Calculate zeta bounds
        double sqrt_term = sqrt(1 - 4 * Q20 / shat);
        zetamin = 0.5 * (1 - sqrt_term);
        zetamax = 0.5 * (1 + sqrt_term);

        // WW function parameters
        Q2min = me * me * csi * csi / (1 - csi);
        Q2max = (S / 4) * thetac * thetac * (1 - csi) + Q2min;

        // Calculate WW functions
        double log_ratio = log(Q2max / Q2min);
        double me2_term = 2 * me * me * csi * (1.0 / Q2max - 1.0 / Q2min);

        fgl = (alphaem / (2 * PI)) * ((1 + (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);
        DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + 2 * me * me * csi * csi * (1.0 / Q2max - 1.0 / Q2min));
    }
    void setVariablesWithQ2(double xB_val, double Q2_val, double csi_val)
    {
        xB = xB_val;
        Q2 = Q2_val;
        csi = csi_val;
        
        y = Q2 / (xB * S);
        
        shat = (csi - xB) * y * S;
        
        double sqrt_term = sqrt(1 - 4 * Q20 / shat);
        zetamin = 0.5 * (1 - sqrt_term);
        zetamax = 0.5 * (1 + sqrt_term);
        
        Q2min = me * me * csi * csi / (1 - csi);
        Q2max = (S / 4) * thetac * thetac * (1 - csi) + Q2min;
        
        double log_ratio = log(Q2max / Q2min);
        double me2_term = 2 * me * me * csi * (1.0 / Q2max - 1.0 / Q2min);
        
        fgl = (alphaem / (2 * PI)) * ((1 + (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);
        DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + 2 * me * me * csi * csi * (1.0 / Q2max - 1.0 / Q2min));
    }

    // Utility functions for variable bounds
    double xBmin(double Q2_val) const
    {
        return Q2_val / S;
    }

    double xBmax(double Q2_val) const
    {
        return Q2_val / (Q2_val + 4 * Q20);
    }

    double csimin(double Q2_val) const
    {
        return xB * (1 + 4 * Q20 / Q2_val);
    }

    double ymin(double Q2_val) const
    {
        return xB * (1 + 4 * Q20 / Q2_val);
    }

    // K factors for different differentials
    double Kxy() const
    {
        return 1.0 / (xB * y * y * csi * csi * csi * S);
    }

    double KQx() const
    {
        return 1.0 / (csi * csi * csi * Q2 * Q2);
    }

    double KQy() const
    {
        return 1.0 / (Q2 * y * y * csi * csi * csi);
    }

    // A coefficients
    double AU() const
    {
        double term1 = (1 + (1 - y) * (1 - y)) * (xB * xB + (csi - xB) * (csi - xB));
        double term2 = (1 - 2 * zeta * (1 - zeta)) / (zeta * (1 - zeta));
        double term3 = 16 * (1 - y) * xB * (csi - xB);
        return 2 * (term1 * term2 + term3) * fgl;
    }

    double AUcfq() const
    {
        double term1 = -8 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB);
        double term2 = sqrt(xB * (csi - xB)) * (1 - 2 * zeta) / sqrt(zeta * (1 - zeta));
        return term1 * term2 * fgl;
    }

    double AUc2fq() const
    {
        return 16 * (1 - y) * xB * (csi - xB) * fgl;
    }

    double AL() const
    {
        double term1 = -2 * y * (2 - y) * csi * (csi - 2 * xB);
        double term2 = (1 - 2 * zeta * (1 - zeta)) / (zeta * (1 - zeta));
        return term1 * term2 * DLfgl;
    }

    double ALcfq() const
    {
        double term1 = 8 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = (1 - 2 * zeta) / sqrt(zeta * (1 - zeta));
        return term1 * term2 * DLfgl;
    }

    // B coefficients
    double BUcf12() const
    {
        double term1 = (1 + (1 - y) * (1 - y)) * (xB * xB + (csi - xB) * (csi - xB));
        double term2 = 8 * (1 - y) * xB * (csi - xB);
        return (term1 - term2) * fgl;
    }

    double BUcfqmf12() const
    {
        double term1 = -2 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB);
        double term2 = sqrt(xB * (csi - xB)) * sqrt(zeta / (1 - zeta));
        return term1 * term2 * fgl;
    }

    double BUcfqpf12() const
    {
        double term1 = 2 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB);
        double term2 = sqrt(xB * (csi - xB)) * sqrt((1 - zeta) / zeta);
        return term1 * term2 * fgl;
    }

    double BUc2fqmf12() const
    {
        return 2 * (1 - y) * xB * (csi - xB) * zeta / (1 - zeta) * fgl;
    }

    double BUc2fqpf12() const
    {
        return 2 * (1 - y) * xB * (csi - xB) * (1 - zeta) / zeta * fgl;
    }

    double BLcf12() const
    {
        return -y * (2 - y) * csi * (csi - 2 * xB) * DLfgl;
    }

    double BLcfqmf12() const
    {
        double term1 = 2 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = sqrt(zeta / (1 - zeta));
        return term1 * term2 * DLfgl;
    }

    double BLcfqpf12() const
    {
        double term1 = -2 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = sqrt((1 - zeta) / zeta);
        return term1 * term2 * DLfgl;
    }

    // Integrated coefficients over zeta range
    double AUzi() const
    {
        double term1 = (1 + (1 - y) * (1 - y)) * (xB * xB + (csi - xB) * (csi - xB));
        double term2 = -2 * (zetamax - zetamin) + 2 * log(zetamax / zetamin);
        double term3 = 16 * (1 - y) * xB * (csi - xB) * (zetamax - zetamin);
        return 2 * (term1 * term2 + term3) * fgl;
    }

    double AUcfqzi() const
    {
        double term1 = -8 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB);
        double term2 = sqrt(xB * (csi - xB)) * (2 * sqrt(zetamax * zetamin) - 1);
        return term1 * term2 * fgl; // half-range
    }

    double AUc2fqzi() const
    {
        return 16 * (1 - y) * xB * (csi - xB) * (zetamax - zetamin) * fgl;
    }

    double ALzi() const
    {
        double term1 = -2 * y * (2 - y) * csi * (csi - 2 * xB);
        double term2 = -2 * (zetamax - zetamin) + 2 * log(zetamax / zetamin);
        return term1 * term2 * DLfgl;
    }

    double ALcfqzi() const
    {
        double term1 = 8 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = 2 * sqrt(zetamax * zetamin) - 1;
        return term1 * term2 * DLfgl; // half-range
    }

    double BUcf12zi() const
    {
        double term1 = (1 + (1 - y) * (1 - y)) * (xB * xB + (csi - xB) * (csi - xB));
        double term2 = 8 * (1 - y) * xB * (csi - xB);
        return (term1 - term2) * (zetamax - zetamin) * fgl;
    }

    double BUcfqmf12zi() const
    {
        double term1 = -2 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB) * sqrt(xB * (csi - xB));
        double term2 = -atan(sqrt(zetamin / zetamax)) + atan(sqrt(zetamax / zetamin));
        return term1 * term2 * fgl;
    }

    double BUcfqpf12zi() const
    {
        double term1 = 2 * (2 - y) * sqrt(1 - y) * (csi - 2 * xB) * sqrt(xB * (csi - xB));
        double term2 = -atan(sqrt(zetamin / zetamax)) + atan(sqrt(zetamax / zetamin));
        return term1 * term2 * fgl;
    }

    double BUc2fqmf12zi() const
    {
        double term1 = 2 * (1 - y) * xB * (csi - xB);
        double term2 = -(zetamax - zetamin) + log(zetamax / zetamin);
        return term1 * term2 * fgl;
    }

    double BUc2fqpf12zi() const 
    {
        double term1 = 2 * (1 - y) * xB * (csi - xB);
        double term2 = -(zetamax - zetamin) + log(zetamax / zetamin);
        return term1 * term2 * fgl;
    }

    double BLcf12zi() const
    {
        return -y * (2 - y) * csi * (csi - 2 * xB) * (zetamax - zetamin) * DLfgl;
    }

    double BLcfqmf12zi() const
    {
        double term1 = 2 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = -atan(sqrt(zetamin / zetamax)) + atan(sqrt(zetamax / zetamin));
        return term1 * term2 * DLfgl;
    }

    double BLcfqpf12zi() const
    {
        double term1 = -2 * y * sqrt(1 - y) * csi * sqrt(xB * (csi - xB));
        double term2 = -atan(sqrt(zetamin / zetamax)) + atan(sqrt(zetamax / zetamin));
        return term1 * term2 * DLfgl;
    }

    // Getter functions for current variable values
    double getXB() const { return xB; }
    double getY() const { return y; }
    double getCsi() const { return csi; }
    double getQ2() const { return Q2; }
    double getZetamin() const { return zetamin; }
    double getZetamax() const { return zetamax; }
    double getFgl() const { return fgl; }
    double getDLfgl() const { return DLfgl; }

    // Static getters for constants
    static double getAlphaem() { return alphaem; }
    static double getThetac() { return thetac; }
    static double getMe() { return me; }
    static double getSqrtS() { return SqrtS; }
    static double getS() { return S; }
    static double getQ20() { return Q20; }
    static double getQ2M() { return Q2M; }
    static double getCost() { return Cost; }
};

static double fixed_xB ; // Global variable to hold fixed xB value
static double fixed_y ; // Global variable to hold fixed y value
static double fixed_Q2 ; // Global variable to hold fixed Q2 value

//Cut function to check if the input values are within the allowed ranges
bool cut(const PhysicsCalculator& pscut)
{
    double S = PhysicsCalculator::getS();
    double Q20 = PhysicsCalculator::getQ20();
    double Q2M = PhysicsCalculator::getQ2M();

    double xB = pscut.getXB();
    double y = pscut.getY();
    double csi = pscut.getCsi();
    double Q2 = pscut.getQ2();

    if (Q2 < Q20 || Q2 > Q2M)
        return false;
    if (xB < pscut.xBmin(Q2) || xB > pscut.xBmax(Q2))
        return false;
    if (csi < pscut.csimin(Q2) || csi > 1.0)
        return false;
    if (y < pscut.ymin(Q2) || y > 1.0)
        return false;

    return true;
}

//example functions to calculate azimuthal moments NOT used in the main program
void Values(const std::vector<double> &x, double results[])
{   
    PhysicsCalculator pc;

    pc.setVariables(x[0], x[1], x[2]);
    if(cut(pc) == false)
    {
        std::cerr << "Cut condition failed for input values: "
                  << "xB = " << x[0] << ", y = " << x[1] << ", csi = " << x[2] << std::endl;
    }
    else
    {
        std::cout << "Cut condition passed for input values: "
                  << "xB = " << x[0] << ", y = " << x[1] << ", csi = " << x[2] << std::endl;
    }
    if(cut(pc)==true){
    results[0]  = pc.AUzi();
    results[1]  = pc.AUcfqzi();
    results[2]  = pc.AUc2fqzi();
    results[3]  = pc.ALzi();
    results[4]  = pc.ALcfqzi();
    results[5]  = pc.BUcf12zi();
    results[6]  = pc.BUcfqmf12zi();
    results[7]  = pc.BUcfqpf12zi();
    results[8]  = pc.BUc2fqmf12zi();
    results[9]  = pc.BUc2fqpf12zi();
    results[10] = pc.BLcf12zi();
    results[11] = pc.BLcfqmf12zi();
    results[12] = pc.BLcfqpf12zi();
    }
    else
    {
        std::fill(results, results + 13, 0.0); // Fill results with zeros if cut fails
    }

    // Print the results for debugging
    std::cout << "Results of Azimuthal Moments Calculation:" << std::endl
                << "xB: " << pc.getXB() << std::endl
                << "y: " << pc.getY() << std::endl
                << "csi: " << pc.getCsi() << std::endl
                << "AUzi: " << results[0] << std::endl
                << "AUcfqzi: " << results[1] << std::endl
                << "AUc2fqzi: " << results[2] << std::endl
                << "ALzi: " << results[3] << std::endl
                << "ALcfqzi: " << results[4] << std::endl
                << "BUcf12zi: " << results[5] << std::endl
                << "BUcfqmf12zi: " << results[6] << std::endl
                << "BUcfqpf12zi: " << results[7] << std::endl
                << "BUc2fqmf12zi: " << results[8] << std::endl
                << "BUc2fqpf12zi: " << results[9] << std::endl
                << "BLcf12zi: " << results[10] << std::endl
                << "BLcfqmf12zi: " << results[11] << std::endl
                << "BLcfqpf12zi: " << results[12] << std::endl;
}

void ValuesfixedxB(const std::vector<double> &x, double results[])
{
    // x[0] = y, x[1] = csi; fixed_xB is provided externally
    PhysicsCalculator pc;

    pc.setVariables(fixed_xB, x[0], x[1]);
    if (cut(pc))
    {
        results[0]  = pc.AUzi() *pc.Kxy();
        results[1]  = pc.AUcfqzi()*pc.Kxy();
        results[2]  = pc.AUc2fqzi()*pc.Kxy();
        results[3]  = pc.ALzi()*pc.Kxy();
        results[4]  = pc.ALcfqzi()*pc.Kxy();
        results[5]  = pc.BUcf12zi()*pc.Kxy();
        results[6]  = pc.BUcfqmf12zi()*pc.Kxy();
        results[7]  = pc.BUcfqpf12zi()*pc.Kxy();
        results[8]  = pc.BUc2fqmf12zi()*pc.Kxy();
        results[9]  = pc.BUc2fqpf12zi()*pc.Kxy();
        results[10] = pc.BLcf12zi()*pc.Kxy();
        results[11] = pc.BLcfqmf12zi()*pc.Kxy();
        results[12] = pc.BLcfqpf12zi()*pc.Kxy();
    }
    else
    {
        std::fill(results, results + 13, 0.0);
    }
}
void Valuesfixedy(const std::vector<double> &x, double results[])
{
    // x[0] = xB, x[1] = csi; fixed_y is provided externally
    PhysicsCalculator pc;

    pc.setVariables(x[0], fixed_y, x[1]);
    if (cut(pc))
    {
        results[0]  = pc.AUzi()*pc.Kxy();
        results[1]  = pc.AUcfqzi()*pc.Kxy(); 
        results[2]  = pc.AUc2fqzi()*pc.Kxy();
        results[3]  = pc.ALzi()*pc.Kxy();
        results[4]  = pc.ALcfqzi()*pc.Kxy();
        results[5]  = pc.BUcf12zi()*pc.Kxy();
        results[6]  = pc.BUcfqmf12zi()*pc.Kxy();
        results[7]  = pc.BUcfqpf12zi()*pc.Kxy();
        results[8]  = pc.BUc2fqmf12zi()*pc.Kxy();
        results[9]  = pc.BUc2fqpf12zi()*pc.Kxy();
        results[10] = pc.BLcf12zi()*pc.Kxy();
        results[11] = pc.BLcfqmf12zi()*pc.Kxy();
        results[12] = pc.BLcfqpf12zi()*pc.Kxy();
    }
    else
    {
        std::fill(results, results + 13, 0.0);
    }
}

void ValuesfixedQ2(const std::vector<double> &x, double results[])
{
    // x[0] = xB, x[1] = csi; fixed Q2 is provided externally
    PhysicsCalculator pc;
    pc.setVariablesWithQ2(x[0], fixed_Q2, x[1]);
    if (cut(pc))
    {
        results[0]  = pc.AUzi()*pc.KQx();
        results[1]  = pc.AUcfqzi()*pc.KQx();
        results[2]  = pc.AUc2fqzi()*pc.KQx();
        results[3]  = pc.ALzi()*pc.KQx();
        results[4]  = pc.ALcfqzi()*pc.KQx();
        results[5]  = pc.BUcf12zi()*pc.KQx();
        results[6]  = pc.BUcfqmf12zi()*pc.KQx();
        results[7]  = pc.BUcfqpf12zi()*pc.KQx();
        results[8]  = pc.BUc2fqmf12zi()*pc.KQx();
        results[9]  = pc.BUc2fqpf12zi()*pc.KQx();
        results[10] = pc.BLcf12zi()*pc.KQx();
        results[11] = pc.BLcfqmf12zi()*pc.KQx();
        results[12] = pc.BLcfqpf12zi()*pc.KQx();
    }
    else
    {
        std::fill(results, results + 13, 0.0);
    }
}

int integrand_fixedxB(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata)
{
    std::vector<double> xv;
    for (int i = 0; i < *ndim; i++)
    {
        xv.push_back(x[i]);
    }

    ValuesfixedxB(xv, f);

    return 0;
}

int integrand_fixedy(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata)
{
    std::vector<double> xv;
    for (int i = 0; i < *ndim; i++)
    {
        xv.push_back(x[i]);
    }

    Valuesfixedy(xv, f);

    return 0;
}
int integrand_fixedQ2(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata)
{
    std::vector<double> xv;
    for (int i = 0; i < *ndim; i++)
    {
        xv.push_back(x[i]);
    }

    ValuesfixedQ2(xv, f);

    return 0;
}


struct IntegrationResults {
    int maxeval;
    int nstart;
    int neval;
    int fail;
    double integral[13];
    double error[13];
    double prob[13];
    double ratio[13];
    double ratioerror[13];
    double elapsed_seconds;
};

int main() {
    // Fixed integration settings
    const int maxeval = static_cast<int>(1e6);
    const int nstart  = static_cast<int>(1e5);

    // Fixed values for xB, y, and Q2
    std::vector<double> xB_values;
    for(double i = 0; i <=0.2; i += 0.001) xB_values.push_back(i);
    std::vector<double> y_values;
    for(double i = 0; i <=1; i += 0.01) y_values.push_back(i);
    std::vector<double> Q2_values;
    for(double i = 0; i <= 1; i += 0.001) Q2_values.push_back(i * (PhysicsCalculator::getQ2M() - PhysicsCalculator::getQ20()) + PhysicsCalculator::getQ20()); // Assuming Q2 ranges from 0 to Q2M
    
    //set parameters for the integration
    int ndim = 2, ncomp = 13, nvec = 1;
    double epsrel = 1e-4, epsabs = 1e-9;
    int flags = 2, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;
    
    std::ofstream outFile("Prova0.1.txt");  
    outFile << std::setw(10) << "xB<0.2" 
            << std::setw(15) << "<dsig|+1,0>" << std::setw(15) << "err[1]"
            << std::setw(15) << "<dsig|+2,0>" << std::setw(15) << "err[2]"
            << std::setw(15) << "<ALL|0,0>" << std::setw(15) << "err[3]"
            << std::setw(15) << "<ALL|+1,0>" << std::setw(15) << "err[4]"
            << std::setw(15) << "<dsig|0,+1> " << std::setw(15) << "err[5]"
            << std::setw(15) << "<dsig|+1,-1>" << std::setw(15) << "err[6]"
            << std::setw(15) << "dsig|+1,+1>" << std::setw(15) << "err[7]"
            << std::setw(15) << "<dsig|+2,+1> " << std::setw(15) << "err[8]"
            << std::setw(15) << "<dsig|+2,-1> " << std::setw(15) << "err[9]"
            << std::setw(15) << "<ALL|0,+1>" << std::setw(15) << "err[10]"
            << std::setw(15) << "<ALL|+1,-1>" << std::setw(15) << "err[11]"
            << std::setw(15) << "<ALL|+1,+1>" << std::setw(15) << "err[12]"        
            << std::setw(15) << "time[s]"
            << std::setw(10) << "neval" << std::setw(10) << "fail"
            << std::endl;

    for (double xB_val : xB_values) {
        
        fixed_xB = xB_val; 
        std::cout << "Running scan for fixed xB = " << fixed_xB << std::endl;
        
        IntegrationResults res;
        res.maxeval = maxeval;
        res.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();

        Vegas(ndim, ncomp, integrand_fixedxB, nullptr,
                nvec, epsrel, epsabs,
                flags, seed, mineval, maxeval,
                nstart, nincrease, nbatch,
                gridno, statefile, spin,
                &res.neval, &res.fail,
                res.integral, res.error, res.prob);

        auto t1 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
        
        // Calculate ratios and ratio errors 
        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];
            res.ratioerror[i] = abs(res.ratio[i] * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2))
            );
        }
        outFile << std::setw(10) << xB_val  
                << std::setw(15) << res.integral[0] << std::setw(15) << res.error[0]
                << std::setw(15) << res.integral[1] << std::setw(15) << res.error[1]
                << std::setw(15) << res.integral[2] << std::setw(15) << res.error[2]
                << std::setw(15) << res.integral[3] << std::setw(15) << res.error[3]
                << std::setw(15) << res.integral[4] << std::setw(15) << res.error[4]
                << std::setw(15) << res.integral[5] << std::setw(15) << res.error[5]
                << std::setw(15) << res.integral[6] << std::setw(15) << res.error[6]
                << std::setw(15) << res.integral[7] << std::setw(15) << res.error[7]
                << std::setw(15) << res.integral[8] << std::setw(15) << res.error[8]
                << std::setw(15) << res.integral[9] << std::setw(15) << res.error[9]
                << std::setw(15) << res.integral[10] << std::setw(15) << res.error[10]
                << std::setw(15) << res.integral[11] << std::setw(15) << res.error[11]
                << std::setw(15) << res.integral[12] << std::setw(15) << res.error[12]
                << std::setw(15) << res.elapsed_seconds
                << std::setw(10) << res.neval << std::setw(10) << res.fail
                << std::endl;

        std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    }

    outFile.close();

    std::ofstream outFile1("Prova1.txt");  
    outFile1 << std::setw(10) << "y" 
            << std::setw(15) << "<dsig|+1,0>" << std::setw(15) << "err[1]"
            << std::setw(15) << "<dsig|+2,0>" << std::setw(15) << "err[2]"
            << std::setw(15) << "<ALL|0,0>" << std::setw(15) << "err[3]"
            << std::setw(15) << "<ALL|+1,0>" << std::setw(15) << "err[4]"
            << std::setw(15) << "<dsig|0,+1> " << std::setw(15) << "err[5]"
            << std::setw(15) << "<dsig|+1,-1>" << std::setw(15) << "err[6]"
            << std::setw(15) << "dsig|+1,+1>" << std::setw(15) << "err[7]"
            << std::setw(15) << "<dsig|+2,+1> " << std::setw(15) << "err[8]"
            << std::setw(15) << "<dsig|+2,-1> " << std::setw(15) << "err[9]"
            << std::setw(15) << "<ALL|0,+1>" << std::setw(15) << "err[10]"
            << std::setw(15) << "<ALL|+1,-1>" << std::setw(15) << "err[11]"
            << std::setw(15) << "<ALL|+1,+1>" << std::setw(15) << "err[12]"        
            << std::setw(15) << "time[s]"
            << std::setw(10) << "neval" << std::setw(10) << "fail"
            << std::endl;

    for (double y_val : y_values) {
        
        fixed_y = y_val; 
        std::cout << "Running scan for fixed y = " << fixed_y << std::endl;
        
        IntegrationResults res;
        res.maxeval = maxeval;
        res.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();

        Vegas(ndim, ncomp, integrand_fixedy, nullptr,
                nvec, epsrel, epsabs,
                flags, seed, mineval, maxeval,
                nstart, nincrease, nbatch,
                gridno, statefile, spin,
                &res.neval, &res.fail,
                res.integral, res.error, res.prob);

        auto t1 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
        // Calculate ratios and ratio errors
        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];
            res.ratioerror[i] =abs( res.ratio[i] * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2))
            );
        }
        
        outFile1 << std::setw(10) << y_val 
                << std::setw(15) << 2 * res.ratio[1] << std::setw(15) << 2 * res.ratioerror[1]
                << std::setw(15) << res.ratio[2] << std::setw(15) << res.ratioerror[2]
                << std::setw(15) << res.ratio[3] << std::setw(15) << res.ratioerror[3]
                << std::setw(15) << 2 * res.ratio[4] << std::setw(15) << 2 * res.ratioerror[4]
                << std::setw(15) << res.ratio[5] << std::setw(15) << res.ratioerror[5]
                << std::setw(15) << res.ratio[6] << std::setw(15) << res.ratioerror[6]
                << std::setw(15) << res.ratio[7] << std::setw(15) << res.ratioerror[7]
                << std::setw(15) << res.ratio[8] << std::setw(15) << res.ratioerror[8]
                << std::setw(15) << res.ratio[9] << std::setw(15) << res.ratioerror[9]
                << std::setw(15) << res.ratio[10] << std::setw(15) << res.ratioerror[10]
                << std::setw(15) << res.ratio[11] << std::setw(15) << res.ratioerror[11]
                << std::setw(15) << res.ratio[12] << std::setw(15) << res.ratioerror[12]
                << std::setw(15) << res.elapsed_seconds
                << std::setw(10) << res.neval << std::setw(10) << res.fail
                << std::endl;

        std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    }

    outFile1.close();      
    std::ofstream outFile2("Prova2.txt");  
    outFile2 << std::setw(10) << "Q2" 
            << std::setw(15) << "<dsig|+1,0>" << std::setw(15) << "err[1]"
            << std::setw(15) << "<dsig|+2,0>" << std::setw(15) << "err[2]"
            << std::setw(15) << "<ALL|0,0>" << std::setw(15) << "err[3]"
            << std::setw(15) << "<ALL|+1,0>" << std::setw(15) << "err[4]"
            << std::setw(15) << "<dsig|0,+1> " << std::setw(15) << "err[5]"
            << std::setw(15) << "<dsig|+1,-1>" << std::setw(15) << "err[6]"
            << std::setw(15) << "dsig|+1,+1>" << std::setw(15) << "err[7]"
            << std::setw(15) << "<dsig|+2,+1> " << std::setw(15) << "err[8]"
            << std::setw(15) << "<dsig|+2,-1> " << std::setw(15) << "err[9]"
            << std::setw(15) << "<ALL|0,+1>" << std::setw(15) << "err[10]"
            << std::setw(15) << "<ALL|+1,-1>" << std::setw(15) << "err[11]"
            << std::setw(15) << "<ALL|+1,+1>" << std::setw(15) << "err[12]"        
            << std::setw(15) << "time[s]"
            << std::setw(10) << "neval" << std::setw(10) << "fail"
            << std::endl;

    for (double Q2_val : Q2_values) {
        
        fixed_Q2 = Q2_val; 
        std::cout << "Running scan for fixed Q2 = " << fixed_Q2 << std::endl;
        
        IntegrationResults res;
        res.maxeval = maxeval;
        res.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();

        Vegas(ndim, ncomp, integrand_fixedQ2, nullptr,
                nvec, epsrel, epsabs,
                flags, seed, mineval, maxeval,
                nstart, nincrease, nbatch,
                gridno, statefile, spin,
                &res.neval, &res.fail,
                res.integral, res.error, res.prob);

        auto t1 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
        
        //Calculate ratios and ratio errors
        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];
            res.ratioerror[i] = abs(res.ratio[i]) * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2)
            );
        }
        
        outFile2<< std::setw(10) << Q2_val 
                << std::setw(15) << 2 * res.ratio[1] << std::setw(15) << 2 * res.ratioerror[1]
                << std::setw(15) << res.ratio[2] << std::setw(15) << res.ratioerror[2]
                << std::setw(15) << res.ratio[3] << std::setw(15) << res.ratioerror[3]
                << std::setw(15) << 2 * res.ratio[4] << std::setw(15) << 2 * res.ratioerror[4]
                << std::setw(15) << res.ratio[5] << std::setw(15) << res.ratioerror[5]
                << std::setw(15) << res.ratio[6] << std::setw(15) << res.ratioerror[6]
                << std::setw(15) << res.ratio[7] << std::setw(15) << res.ratioerror[7]
                << std::setw(15) << res.ratio[8] << std::setw(15) << res.ratioerror[8]
                << std::setw(15) << res.ratio[9] << std::setw(15) << res.ratioerror[9]
                << std::setw(15) << res.ratio[10] << std::setw(15) << res.ratioerror[10]
                << std::setw(15) << res.ratio[11] << std::setw(15) << res.ratioerror[11]
                << std::setw(15) << res.ratio[12] << std::setw(15) << res.ratioerror[12]
                << std::setw(15) << res.elapsed_seconds
                << std::setw(10) << res.neval << std::setw(10) << res.fail
                << std::endl;

        std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    }

    outFile2.close();
    return 0;
}