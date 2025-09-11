#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include <array>
#include <cuba.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstring>
#include <chrono>
#include <iomanip>

#include "FragFunct.h"
#include "COL.h"

#include "constants.h"

std::vector<double> pars = {0.5316, -.7290, .8246, 3.2293, 0.8568, -.3539, 1.7524, 0.0, 0.2551};
std::vector<double> COL_z(13), COL_ppz1(13), COL_ppz2(13), COL_pmz1(13), COL_pmz2(13);
std::vector<double> FF(13), FF_ppz1(13), FF_ppz2(13), FF_pmz1(13), FF_pmz2(13),\
    FF_ppz1_fixedQ2(13), FF_ppz2_fixedQ2(13), FF_pmz1_fixedQ2(13), FF_pmz2_fixedQ2(13),\
    FF_evo(13);
        
double f[13], h1[13]; //for HOPPET and transversity
    
const double charges[13] = {-2./3., 1./3., -2./3., 1./3., -2./3., 1./3., 0, -1./3., 2./3., -1./3., 2./3., -1./3., 2./3.};

double  numU = 0.0, denU = 0.0,\
        numL = 0.0, denL = 0.0,\
        numC = 0.0, denC = 0.0,\
        numFact = 0.0, denFact = 0.0,\
        z1 = 0.0, z2 = 0.0, pperp2 = 0.12;
        
double MC2 = pars.back();
double pperp2Col = MC2 * pperp2 / (MC2 + pperp2);
double prefact = (pi * EulerConst / 2) * (pow(pperp2Col, 3) / (pperp2 * pperp2 * MC2));

int charge = 1, hadron = 1;

FRAG::FF myFF("DEHSS","NLO");
COL::COLLINS myCol;

void Collins_FF(int & hadron, int & charge, double & z1, double & z2, double & Q2){

    //call to the Fragmentation Functions

//    cout << hadron << "\t" << charge << "\t" << z1 << "\t" << z2 << "\t" << Q2 << endl;
      
    myFF.FF_eval(hadron, charge, z1, Q2);
    for(int i = 0; i < FF_ppz1.size(); i++) FF_ppz1[i] = myFF.theFF[i] / z1;
        
    myFF.FF_eval(hadron, charge, z2, Q2);
    for(int i = 0; i < FF_ppz2.size(); i++) FF_ppz2[i] = myFF.theFF[i] / z2;
        
    charge*=-1; //to call pi- FFs
        
    myFF.FF_eval(hadron, charge, z1, Q2);
    for(int i = 0; i < FF_pmz1.size(); i++) FF_pmz1[i] = myFF.theFF[i] / z1;
        
    myFF.FF_eval(hadron, charge, z2, Q2);
    for(int i = 0; i < FF_pmz2.size(); i++) FF_pmz2[i] = myFF.theFF[i] / z2;

    // for(int i = 0; i < FF_pmz2.size(); i++) cout << FF_pmz2[i] << "\t";
    // cout << endl;
       
}


void Collins_epem_loop(const std::vector<double> &COL_ppz1_in, const std::vector<double> &COL_ppz2_in, const std::vector<double> &COL_pmz1_in, const std::vector<double> &COL_pmz2_in, const std::vector<double> &FF_ppz1_in, const std::vector<double> &FF_ppz2_in, const std::vector<double> &FF_pmz1_in, const std::vector<double> &FF_pmz2_in){
    
    for(int i = 3; i <= 9; i++){
            
        if(i == 3){ //sb
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i+6] + COL_pmz1_in[i] * COL_ppz2_in[i+6]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i+6] + FF_pmz1_in[i] * FF_ppz2_in[i+6]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i+6] + COL_pmz1_in[i] * COL_pmz2_in[i+6]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i+6] + FF_pmz1_in[i] * FF_pmz2_in[i+6]);
        }                                
                                         
        if(i == 4){ //ub                 
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i+4] + COL_pmz1_in[i] * COL_ppz2_in[i+4]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i+4] + FF_pmz1_in[i] * FF_ppz2_in[i+4]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i+4] + COL_pmz1_in[i] * COL_pmz2_in[i+4]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i+4] + FF_pmz1_in[i] * FF_pmz2_in[i+4]);
        }                                
                                         
        if(i == 5){ //db                 
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i+2] + COL_pmz1_in[i] * COL_ppz2_in[i+2]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i+2] + FF_pmz1_in[i] * FF_ppz2_in[i+2]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i+2] + COL_pmz1_in[i] * COL_pmz2_in[i+2]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i+2] + FF_pmz1_in[i] * FF_pmz2_in[i+2]);
        }                                
                                         
        if(i == 6){ //g                  
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i] + COL_pmz1_in[i] * COL_ppz2_in[i]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i] + FF_pmz1_in[i] * FF_ppz2_in[i]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i] + COL_pmz1_in[i] * COL_pmz2_in[i]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i] + FF_pmz1_in[i] * FF_pmz2_in[i]);
        }                                
                                         
        if(i == 7){ //d                  
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i-2] + COL_pmz1_in[i] * COL_ppz2_in[i-2]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i-2] + FF_pmz1_in[i] * FF_ppz2_in[i-2]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i-2] + COL_pmz1_in[i] * COL_pmz2_in[i-2]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i-2] + FF_pmz1_in[i] * FF_pmz2_in[i-2]);            
        }                                
                                         
        if(i == 8){ //u                  
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i-4] + COL_pmz1_in[i] * COL_ppz2_in[i-4]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i-4] + FF_pmz1_in[i] * FF_ppz2_in[i-4]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i-4] + COL_pmz1_in[i] * COL_pmz2_in[i-4]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i-4] + FF_pmz1_in[i] * FF_pmz2_in[i-4]);
        }                                
                                         
        if(i == 9){ //s                  
            numU += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i-6] + COL_pmz1_in[i] * COL_ppz2_in[i-6]);
            denU += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i-6] + FF_pmz1_in[i] * FF_ppz2_in[i-6]);
            numL += pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i-6] + COL_pmz1_in[i] * COL_pmz2_in[i-6]);
            denL += pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i-6] + FF_pmz1_in[i] * FF_pmz2_in[i-6]);
        }  
    }
    
}



class PhysicsCalculator
{
private:
    // Constants
    static constexpr double alphaem = 1.0 / 137.036;
    static constexpr double thetac = 3.0e-2; // rad [value for FCC-ee]; search for adequate values for other colliders
    static constexpr double me = 5.11e-4;    // in GeV
    static constexpr double SqrtS = 92.0;    // GeV [value for FCC-ee: 92.]
    static constexpr double S = SqrtS * SqrtS;
    static constexpr double PI = 3.14159265358979323846;

    // Variable extremes
    static constexpr double Q20 = 10.0; // GeV
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

    // A coefficients // Il extra 2 term is for the denominator in the azimuthal moments... ->> PERCHé??
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

    // Coefficienti AU,AL,BU,BL integrati sul range permesso di zeta *)
    // A eccezione di quelli che integrano a zero tra zetamin e zetamax per simmetria, integrati solo da zeta = 1/2 a zetamax *)
    // per cui bisogna ricordarsi che consistentemente vanno divisi per AU/2 *)
    // Nota: usato 1 - zetamax = zetamin e viceversa *)
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
        double term2 = sqrt(xB * (csi - xB)) * ( 2 * sqrt(zetamax * zetamin) - 1); 
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
       
        double Q2 = fixed_Q2;

        Collins_FF(hadron, charge, z1, z2, Q2);


        //call for the Collins functions
        if(myCol.evo == "DGLAP" || myCol.evo == "none"){ 
            
            myCol.eval(z1, Q2, charge, FF_ppz1, pars);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, Q2, charge, FF_ppz2, pars);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = myCol.COL_z[i];
            
            charge *= -1; //to call pi- Collins
                
            myCol.eval(z1, Q2, charge, FF_pmz1, pars);
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, Q2, charge, FF_pmz2, pars);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = myCol.COL_z[i];
        }
        
        if(myCol.evo == "CT3"){
        
            //call for the Collins functions
            hoppetEvalcf(z1, sqrt(Q2), f);
//             COL_ppz1 = CollinsHoppetEval(z1, int_charge, f);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = f[i] / z1;

            hoppetEvalcf(z2, sqrt(Q2), f);
//             COL_ppz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = f[i] / z2;
            
            charge *= -1; //to call pi- Collins
                
            hoppetEvalcff(z1, sqrt(Q2), f);      //modified to use proper hoppet calls for pip and pim
//             COL_pmz1 = CollinsHoppetEval(z1, int_charge, f);    
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = f[i] / z1;

            hoppetEvalcff(z2, sqrt(Q2), f);
//             COL_pmz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = f[i] / z2;

        }
        
        //calling the loop to calculate numerator and denominator of A0
        Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1, FF_ppz2, FF_pmz1, FF_pmz2);

        results[0]  = pc.AUzi()*pc.KQx();// * denU; // moltiplicare per D1D1
        results[1]  = pc.AUcfqzi()*pc.KQx();// * denU; // moltiplicare per D1D1
        results[2]  = pc.AUc2fqzi()*pc.KQx();// * denU; // moltiplicare per D1D1
        results[3]  = pc.ALzi()*pc.KQx();// * denU; // moltiplicare per D1D1
        results[4]  = pc.ALcfqzi()*pc.KQx();// * denU; // moltiplicare per D1D1
        results[5]  = pc.BUcf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[6]  = pc.BUcfqmf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[7]  = pc.BUcfqpf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[8]  = pc.BUc2fqmf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[9]  = pc.BUc2fqpf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[10] = pc.BLcf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[11] = pc.BLcfqmf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
        results[12] = pc.BLcfqpf12zi()*pc.KQx();// * numU * prefact; // moltiplicare per H1p1H1p
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

int main(int argc,char *argv[]) {
    // Fixed integration settings
    const int maxeval = static_cast<int>(1e7);
    const int nstart  = static_cast<int>(1e6);

    //TO DO: input arguments here
    stringstream modelstr, widthsstr, evostr, epemstr;
    modelstr << argv[2];
    widthsstr << argv[4];
    evostr << argv[6];
    epemstr << argv[8];
    z1 = atof(argv[10]);
    z2 = atof(argv[12]);
    
    std::string model = modelstr.str();
    std::string widths = widthsstr.str();
    std::string evo = evostr.str();
    std::string epem_sign = epemstr.str();

    myFF.use_LHAPDF(true);
    cout << myFF.FFset << "\t" << myFF.FForder << endl;
    if(myFF.useLHAPDF){

        if(myFF.FFset == "NNFF10"){

            if(myFF.FForder == "LO") myFF.set_LHAPDF_FFset(NNFF10setsLO);
            if(myFF.FForder == "NLO") myFF.set_LHAPDF_FFset(NNFF10setsNLO);
        }
        if(myFF.FFset == "DEHSS") myFF.set_LHAPDF_FFset(DEHSSsetsNLO);

    }


    myCol.set_model(model);
    myCol.set_widths(widths);
    myCol.set_evolution(evo);
    
//     cout << MC2 << "\t" << pperp2 << "\t"<< pperp2Col  << "\t" << prefact <<endl;

    // Fixed values for xB, y, and Q2
    std::vector<double> xB_values;
    for(double i = 0; i <=1; i += 0.001) xB_values.push_back(i);
    std::vector<double> y_values;
    for(double i = 0; i <=1; i += 0.001) y_values.push_back(i);
    std::vector<double> Q2_values;
    for(double i = 0; i <= 1; i += 0.001) Q2_values.push_back(i * (PhysicsCalculator::getQ2M() - PhysicsCalculator::getQ20()) + PhysicsCalculator::getQ20()); // Assuming Q2 ranges from 0 to Q2M
    
    //set parameters for the integration
    int ndim = 2, ncomp = 13, nvec = 1;
    double epsrel = 1e-6, epsabs = 1e-12;
    int flags = 2, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;
    
    std::ofstream outFile("Fixed_xB.txt");  
    outFile << std::setw(10) << "xB" 
            << std::setw(15) << "<dsig|+1,0>" << std::setw(15) << "err[1]"
            << std::setw(15) << "<dsig|+2,0>" << std::setw(15) << "err[2]"
            << std::setw(15) << "<ALL|0,0>" << std::setw(15) << "err[3]"
            << std::setw(15) << "<ALL|+1,0>" << std::setw(15) << "err[4]"
            << std::setw(15) << "<dsig|0,+1> " << std::setw(15) << "err[5]"
            << std::setw(15) << "<dsig|+1,-1>" << std::setw(15) << "err[6]"
            << std::setw(15) << "<dsig|+1,+1>" << std::setw(15) << "err[7]"
            << std::setw(15) << "<dsig|+2,+1> " << std::setw(15) << "err[8]"
            << std::setw(15) << "<dsig|+2,-1> " << std::setw(15) << "err[9]"
            << std::setw(15) << "<ALL|0,+1>" << std::setw(15) << "err[10]"
            << std::setw(15) << "<ALL|+1,-1>" << std::setw(15) << "err[11]"
            << std::setw(15) << "<ALL|+1,+1>" << std::setw(15) << "err[12]"
            << std::setw(15) << "AUzi" << std::setw(15) << "err[13]"
            << std::setw(15) << "AUcfqzi" << std::setw(15) << "err[14]"
            << std::setw(15) << "AUc2fqzi" << std::setw(15) << "err[15]"
            << std::setw(15) << "ALzi" << std::setw(15) << "err[16]"
            << std::setw(15) << "ALcfqzi" << std::setw(15) << "err[17]"
            << std::setw(15) << "BUcf12zi" << std::setw(15) << "err[18]"
            << std::setw(15) << "BUcfqmf12zi" << std::setw(15) << "err[19]"
            << std::setw(15) << "BUcfqpf12zi" << std::setw(15) << "err[20]"
            << std::setw(15) << "BUc2fqmf12zi" << std::setw(15) << "err[21]"
            << std::setw(15) << "BUc2fqpf12zi" << std::setw(15) << "err[22]"
            << std::setw(15) << "BLcf12zi" << std::setw(15) << "err[23]"
            << std::setw(15) << "BLcfqmf12zi" << std::setw(15) << "err[24]"
            << std::setw(15) << "BLcfqpf12zi" << std::setw(15) << "err[25]"
            << std::setw(15) << "time[s]"
            << std::setw(10) << "neval" << std::setw(10) << "fail"
            << std::endl;

    // for (double xB_val : xB_values) {
        
    //     fixed_xB = xB_val; 
    //     std::cout << "Running scan for fixed xB = " << fixed_xB << std::endl;
        
    //     IntegrationResults res;
    //     res.maxeval = maxeval;
    //     res.nstart  = nstart;

    //     auto t0 = std::chrono::high_resolution_clock::now();

    //     Vegas(ndim, ncomp, integrand_fixedxB, nullptr,
    //             nvec, epsrel, epsabs,
    //             flags, seed, mineval, maxeval,
    //             nstart, nincrease, nbatch,
    //             gridno, statefile, spin,
    //             &res.neval, &res.fail,
    //             res.integral, res.error, res.prob);

    //     auto t1 = std::chrono::high_resolution_clock::now();
    //     res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
        
    //     // Calculate ratios and ratio errors 
    //     for (int i = 0; i < ncomp; ++i) {
    //         res.ratio[i] = res.integral[i] / res.integral[0];
    //         res.ratioerror[i] = abs(res.ratio[i] * std::sqrt(
    //             std::pow(res.error[i] / res.integral[i], 2) +
    //             std::pow(res.error[0] / res.integral[0], 2))
    //         );
    //     }
    //     outFile << std::setw(10) << xB_val  
    //             << std::setw(15) << 2 * res.ratio[1] << std::setw(15) << 2 * res.ratioerror[1]
    //             << std::setw(15) << res.ratio[2] << std::setw(15) << res.ratioerror[2]
    //             << std::setw(15) << res.ratio[3] << std::setw(15) << res.ratioerror[3]
    //             << std::setw(15) << 2 * res.ratio[4] << std::setw(15) << 2 * res.ratioerror[4]
    //             << std::setw(15) << res.ratio[5] << std::setw(15) << res.ratioerror[5]
    //             << std::setw(15) << res.ratio[6] << std::setw(15) << res.ratioerror[6]
    //             << std::setw(15) << res.ratio[7] << std::setw(15) << res.ratioerror[7]
    //             << std::setw(15) << res.ratio[8] << std::setw(15) << res.ratioerror[8]
    //             << std::setw(15) << res.ratio[9] << std::setw(15) << res.ratioerror[9]
    //             << std::setw(15) << res.ratio[10] << std::setw(15) << res.ratioerror[10]
    //             << std::setw(15) << res.ratio[11] << std::setw(15) << res.ratioerror[11]
    //             << std::setw(15) << res.ratio[12] << std::setw(15) << res.ratioerror[12]
    //             << std::setw(15) << res.integral[0] << std::setw(15) << res.error[0]
    //             << std::setw(15) << res.integral[1] << std::setw(15) << res.error[1]
    //             << std::setw(15) << res.integral[2] << std::setw(15) << res.error[2]
    //             << std::setw(15) << res.integral[3] << std::setw(15) << res.error[3]
    //             << std::setw(15) << res.integral[4] << std::setw(15) << res.error[4]
    //             << std::setw(15) << res.integral[5] << std::setw(15) << res.error[5]
    //             << std::setw(15) << res.integral[6] << std::setw(15) << res.error[6]
    //             << std::setw(15) << res.integral[7] << std::setw(15) << res.error[7]
    //             << std::setw(15) << res.integral[8] << std::setw(15) << res.error[8]
    //             << std::setw(15) << res.integral[9] << std::setw(15) << res.error[9]
    //             << std::setw(15) << res.integral[10] << std::setw(15) << res.error[10]
    //             << std::setw(15) << res.integral[11] << std::setw(15) << res.error[11]
    //             << std::setw(15) << res.integral[12] << std::setw(15) << res.error[12]
    //             << std::setw(15) << res.elapsed_seconds
    //             << std::setw(10) << res.neval << std::setw(10) << res.fail
    //             << std::endl;

    //     std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    // }

    // outFile.close();

    // std::ofstream outFile1("Fixed_y.txt");  
    // outFile1 << std::setw(10) << "y" 
    //         << std::setw(15) << "<dsig|+1,0>" << std::setw(15) << "err[1]"
    //         << std::setw(15) << "<dsig|+2,0>" << std::setw(15) << "err[2]"
    //         << std::setw(15) << "<ALL|0,0>" << std::setw(15) << "err[3]"
    //         << std::setw(15) << "<ALL|+1,0>" << std::setw(15) << "err[4]"
    //         << std::setw(15) << "<dsig|0,+1> " << std::setw(15) << "err[5]"
    //         << std::setw(15) << "<dsig|+1,-1>" << std::setw(15) << "err[6]"
    //         << std::setw(15) << "dsig|+1,+1>" << std::setw(15) << "err[7]"
    //         << std::setw(15) << "<dsig|+2,+1> " << std::setw(15) << "err[8]"
    //         << std::setw(15) << "<dsig|+2,-1> " << std::setw(15) << "err[9]"
    //         << std::setw(15) << "<ALL|0,+1>" << std::setw(15) << "err[10]"
    //         << std::setw(15) << "<ALL|+1,-1>" << std::setw(15) << "err[11]"
    //         << std::setw(15) << "<ALL|+1,+1>" << std::setw(15) << "err[12]"
    //         << std::setw(15) << "AUzi" << std::setw(15) << "err[13]"
    //         << std::setw(15) << "AUcfqzi" << std::setw(15) << "err[14]"
    //         << std::setw(15) << "AUc2fqzi" << std::setw(15) << "err[15]"
    //         << std::setw(15) << "ALzi" << std::setw(15) << "err[16]"
    //         << std::setw(15) << "ALcfqzi" << std::setw(15) << "err[17]"
    //         << std::setw(15) << "BUcf12zi" << std::setw(15) << "err[18]"
    //         << std::setw(15) << "BUcfqmf12zi" << std::setw(15) << "err[19]"
    //         << std::setw(15) << "BUcfqpf12zi" << std::setw(15) << "err[20]"
    //         << std::setw(15) << "BUc2fqmf12zi" << std::setw(15) << "err[21]"
    //         << std::setw(15) << "BUc2fqpf12zi" << std::setw(15) << "err[22]"
    //         << std::setw(15) << "BLcf12zi" << std::setw(15) << "err[23]"
    //         << std::setw(15) << "BLcfqmf12zi" << std::setw(15) << "err[24]"
    //         << std::setw(15) << "BLcfqpf12zi" << std::setw(15) << "err[25]"        
    //         << std::setw(15) << "time[s]"
    //         << std::setw(10) << "neval" << std::setw(10) << "fail"
    //         << std::endl;

    // for (double y_val : y_values) {
        
    //     fixed_y = y_val; 
    //     std::cout << "Running scan for fixed y = " << fixed_y << std::endl;
        
    //     IntegrationResults res;
    //     res.maxeval = maxeval;
    //     res.nstart  = nstart;

    //     auto t0 = std::chrono::high_resolution_clock::now();

    //     Vegas(ndim, ncomp, integrand_fixedy, nullptr,
    //             nvec, epsrel, epsabs,
    //             flags, seed, mineval, maxeval,
    //             nstart, nincrease, nbatch,
    //             gridno, statefile, spin,
    //             &res.neval, &res.fail,
    //             res.integral, res.error, res.prob);

    //     auto t1 = std::chrono::high_resolution_clock::now();
    //     res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    //     // Calculate ratios and ratio errors
    //     for (int i = 0; i < ncomp; ++i) {
    //         res.ratio[i] = res.integral[i] / res.integral[0];
    //         res.ratioerror[i] =abs( res.ratio[i] * std::sqrt(
    //             std::pow(res.error[i] / res.integral[i], 2) +
    //             std::pow(res.error[0] / res.integral[0], 2))
    //         );
    //     }
        
    //     outFile1 << std::setw(10) << y_val 
    //             << std::setw(15) << 2 * res.ratio[1] << std::setw(15) << 2 * res.ratioerror[1]
    //             << std::setw(15) << res.ratio[2] << std::setw(15) << res.ratioerror[2]
    //             << std::setw(15) << res.ratio[3] << std::setw(15) << res.ratioerror[3]
    //             << std::setw(15) << 2 * res.ratio[4] << std::setw(15) << 2 * res.ratioerror[4]
    //             << std::setw(15) << res.ratio[5] << std::setw(15) << res.ratioerror[5]
    //             << std::setw(15) << res.ratio[6] << std::setw(15) << res.ratioerror[6]
    //             << std::setw(15) << res.ratio[7] << std::setw(15) << res.ratioerror[7]
    //             << std::setw(15) << res.ratio[8] << std::setw(15) << res.ratioerror[8]
    //             << std::setw(15) << res.ratio[9] << std::setw(15) << res.ratioerror[9]
    //             << std::setw(15) << res.ratio[10] << std::setw(15) << res.ratioerror[10]
    //             << std::setw(15) << res.ratio[11] << std::setw(15) << res.ratioerror[11]
    //             << std::setw(15) << res.ratio[12] << std::setw(15) << res.ratioerror[12]
    //             << std::setw(15) << res.integral[0] << std::setw(15) << res.error[0]
    //             << std::setw(15) << res.integral[1] << std::setw(15) << res.error[1]
    //             << std::setw(15) << res.integral[2] << std::setw(15) << res.error[2]
    //             << std::setw(15) << res.integral[3] << std::setw(15) << res.error[3]
    //             << std::setw(15) << res.integral[4] << std::setw(15) << res.error[4]
    //             << std::setw(15) << res.integral[5] << std::setw(15) << res.error[5]
    //             << std::setw(15) << res.integral[6] << std::setw(15) << res.error[6]
    //             << std::setw(15) << res.integral[7] << std::setw(15) << res.error[7]
    //             << std::setw(15) << res.integral[8] << std::setw(15) << res.error[8]
    //             << std::setw(15) << res.integral[9] << std::setw(15) << res.error[9]
    //             << std::setw(15) << res.integral[10] << std::setw(15) << res.error[10]
    //             << std::setw(15) << res.integral[11] << std::setw(15) << res.error[11]
    //             << std::setw(15) << res.integral[12] << std::setw(15) << res.error[12]
    //             << std::setw(15) << res.elapsed_seconds
    //             << std::setw(10) << res.neval << std::setw(10) << res.fail
    //             << std::endl;

    //     std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    // }

    // outFile1.close();      
    std::ofstream outFile2("Fixed_Q2.txt");  
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
            << std::setw(15) << "AUzi" << std::setw(15) << "err[13]"
            << std::setw(15) << "AUcfqzi" << std::setw(15) << "err[14]"
            << std::setw(15) << "AUc2fqzi" << std::setw(15) << "err[15]"
            << std::setw(15) << "ALzi" << std::setw(15) << "err[16]"
            << std::setw(15) << "ALcfqzi" << std::setw(15) << "err[17]"
            << std::setw(15) << "BUcf12zi" << std::setw(15) << "err[18]"
            << std::setw(15) << "BUcfqmf12zi" << std::setw(15) << "err[19]"
            << std::setw(15) << "BUcfqpf12zi" << std::setw(15) << "err[20]"
            << std::setw(15) << "BUc2fqmf12zi" << std::setw(15) << "err[21]"
            << std::setw(15) << "BUc2fqpf12zi" << std::setw(15) << "err[22]"
            << std::setw(15) << "BLcf12zi" << std::setw(15) << "err[23]"
            << std::setw(15) << "BLcfqmf12zi" << std::setw(15) << "err[24]"
            << std::setw(15) << "BLcfqpf12zi" << std::setw(15) << "err[25]"        
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
            
            if (i >= 5){
                
                if(epem_sign == "UL") res.ratio[i] *= prefact * ( numU / denU - numL / denL );
                if(epem_sign == "UC") res.ratio[i] *= prefact * ( numU / denU - (numU + numL) / (denU + denL) );
            }
            
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

    outFile2.close();
    return 0;
}
