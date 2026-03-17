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
#include <numeric>

#include "FragFunct.h"
#include "COL.h"
#include "photon_flux.h"
#include <rapidcsv.h>

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
        sqrts = 0.0, thetac = 0.0, Q20 = 0.0,\
        z1 = 0.0, z2 = 0.0, pperp2 = 0.12;

double MC2 = pars.back();
double pperp2Col = MC2 * pperp2 / (MC2 + pperp2);
double prefact = (pi * EulerConst / 2) * (pow(pperp2Col, 3) / (pperp2 * pperp2 * MC2));

int charge = 1, hadron = 1;

FRAG::FF myFF("DEHSS","NLO");
COL::COLLINS myCol;
EPA::EPA_flux flux;

rapidcsv::Document * CSV = new rapidcsv::Document;


double computeMean(const std::vector<double>& values, const std::vector<double>& probabilities) {
    // Ensure that values and probabilities have the same size
    if (values.size() != probabilities.size()) {
        std::cerr << "Error: Size mismatch between values and probabilities." << std::endl;
        return 0.0; // You might want to handle this error differently based on your requirements
    }

    // Find the median value
    double meanValue = 0.0;
    for (size_t i = 0; i < values.size(); ++i) meanValue += probabilities[i] * values[i];

    return meanValue;
}

double computeMedian(const std::vector<double>& values, const std::vector<double>& probabilities) {
    // Ensure that values and probabilities have the same size
    if (values.size() != probabilities.size()) {
        std::cerr << "Error: Size mismatch between values and probabilities." << std::endl;
        return 0.0; // You might want to handle this error differently based on your requirements
    }

    // Create a vector of pairs (value, probability)
    std::vector<std::pair<double, double>> valueProbPairs;
    for (size_t i = 0; i < values.size(); ++i) {
        valueProbPairs.push_back({values[i], probabilities[i]});
    }

    // Sort the pairs by values
    std::sort(valueProbPairs.begin(), valueProbPairs.end(),
              [](const auto& lhs, const auto& rhs) {
                  return lhs.first < rhs.first;
              });

    // Calculate the cumulative probabilities
    std::vector<double> cumulativeProbabilities(valueProbPairs.size(), 0.0);
    cumulativeProbabilities[0] = valueProbPairs[0].second;

    for (size_t i = 1; i < valueProbPairs.size(); ++i) {
        cumulativeProbabilities[i] = cumulativeProbabilities[i - 1] + valueProbPairs[i].second;
    }

    // Find the median value
    double medianValue = 0.0;
    for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
        if (cumulativeProbabilities[i] >= 0.5) {
            medianValue = valueProbPairs[i].first;
            break;
        }
    }

    return medianValue;
}

std::pair<double, double> calculateDiscreteInterval(
    const std::vector<double>& values,
    const std::vector<double>& probabilities,
    double confidenceLevel
) {
    if (confidenceLevel <= 0.0 || confidenceLevel >= 1.0) {
        std::cerr << "Error: Confidence level must be in (0, 1)." << std::endl;
        return {0.0, 0.0};
    }

    // Normalize probabilities to ensure they sum to 1
    std::vector<double> normalizedProbabilities = probabilities;
    double sum = std::accumulate(normalizedProbabilities.begin(), normalizedProbabilities.end(), 0.0);
    for (auto& prob : normalizedProbabilities) {
        prob /= sum;
    }

    // Create a vector of pairs (value, probability)
    std::vector<std::pair<double, double>> valueProbPairs;
    for (size_t i = 0; i < values.size(); ++i) {
        valueProbPairs.push_back({values[i], normalizedProbabilities[i]});
    }

    // Sort the pairs by values
    std::sort(valueProbPairs.begin(), valueProbPairs.end(),
              [](const auto& lhs, const auto& rhs) {
                  return lhs.first < rhs.first;
              });

    // Calculate the cumulative probabilities
    std::vector<double> cumulativeProbabilities(valueProbPairs.size(), 0.0);
    cumulativeProbabilities[0] = valueProbPairs[0].second;

    for (size_t i = 1; i < valueProbPairs.size(); ++i) {
        cumulativeProbabilities[i] = cumulativeProbabilities[i - 1] + valueProbPairs[i].second;
    }

    // Find the lower bound
    double lowerBound = 0.0;
    for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
        if (cumulativeProbabilities[i] >= (1.0 - confidenceLevel) / 2.0) {
            lowerBound = valueProbPairs[i].first;
            break;
        }
    }

    // Find the upper bound
    double upperBound = 0.0;
    for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
        if (cumulativeProbabilities[i] > 1.0 - (1.0 - confidenceLevel) / 2.0) {
            upperBound = valueProbPairs[i].first;
            break;
        }
    }

    return {lowerBound, upperBound};
}



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
    static double _thetac; // rad [value for FCC-ee]; search for adequate values for other colliders
    static constexpr double me = 5.11e-4;    // in GeV
    static constexpr double mmu = 1.056e-1;    // in GeV
    static double _SqrtS;    // GeV [value for FCC-ee: 92.]
    static double _S;
    static constexpr double PI = 3.14159265358979323846;

    // Variable extremes
    // static constexpr double Q20 = 10.0; // GeV
    // static constexpr double Q2M = S - 4 * Q20;
    static double _Q20; // GeV
    static double _Q2M;
    static constexpr double csimax = 1.0;
    static constexpr double ymax = 1.0;

    // Cost factor
    static constexpr double Cost = 3 * alphaem * alphaem * alphaem / (4 * PI);

    // Member variables for current calculation
    double xB, y, csi, Q2;//, SqrtS, Q20, Q2M, S, thetac;
    double zeta, zetamin, zetamax;
    double shat;
    double Q2min, Q2max;
    double fgl, DLfgl;


public:

    double SqrtS, S, Q20, Q2M, thetac;

    // Constructor
    PhysicsCalculator() : xB(0), y(0), csi(0), Q2(0), zeta(0), zetamin(0), zetamax(0), shat(0), Q2min(0), Q2max(0), fgl(0), DLfgl(0) {}

    void set_SqrtS(double val){

        SqrtS = val;
        S = SqrtS * SqrtS;

    }

    void set_Q20(double val){

        Q20 = val;
        Q2M = S - 4 * Q20;

    }

    void setExperiment(double sqrts_val, double Q20_val, double thetac_val){

        set_SqrtS(sqrts_val);
        set_Q20(Q20_val);
        thetac = thetac_val;

    }

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


        ww = flux.eval(csi);
        flux.set_polarisation(-1);
        DLfgl = flux.eval(csi);
        flux.set_polarisation(1);

        // // WW function parameters
        // Q2min = me * me * csi * csi / (1 - csi);
        // Q2max = (S / 4) * thetac * thetac * (1 - csi) + Q2min;
        //
        // // Calculate WW functions
        // double log_ratio = log(Q2max / Q2min);
        // double me2_term = 2 * me * me * csi * (1.0 / Q2max - 1.0 / Q2min);
        //
        // fgl = (alphaem / (2 * PI)) * ((1 + (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);
        // // DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + 2 * me * me * csi * csi * (1.0 / Q2max - 1.0 / Q2min));
        // DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);
        //
        // std::cout << "Check WW:   " << fgl << "  "<< ww << " ... " << DLfgl << "  "<< Dww << std::endl;
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

        fgl = flux.eval(csi);
        flux.set_polarisation(-1);
        Dfgl = flux.eval(csi);
        flux.set_polarisation(1);

        // Q2min = me * me * csi * csi / (1 - csi);
        // Q2max = (S / 4) * thetac * thetac * (1 - csi) + Q2min;
        //
        // double log_ratio = log(Q2max / Q2min);
        // double me2_term = 2 * me * me * csi * (1.0 / Q2max - 1.0 / Q2min);
        //
        // fgl = (alphaem / (2 * PI)) * ((1 + (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);
        // // DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + 2 * me * me * csi * csi * (1.0 / Q2max - 1.0 / Q2min));
        // DLfgl = (alphaem / (2 * PI)) * ((1 - (1 - csi) * (1 - csi)) / csi * log_ratio + me2_term);

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

    // Getter functions for current variable values1
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
    // static
    double getThetac() { return thetac; }
    static double getMe() { return me; }
    // static
    double getSqrtS() { return SqrtS; }
    // static
    double getS() { return S; }
    // static
    double getQ20() { return Q20; }
    // static
    double getQ2M() { return Q2M; }
    static double getCost() { return Cost; }
};

static double fixed_xB ; // Global variable to hold fixed xB value
static double fixed_y ; // Global variable to hold fixed y value
static double fixed_Q2 ; // Global variable to hold fixed Q2 value

//Cut function to check if the input values are within the allowed ranges
bool cut(const PhysicsCalculator& pscut)
{
    // double S = PhysicsCalculator::getS();
    // double Q20 = PhysicsCalculator::getQ20();
    // double Q2M = PhysicsCalculator::getQ2M();
    double S = pscut.S;
    double Q20 = pscut.Q20;
    double Q2M = pscut.Q2M;

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
    pc.setExperiment(sqrts, Q20, thetac);

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
    pc.setExperiment(sqrts, Q20, thetac);

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

void ValuesfixedQ2(const std::vector<double> &x, double results[], void *userdata)
{
    // x[0] = xB, x[1] = csi; fixed Q2 is provided externally

    double sqrts = *static_cast<double **>(userdata)[0];
    double Q20 = *static_cast<double **>(userdata)[1];
    double thetac = *static_cast<double **>(userdata)[2];
    // double z1_min = *static_cast<double **>(userdata)[2];
    // double z1_max = *static_cast<double **>(userdata)[2];

    // std::cout << "In ValuesfixedQ2: "<< sqrts << "\t" << Q20 << "\t" << thetac << std::endl;

    PhysicsCalculator pc;
    pc.setExperiment(sqrts, Q20, thetac);

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

    ValuesfixedQ2(xv, f, userdata);

    return 0;
}

int integrand_Collins(const int *ndim, const double x[], const int *ncomp, double ff[], void *userdata)
{

#define f0 ff[0]
#define f1 ff[1]
#define f2 ff[2]
#define f3 ff[3]

    // std::cout<< "Integrating Collins functions on z1" << std::endl;

    double Q2 = *static_cast<double **>(userdata)[3];
    double z1_min = *static_cast<double **>(userdata)[4];
    double z1_max = *static_cast<double **>(userdata)[5];
    double z2     = *static_cast<double **>(userdata)[6];

    // std::cout << "In integrand: z1_min = "<< z1_min << "  z1_max = "<< z1_max << "  z2 = " << z2 << "   Q2 = "<< Q2 <<std::endl;

    double z1 = z1_min + x[0] * (z1_max - z1_min);

    double f[13];

    // Collins_PT(z1, Q2, pars, int_fname); //note: av_z1 is not used here, just needed to call the function to fix FC_ee factor

    Collins_FF(hadron, charge, z1, z2, Q2);

//     cout<<"In integrand: int_had = "<<int_hadron<<"\t int_ch = "<<charge<<endl;


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
//         COL_ppz1 = CollinsHoppetEval(z1, charge, f);
        for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = f[i] / z1;

        hoppetEvalcf(z2, sqrt(Q2), f);
//         COL_ppz2 = CollinsHoppetEval(z2, charge, f);
        for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = f[i] / z2;

        charge *= -1; //to call pi- Collins

        hoppetEvalcff(z1, sqrt(Q2), f);      //modified to use proper hoppet calls for pip and pim
//         COL_pmz1 = CollinsHoppetEval(z1, charge, f);
        for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = f[i] / z1;

        hoppetEvalcff(z2, sqrt(Q2), f);
//         COL_pmz2 = CollinsHoppetEval(z2, charge, f);
        for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = f[i] / z2;

    }

    //calling the loop to calculate numerator and denominator of A0
    Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1, FF_ppz2, FF_pmz1, FF_pmz2);

    f0 = denU;
    f1 = numU;
    f2 = denL;
    f3 = numL;

    return 0;

}



struct IntegrationResults {
    int maxeval;
    int nstart;
    int neval;
    int fail;
    int nregions;
    double integral[13];
    double error[13];
    double prob[13];
    double ratio[13];
    double ratioerror[13];
    double elapsed_seconds;
};


PhysicsCalculator pc;


int main(int argc,char *argv[]) {
    // Fixed integration settings
    const int maxeval = static_cast<int>(1e7);
    const int nstart  = static_cast<int>(1e6);

    //TO DO: input arguments here
    stringstream modelstr, widthsstr, evostr, MCnamestring, EPAnamestring;
    modelstr << argv[2];
    widthsstr << argv[4];
    evostr << argv[6];
    double z1_min = atof(argv[8]);
    double z1_max = atof(argv[9]);
    double z2_min = atof(argv[11]);
    double z2_max = atof(argv[12]);
    sqrts = atof(argv[14]);
    Q20 = atof(argv[16]);
    double Q2val = atof(argv[18]);
    thetac = atof(argv[20]);
    EPAnamestring << argv[22];
    MCnamestring << argv[24];

    std::string MCname = MCnamestring.str(), EPAname = EPAnamestring.str();

    CSV->Load(MCname);
    int Nrep = CSV->GetRowCount();

    std::vector<double> z2_values, params;
    for (double i = z2_min; i <= z2_max + 0.05; i+=0.05) z2_values.push_back(i);
    // for (int i = 0; i < z2_values.size(); i++) std::cout << z2_values[i] <<"\t";
    // std::cout<<std::endl;
    // .push_back(i * (z2_max - z2_min) + z2_min);

    flux.set_source(EPAname);
    if (flux._leptonMasses.find(EPAname) != flux._leptonMasses.end()){

        std::cout << EPAname << std::endl;
        flux.set_thetac(thetac);
        flux.set_s(sqrts * sqrts);
    }

    // std::cout << sqrts << "\t" << Q20 << "\t" << Q2val << "\t" << thetac << std::endl;


    std::string model = modelstr.str();
    std::string widths = widthsstr.str();
    std::string evo = evostr.str();

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
    // std::vector<double> xB_values;
    // for(double i = 0; i <=1; i += 0.001) xB_values.push_back(i);
    // std::vector<double> y_values;
    // for(double i = 0; i <=1; i += 0.001) y_values.push_back(i);
    // std::vector<double> Q2_values = {3.0, 4.0, 9.0, 20.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0};
    // for(double i = 0; i <= 1; i += 0.001) Q2_values.push_back(i * (PhysicsCalculator::getQ2M() - PhysicsCalculator::getQ20()) + PhysicsCalculator::getQ20()); // Assuming Q2 ranges from 0 to Q2M
    // for (double i = 0; i <= 0)


    //set parameters for the integration
    int ndim = 2, ncomp = 13, nvec = 1, verbose = 0, last = 4, key = 13;
    double epsrel = 1e-6, epsabs = 1e-12;
    int flags = 2, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;

    // std::ofstream outFile("Fixed_xB.txt");
    // outFile << "," << "xB"
    //         << "," << "<dsig|+1,0>" << "," << "err[1]"
    //         << "," << "<dsig|+2,0>" << "," << "err[2]"
    //         << "," << "<ALL|0,0>" << "," << "err[3]"
    //         << "," << "<ALL|+1,0>" << "," << "err[4]"
    //         << "," << "<dsig|0,+1> " << "," << "err[5]"
    //         << "," << "<dsig|+1,-1>" << "," << "err[6]"
    //         << "," << "<dsig|+1,+1>" << "," << "err[7]"
    //         << "," << "<dsig|+2,+1> " << "," << "err[8]"
    //         << "," << "<dsig|+2,-1> " << "," << "err[9]"
    //         << "," << "<ALL|0,+1>" << "," << "err[10]"
    //         << "," << "<ALL|+1,-1>" << "," << "err[11]"
    //         << "," << "<ALL|+1,+1>" << "," << "err[12]"
    //         << "," << "AUzi" << "," << "err[13]"
    //         << "," << "AUcfqzi" << "," << "err[14]"
    //         << "," << "AUc2fqzi" << "," << "err[15]"
    //         << "," << "ALzi" << "," << "err[16]"
    //         << "," << "ALcfqzi" << "," << "err[17]"
    //         << "," << "BUcf12zi" << "," << "err[18]"
    //         << "," << "BUcfqmf12zi" << "," << "err[19]"
    //         << "," << "BUcfqpf12zi" << "," << "err[20]"
    //         << "," << "BUc2fqmf12zi" << "," << "err[21]"
    //         << "," << "BUc2fqpf12zi" << "," << "err[22]"
    //         << "," << "BLcf12zi" << "," << "err[23]"
    //         << "," << "BLcfqmf12zi" << "," << "err[24]"
    //         << "," << "BLcfqpf12zi" << "," << "err[25]"
    //         << "," << "time[s]"
    //         << "," << "neval" << "," << "fail"
    //         << std::endl;

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
    //     outFile << "," << xB_val
    //             << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
    //             << "," << res.ratio[2] << "," << res.ratioerror[2]
    //             << "," << res.ratio[3] << "," << res.ratioerror[3]
    //             << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
    //             << "," << res.ratio[5] << "," << res.ratioerror[5]
    //             << "," << res.ratio[6] << "," << res.ratioerror[6]
    //             << "," << res.ratio[7] << "," << res.ratioerror[7]
    //             << "," << res.ratio[8] << "," << res.ratioerror[8]
    //             << "," << res.ratio[9] << "," << res.ratioerror[9]
    //             << "," << res.ratio[10] << "," << res.ratioerror[10]
    //             << "," << res.ratio[11] << "," << res.ratioerror[11]
    //             << "," << res.ratio[12] << "," << res.ratioerror[12]
    //             << "," << res.integral[0] << "," << res.error[0]
    //             << "," << res.integral[1] << "," << res.error[1]
    //             << "," << res.integral[2] << "," << res.error[2]
    //             << "," << res.integral[3] << "," << res.error[3]
    //             << "," << res.integral[4] << "," << res.error[4]
    //             << "," << res.integral[5] << "," << res.error[5]
    //             << "," << res.integral[6] << "," << res.error[6]
    //             << "," << res.integral[7] << "," << res.error[7]
    //             << "," << res.integral[8] << "," << res.error[8]
    //             << "," << res.integral[9] << "," << res.error[9]
    //             << "," << res.integral[10] << "," << res.error[10]
    //             << "," << res.integral[11] << "," << res.error[11]
    //             << "," << res.integral[12] << "," << res.error[12]
    //             << "," << res.elapsed_seconds
    //             << "," << res.neval << "," << res.fail
    //             << std::endl;

    //     std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    // }

    // outFile.close();

    // std::ofstream outFile1("Fixed_y.txt");
    // outFile1 << "," << "y"
    //         << "," << "<dsig|+1,0>" << "," << "err[1]"
    //         << "," << "<dsig|+2,0>" << "," << "err[2]"
    //         << "," << "<ALL|0,0>" << "," << "err[3]"
    //         << "," << "<ALL|+1,0>" << "," << "err[4]"
    //         << "," << "<dsig|0,+1> " << "," << "err[5]"
    //         << "," << "<dsig|+1,-1>" << "," << "err[6]"
    //         << "," << "dsig|+1,+1>" << "," << "err[7]"
    //         << "," << "<dsig|+2,+1> " << "," << "err[8]"
    //         << "," << "<dsig|+2,-1> " << "," << "err[9]"
    //         << "," << "<ALL|0,+1>" << "," << "err[10]"
    //         << "," << "<ALL|+1,-1>" << "," << "err[11]"
    //         << "," << "<ALL|+1,+1>" << "," << "err[12]"
    //         << "," << "AUzi" << "," << "err[13]"
    //         << "," << "AUcfqzi" << "," << "err[14]"
    //         << "," << "AUc2fqzi" << "," << "err[15]"
    //         << "," << "ALzi" << "," << "err[16]"
    //         << "," << "ALcfqzi" << "," << "err[17]"
    //         << "," << "BUcf12zi" << "," << "err[18]"
    //         << "," << "BUcfqmf12zi" << "," << "err[19]"
    //         << "," << "BUcfqpf12zi" << "," << "err[20]"
    //         << "," << "BUc2fqmf12zi" << "," << "err[21]"
    //         << "," << "BUc2fqpf12zi" << "," << "err[22]"
    //         << "," << "BLcf12zi" << "," << "err[23]"
    //         << "," << "BLcfqmf12zi" << "," << "err[24]"
    //         << "," << "BLcfqpf12zi" << "," << "err[25]"
    //         << "," << "time[s]"
    //         << "," << "neval" << "," << "fail"
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

    //     outFile1 << "," << y_val
    //             << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
    //             << "," << res.ratio[2] << "," << res.ratioerror[2]
    //             << "," << res.ratio[3] << "," << res.ratioerror[3]
    //             << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
    //             << "," << res.ratio[5] << "," << res.ratioerror[5]
    //             << "," << res.ratio[6] << "," << res.ratioerror[6]
    //             << "," << res.ratio[7] << "," << res.ratioerror[7]
    //             << "," << res.ratio[8] << "," << res.ratioerror[8]
    //             << "," << res.ratio[9] << "," << res.ratioerror[9]
    //             << "," << res.ratio[10] << "," << res.ratioerror[10]
    //             << "," << res.ratio[11] << "," << res.ratioerror[11]
    //             << "," << res.ratio[12] << "," << res.ratioerror[12]
    //             << "," << res.integral[0] << "," << res.error[0]
    //             << "," << res.integral[1] << "," << res.error[1]
    //             << "," << res.integral[2] << "," << res.error[2]
    //             << "," << res.integral[3] << "," << res.error[3]
    //             << "," << res.integral[4] << "," << res.error[4]
    //             << "," << res.integral[5] << "," << res.error[5]
    //             << "," << res.integral[6] << "," << res.error[6]
    //             << "," << res.integral[7] << "," << res.error[7]
    //             << "," << res.integral[8] << "," << res.error[8]
    //             << "," << res.integral[9] << "," << res.error[9]
    //             << "," << res.integral[10] << "," << res.error[10]
    //             << "," << res.integral[11] << "," << res.error[11]
    //             << "," << res.integral[12] << "," << res.error[12]
    //             << "," << res.elapsed_seconds
    //             << "," << res.neval << "," << res.fail
    //             << std::endl;

    //     std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    // }

    // outFile1.close();
    std::ofstream outFile_UL("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) + "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" + LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UL.txt");
    std::ofstream outFile_UC("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) + "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" + LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UC.txt");
    std::ofstream outFile_sep("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) + "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" + LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_separated.txt");
    outFile_UL<< "z2"
            << "," << "<dsig|1;0>" << "," << "err[1]"
            << "," << "<dsig|2;0>" << "," << "err[2]"
            << "," << "<ALL|0;0>" << "," << "err[3]"
            << "," << "<ALL|1;0>" << "," << "err[4]"
            << "," << "<dsig|0;1>" << "," << "err[5]"
            << "," << "<dsig|1;-1>" << "," << "err[6]"
            << "," << "<dsig|1;1>" << "," << "err[7]"
            << "," << "<dsig|2;1>" << "," << "err[8]"
            << "," << "<dsig|2;-1>" << "," << "err[9]"
            << "," << "<ALL|0;1>" << "," << "err[10]"
            << "," << "<ALL|1;-1>" << "," << "err[11]"
            << "," << "<ALL|1;1>" << "," << "err[12]"
            << "," << "AUzi" << "," << "err[13]"
            << "," << "AUcfqzi" << "," << "err[14]"
            << "," << "AUc2fqzi" << "," << "err[15]"
            << "," << "ALzi" << "," << "err[16]"
            << "," << "ALcfqzi" << "," << "err[17]"
            << "," << "BUcf12zi" << "," << "err[18]"
            << "," << "BUcfqmf12zi" << "," << "err[19]"
            << "," << "BUcfqpf12zi" << "," << "err[20]"
            << "," << "BUc2fqmf12zi" << "," << "err[21]"
            << "," << "BUc2fqpf12zi" << "," << "err[22]"
            << "," << "BLcf12zi" << "," << "err[23]"
            << "," << "BLcfqmf12zi" << "," << "err[24]"
            << "," << "BLcfqpf12zi" << "," << "err[25]"
            << "," << "time[s]"
            << "," << "neval" << "," << "fail"
            << std::endl;

    outFile_UC << "z2"
            << "," << "<dsig|1;0>" << "," << "err[1]"
            << "," << "<dsig|2;0>" << "," << "err[2]"
            << "," << "<ALL|0;0>" << "," << "err[3]"
            << "," << "<ALL|1;0>" << "," << "err[4]"
            << "," << "<dsig|0;1>" << "," << "err[5]"
            << "," << "<dsig|1;-1>" << "," << "err[6]"
            << "," << "<dsig|1;1>" << "," << "err[7]"
            << "," << "<dsig|2;1>" << "," << "err[8]"
            << "," << "<dsig|2;-1>" << "," << "err[9]"
            << "," << "<ALL|0;1>" << "," << "err[10]"
            << "," << "<ALL|1;-1>" << "," << "err[11]"
            << "," << "<ALL|1;1>" << "," << "err[12]"
            << "," << "AUzi" << "," << "err[13]"
            << "," << "AUcfqzi" << "," << "err[14]"
            << "," << "AUc2fqzi" << "," << "err[15]"
            << "," << "ALzi" << "," << "err[16]"
            << "," << "ALcfqzi" << "," << "err[17]"
            << "," << "BUcf12zi" << "," << "err[18]"
            << "," << "BUcfqmf12zi" << "," << "err[19]"
            << "," << "BUcfqpf12zi" << "," << "err[20]"
            << "," << "BUc2fqmf12zi" << "," << "err[21]"
            << "," << "BUc2fqpf12zi" << "," << "err[22]"
            << "," << "BLcf12zi" << "," << "err[23]"
            << "," << "BLcfqmf12zi" << "," << "err[24]"
            << "," << "BLcfqpf12zi" << "," << "err[25]"
            << "," << "time[s]"
            << "," << "neval" << "," << "fail"
            << std::endl;

        outFile_sep << "z2"
            << "," << "<dsig|1;0>U" << "," << "err[1]"
            << "," << "<dsig|1;0>L" << "," << "err[2]"
            << "," << "<dsig|2;0>U" << "," << "err[3]"
            << "," << "<dsig|2;0>L" << "," << "err[4]"
            << "," << "<ALL|0;0>U" << "," << "err[5]"
            << "," << "<ALL|0;0>L" << "," << "err[6]"
            << "," << "<ALL|1;0>U" << "," << "err[7]"
            << "," << "<ALL|1;0>L" << "," << "err[8]"
            << "," << "<dsig|0;1>U" << "," << "err[9]"
            << "," << "<dsig|0;1>L" << "," << "err[10]"
            << "," << "<dsig|1;-1>U" << "," << "err[11]"
            << "," << "<dsig|1;-1>L" << "," << "err[12]"
            << "," << "<dsig|1;1>U" << "," << "err[13]"
            << "," << "<dsig|1;1>L" << "," << "err[14]"
            << "," << "<dsig|2;1>U" << "," << "err[15]"
            << "," << "<dsig|2;1>L" << "," << "err[16]"
            << "," << "<dsig|2;-1>U" << "," << "err[17]"
            << "," << "<dsig|2;-1>L" << "," << "err[18]"
            << "," << "<ALL|0;1>U" << "," << "err[19]"
            << "," << "<ALL|0;1>L" << "," << "err[20]"
            << "," << "<ALL|1;-1>U" << "," << "err[21]"
            << "," << "<ALL|1;-1>L" << "," << "err[22]"
            << "," << "<ALL|1;1>U" << "," << "err[23]"
            << "," << "<ALL|1;1>L" << "," << "err[24]"
            << "," << "time[s]"
            << "," << "neval" << "," << "fail"
            << std::endl;

    for (double z2 : z2_values) {

        fixed_Q2 = Q2val;
        std::cout << "Running scan for fixed Q2 = " << Q2val << "  z2 = " << z2 <<std::endl;

        // for
        void *USERDATA[] = {&sqrts, &Q20, &thetac, &Q2val, &z1_min, &z1_max, &z2};

        IntegrationResults col;
        col.maxeval = maxeval;
        col.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();

        Cuhre(ndim, ncomp, integrand_Collins, USERDATA, nvec,
            epsrel, epsabs, verbose | last,
            mineval, maxeval, key,
            statefile, spin,
            &col.nregions, &col.neval, &col.fail,
            col.integral, col.error, col.prob);

        double ratio_UL = col.integral[1] / col.integral[0] - col.integral[3] / col.integral[2];
        double ratio_UC = col.integral[1] / col.integral[0] - (col.integral[1] + col.integral[3]) / (col.integral[0] + col.integral[2]);
        double ratio_U = col.integral[1] / col.integral[0];
        double ratio_L = col.integral[3] / col.integral[2];

        std::cout << "ratio_UL = " << ratio_UL << " ratio_UC = " << ratio_UC << std::endl;
        std::cout << "ratio_U = " << ratio_U << " ratio_L = " << ratio_L << std::endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        col.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
        std::cout << " --> time=" << col.elapsed_seconds << "s, neval=" << col.neval << std::endl;


        auto t2 = std::chrono::high_resolution_clock::now();

        IntegrationResults res;

        res.maxeval = maxeval;
        res.nstart  = nstart;
        Vegas(ndim, ncomp, integrand_fixedQ2, USERDATA,
                nvec, epsrel, epsabs,
                flags, seed, mineval, maxeval,
                nstart, nincrease, nbatch,
                gridno, statefile, spin,
                &res.neval, &res.fail,
                res.integral, res.error, res.prob);


        auto t3 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t3 - t2).count();

        //Calculate ratios and ratio errors
        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];

            if (i >= 5){
                res.ratio[i] *= prefact * ratio_UL;
            }

            res.ratioerror[i] = abs(res.ratio[i]) * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2)
            );
        }

        outFile_UL<< "," << z2
                << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
                << "," << res.ratio[2] << "," << res.ratioerror[2]
                << "," << res.ratio[3] << "," << res.ratioerror[3]
                << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
                << "," << res.ratio[5] << "," << res.ratioerror[5]
                << "," << res.ratio[6] << "," << res.ratioerror[6]
                << "," << res.ratio[7] << "," << res.ratioerror[7]
                << "," << res.ratio[8] << "," << res.ratioerror[8]
                << "," << res.ratio[9] << "," << res.ratioerror[9]
                << "," << res.ratio[10] << "," << res.ratioerror[10]
                << "," << res.ratio[11] << "," << res.ratioerror[11]
                << "," << res.ratio[12] << "," << res.ratioerror[12]
                << "," << res.integral[0] << "," << res.error[0]
                << "," << res.integral[1] << "," << res.error[1]
                << "," << res.integral[2] << "," << res.error[2]
                << "," << res.integral[3] << "," << res.error[3]
                << "," << res.integral[4] << "," << res.error[4]
                << "," << res.integral[5] << "," << res.error[5]
                << "," << res.integral[6] << "," << res.error[6]
                << "," << res.integral[7] << "," << res.error[7]
                << "," << res.integral[8] << "," << res.error[8]
                << "," << res.integral[9] << "," << res.error[9]
                << "," << res.integral[10] << "," << res.error[10]
                << "," << res.integral[11] << "," << res.error[11]
                << "," << res.integral[12] << "," << res.error[12]
                << "," << res.elapsed_seconds
                << "," << res.neval << "," << res.fail
                << std::endl;

        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];

            if (i >= 5){
                res.ratio[i] *= prefact * ratio_UC;
            }

            res.ratioerror[i] = abs(res.ratio[i]) * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2)
            );
        }
       outFile_UC<< "," << z2
                << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
                << "," << res.ratio[2] << "," << res.ratioerror[2]
                << "," << res.ratio[3] << "," << res.ratioerror[3]
                << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
                << "," << res.ratio[5] << "," << res.ratioerror[5]
                << "," << res.ratio[6] << "," << res.ratioerror[6]
                << "," << res.ratio[7] << "," << res.ratioerror[7]
                << "," << res.ratio[8] << "," << res.ratioerror[8]
                << "," << res.ratio[9] << "," << res.ratioerror[9]
                << "," << res.ratio[10] << "," << res.ratioerror[10]
                << "," << res.ratio[11] << "," << res.ratioerror[11]
                << "," << res.ratio[12] << "," << res.ratioerror[12]
                << "," << res.integral[0] << "," << res.error[0]
                << "," << res.integral[1] << "," << res.error[1]
                << "," << res.integral[2] << "," << res.error[2]
                << "," << res.integral[3] << "," << res.error[3]
                << "," << res.integral[4] << "," << res.error[4]
                << "," << res.integral[5] << "," << res.error[5]
                << "," << res.integral[6] << "," << res.error[6]
                << "," << res.integral[7] << "," << res.error[7]
                << "," << res.integral[8] << "," << res.error[8]
                << "," << res.integral[9] << "," << res.error[9]
                << "," << res.integral[10] << "," << res.error[10]
                << "," << res.integral[11] << "," << res.error[11]
                << "," << res.integral[12] << "," << res.error[12]
                << "," << res.elapsed_seconds
                << "," << res.neval << "," << res.fail
                << std::endl;

        for (int i = 0; i < ncomp; ++i) {
            res.ratio[i] = res.integral[i] / res.integral[0];

            if (i >= 5){
                res.ratio[i] *= prefact;// * ratio_U;
            }

            res.ratioerror[i] = abs(res.ratio[i]) * std::sqrt(
                std::pow(res.error[i] / res.integral[i], 2) +
                std::pow(res.error[0] / res.integral[0], 2)
            );
        }
       outFile_sep<< "," << z2
                << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
                << "," << 2 * res.ratio[1] << "," << 2 * res.ratioerror[1]
                << "," << res.ratio[2] << "," << res.ratioerror[2]
                << "," << res.ratio[2] << "," << res.ratioerror[2]
                << "," << res.ratio[3] << "," << res.ratioerror[3]
                << "," << res.ratio[3] << "," << res.ratioerror[3]
                << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
                << "," << 2 * res.ratio[4] << "," << 2 * res.ratioerror[4]
                << "," << res.ratio[5] * ratio_U << "," << res.ratioerror[5]
                << "," << res.ratio[5] * ratio_L << "," << res.ratioerror[5]
                << "," << res.ratio[6] * ratio_U  << "," << res.ratioerror[6]
                << "," << res.ratio[6] * ratio_L  << "," << res.ratioerror[6]
                << "," << res.ratio[7] * ratio_U  << "," << res.ratioerror[7]
                << "," << res.ratio[7] * ratio_L  << "," << res.ratioerror[7]
                << "," << res.ratio[8] * ratio_U  << "," << res.ratioerror[8]
                << "," << res.ratio[8] * ratio_L  << "," << res.ratioerror[8]
                << "," << res.ratio[9] * ratio_U  << "," << res.ratioerror[9]
                << "," << res.ratio[9] * ratio_L  << "," << res.ratioerror[9]
                << "," << res.ratio[10] * ratio_U  << "," << res.ratioerror[10]
                << "," << res.ratio[10] * ratio_L  << "," << res.ratioerror[10]
                << "," << res.ratio[11] * ratio_U  << "," << res.ratioerror[11]
                << "," << res.ratio[11] * ratio_L  << "," << res.ratioerror[11]
                << "," << res.ratio[12] * ratio_U  << "," << res.ratioerror[12]
                << "," << res.ratio[12] * ratio_L  << "," << res.ratioerror[12]
                << "," << res.elapsed_seconds
                << "," << res.neval << "," << res.fail
                << std::endl;




        std::cout << " --> time=" << res.elapsed_seconds << "s, neval=" << res.neval << std::endl;
    }

    outFile_UL.close();
    outFile_UC.close();
    outFile_sep.close();

    return 0;
}
