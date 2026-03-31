#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <cuba.h>

#include "FragFunct.h"
#include "COL.h"
#include <rapidcsv.h>

#include "constants.h"

namespace
{
constexpr int kNumComponents = 4;
constexpr int kUserDataSqrts = 0;
constexpr int kUserDataQ20 = 1;
constexpr int kUserDataThetac = 2;
constexpr int kUserDataFixedValue = 3;
constexpr int kUserDataZ1Min = 4;
constexpr int kUserDataZ1Max = 5;
constexpr int kUserDataZ2 = 6;

double userdata_value(void* userdata, int index)
{
    return *static_cast<double **>(userdata)[index];
}

double safe_ratio(double numerator, double denominator)
{
    constexpr double kTolerance = 1e-60;
    if (std::abs(denominator) < kTolerance) {
        return 0.0;
    }
    return numerator / denominator;
}

double safe_ratio_error(double numerator, double numerator_error, double denominator, double denominator_error)
{
    constexpr double kTolerance = 1e-60;
    if (std::abs(numerator) < kTolerance || std::abs(denominator) < kTolerance) {
        return 0.0;
    }

    const double ratio = numerator / denominator;
    const double relative_numerator = numerator_error / numerator;
    const double relative_denominator = denominator_error / denominator;

    return std::abs(ratio) * std::sqrt(relative_numerator * relative_numerator +
                                       relative_denominator * relative_denominator);
}
}

std::vector<double> pars = {0.5316, -.7290, .8246, 3.2293, 0.8568, -.3539, 1.7524, 0.0, 0.2551};
std::vector<double> COL_z(13), COL_ppz1(13), COL_ppz2(13), COL_pmz1(13), COL_pmz2(13);
std::vector<double> FF(13), FF_ppz1(13), FF_ppz2(13), FF_pmz1(13), FF_pmz2(13),\
    FF_ppz1_fixedQ2(13), FF_ppz2_fixedQ2(13), FF_pmz1_fixedQ2(13), FF_pmz2_fixedQ2(13),\
    FF_evo(13);

double f[13], h1[13]; //for HOPPET and transversity

const double charges[13] = {-2./3., 1./3., -2./3., 1./3., -2./3., 1./3., 0, -1./3., 2./3., -1./3., 2./3., -1./3., 2./3.};

double  numU = 0.0, denU = 0.0,\
        numL = 0.0, denL = 0.0,\
        sqrts = 0.0, thetac = 0.0, Q20 = 0.0,\
        z1 = 0.0, z2 = 0.0, pperp2 = 0.12;

double MC2 = pars.back();
double pperp2Col = MC2 * pperp2 / (MC2 + pperp2);
double prefact = (pi * EulerConst / 2) * (pow(pperp2Col, 3) / (pperp2 * pperp2 * MC2));

constexpr int kPionHadron = 1;
constexpr int kPositiveCharge = 1;
constexpr int kNegativeCharge = -1;
int charge 1, hadron = 1;

FRAG::FF myFF("DEHSS","NLO");
COL::COLLINS myCol;

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



void Collins_FF(int & hadron, double & z1, double & z2, double & hard_scale_sq)
{
    myFF.FF_eval(hadron, charge, z1, hard_scale_sq);
    for (int i = 0; i < FF_ppz1.size(); ++i) FF_ppz1[i] = myFF.theFF[i] / z1;

    myFF.FF_eval(hadron, charge, z2, hard_scale_sq);
    for (int i = 0; i < FF_ppz2.size(); ++i) FF_ppz2[i] = myFF.theFF[i] / z2;

    charge *= -1; //to call pi- FFs

    myFF.FF_eval(hadron, charge, z1, hard_scale_sq);
    for (int i = 0; i < FF_pmz1.size(); ++i) FF_pmz1[i] = myFF.theFF[i] / z1;

    myFF.FF_eval(hadron, charge, z2, hard_scale_sq);
    for (int i = 0; i < FF_pmz2.size(); ++i) FF_pmz2[i] = myFF.theFF[i] / z2;
}


void Collins_epem_loop(const std::vector<double>& COL_ppz1_in,
                       const std::vector<double>& COL_ppz2_in,
                       const std::vector<double>& COL_pmz1_in,
                       const std::vector<double>& COL_pmz2_in,
                       const std::vector<double>& FF_ppz1_in,
                       const std::vector<double>& FF_ppz2_in,
                       const std::vector<double>& FF_pmz1_in,
                       const std::vector<double>& FF_pmz2_in)
{
    denU = 0.0;
    numU = 0.0;
    denL = 0.0;
    numL = 0.0;
    for (int i = 3; i <= 9; ++i) {
        if (i == 3) { // sb
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i + 6] + COL_pmz1_in[i] * COL_ppz2_in[i + 6]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i + 6] + FF_pmz1_in[i] * FF_ppz2_in[i + 6]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i + 6] + COL_pmz1_in[i] * COL_pmz2_in[i + 6]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i + 6] + FF_pmz1_in[i] * FF_pmz2_in[i + 6]);
        }

        if (i == 4) { // ub
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i + 4] + COL_pmz1_in[i] * COL_ppz2_in[i + 4]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i + 4] + FF_pmz1_in[i] * FF_ppz2_in[i + 4]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i + 4] + COL_pmz1_in[i] * COL_pmz2_in[i + 4]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i + 4] + FF_pmz1_in[i] * FF_pmz2_in[i + 4]);
        }

        if (i == 5) { // db
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i + 2] + COL_pmz1_in[i] * COL_ppz2_in[i + 2]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i + 2] + FF_pmz1_in[i] * FF_ppz2_in[i + 2]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i + 2] + COL_pmz1_in[i] * COL_pmz2_in[i + 2]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i + 2] + FF_pmz1_in[i] * FF_pmz2_in[i + 2]);
        }

        if (i == 6) { // g
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i] + COL_pmz1_in[i] * COL_ppz2_in[i]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i] + FF_pmz1_in[i] * FF_ppz2_in[i]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i] + COL_pmz1_in[i] * COL_pmz2_in[i]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i] + FF_pmz1_in[i] * FF_pmz2_in[i]);
        }

        if (i == 7) { // d
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i - 2] + COL_pmz1_in[i] * COL_ppz2_in[i - 2]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i - 2] + FF_pmz1_in[i] * FF_ppz2_in[i - 2]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i - 2] + COL_pmz1_in[i] * COL_pmz2_in[i - 2]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i - 2] + FF_pmz1_in[i] * FF_pmz2_in[i - 2]);
        }

        if (i == 8) { // u
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i - 4] + COL_pmz1_in[i] * COL_ppz2_in[i - 4]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i - 4] + FF_pmz1_in[i] * FF_ppz2_in[i - 4]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i - 4] + COL_pmz1_in[i] * COL_pmz2_in[i - 4]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i - 4] + FF_pmz1_in[i] * FF_pmz2_in[i - 4]);
        }

        if (i == 9) { // s
            numU += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_pmz2_in[i - 6] + COL_pmz1_in[i] * COL_ppz2_in[i - 6]);
            denU += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_pmz2_in[i - 6] + FF_pmz1_in[i] * FF_ppz2_in[i - 6]);
            numL += std::pow(charges[i], 4) * (COL_ppz1_in[i] * COL_ppz2_in[i - 6] + COL_pmz1_in[i] * COL_pmz2_in[i - 6]);
            denL += std::pow(charges[i], 4) * (FF_ppz1_in[i] * FF_ppz2_in[i - 6] + FF_pmz1_in[i] * FF_pmz2_in[i - 6]);
        }
    }
}

std::vector<double> integration_point_to_vector(const int* ndim, const double x[])
{
    return std::vector<double>(x, x + *ndim);
}



class PhysicsCalculator
{
private:
    static constexpr double alphaem = 1.0 / 137.036;
    static constexpr double me = 5.11e-4;
    static constexpr double PI = 3.14159265358979323846;

    double etaq = 0.0;
    double etabarq = 0.0;
    double kT = 0.0;
    double csi1 = 0.0;
    double csi2 = 0.0;
    double Q2min_csi1 = 0.0;
    double Q2max_csi1 = 0.0;
    double Q2min_csi2 = 0.0;
    double Q2max_csi2 = 0.0;
    double fgl1 = 0.0;
    double DLfgl1 = 0.0;
    double fgl2 = 0.0;
    double DLfgl2 = 0.0;
    double jacob = 0.0;
    bool valid_kinematics = false;

public:
    double SqrtS = 0.0;
    double S = 0.0;
    double kT0 = 0.0;
    double kTM = 0.0;
    double thetac = 0.0;

    void set_SqrtS(double val)
    {
        SqrtS = val;
        S = SqrtS * SqrtS;
        kTM = SqrtS / 2.0;
    }

    void set_kT0(double val)
    {
        kT0 = val;
    }

    void setExperiment(double sqrts_val, double kT0_val, double thetac_val)
    {
        set_SqrtS(sqrts_val);
        set_kT0(kT0_val);
        thetac = thetac_val;
    }

    void setVariables(double x1_val, double x2_val, double kT_val)
    {
        valid_kinematics = false;
        kT = kT_val;
        jacob = 0.0;
        if (SqrtS <= 0.0 || thetac <= 0.0 || kT <= 0.0 || kT < kT0 || kT >= kTM) {
            return;
        }

        const double sqrt_arg = S / (kT * kT) - 4.0;
        if (sqrt_arg <= 0.0) {
            return;
        }

        const double sqrt_term = std::sqrt(sqrt_arg);
        const double etaqmax = std::log(0.5 * (SqrtS / kT + sqrt_term));
        const double etaqmin = std::log(0.5 * (SqrtS / kT - sqrt_term));
        etaq = (etaqmax - etaqmin) * x1_val + etaqmin;

        const double etabarqmax_arg = SqrtS / kT - std::exp(etaq);
        const double etabarqmin_arg = SqrtS / kT - std::exp(-etaq);
        if (etabarqmax_arg <= 0.0 || etabarqmin_arg <= 0.0) {
            return;
        }

        const double etabarqmax = std::log(etabarqmax_arg);
        const double etabarqmin = -std::log(etabarqmin_arg);
        jacob = (etaqmax - etaqmin) * (etabarqmax - etabarqmin);
        etabarq = (etabarqmax - etabarqmin) * x2_val + etabarqmin;

        csi1 = (kT / SqrtS) * (std::exp(etaq) + std::exp(etabarq));
        csi2 = (kT / SqrtS) * (std::exp(-etaq) + std::exp(-etabarq));
        if (!(csi1 > 0.0 && csi1 < 1.0 && csi2 > 0.0 && csi2 < 1.0)) {
            return;
        }

        Q2min_csi1 = me * me * csi1 * csi1 / (1.0 - csi1);
        Q2max_csi1 = (S / 4.0) * thetac * thetac * (1.0 - csi1) + Q2min_csi1;
        Q2min_csi2 = me * me * csi2 * csi2 / (1.0 - csi2);
        Q2max_csi2 = (S / 4.0) * thetac * thetac * (1.0 - csi2) + Q2min_csi2;
        if (!(Q2min_csi1 > 0.0 && Q2max_csi1 > Q2min_csi1 &&
              Q2min_csi2 > 0.0 && Q2max_csi2 > Q2min_csi2)) {
            return;
        }

        const double log_ratio_csi1 = std::log(Q2max_csi1 / Q2min_csi1);
        const double me2_term_csi1 = 2.0 * me * me * csi1 * (1.0 / Q2max_csi1 - 1.0 / Q2min_csi1);
        const double log_ratio_csi2 = std::log(Q2max_csi2 / Q2min_csi2);
        const double me2_term_csi2 = 2.0 * me * me * csi2 * (1.0 / Q2max_csi2 - 1.0 / Q2min_csi2);

        fgl1 = (alphaem / (2.0 * PI)) * ((1.0 + (1.0 - csi1) * (1.0 - csi1)) / csi1 * log_ratio_csi1 + me2_term_csi1);
        DLfgl1 = (alphaem / (2.0 * PI)) * ((1.0 - (1.0 - csi1) * (1.0 - csi1)) / csi1 * log_ratio_csi1 +
                                           2.0 * me * me * csi1 * csi1 * (1.0 / Q2max_csi1 - 1.0 / Q2min_csi1));
        fgl2 = (alphaem / (2.0 * PI)) * ((1.0 + (1.0 - csi2) * (1.0 - csi2)) / csi2 * log_ratio_csi2 + me2_term_csi2);
        DLfgl2 = (alphaem / (2.0 * PI)) * ((1.0 - (1.0 - csi2) * (1.0 - csi2)) / csi2 * log_ratio_csi2 +
                                           2.0 * me * me * csi2 * csi2 * (1.0 / Q2max_csi2 - 1.0 / Q2min_csi2));

        valid_kinematics = std::isfinite(fgl1) && std::isfinite(DLfgl1) &&
                           std::isfinite(fgl2) && std::isfinite(DLfgl2);
        std::cout << "Set variables: etaq = " << etaq << ", etabarq = " << etabarq << ", kT = " << kT
                  << ", csi1 = " << csi1 << ", csi2 = " << csi2
                  << ", Q2min_csi1 = " << Q2min_csi1 << ", Q2max_csi1 = " << Q2max_csi1
                  << ", Q2min_csi2 = " << Q2min_csi2 << ", Q2max_csi2 = " << Q2max_csi2
                  << ", fgl1 = " << fgl1 << ", DLfgl1 = " << DLfgl1
                  << ", fgl2 = " << fgl2 << ", DLfgl2 = " << DLfgl2
                  << ", jacob = " << jacob
                  << ", valid_kinematics = " << valid_kinematics
                  << std::endl;
    }

    bool isValid() const { return valid_kinematics; }
    double K() const { return valid_kinematics ? 1.0 / (S * kT * kT * (1.0 + std::cosh(etaq - etabarq))) : 0.0; }
    double AU() const { return valid_kinematics ? std::cosh(etaq - etabarq) * fgl1 * fgl2 : 0.0; }
    double BU() const { return valid_kinematics ? fgl1 * fgl2 / 4.0 : 0.0; }
    double AL() const { return valid_kinematics ? std::cosh(etaq - etabarq) * DLfgl1 * DLfgl2 : 0.0; }
    double BL() const { return valid_kinematics ? DLfgl1 * DLfgl2 / 4.0 : 0.0; }

    double getEtaq() const { return etaq; }
    double getEtabarq() const { return etabarq; }
    double getKT() const { return kT; }
    double getCsi1() const { return csi1; }
    double getCsi2() const { return csi2; }
    double getQ2minCsi1() const { return Q2min_csi1; }
    double getQ2maxCsi1() const { return Q2max_csi1; }
    double getQ2minCsi2() const { return Q2min_csi2; }
    double getQ2maxCsi2() const { return Q2max_csi2; }
    double getJacob() const { return jacob; }
};

void fill_azimuthal_results(const PhysicsCalculator& pc, double norm, double results[])
{
    std::fill(results, results + kNumComponents, 0.0);
    const double total_norm = norm * pc.getJacob();
    results[0]  = pc.AU() * total_norm;
    results[1]  = pc.BU() * total_norm;
    results[2]  = pc.AL() * total_norm;
    results[3]  = pc.BL() * total_norm;
}

bool cut(const PhysicsCalculator& pscut)
{
    return pscut.isValid();
}

void Values(const std::vector<double> &x, double results[])
{
    std::fill(results, results + kNumComponents, 0.0);
    if (x.size() < 3) {
        return;
    }

    PhysicsCalculator pc;
    pc.setExperiment(sqrts, Q20, thetac);
    pc.setVariables(x[0], x[1], x[2]);

    if (cut(pc)) {
        fill_azimuthal_results(pc, 1.0, results);
    }
}


void ValuesfixedkT(const std::vector<double> &x, double results[], void *userdata)
{
    // x[0] PROP ETAQ, x[1] PROP ETABARQ; fixed kT is passed via userdata
    const double sqrts = userdata_value(userdata, kUserDataSqrts);
    const double kT0 = userdata_value(userdata, kUserDataQ20);
    const double thetac = userdata_value(userdata, kUserDataThetac);
    const double fixed_kT = userdata_value(userdata, kUserDataFixedValue);

    PhysicsCalculator pc;
    pc.setExperiment(sqrts, kT0, thetac);

    pc.setVariables(x[0], x[1], fixed_kT);

    if (cut(pc))
    {
        fill_azimuthal_results(pc, pc.K(), results);
    } 
    else
    {
        std::fill(results, results + kNumComponents, 0.0);
    }
}

int integrand_fixedkT(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata)
{
    (void)ncomp;
    
    ValuesfixedkT(integration_point_to_vector(ndim, x), f, userdata);

    return 0;
}

int integrand_Collins(const int *ndim, const double x[], const int *ncomp, double ff[], void *userdata)
{
    (void)ndim;
    (void)ncomp;

#define f0 ff[0]
#define f1 ff[1]
#define f2 ff[2]
#define f3 ff[3]

    std::fill(ff, ff + kNumComponents, 0.0);

    const double hard_scale = userdata_value(userdata, kUserDataFixedValue);
    const double hard_scale_sq = hard_scale * hard_scale;
    const double z1_min = userdata_value(userdata, kUserDataZ1Min);
    const double z1_max = userdata_value(userdata, kUserDataZ1Max);
    const double z2 = userdata_value(userdata, kUserDataZ2);

    const double z1 = z1_min + x[0] * (z1_max - z1_min);
    double evolved_collins[13];


    Collins_FF(hadron, z1, z2, hard_scale_sq);

    if (myCol.evo == "DGLAP" || myCol.evo == "none") {
        myCol.eval(z1, hard_scale_sq, charge, FF_ppz1, pars);
        for (int i = 0; i < COL_ppz1.size(); ++i) COL_ppz1[i] = myCol.COL_z[i];

        myCol.eval(z2, hard_scale_sq, charge, FF_ppz2, pars);
        for (int i = 0; i < COL_ppz2.size(); ++i) COL_ppz2[i] = myCol.COL_z[i];

        charge *= -1; //to call pi- Collins FFs

        myCol.eval(z1, hard_scale_sq, charge, FF_pmz1, pars);
        for (int i = 0; i < COL_pmz1.size(); ++i) COL_pmz1[i] = myCol.COL_z[i];

        myCol.eval(z2, hard_scale_sq, charge, FF_pmz2, pars);
        for (int i = 0; i < COL_pmz2.size(); ++i) COL_pmz2[i] = myCol.COL_z[i];
    }

    if (myCol.evo == "CT3") {
        hoppetEvalcf(z1, hard_scale, evolved_collins);
        for (int i = 0; i < COL_ppz1.size(); ++i) COL_ppz1[i] = evolved_collins[i] / z1;

        hoppetEvalcf(z2, hard_scale, evolved_collins);
        for (int i = 0; i < COL_ppz2.size(); ++i) COL_ppz2[i] = evolved_collins[i] / z2;

        hoppetEvalcff(z1, hard_scale, evolved_collins);
        for (int i = 0; i < COL_pmz1.size(); ++i) COL_pmz1[i] = evolved_collins[i] / z1;

        hoppetEvalcff(z2, hard_scale, evolved_collins);
        for (int i = 0; i < COL_pmz2.size(); ++i) COL_pmz2[i] = evolved_collins[i] / z2;
    }

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
    double integral[kNumComponents];
    double error[kNumComponents];
    double prob[kNumComponents];
    double ratio[kNumComponents];
    double ratioerror[kNumComponents];
    double elapsed_seconds;
};


int main(int argc, char *argv[])
{
    (void)argc;

    const int maxeval = static_cast<int>(1e7);
    const int nstart  = static_cast<int>(1e6);

    stringstream modelstr, widthsstr, evostr, MCnamestring;
    modelstr << argv[2];
    widthsstr << argv[4];
    evostr << argv[6];

    double z1_min = atof(argv[8]);
    double z1_max = atof(argv[9]);
    double z2_min = atof(argv[11]);
    double z2_max = atof(argv[12]);
    sqrts = atof(argv[14]);
    Q20 = atof(argv[16]); // kept for input compatibility, interpreted here as kT_min
    double fixed_kT = atof(argv[18]);
    thetac = atof(argv[20]);
    MCnamestring << argv[22];

    const std::string MCname = MCnamestring.str();
    CSV->Load(MCname);

    std::vector<double> z2_values;
    for (double value = z2_min; value <= z2_max + 0.05; value += 0.05) {
        z2_values.push_back(value);
    }

    const std::string model = modelstr.str();
    const std::string widths = widthsstr.str();
    const std::string evo = evostr.str();

    myFF.use_LHAPDF(true);
    cout << myFF.FFset << "\t" << myFF.FForder << endl;
    if (myFF.useLHAPDF) {
        if (myFF.FFset == "NNFF10") {
            if (myFF.FForder == "LO") myFF.set_LHAPDF_FFset(NNFF10setsLO);
            if (myFF.FForder == "NLO") myFF.set_LHAPDF_FFset(NNFF10setsNLO);
        }
        if (myFF.FFset == "DEHSS") myFF.set_LHAPDF_FFset(DEHSSsetsNLO);
    }

    myCol.set_model(model);
    myCol.set_widths(widths);
    myCol.set_evolution(evo);

    int ndim = 2, ncomp = kNumComponents, nvec = 1, verbose = 0, last = 4, key = 13;
    double epsrel = 1e-6, epsabs = 1e-12;
    int flags = 2, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;

    std::ofstream outFile_UL("Fixed_kT_" + LHAPDF::to_str(fixed_kT) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UL.txt");
    std::ofstream outFile_UC("Fixed_kT_" + LHAPDF::to_str(fixed_kT) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UC.txt");
    std::ofstream outFile_sep("Fixed_kT_" + LHAPDF::to_str(fixed_kT) + "_Vs_" + LHAPDF::to_str(sqrts) +
                              "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                              LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_separated.txt");

    outFile_UL << "z2"
               << ",AL_over_AU" << ",err_AL_over_AU"
               << ",BU_over_AU_times_Collins_UL" << ",err_BU_over_AU_times_Collins_UL"
               << ",BL_over_AU_times_Collins_UL" << ",err_BL_over_AU_times_Collins_UL"
               << ",AU" << ",err_AU"
               << ",BU" << ",err_BU"
               << ",AL" << ",err_AL"
               << ",BL" << ",err_BL"
               << ",Collins_UL" << ",Collins_U" << ",Collins_L"
               << ",time[s]" << ",neval" << ",fail"
               << std::endl;

    outFile_UC << "z2"
               << ",AL_over_AU" << ",err_AL_over_AU"
               << ",BU_over_AU_times_Collins_UC" << ",err_BU_over_AU_times_Collins_UC"
               << ",BL_over_AU_times_Collins_UC" << ",err_BL_over_AU_times_Collins_UC"
               << ",AU" << ",err_AU"
               << ",BU" << ",err_BU"
               << ",AL" << ",err_AL"
               << ",BL" << ",err_BL"
               << ",Collins_UC" << ",Collins_U" << ",Collins_L"
               << ",time[s]" << ",neval" << ",fail"
               << std::endl;

    outFile_sep << "z2"
                << ",AL_over_AU" << ",err_AL_over_AU"
                << ",BU_over_AU_times_Collins_U" << ",err_BU_over_AU_times_Collins_U"
                << ",BU_over_AU_times_Collins_L" << ",err_BU_over_AU_times_Collins_L"
                << ",BL_over_AU_times_Collins_U" << ",err_BL_over_AU_times_Collins_U"
                << ",BL_over_AU_times_Collins_L" << ",err_BL_over_AU_times_Collins_L"
                << ",Collins_U" << ",Collins_L"
                << ",time[s]" << ",neval" << ",fail"
                << std::endl;

    for (double z2_value : z2_values) {
        std::cout << "Running scan for fixed kT = " << fixed_kT << "  z2 = " << z2_value << std::endl;

        void *USERDATA[] = {&sqrts, &Q20, &thetac, &fixed_kT, &z1_min, &z1_max, &z2_value};

        IntegrationResults col{};
        col.maxeval = maxeval;
        col.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();
        Cuhre(ndim, ncomp, integrand_Collins, USERDATA, nvec,
              epsrel, epsabs, verbose | last,
              mineval, maxeval, key,
              statefile, spin,
              &col.nregions, &col.neval, &col.fail,
              col.integral, col.error, col.prob);
        auto t1 = std::chrono::high_resolution_clock::now();
        col.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();

        const double ratio_U = safe_ratio(col.integral[1], col.integral[0]);
        const double ratio_L = safe_ratio(col.integral[3], col.integral[2]);
        const double ratio_UL = ratio_U - ratio_L;
        const double ratio_UC = ratio_U - safe_ratio(col.integral[1] + col.integral[3], col.integral[0] + col.integral[2]);

        std::cout << "Collins ratios: UL = " << ratio_UL
                  << "  UC = " << ratio_UC
                  << "  U = " << ratio_U
                  << "  L = " << ratio_L << std::endl;
        std::cout << " --> Collins integration time=" << col.elapsed_seconds
                  << "s, neval=" << col.neval << std::endl;

        IntegrationResults res{};
        res.maxeval = maxeval;
        res.nstart  = nstart;

        auto t2 = std::chrono::high_resolution_clock::now();
        Vegas(ndim, ncomp, integrand_fixedkT, USERDATA,
              nvec, epsrel, epsabs,
              flags, seed, mineval, maxeval,
              nstart, nincrease, nbatch,
              gridno, statefile, spin,
              &res.neval, &res.fail,
              res.integral, res.error, res.prob);
        auto t3 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t3 - t2).count();

        const double al_over_au = safe_ratio(res.integral[2], res.integral[0]);
        const double al_over_au_error = safe_ratio_error(res.integral[2], res.error[2], res.integral[0], res.error[0]);
        const double bu_over_au = safe_ratio(res.integral[1], res.integral[0]);
        const double bu_over_au_error = safe_ratio_error(res.integral[1], res.error[1], res.integral[0], res.error[0]);
        const double bl_over_au = safe_ratio(res.integral[3], res.integral[0]);
        const double bl_over_au_error = safe_ratio_error(res.integral[3], res.error[3], res.integral[0], res.error[0]);

        const double scale_ul = prefact * ratio_UL;
        const double scale_uc = prefact * ratio_UC;
        const double scale_u = prefact * ratio_U;
        const double scale_l = prefact * ratio_L;

        outFile_UL << "," << z2_value
                   << "," << al_over_au << "," << al_over_au_error
                   << "," << bu_over_au * scale_ul << "," << bu_over_au_error * std::abs(scale_ul)
                   << "," << bl_over_au * scale_ul << "," << bl_over_au_error * std::abs(scale_ul)
                   << "," << res.integral[0] << "," << res.error[0]
                   << "," << res.integral[1] << "," << res.error[1]
                   << "," << res.integral[2] << "," << res.error[2]
                   << "," << res.integral[3] << "," << res.error[3]
                   << "," << ratio_UL << "," << ratio_U << "," << ratio_L
                   << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                   << std::endl;

        outFile_UC << "," << z2_value
                   << "," << al_over_au << "," << al_over_au_error
                   << "," << bu_over_au * scale_uc << "," << bu_over_au_error * std::abs(scale_uc)
                   << "," << bl_over_au * scale_uc << "," << bl_over_au_error * std::abs(scale_uc)
                   << "," << res.integral[0] << "," << res.error[0]
                   << "," << res.integral[1] << "," << res.error[1]
                   << "," << res.integral[2] << "," << res.error[2]
                   << "," << res.integral[3] << "," << res.error[3]
                   << "," << ratio_UC << "," << ratio_U << "," << ratio_L
                   << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                   << std::endl;

        outFile_sep << "," << z2_value
                    << "," << al_over_au << "," << al_over_au_error
                    << "," << bu_over_au * scale_u << "," << bu_over_au_error * std::abs(scale_u)
                    << "," << bu_over_au * scale_l << "," << bu_over_au_error * std::abs(scale_l)
                    << "," << bl_over_au * scale_u << "," << bl_over_au_error * std::abs(scale_u)
                    << "," << bl_over_au * scale_l << "," << bl_over_au_error * std::abs(scale_l)
                    << "," << ratio_U << "," << ratio_L
                    << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                    << std::endl;

        std::cout << " --> kT integration time=" << res.elapsed_seconds
                  << "s, neval=" << res.neval << std::endl;
    }

    outFile_UL.close();
    outFile_UC.close();
    outFile_sep.close();

    return 0;
}
