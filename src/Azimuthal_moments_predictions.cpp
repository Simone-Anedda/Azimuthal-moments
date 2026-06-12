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
#include "PhysicsCalculator.h"
#include "photon_flux.h"
#include <rapidcsv.h>

#include "constants.h"

namespace
{
constexpr int kNumComponents = 13;
constexpr int kCollinsChannels = 8;
constexpr int kUserDataSqrts = 0;
constexpr int kUserDataQ20 = 1;
constexpr int kUserDataThetac = 2;
constexpr int kUserDataFixedValue = 3;
constexpr int kUserDataZ1Min = 4;
constexpr int kUserDataZ1Max = 5;
constexpr int kUserDataZ2 = 6;
// index 7 is reserved for std::vector<double>** pars


double userdata_value(void* userdata, int index)
{
    return *static_cast<double **>(userdata)[index];
}

double safe_ratio(double numerator, double denominator)
{
    return numerator / denominator;
}

double safe_ratio_error(double numerator, double numerator_error, double denominator, double denominator_error)
{
    double ratio = numerator / denominator;
    double relative_numerator = numerator_error / numerator;
    double relative_denominator = denominator_error / denominator;

    return std::abs(ratio) * std::sqrt(relative_numerator * relative_numerator +
                                       relative_denominator * relative_denominator);
}
}

std::vector<double> mean_pars = {0.5316, -.7290, .8246, 3.2293, 0.8568, -.3539, 1.7524, 0.0, 0.2551};
std::vector<double> COL_z(13), COL_ppz1(13), COL_ppz2(13), COL_pmz1(13), COL_pmz2(13);
std::vector<double> FF(13), FF_ppz1(13), FF_ppz2(13), FF_pmz1(13), FF_pmz2(13),\
    FF_ppz1_fixedQ2(13), FF_ppz2_fixedQ2(13), FF_pmz1_fixedQ2(13), FF_pmz2_fixedQ2(13),\
    FF_evo(13);

const double charges[13] = {-2./3., 1./3., -2./3., 1./3., -2./3., 1./3., 0, -1./3., 2./3., -1./3., 2./3., -1./3., 2./3.};

double sqrts = 0.0, thetac = 0.0, Q20 = 0.0, pperp2 = 0.12;
double Coll_piPpiM = 0.0, Coll_piMpiP = 0.0, Coll_piPpiP = 0.0, Coll_piMpiM = 0.0;
double Unp_piPpiM = 0.0, Unp_piMpiP = 0.0, Unp_piPpiP = 0.0, Unp_piMpiM = 0.0;
double mean_MC2 = mean_pars.back();
double mean_pperp2Col = mean_MC2 * pperp2 / (mean_MC2 + pperp2);
double mean_prefact = (pi * EulerConst / 2) * (pow(mean_pperp2Col, 3) / (pperp2 * pperp2 * mean_MC2));

int charge = 1, hadron = 1;

FRAG::FF myFF("DEHSS","NLO");
COL::COLLINS myCol;
EPA::EPA_flux * flux = new EPA::EPA_flux();

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

double get_prefact(double & MC2){

    double pperp2Col = MC2 * pperp2 / (MC2 + pperp2);
    return (pi * EulerConst / 2) * (pow(pperp2Col, 3) / (pperp2 * pperp2 * MC2));

}


void Collins_FF(int & hadron, int & charge, double & z1, double & z2, double & Q2)
{
    myFF.FF_eval(hadron, charge, z1, Q2);
    for (int i = 0; i < FF_ppz1.size(); ++i) FF_ppz1[i] = myFF.theFF[i] / z1;

    myFF.FF_eval(hadron, charge, z2, Q2);
    for (int i = 0; i < FF_ppz2.size(); ++i) FF_ppz2[i] = myFF.theFF[i] / z2;

    charge *= -1; //to call pi- FFs

    myFF.FF_eval(hadron, charge, z1, Q2);
    for (int i = 0; i < FF_pmz1.size(); ++i) FF_pmz1[i] = myFF.theFF[i] / z1;

    myFF.FF_eval(hadron, charge, z2, Q2);
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
    Coll_piMpiM = 0.0;
    Coll_piPpiM = 0.0;
    Coll_piMpiP = 0.0;
    Coll_piPpiP = 0.0;
    Unp_piMpiM = 0.0;
    Unp_piPpiM = 0.0;
    Unp_piMpiP = 0.0;
    Unp_piPpiP = 0.0;

    for (int i = 3; i <= 9; ++i) {
        if (i == 3) { // sb
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i + 6];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i + 6];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i + 6];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i + 6];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i + 6];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i + 6];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i + 6];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i + 6];
        }

        if (i == 4) { // ub
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i + 4];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i + 4];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i + 4];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i + 4];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i + 4];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i + 4];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i + 4];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i + 4];
        }

        if (i == 5) { // db
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i + 2];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i + 2];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i + 2];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i + 2];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i + 2];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i + 2];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i + 2];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i + 2];
        }

        if (i == 6) { // g
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i];
        }

        if (i == 7) { // d
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i - 2];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i - 2];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i - 2];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i - 2];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i - 2];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i - 2];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i - 2];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i - 2];
        }

        if (i == 8) { // u
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i - 4];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i - 4];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i - 4];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i - 4];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i - 4];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i - 4];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i - 4];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i - 4];
        }

        if (i == 9) { // s
            Coll_piPpiM += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_pmz2_in[i - 6];
            Coll_piMpiP += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_ppz2_in[i - 6];
            Coll_piPpiP += std::pow(charges[i], 4) * COL_ppz1_in[i] * COL_ppz2_in[i - 6];
            Coll_piMpiM += std::pow(charges[i], 4) * COL_pmz1_in[i] * COL_pmz2_in[i - 6];
            Unp_piPpiM += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_pmz2_in[i - 6];
            Unp_piMpiP += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_ppz2_in[i - 6];
            Unp_piPpiP += std::pow(charges[i], 4) * FF_ppz1_in[i] * FF_ppz2_in[i - 6];
            Unp_piMpiM += std::pow(charges[i], 4) * FF_pmz1_in[i] * FF_pmz2_in[i - 6];
        }
    }
}


void fill_azimuthal_results(const YStarYKinematics& kin, double norm, double results[])
{
    std::fill(results, results + kNumComponents, 0.0);
    results[0]  = kin.AUzi * norm;
    results[1]  = kin.AUcfqzi * norm;
    results[2]  = kin.AUc2fqzi * norm;
    results[3]  = kin.ALzi * norm;
    results[4]  = kin.ALcfqzi * norm;
    results[5]  = kin.BUcf12zi * norm;
    results[6]  = kin.BUcfqmf12zi * norm;
    results[7]  = kin.BUcfqpf12zi * norm;
    results[8]  = kin.BUc2fqmf12zi * norm;
    results[9]  = kin.BUc2fqpf12zi * norm;
    results[10] = kin.BLcf12zi * norm;
    results[11] = kin.BLcfqmf12zi * norm;
    results[12] = kin.BLcfqpf12zi * norm;
}

bool cut(const YStarYKinematics& kin)
{
    return kin.valid;
}


int integrand(const int *ndim, const double x[], const int *ncomp, double ff[], void *userdata)
{
    (void)ndim;
    (void)ncomp;

#define f0  ff[0]  // Coll_piPpiM
#define f1  ff[1]  // Coll_piMpiP
#define f2  ff[2]  // Coll_piPpiP
#define f3  ff[3]  // Coll_piMpiM
#define f4  ff[4]  // Unp_piPpiM
#define f5  ff[5]  // Unp_piMpiP
#define f6  ff[6]  // Unp_piPpiP
#define f7  ff[7]  // Unp_piMpiM
    // ff[8..20] hold the 13 azimuthal-moment components

    double sqrts  = userdata_value(userdata, kUserDataSqrts);
    double Q20    = userdata_value(userdata, kUserDataQ20);
    double thetac = userdata_value(userdata, kUserDataThetac);
    double Q2     = userdata_value(userdata, kUserDataFixedValue);
    double z1_min = userdata_value(userdata, kUserDataZ1Min);
    double z1_max = userdata_value(userdata, kUserDataZ1Max);
    double z2     = userdata_value(userdata, kUserDataZ2);
    auto pars = *static_cast<std::vector<double>**>(userdata)[7];

    double z1  = z1_min + x[0] * (z1_max - z1_min);
    double xB  = x[1];
    double csi = x[2];

    // Azimuthal-moments part: depends on (xB, csi) with fixed Q2
    YStarYKinematics kin = PhysicsCalculator::computeYStarYQ2(sqrts, Q20, thetac, xB, Q2, csi, *flux);
    double azim[kNumComponents];
    if (cut(kin)) {
        fill_azimuthal_results(kin, kin.KQx, azim);
    } else {
        std::fill(azim, azim + kNumComponents, 0.0);
    }
    for (int i = 0; i < kNumComponents; ++i) ff[8 + i] = azim[i];

    // Collins/FF part: depends on (z1, z2, Q2)
    double f[13];
    Collins_FF(hadron, charge, z1, z2, Q2);

    if (myCol.evo == "DGLAP" || myCol.evo == "none") {
        myCol.eval(z1, Q2, charge, FF_ppz1, pars);
        for (int i = 0; i < COL_ppz1.size(); ++i) COL_ppz1[i] = myCol.COL_z[i];

        myCol.eval(z2, Q2, charge, FF_ppz2, pars);
        for (int i = 0; i < COL_ppz2.size(); ++i) COL_ppz2[i] = myCol.COL_z[i];

        charge *= -1; //to call pi- Collins

        myCol.eval(z1, Q2, charge, FF_pmz1, pars);
        for (int i = 0; i < COL_pmz1.size(); ++i) COL_pmz1[i] = myCol.COL_z[i];

        myCol.eval(z2, Q2, charge, FF_pmz2, pars);
        for (int i = 0; i < COL_pmz2.size(); ++i) COL_pmz2[i] = myCol.COL_z[i];
    }

    if (myCol.evo == "CT3") {
        hoppetEvalcf(z1, sqrt(Q2), f);
        for (int i = 0; i < COL_ppz1.size(); ++i) COL_ppz1[i] = f[i] / z1;

        hoppetEvalcf(z2, sqrt(Q2), f);
        for (int i = 0; i < COL_ppz2.size(); ++i) COL_ppz2[i] = f[i] / z2;

        hoppetEvalcff(z1, sqrt(Q2), f);
        for (int i = 0; i < COL_pmz1.size(); ++i) COL_pmz1[i] = f[i] / z1;

        hoppetEvalcff(z2, sqrt(Q2), f);
        for (int i = 0; i < COL_pmz2.size(); ++i) COL_pmz2[i] = f[i] / z2;
    }

    Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1, FF_ppz2, FF_pmz1, FF_pmz2);

    f0 = Coll_piPpiM;
    f1 = Coll_piMpiP;
    f2 = Coll_piPpiP;
    f3 = Coll_piMpiM;
    f4 = Unp_piPpiM;
    f5 = Unp_piMpiP;
    f6 = Unp_piPpiP;
    f7 = Unp_piMpiM;

    return 0;
}


struct IntegrationResults {
    int maxeval;
    int nstart;
    int neval;
    int fail;
    int nregions;
    double integral[kNumComponents + kCollinsChannels];
    double error[kNumComponents + kCollinsChannels];
    double prob[kNumComponents + kCollinsChannels];
    double ratio[kNumComponents];
    double ratioerror[kNumComponents];
    double elapsed_seconds;
};


int main(int argc, char *argv[]) {
    (void)argc;

    const int maxeval = static_cast<int>(1e8);
    const int nstart  = static_cast<int>(1e7);

    stringstream modelstr, widthsstr, evostr, MCnamestring, EPAnamestr;
    modelstr     << argv[2];
    widthsstr    << argv[4];
    evostr       << argv[6];

    double z1_min = atof(argv[8]);
    double z1_max = atof(argv[9]);
    double z2_min = atof(argv[11]);
    double z2_max = atof(argv[12]);
    sqrts         = atof(argv[14]);
    Q20           = atof(argv[16]);
    double Q2val  = atof(argv[18]);
    thetac        = atof(argv[20]);
    EPAnamestr    << argv[22];
    MCnamestring  << argv[26];
    int iset_min  = atoi(argv[28]);
    int iset_max  = atoi(argv[30]);

    const std::string MCname = MCnamestring.str(), EPAname = EPAnamestr.str();
    CSV->Load(MCname);
    int Nset = CSV->GetRowCount();

    if (iset_min < 1 || iset_max >= Nset) {

        std::cerr << "iset out of range: iset_min = " << iset_min << ", iset_max = " << iset_max << ", Nset = " << Nset << std::endl;
        return 0;
    }

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

    cout << EPAname << endl;
    flux -> set_source(EPAname);
    flux -> set_thetac(thetac);
    flux -> set_s(sqrts * sqrts);

    int ndim = 3, ncomp = kNumComponents + kCollinsChannels, nvec = 1, verbose = 0, last = 4, key = 13;
    double epsrel = 1e-4, epsabs = 1e-5;
    int flags = 0, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;

    std::ofstream outFile_UL("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) +
                             "_iset_" + LHAPDF::to_str(iset_min) + "_" + LHAPDF::to_str(iset_max) + "_UL.txt");
    std::ofstream outFile_UC("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) +
                             "_iset_" + LHAPDF::to_str(iset_min) + "_" + LHAPDF::to_str(iset_max) + "_UC.txt");
    std::ofstream outFile_sep("Fixed_Q2_" + LHAPDF::to_str(Q2val) + "_Vs_" + LHAPDF::to_str(sqrts) +
                              "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                              LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) +
                              "_iset_" + LHAPDF::to_str(iset_min) + "_" + LHAPDF::to_str(iset_max) + "_separated.txt");

    outFile_UL << "z2"
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

    int min, max;
    if (iset_min == 0) min = 0;
    else if (iset_min > 0) min = iset_min;

    if (iset_max < Nset) max = iset_max + 1;
    else max = Nset;

    for (int iset = min; iset < max; iset++) {

        std::vector<double> pars = CSV->GetRow<double>(iset);

        double prefact = get_prefact(pars.back());

        for (double z2 : z2_values) {

            void *USERDATA[] = {&sqrts, &Q20, &thetac, &Q2val, &z1_min, &z1_max, &z2, &pars};

            std::cout << "Running scan for iset = " << iset << ", Q2 = " << Q2val << ", z2 = " << z2 << std::endl;

            IntegrationResults res{};
            res.maxeval = maxeval;
            res.nstart  = nstart;

            auto t0 = std::chrono::high_resolution_clock::now();

            Vegas(ndim, ncomp, integrand, USERDATA,
                  nvec, epsrel, epsabs,
                  flags, seed, mineval, maxeval,
                  nstart, nincrease, nbatch,
                  gridno, statefile, spin,
                  &res.neval, &res.fail,
                  res.integral, res.error, res.prob);

            auto t1 = std::chrono::high_resolution_clock::now();
            res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();

            // Collins ratios from components 0..7
            double numU = res.integral[0] + res.integral[1];
            double numL = res.integral[2] + res.integral[3];
            double denU = res.integral[4] + res.integral[5];
            double denL = res.integral[6] + res.integral[7];

            double ratio_U  = safe_ratio(numU, denU);
            double ratio_L  = safe_ratio(numL, denL);
            double ratio_UL = ratio_U - ratio_L;
            double ratio_UC = ratio_U - safe_ratio(numU + numL, denU + denL);

            std::cout << "Collins ratios: UL = " << ratio_UL
                      << "  UC = " << ratio_UC
                      << "  U = " << ratio_U
                      << "  L = " << ratio_L << std::endl;
            std::cout << " --> Integration time = " << res.elapsed_seconds
                      << "s, neval = " << res.neval << std::endl;

            // Azimuthal moments are stored in res.integral[8..20]; their reference
            // denominator is res.integral[8] (AUzi), which plays the role that
            // res.integral[0] played in the previous (separate) integrand_fixedQ2.

            // --- outFile_UL: scaling by prefact * ratio_UL ---
            for (int i = 0; i < kNumComponents; ++i) {
                res.ratio[i] = res.integral[kCollinsChannels + i] / res.integral[kCollinsChannels];
                if (i >= 5) {
                    res.ratio[i] *= prefact * ratio_UL;
                }
                res.ratioerror[i] = std::abs(res.ratio[i]) * std::sqrt(
                    std::pow(res.error[kCollinsChannels + i] / res.integral[kCollinsChannels + i], 2) +
                    std::pow(res.error[kCollinsChannels] / res.integral[kCollinsChannels], 2)
                );
            }

            outFile_UL << z2
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
                       << "," << res.integral[kCollinsChannels + 0]  << "," << res.error[kCollinsChannels + 0]
                       << "," << res.integral[kCollinsChannels + 1]  << "," << res.error[kCollinsChannels + 1]
                       << "," << res.integral[kCollinsChannels + 2]  << "," << res.error[kCollinsChannels + 2]
                       << "," << res.integral[kCollinsChannels + 3]  << "," << res.error[kCollinsChannels + 3]
                       << "," << res.integral[kCollinsChannels + 4]  << "," << res.error[kCollinsChannels + 4]
                       << "," << res.integral[kCollinsChannels + 5]  << "," << res.error[kCollinsChannels + 5]
                       << "," << res.integral[kCollinsChannels + 6]  << "," << res.error[kCollinsChannels + 6]
                       << "," << res.integral[kCollinsChannels + 7]  << "," << res.error[kCollinsChannels + 7]
                       << "," << res.integral[kCollinsChannels + 8]  << "," << res.error[kCollinsChannels + 8]
                       << "," << res.integral[kCollinsChannels + 9]  << "," << res.error[kCollinsChannels + 9]
                       << "," << res.integral[kCollinsChannels + 10] << "," << res.error[kCollinsChannels + 10]
                       << "," << res.integral[kCollinsChannels + 11] << "," << res.error[kCollinsChannels + 11]
                       << "," << res.integral[kCollinsChannels + 12] << "," << res.error[kCollinsChannels + 12]
                       << "," << res.elapsed_seconds
                       << "," << res.neval << "," << res.fail
                       << std::endl;

            // --- outFile_UC: scaling by prefact * ratio_UC ---
            for (int i = 0; i < kNumComponents; ++i) {
                res.ratio[i] = res.integral[kCollinsChannels + i] / res.integral[kCollinsChannels];
                if (i >= 5) {
                    res.ratio[i] *= prefact * ratio_UC;
                }
                res.ratioerror[i] = std::abs(res.ratio[i]) * std::sqrt(
                    std::pow(res.error[kCollinsChannels + i] / res.integral[kCollinsChannels + i], 2) +
                    std::pow(res.error[kCollinsChannels] / res.integral[kCollinsChannels], 2)
                );
            }

            outFile_UC << z2
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
                       << "," << res.integral[kCollinsChannels + 0]  << "," << res.error[kCollinsChannels + 0]
                       << "," << res.integral[kCollinsChannels + 1]  << "," << res.error[kCollinsChannels + 1]
                       << "," << res.integral[kCollinsChannels + 2]  << "," << res.error[kCollinsChannels + 2]
                       << "," << res.integral[kCollinsChannels + 3]  << "," << res.error[kCollinsChannels + 3]
                       << "," << res.integral[kCollinsChannels + 4]  << "," << res.error[kCollinsChannels + 4]
                       << "," << res.integral[kCollinsChannels + 5]  << "," << res.error[kCollinsChannels + 5]
                       << "," << res.integral[kCollinsChannels + 6]  << "," << res.error[kCollinsChannels + 6]
                       << "," << res.integral[kCollinsChannels + 7]  << "," << res.error[kCollinsChannels + 7]
                       << "," << res.integral[kCollinsChannels + 8]  << "," << res.error[kCollinsChannels + 8]
                       << "," << res.integral[kCollinsChannels + 9]  << "," << res.error[kCollinsChannels + 9]
                       << "," << res.integral[kCollinsChannels + 10] << "," << res.error[kCollinsChannels + 10]
                       << "," << res.integral[kCollinsChannels + 11] << "," << res.error[kCollinsChannels + 11]
                       << "," << res.integral[kCollinsChannels + 12] << "," << res.error[kCollinsChannels + 12]
                       << "," << res.elapsed_seconds
                       << "," << res.neval << "," << res.fail
                       << std::endl;

            // --- outFile_sep: scaling by prefact, U and L applied per column ---
            for (int i = 0; i < kNumComponents; ++i) {
                res.ratio[i] = res.integral[kCollinsChannels + i] / res.integral[kCollinsChannels];
                if (i >= 5) {
                    res.ratio[i] *= prefact;
                }
                res.ratioerror[i] = std::abs(res.ratio[i]) * std::sqrt(
                    std::pow(res.error[kCollinsChannels + i] / res.integral[kCollinsChannels + i], 2) +
                    std::pow(res.error[kCollinsChannels] / res.integral[kCollinsChannels], 2)
                );
            }

            outFile_sep << z2
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
                        << "," << res.ratio[6] * ratio_U << "," << res.ratioerror[6]
                        << "," << res.ratio[6] * ratio_L << "," << res.ratioerror[6]
                        << "," << res.ratio[7] * ratio_U << "," << res.ratioerror[7]
                        << "," << res.ratio[7] * ratio_L << "," << res.ratioerror[7]
                        << "," << res.ratio[8] * ratio_U << "," << res.ratioerror[8]
                        << "," << res.ratio[8] * ratio_L << "," << res.ratioerror[8]
                        << "," << res.ratio[9] * ratio_U << "," << res.ratioerror[9]
                        << "," << res.ratio[9] * ratio_L << "," << res.ratioerror[9]
                        << "," << res.ratio[10] * ratio_U << "," << res.ratioerror[10]
                        << "," << res.ratio[10] * ratio_L << "," << res.ratioerror[10]
                        << "," << res.ratio[11] * ratio_U << "," << res.ratioerror[11]
                        << "," << res.ratio[11] * ratio_L << "," << res.ratioerror[11]
                        << "," << res.ratio[12] * ratio_U << "," << res.ratioerror[12]
                        << "," << res.ratio[12] * ratio_L << "," << res.ratioerror[12]
                        << "," << res.elapsed_seconds
                        << "," << res.neval << "," << res.fail
                        << std::endl;

            std::cout << " --> time = " << res.elapsed_seconds << "s, neval = " << res.neval << std::endl;
        }

    }

    outFile_UL.close();
    outFile_UC.close();
    outFile_sep.close();

    return 0;
}
