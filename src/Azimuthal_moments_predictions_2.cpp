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
constexpr int kNumComponents = 3;
constexpr int kUserDataSqrts = 0;
constexpr int kUserDataQ20 = 1;
constexpr int kUserDataThetac = 2;
constexpr int kUserDataKTMax = 3;
constexpr int kUserDataKTMin = 4;
constexpr int kUserDataZ1Min = 5;
constexpr int kUserDataZ1Max = 6;
constexpr int kUserDataZ2 = 7;

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

std::vector<double> pars = {0.5316, -.7290, .8246, 3.2293, 0.8568, -.3539, 1.7524, 0.0, 0.2551};
std::vector<double> COL_z(13), COL_ppz1(13), COL_ppz2(13), COL_pmz1(13), COL_pmz2(13);
std::vector<double> FF(13), FF_ppz1(13), FF_ppz2(13), FF_pmz1(13), FF_pmz2(13),\
    FF_ppz1_fixedQ2(13), FF_ppz2_fixedQ2(13), FF_pmz1_fixedQ2(13), FF_pmz2_fixedQ2(13),\
    FF_evo(13);

double f[13], h1[13]; //for HOPPET and transversity

const double charges[13] = {-2./3., 1./3., -2./3., 1./3., -2./3., 1./3., 0, -1./3., 2./3., -1./3., 2./3., -1./3., 2./3.};

double sqrts = 0.0, thetac = 0.0, Q20 = 0.0,\
        z1 = 0.0, z2 = 0.0, pperp2 = 0.12;
double Coll_piPpiM = 0.0, Coll_piMpiP = 0.0, Coll_piPpiP = 0.0, Coll_piMpiM = 0.0;
double Unp_piPpiM = 0.0, Unp_piMpiP = 0.0, Unp_piPpiP = 0.0, Unp_piMpiM = 0.0;
double MC2 = pars.back();
double pperp2Col = MC2 * pperp2 / (MC2 + pperp2);
double prefact = (pi * EulerConst / 2) * (pow(pperp2Col, 3) / (pperp2 * pperp2 * MC2));

int charge = 1, hadron = 1;

FRAG::FF myFF("DEHSS","NLO");
COL::COLLINS myCol;
EPA::EPA_flux *flux1 = new EPA::EPA_flux();
EPA::EPA_flux *flux2 = new EPA::EPA_flux();

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



void Collins_FF( int & hadron, int & charge , double & z1, double & z2, double & Q2)
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


bool cut(YYKinematics& kin)
{
    return kin.valid;
}

int integrand(const int *ndim, const double x[], const int *ncomp, double ff[], void *userdata)
{
    (void)ndim;
    (void)ncomp;

#define f0 ff[0] // Coll_piPpiM
#define f1 ff[1] // Coll_piMpiP
#define f2 ff[2] // Coll_piPpiP
#define f3 ff[3] // Coll_piMpiM
#define f4 ff[4] // Unp_piPpiM
#define f5 ff[5] // Unp_piMpiP
#define f6 ff[6] // Unp_piPpiP
#define f7 ff[7] // Unp_piMpiM
#define f8 ff[8] // AU
#define f9 ff[9] // BU
#define f10 ff[10] // BL


    //double Q        = userdata_value(userdata, kUserDataFixedValue);
    //double Q2       = Q * Q;
    double z1_min   = userdata_value(userdata, kUserDataZ1Min);
    double z1_max   = userdata_value(userdata, kUserDataZ1Max);
    double z2       = userdata_value(userdata, kUserDataZ2);
    double sqrts    = userdata_value(userdata, kUserDataSqrts);
    double kT0      = userdata_value(userdata, kUserDataQ20);
    double kT_max   = userdata_value(userdata, kUserDataKTMax);
    double kT_min   = userdata_value(userdata, kUserDataKTMin);
    double thetac   = userdata_value(userdata, kUserDataThetac);
 
    double z1 = z1_min + x[0] * (z1_max - z1_min);

    double kT;
    double jabob_kT;
    if (kT_max==kT_min)
    {
        kT = kT_min;
        jabob_kT = 1.0;
    }
    else
    {
        kT = kT_min + x[1] * (kT_max - kT_min);
        jabob_kT = kT_max - kT_min;
    }
    double Q  = kT;
    double Q2 = Q * Q;

    double f[13];

    const YYKinematics kin = PhysicsCalculator::computeYY(sqrts, kT0, thetac, x[3], x[4], kT, *flux1, flux2);

    if (cut(kin))
    {
        f8 = kin.AU * jabob_kT;
        f9 = kin.BU * jabob_kT;
        f10 = kin.BL * jabob_kT;
    } 
    else
    {
        f8 = 0.0;
        f9 = 0.0;
        f10 = 0.0;
    }



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
        hoppetEvalcf(z1, Q, f);
        for (int i = 0; i < COL_ppz1.size(); ++i) COL_ppz1[i] = f[i] / z1;

        hoppetEvalcf(z2, Q, f);
        for (int i = 0; i < COL_ppz2.size(); ++i) COL_ppz2[i] = f[i] / z2;

        hoppetEvalcff(z1, Q, f);
        for (int i = 0; i < COL_pmz1.size(); ++i) COL_pmz1[i] = f[i] / z1;

        hoppetEvalcff(z2, Q, f );
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

    const int maxeval = static_cast<int>(1e8);
    const int nstart  = static_cast<int>(1e7);

    stringstream modelstr, widthsstr, evostr, MCnamestring, EPA1namestr, EPA2namestr;
    modelstr    << argv[2];
    widthsstr   << argv[4];
    evostr      << argv[6];

    double z1_min = atof(argv[8]);
    double z1_max = atof(argv[9]);
    double z2_min = atof(argv[11]);
    double z2_max = atof(argv[12]);
    sqrts         = atof(argv[14]);
    Q20           = atof(argv[16]); // kept for input compatibility, interpreted here as kT_min
    double kT_min = atof(argv[18]);
    double kT_max = atof(argv[19]);
    thetac        = atof(argv[21]);
    EPA1namestr   << argv[23];
    EPA2namestr   << argv[25];
    MCnamestring  << argv[27];

    const std::string MCname = MCnamestring.str(), EPA1name = EPA1namestr.str(), EPA2name = EPA2namestr.str();
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

    flux1 -> set_source(EPA1name);
    flux1 -> set_thetac(thetac);
    flux1 -> set_s(sqrts * sqrts);

    if (EPA1name == EPA2name){

        cout << "Set EPA flux for " << EPA1name << endl;
        delete flux2;
        flux2 = nullptr;
    }
    else{

        cout << "Setting EPA fluxes for " << EPA1name << " and " << EPA2name << endl;
        flux2 -> set_source(EPA2name);
        flux2 -> set_thetac(thetac);
        flux2 -> set_s(sqrts * sqrts);
    }

    int ndim = 4, ncomp = kNumComponents + 8, ncompColl = 8, nvec = 1, verbose = 0, last = 4, key = 13;
    double epsrel = 1e-6, epsabs = 1e-12;
    int flags = 2, seed = 0, mineval = 0, nincrease = 0, nbatch = 1000, gridno = 0;
    char statefile[64] = "";
    void* spin = nullptr;

    std::ofstream outFile_UL("kT_max_" + LHAPDF::to_str(kT_max) + "_kT_min_" + LHAPDF::to_str(kT_min) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UL.txt");
    std::ofstream outFile_UC("kT_max_" + LHAPDF::to_str(kT_max) + "_kT_min_" + LHAPDF::to_str(kT_min) + "_Vs_" + LHAPDF::to_str(sqrts) +
                             "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                             LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_UC.txt");
    std::ofstream outFile_U_L("kT_max_" + LHAPDF::to_str(kT_max) + "_kT_min_" + LHAPDF::to_str(kT_min) + "_Vs_" + LHAPDF::to_str(sqrts) +
                              "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                              LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_U_L.txt");                         
    std::ofstream outFile_sep("kT_max_" + LHAPDF::to_str(kT_max) + "_kT_min_" + LHAPDF::to_str(kT_min) + "_Vs_" + LHAPDF::to_str(sqrts) +
                              "_thetac_" + LHAPDF::to_str(thetac) + "_z1_" +
                              LHAPDF::to_str(z1_min) + "_" + LHAPDF::to_str(z1_max) + "_separated.txt");

    outFile_UL << "z2"
               << ",<dsig|c12>_UL" <<",err_<dsig|c12>_UL"
               << ",<ALL|c12>_UL" << ",err_<ALL|c12>_UL"
               << ",BU_o_AU" << ",err_BU_o_AU"
               << ",BL_o_AU" << ",err_BL_o_AU"
               << ",Collins_UL" << ",times[s]" << ",neval" << ",fail"
               << std::endl;

    outFile_UC << "z2"
               << ",<dsig|c12>_UC" <<",err_<dsig|c12>_UC"
               << ",<ALL|c12>_UC" << ",err_<ALL|c12>_UC"
               << ",BU_o_AU" << ",err_BU_o_AU"
               << ",BL_o_AU" << ",err_BL_o_AU"
               << ",Collins_UC" << ",times[s]" << ",neval" << ",fail"
               << std::endl;

    outFile_U_L << "z2"
                << ",<dsig|c12>_U" <<",err_<dsig|c12>_U"
                << ",<ALL|c12>_U" << ",err_<ALL|c12>_U"
                << ",<dsig|c12>_L" <<",err_<dsig|c12>_L"
                << ",<ALL|c12>_L" << ",err_<ALL|c12>_L"
                << ",BU_o_AU" << ",err_BU_o_AU"
                << ",BL_o_AU" << ",err_BL_o_AU"
                << ",Collins_U" << ",Collins_L"
                << ",time[s]" << ",neval" << ",fail"
                << std::endl;

    outFile_sep << "z2"
                << ",<dsig|c12>_piPpiM" <<",err_<dsig|c12>_piPpiM"
                << ",<dsig|c12>_piMpiP" <<",err_<dsig|c12>_piMpiP"
                << ",<dsig|c12>_piPpiP" <<",err_<dsig|c12>_piPpiP"
                << ",<dsig|c12>_piMpiM" <<",err_<dsig|c12>_piMpiM"
                << ",<ALL|c12>_piPpiM" << ",err_<ALL|c12>_piPpiM"
                << ",<ALL|c12>_piMpiP" << ",err_<ALL|c12>_piMpiP"
                << ",<ALL|c12>_piPpiP" << ",err_<ALL|c12>_piPpiP"
                << ",<ALL|c12>_piMpiM" << ",err_<ALL|c12>_piMpiM"
                << ",BU_o_AU" << ",err_BU_o_AU"
                << ",BL_o_AU" << ",err_BL_o_AU"
                << ",Collins_piPpiM" << ",Collins_piMpiP" << ",Collins_piPpiP" << ",Collins_piMpiM"
                << ",time[s]" << ",neval" << ",fail"
                << std::endl;

    for (double z2_value : z2_values) {
        void *USERDATA[] = {&sqrts, &Q20, &thetac, &kT_min, &kT_max, &z1_min, &z1_max, &z2_value};
        
        std::cout << "Running scan for kT =[" << kT_min << "," << kT_max << "], z2 = " << z2_value << std::endl;
        
        IntegrationResults res{};
        res.maxeval = maxeval;
        res.nstart  = nstart;

        auto t0 = std::chrono::high_resolution_clock::now();
        Cuhre(ndim, ncomp, integrand, USERDATA, nvec,
              epsrel, epsabs, verbose | last,
              mineval, maxeval, key,
              statefile, spin,
              &res.nregions, &res.neval, &res.fail,
              res.integral, res.error, res.prob);
        auto t1 = std::chrono::high_resolution_clock::now();
        res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();

        double numU = res.integral[0] + res.integral[1];
        double numL = res.integral[2] + res.integral[3];
        double denU = res.integral[4] + res.integral[5];
        double denL = res.integral[6] + res.integral[7];

        double ratio_piPpiM = safe_ratio(res.integral[0], res.integral[4]);
        double ratio_piMpiP = safe_ratio(res.integral[1], res.integral[5]);
        double ratio_piPpiP = safe_ratio(res.integral[2], res.integral[6]);
        double ratio_piMpiM = safe_ratio(res.integral[3], res.integral[7]);
        
        double ratio_U  = safe_ratio( numU, denU );
        double ratio_L  = safe_ratio( numL, denL );
        double ratio_UL = ratio_U - ratio_L;
        double ratio_UC = ratio_U - safe_ratio( numU + numL, denU + denL );

        double BU_o_AU = safe_ratio(res.integral[9], res.integral[8]);
        double BU_o_AU_error = safe_ratio_error(res.integral[9], res.error[9], res.integral[8], res.error[8]);
        double BL_o_AU = safe_ratio(res.integral[10], res.integral[8]);
        double BL_o_AU_error = safe_ratio_error(res.integral[10], res.error[10], res.integral[8], res.error[8]);

        std::cout << "Collins ratios: UL = " << ratio_UL
                  << "  UC = " << ratio_UC
                  << "  U = " << ratio_U
                  << "  L = " << ratio_L << std::endl;
        std::cout << " --> Integration time=" << res.elapsed_seconds
                  << "s, neval=" << res.neval << std::endl;

        double scale_UL = prefact * ratio_UL;
        double scale_UC = prefact * ratio_UC;
        double scale_U = prefact * ratio_U;
        double scale_L = prefact * ratio_L;
        double scale_piPpiM = prefact * ratio_piPpiM;
        double scale_piMpiP = prefact * ratio_piMpiP;
        double scale_piPpiP = prefact * ratio_piPpiP;
        double scale_piMpiM = prefact * ratio_piMpiM;

        outFile_UL << "," << z2_value
                   << "," << BU_o_AU * scale_UL << "," << BU_o_AU_error * std::abs(scale_UL)
                   << "," << BL_o_AU * scale_UL << "," << BL_o_AU_error * std::abs(scale_UL)
                   << "," << BU_o_AU << "," << BU_o_AU_error
                   << "," << BL_o_AU << "," << BL_o_AU_error
                   << "," << ratio_UL<< "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                   << std::endl;

        outFile_UC << "," << z2_value
                   << "," << BU_o_AU * scale_UC << "," << BU_o_AU_error * std::abs(scale_UC)
                   << "," << BL_o_AU * scale_UC << "," << BL_o_AU_error * std::abs(scale_UC)
                   << "," << BU_o_AU << "," << BU_o_AU_error
                   << "," << BL_o_AU << "," << BL_o_AU_error
                   << "," << ratio_UC << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                   << std::endl;

        outFile_U_L << "," << z2_value
                    << "," << BU_o_AU * scale_U << "," << BU_o_AU_error * std::abs(scale_U)
                    << "," << BL_o_AU * scale_U << "," << BL_o_AU_error * std::abs(scale_U)
                    << "," << BU_o_AU * scale_L << "," << BU_o_AU_error * std::abs(scale_L)
                    << "," << BL_o_AU * scale_L << "," << BL_o_AU_error * std::abs(scale_L)
                    << "," << BU_o_AU << "," << BU_o_AU_error
                    << "," << BL_o_AU << "," << BL_o_AU_error
                    << "," << ratio_U << "," << ratio_L
                    << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                    << std::endl;

        outFile_sep << "," << z2_value
                    << "," << BU_o_AU * scale_piPpiM << "," << BU_o_AU_error * std::abs(scale_piPpiM)
                    << "," << BL_o_AU * scale_piPpiM << "," << BL_o_AU_error * std::abs(scale_piPpiM)
                    << "," << BU_o_AU * scale_piMpiP << "," << BU_o_AU_error * std::abs(scale_piMpiP)
                    << "," << BL_o_AU * scale_piMpiP << "," << BL_o_AU_error * std::abs(scale_piMpiP)
                    << "," << BU_o_AU * scale_piPpiP << "," << BU_o_AU_error * std::abs(scale_piPpiP)
                    << "," << BL_o_AU * scale_piPpiP << "," << BL_o_AU_error * std::abs(scale_piPpiP)
                    << "," << BU_o_AU * scale_piMpiM << "," << BU_o_AU_error * std::abs(scale_piMpiM)
                    << "," << BL_o_AU * scale_piMpiM << "," << BL_o_AU_error * std::abs(scale_piMpiM)
                    << "," << ratio_piPpiM << "," << ratio_piMpiP << "," << ratio_piPpiP << "," << ratio_piMpiM
                    << "," << res.elapsed_seconds << "," << res.neval << "," << res.fail
                    << std::endl;

        std::cout << " --> kT integration time=" << res.elapsed_seconds
                  << "s, neval=" << res.neval << std::endl;
    }

    outFile_UL.close();
    outFile_UC.close();
    outFile_U_L.close();
    outFile_sep.close();

    return 0;
}
