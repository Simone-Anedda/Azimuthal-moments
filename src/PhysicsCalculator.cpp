#include "PhysicsCalculator.h"

#include <cmath>
#include <initializer_list>

#include "photon_flux.h"

namespace
{
constexpr double kAlphaEm = 1.0 / 137.036;
constexpr double kElectronMass = 5.11e-4;
constexpr double kPiValue = 3.14159265358979323846;

bool are_finite(std::initializer_list<double> values)
{
    for (const double value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

void evaluate_flux(EPA::EPA_flux& flux, double csi, double& unpolarised, double& polarised)
{
    unpolarised = flux.eval(csi);
    flux.set_polarisation(-1);
    polarised = flux.eval(csi);
    flux.set_polarisation(1);
}
}

YYKinematics PhysicsCalculator::computeYY(double sqrts,
                                          double kT0,
                                          double x1,
                                          double x2,
                                          double kT,
                                          EPA::EPA_flux& flux1,
                                          EPA::EPA_flux* flux2)
{
    YYKinematics out;
    out.sqrts = sqrts;
    out.S = sqrts * sqrts;
    out.kT0 = kT0;
    out.kTM = sqrts / 2.0;
    out.kT = kT;

    if (out.sqrts <= 0.0 || out.kT <= 0.0 || out.kT < out.kT0 || out.kT >= out.kTM) {
        return out;
    }

    double sqrt_arg = out.S / (out.kT * out.kT) - 4.0;
    if (sqrt_arg <= 0.0) {
        return out;
    }

    double sqrt_term = std::sqrt(sqrt_arg);
    double etaqmax = std::log(0.5 * (out.sqrts / out.kT + sqrt_term));
    double etaqmin = std::log(0.5 * (out.sqrts / out.kT - sqrt_term));
    out.etaq = (etaqmax - etaqmin) * x1 + etaqmin;

    double etaqbarmax_arg = out.sqrts / out.kT - std::exp(out.etaq);
    double etaqbarmin_arg = out.sqrts / out.kT - std::exp(-out.etaq);
    if (etaqbarmax_arg <= 0.0 || etaqbarmin_arg <= 0.0) {
        return out;
    }

    double etaqbarmax =   std::log(etaqbarmax_arg);
    double etaqbarmin = - std::log(etaqbarmin_arg);
    out.jacob = (etaqmax - etaqmin) * (etaqbarmax - etaqbarmin);
    out.etaqbar = (etaqbarmax - etaqbarmin) * x2 + etaqbarmin;

    out.csi1 = (out.kT / out.sqrts) * (std::exp(out.etaq) + std::exp(out.etaqbar));
    out.csi2 = (out.kT / out.sqrts) * (std::exp(-out.etaq) + std::exp(-out.etaqbar));
    if (!(out.csi1 > 0.0 && out.csi1 < 1.0 && out.csi2 > 0.0 && out.csi2 < 1.0)) {
        return out;
    }

    evaluate_flux(flux1, out.csi1, out.fgl1, out.DLfgl1);

    EPA::EPA_flux* second_flux = flux2 != nullptr ? flux2 : &flux1;
    evaluate_flux(*second_flux, out.csi2, out.fgl2, out.DLfgl2);

    if (!are_finite({out.fgl1, out.DLfgl1, out.fgl2, out.DLfgl2, out.jacob})) {
        return out;
    }

    double cosh_diff = std::cosh(out.etaq - out.etaqbar);
    double K = 1.0 / (out.S * out.kT * out.kT * (1.0 + cosh_diff));
    out.AU = K * cosh_diff * out.fgl1 * out.fgl2 * out.jacob;
    out.BU = K * out.fgl1 * out.fgl2 * out.jacob / 4.0;
    out.BL = K * out.DLfgl1 * out.DLfgl2 * out.jacob / 4.0;

    out.valid = are_finite({out.AU, out.BU, out.BL, out.jacob});
    return out;
}

YStarYKinematics PhysicsCalculator::computeYStarY(double sqrts,
                                                  double Q20,
                                                  double thetac,
                                                  double xB,
                                                  double y,
                                                  double csi)
{
    const double S = sqrts * sqrts;
    const double Q2 = xB * y * S;
    return finalizeYStarY(sqrts, Q20, thetac, xB, y, csi, Q2);
}

YStarYKinematics PhysicsCalculator::computeYStarYQ2(double sqrts,
                                                    double Q20,
                                                    double thetac,
                                                    double xB,
                                                    double Q2,
                                                    double csi)
{
    const double S = sqrts * sqrts;
    if (xB <= 0.0 || S <= 0.0) {
        return {};
    }

    const double y = Q2 / (xB * S);
    return finalizeYStarY(sqrts, Q20, thetac, xB, y, csi, Q2);
}

YStarYKinematics PhysicsCalculator::finalizeYStarY(double sqrts,
                                                   double Q20,
                                                   double thetac,
                                                   double xB,
                                                   double y,
                                                   double csi,
                                                   double Q2)
{
    YStarYKinematics out;
    out.sqrts = sqrts;
    out.S = sqrts * sqrts;
    out.Q20 = Q20;
    out.Q2M = out.S - 4.0 * out.Q20;
    out.thetac = thetac;
    out.xB = xB;
    out.y = y;
    out.csi = csi;
    out.Q2 = Q2;

    if (out.sqrts <= 0.0 || out.S <= 0.0 || out.Q20 <= 0.0 || out.Q2M <= 0.0 || out.thetac <= 0.0) {
        return out;
    }

    if (out.xB <= 0.0 || out.y <= 0.0 || out.csi <= 0.0 || out.csi >= 1.0 || out.Q2 <= 0.0) {
        return out;
    }

    out.shat = (out.csi - out.xB) * out.y * out.S;
    if (out.shat <= 0.0) {
        return out;
    }

    const double zeta_sqrt_arg = 1.0 - 4.0 * out.Q20 / out.shat;
    if (zeta_sqrt_arg <= 0.0) {
        return out;
    }

    const double zeta_sqrt = std::sqrt(zeta_sqrt_arg);
    out.zetamin = 0.5 * (1.0 - zeta_sqrt);
    out.zetamax = 0.5 * (1.0 + zeta_sqrt);

    out.Q2min = kElectronMass * kElectronMass * out.csi * out.csi / (1.0 - out.csi);
    out.Q2max = (out.S / 4.0) * out.thetac * out.thetac * (1.0 - out.csi) + out.Q2min;
    if (out.Q2min <= 0.0 || out.Q2max <= out.Q2min) {
        return out;
    }

    const double xBmin = out.Q2 / out.S;
    const double xBmax = out.Q2 / (out.Q2 + 4.0 * out.Q20);
    const double csimin = out.xB * (1.0 + 4.0 * out.Q20 / out.Q2);
    const double ymin = out.xB * (1.0 + 4.0 * out.Q20 / out.Q2);

    if (out.Q2 < out.Q20 || out.Q2 > out.Q2M) {
        return out;
    }
    if (out.xB < xBmin || out.xB > xBmax) {
        return out;
    }
    if (out.csi < csimin || out.csi > 1.0) {
        return out;
    }
    if (out.y < ymin || out.y > 1.0) {
        return out;
    }

    const double log_ratio = std::log(out.Q2max / out.Q2min);
    const double me2_term = 2.0 * kElectronMass * kElectronMass * out.csi * (1.0 / out.Q2max - 1.0 / out.Q2min);

    out.fgl = (kAlphaEm / (2.0 * kPiValue)) *
              ((1.0 + (1.0 - out.csi) * (1.0 - out.csi)) / out.csi * log_ratio + me2_term);
    out.DLfgl = (kAlphaEm / (2.0 * kPiValue)) *
                ((1.0 - (1.0 - out.csi) * (1.0 - out.csi)) / out.csi * log_ratio + me2_term);

    out.Kxy = 1.0 / (out.xB * out.y * out.y * out.csi * out.csi * out.csi * out.S);
    out.KQx = 1.0 / (out.csi * out.csi * out.csi * out.Q2 * out.Q2);
    out.KQy = 1.0 / (out.Q2 * out.y * out.y * out.csi * out.csi * out.csi);

    const double one_minus_y = 1.0 - out.y;
    const double csi_minus_xB = out.csi - out.xB;
    const double xB_csi_term = out.xB * out.xB + csi_minus_xB * csi_minus_xB;
    const double delta_zeta = out.zetamax - out.zetamin;
    const double log_zeta_ratio = std::log(out.zetamax / out.zetamin);
    const double sqrt_xBcsi = std::sqrt(out.xB * csi_minus_xB);
    const double sqrt_zetas = std::sqrt(out.zetamax * out.zetamin);
    const double atan_term = -std::atan(std::sqrt(out.zetamin / out.zetamax)) +
                             std::atan(std::sqrt(out.zetamax / out.zetamin));

    {
        const double term1 = (1.0 + one_minus_y * one_minus_y) * xB_csi_term;
        const double term2 = -2.0 * delta_zeta + 2.0 * log_zeta_ratio;
        const double term3 = 16.0 * one_minus_y * out.xB * csi_minus_xB * delta_zeta;
        out.AUzi = 2.0 * (term1 * term2 + term3) * out.fgl;
    }

    {
        const double term1 = -8.0 * (2.0 - out.y) * std::sqrt(one_minus_y) * (out.csi - 2.0 * out.xB);
        const double term2 = sqrt_xBcsi * (2.0 * sqrt_zetas - 1.0);
        out.AUcfqzi = term1 * term2 * out.fgl;
    }

    out.AUc2fqzi = 16.0 * one_minus_y * out.xB * csi_minus_xB * delta_zeta * out.fgl;

    {
        const double term1 = -2.0 * out.y * (2.0 - out.y) * out.csi * (out.csi - 2.0 * out.xB);
        const double term2 = -2.0 * delta_zeta + 2.0 * log_zeta_ratio;
        out.ALzi = term1 * term2 * out.DLfgl;
    }

    {
        const double term1 = 8.0 * out.y * std::sqrt(one_minus_y) * out.csi * sqrt_xBcsi;
        const double term2 = 2.0 * sqrt_zetas - 1.0;
        out.ALcfqzi = term1 * term2 * out.DLfgl;
    }

    {
        const double term1 = (1.0 + one_minus_y * one_minus_y) * xB_csi_term;
        const double term2 = 8.0 * one_minus_y * out.xB * csi_minus_xB;
        out.BUcf12zi = (term1 - term2) * delta_zeta * out.fgl;
    }

    {
        const double term1 = -2.0 * (2.0 - out.y) * std::sqrt(one_minus_y) *
                             (out.csi - 2.0 * out.xB) * sqrt_xBcsi;
        out.BUcfqmf12zi = term1 * atan_term * out.fgl;
    }

    {
        const double term1 = 2.0 * (2.0 - out.y) * std::sqrt(one_minus_y) *
                             (out.csi - 2.0 * out.xB) * sqrt_xBcsi;
        out.BUcfqpf12zi = term1 * atan_term * out.fgl;
    }

    {
        const double term1 = 2.0 * one_minus_y * out.xB * csi_minus_xB;
        const double term2 = -delta_zeta + log_zeta_ratio;
        out.BUc2fqmf12zi = term1 * term2 * out.fgl;
        out.BUc2fqpf12zi = term1 * term2 * out.fgl;
    }

    out.BLcf12zi = -out.y * (2.0 - out.y) * out.csi * (out.csi - 2.0 * out.xB) * delta_zeta * out.DLfgl;

    {
        const double term1 = 2.0 * out.y * std::sqrt(one_minus_y) * out.csi * sqrt_xBcsi;
        out.BLcfqmf12zi = term1 * atan_term * out.DLfgl;
    }

    {
        const double term1 = -2.0 * out.y * std::sqrt(one_minus_y) * out.csi * sqrt_xBcsi;
        out.BLcfqpf12zi = term1 * atan_term * out.DLfgl;
    }

    out.valid = are_finite({out.fgl, out.DLfgl, out.Kxy, out.KQx, out.KQy,
                            out.AUzi, out.AUcfqzi, out.AUc2fqzi, out.ALzi, out.ALcfqzi,
                            out.BUcf12zi, out.BUcfqmf12zi, out.BUcfqpf12zi,
                            out.BUc2fqmf12zi, out.BUc2fqpf12zi,
                            out.BLcf12zi, out.BLcfqmf12zi, out.BLcfqpf12zi});
    return out;
}
