#pragma once

namespace EPA
{
class EPA_flux;
}

struct YYKinematics
{
    bool valid = false;

    double sqrts = 0.0;
    double S = 0.0;
    double kT0 = 0.0;
    double kTM = 0.0;

    double etaq = 0.0;
    double etaqbar = 0.0;
    double kT = 0.0;
    double csi1 = 0.0;
    double csi2 = 0.0;

    double fgl1 = 0.0;
    double DLfgl1 = 0.0;
    double fgl2 = 0.0;
    double DLfgl2 = 0.0;

    double jacob = 0.0;
    double K = 0.0;
    double AU = 0.0;
    double BU = 0.0;
    double BL = 0.0;
};

struct YStarYKinematics
{
    bool valid = false;

    double sqrts = 0.0;
    double S = 0.0;
    double Q20 = 0.0;
    double Q2M = 0.0;
    double thetac = 0.0;

    double xB = 0.0;
    double y = 0.0;
    double csi = 0.0;
    double Q2 = 0.0;
    double shat = 0.0;

    double zetamin = 0.0;
    double zetamax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double fgl = 0.0;
    double DLfgl = 0.0;

    double Kxy = 0.0;
    double KQx = 0.0;
    double KQy = 0.0;

    double AUzi = 0.0;
    double AUcfqzi = 0.0;
    double AUc2fqzi = 0.0;
    double ALzi = 0.0;
    double ALcfqzi = 0.0;
    double BUcf12zi = 0.0;
    double BUcfqmf12zi = 0.0;
    double BUcfqpf12zi = 0.0;
    double BUc2fqmf12zi = 0.0;
    double BUc2fqpf12zi = 0.0;
    double BLcf12zi = 0.0;
    double BLcfqmf12zi = 0.0;
    double BLcfqpf12zi = 0.0;
};

class PhysicsCalculator
{
public:
    static YYKinematics computeYY(double sqrts,
                                  double kT0,
                                  double thetac,
                                  double x1,
                                  double x2,
                                  double kT,
                                  EPA::EPA_flux& flux1,
                                  EPA::EPA_flux* flux2 = nullptr);

    static YStarYKinematics computeYStarY(double sqrts,
                                          double Q20,
                                          double thetac,
                                          double xB,
                                          double y,
                                          double csi);

    static YStarYKinematics computeYStarYQ2(double sqrts,
                                            double Q20,
                                            double thetac,
                                            double xB,
                                            double Q2,
                                            double csi);

private:
    static YStarYKinematics finalizeYStarY(double sqrts,
                                           double Q20,
                                           double thetac,
                                           double xB,
                                           double y,
                                           double csi,
                                           double Q2);
};
