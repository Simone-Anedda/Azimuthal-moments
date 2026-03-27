//
// Author: Carlo Flore <carlo.flore@unica.it>
//

#include "photon_flux.h"

namespace EPA
{

    using namespace std;

    double EPA_flux::eval(const double &csi)
    {

        if (_Amap.find(_source) != _Amap.end())
        {
            return EPA_nucleus(csi);
        }

        else if (_leptonMasses.find(_source) != _leptonMasses.end())
        {
            return EPA_lepton(csi);
        }

        else return 0.0;
    }

    double EPA_flux::Q2min(const double &x)
    {
        return pow(_ml * x, 2) / (1. - x);
    }

    double EPA_flux::Q2max(const double &x)
    {
        return _s / 4. * pow(_thetac, 2) * (1. - x) + Q2min(x);
    }

    void EPA_flux::set_thetac(const double &thetac)
    {
        _thetac = thetac;
    }

    void EPA_flux::set_s(const double &s)
    {
        _s = s;
    }

    void EPA_flux::set_bmin(const string nucleus)
    {
        _bmin = 2 * _R0 * pow(_Amap.at(nucleus), 1. / 3.) * _fmToGeV;
    }

    void EPA_flux::set_polarisation(const int &pol)
    {
        _pol = pol;
    }

    void EPA_flux::set_source(const string &source)
    {
        _source = source;

        if(_Amap.find(_source) != _Amap.end()){

            set_bmin(_source);

        }

        else if(_leptonMasses.find(_source) != _leptonMasses.end()){

            _ml = _leptonMasses.at(_source);
            std::cout << "ml = " << _ml << std::endl;

        }
    }

    double EPA_flux::EPA_lepton(const double &x)
    {
        double flux = 0.0;

        // photon PDF from electron
        if (x < 1)
        {
            if (Q2min(x) < Q2max(x))
            {
                double pol_factor = (_pol == -1) ? x : 1.0;

                flux = _alphaem / (2. * _pi) *
                       ((1. + _pol * (1. - x) * (1. - x)) / x * log(Q2max(x) / Q2min(x)) +
                        2. * x * pol_factor * _ml * _ml * (1. / Q2max(x) - 1 / Q2min(x)));
            }

            else
                flux = 0.0;
        }

        else
            flux = 0.0;

        return flux;
    }

    double EPA_flux::EPA_nucleus(const double &csi)
    {

        double flux = 0.0;

        // photon PDF from nucleus
        if (csi < 1)
        {

            double xi = csi * _bmin * _MN;

            double Z = _Zmap.at(_source);

            flux = Z * Z * _alphaem / csi *
                   (2 * xi * cyl_bessel_k(0.0, xi) * cyl_bessel_k(1.0, xi) -

                    xi * xi * (pow(cyl_bessel_k(1.0, xi), 2) - pow(cyl_bessel_k(0.0, xi), 2)));
        }

        else
            flux = 0.0;

        return flux;
    }

}
