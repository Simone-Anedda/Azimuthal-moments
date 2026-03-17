//
// Author: Carlo Flore <carlo.flore@unica.it>
//

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>

    
namespace EPA{

    using namespace std;
    
class EPA_flux{

 public:
     
    EPA_flux(){};
    
    EPA_flux(const string source_in){
        
        _source = source_in;
        
        if(_Amap.find(_source) != _Amap.end()){
            
            set_bmin(_source);
        
        }
        
        else if(_leptonMasses.find(_source) != _leptonMasses.end()){
            
            _ml = _leptonMasses.at(_source);
            std::cout << "ml = " << _ml << std::endl;
            
        }
    };
    
    EPA_flux(const string source_in, const int pol_in){
        
        _pol = pol_in;

        _source = source_in;
        
        EPA_flux(_source);
    };
   
    ~EPA_flux(){};
    
    int _pol = 1;

    double _Q2max, _bmin, _ml, _thetac, _s;
    
    std::string _source;
    
    double eval(const double & csi);

    double Q2min(const double & x);

    double Q2max(const double & x);
    
    void set_thetac(const double & thetac);
    
    void set_s(const double & s);
    
    void set_bmin(const string nucleus);
    
    void set_polarisation(const int & pol);

    void set_source(const string & source);
    
    double EPA_lepton(const double & x);
    
    double EPA_nucleus(const double & csi);
        

    std::map<std::string, int> _Zmap = {
            {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
            {"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},
            {"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},{"Zn",30},
            {"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},{"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},
            {"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},{"In",49},{"Sn",50},
            {"Sb",51},{"Te",52},{"I",53},{"Xe",54},{"Cs",55},{"Ba",56},{"La",57},{"Ce",58},{"Pr",59},{"Nd",60},
            {"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},{"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},
            {"Lu",71},{"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},{"Au",79},{"Hg",80},
            {"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},{"Fr",87},{"Ra",88},{"Ac",89},{"Th",90},
            {"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Am",95},{"Cm",96},{"Bk",97},{"Cf",98},{"Es",99},{"Fm",100},
            {"Md",101},{"No",102},{"Lr",103},{"Rf",104},{"Db",105},{"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},
            {"Rg",111},{"Cn",112},{"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118}
        };

    std::map<std::string, int> _Amap = {
            {"H",1},{"He",4},{"Li",7},{"Be",9},{"B",11},{"C",12},{"N",14},{"O",16},{"F",19},{"Ne",20},
            {"Na",23},{"Mg",24},{"Al",27},{"Si",28},{"P",31},{"S",32},{"Cl",35},{"Ar",40},{"K",39},{"Ca",40},
            {"Sc",45},{"Ti",48},{"V",51},{"Cr",52},{"Mn",55},{"Fe",56},{"Co",59},{"Ni",58},{"Cu",63},{"Zn",64},
            {"Ga",69},{"Ge",74},{"As",75},{"Se",80},{"Br",79},{"Kr",84},{"Rb",85},{"Sr",88},{"Y",89},{"Zr",90},
            {"Nb",93},{"Mo",98},{"Tc",98},{"Ru",102},{"Rh",103},{"Pd",106},{"Ag",107},{"Cd",114},{"In",115},{"Sn",120},
            {"Sb",121},{"Te",130},{"I",127},{"Xe",132},{"Cs",133},{"Ba",138},{"La",139},{"Ce",140},{"Pr",141},{"Nd",144},
            {"Pm",145},{"Sm",150},{"Eu",152},{"Gd",158},{"Tb",159},{"Dy",164},{"Ho",165},{"Er",166},{"Tm",169},{"Yb",174},
            {"Lu",175},{"Hf",180},{"Ta",181},{"W",184},{"Re",187},{"Os",192},{"Ir",193},{"Pt",195},{"Au",197},{"Hg",202},
            {"Tl",205},{"Pb",208},{"Bi",209},{"Po",209},{"At",210},{"Rn",222},{"Fr",223},{"Ra",226},{"Ac",227},{"Th",232},
            {"Pa",231},{"U",238},{"Np",237},{"Pu",244},{"Am",243},{"Cm",247},{"Bk",247},{"Cf",251},{"Es",252},{"Fm",257},
            {"Md",258},{"No",259},{"Lr",266},{"Rf",267},{"Db",268},{"Sg",269},{"Bh",270},{"Hs",277},{"Mt",278},{"Ds",281},
            {"Rg",282},{"Cn",285},{"Nh",286},{"Fl",289},{"Mc",290},{"Lv",293},{"Ts",294},{"Og",294}
        };
        
    std::map<std::string, double> _leptonMasses = {
            {"electron", 0.00051099895},  // GeV/c^2
            {"muon",     0.1056583755},   // GeV/c^2
            {"tau",      1.77686}         // GeV/c^2
        };

    private:

        double _R0 = 1.2; //Fermi radius, 1.2 fm

        double _alphaem = 1./137.036;

        double _fmToGeV = 5.068; //fm -> GeV^-1 (ħc = 0.1973 GeV fm)

        double _MN = 0.938; //nucleon mass in GeV

        double _pi = 3.14159265358979323846;

};

}
