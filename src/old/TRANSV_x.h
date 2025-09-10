//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
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
#include <vector>

#include "transv_wrapper.h"
#include "FCN.h"

// #include<LHAPDF/LHAPDF.h>
#include "CollPDF.h"

    
namespace TRANSVX{

    using namespace std;
    
class getTRANSV_x{

 public:
     
    getTRANSV_x(const int & model_in, const int & widths_in, const int & evo_in){
           
        transv_widths = widths_in;
        
        transv_model = model_in;
    
        transv_evo = evo_in;
        
        transv_iset = 0;
        
        TRANSV_x.resize(13);
        
        SB_x.resize(13);
        
        max_TRANSV_x.resize(13);
        
        min_TRANSV_x.resize(13);
        
        factor = 0;
        
        PolyNorm = 0;
        
        hoppet_transv = false;
        
        for(int i = 0; i < max_TRANSV_x.size(); i++){
            max_TRANSV_x[i] = -pow(1,-10);            
            min_TRANSV_x[i] = pow(1,-10);
        }
        
    };
   
    ~getTRANSV_x(){};
         
    int transv_model;
    
    int transv_widths;
    
    int transv_evo;
    
    int transv_iset;
    
    bool hoppet_transv;

    vector<double> TRANSV_x, max_TRANSV_x, min_TRANSV_x;
    
    vector<double> SB_x;
    
    double factor, PolyNorm;
    
    double h1uvMax, h1uvMin, h1dvMax, h1dvMin;
    
    void eval(const double & av_x, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par);
    
    vector<double> eval_SB(const double & av_x, const double & Q2);
    
    void plot(const std::string outname_in, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par);
    
    void plot_SB(const std::string outname_in, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par);
    
    void getPolyNorm(const double & A, const double & B);
    
    double maximum(const double & a, const double & b);
    
    double maximum(const double & a, const double & b, const double & c);
     
};


}