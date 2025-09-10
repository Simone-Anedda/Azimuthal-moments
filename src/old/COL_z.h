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

#include<LHAPDF/LHAPDF.h>
#include "dss_wrapper.h"
#include "FCN.h"
#include "FragFunct.h"

//HOPPET
#include <hoppet_v1.h>
#include <hoppet_v1_collins.h>
#include <hoppet_v1_collins2.h>

    
namespace COLZ{

    using namespace std;
    
class getCOL_z{

 public:
     
    getCOL_z(const int & model_in, const int & widths_in, const int & evo_in){
           
        col_widths = widths_in;
        
        col_model = model_in;
    
        col_evo = evo_in;
        
        COL_z.resize(13);
        
        max_COL_z.resize(13);
        
        min_COL_z.resize(13);
        
        Col_fav = 0;
        
        Col_dis = 0;
        
        N_fav = 0;
        
        N_dis = 0;
        
        factor = 0;
        
        PolyNorm = 0;
        
        hoppet_col = false;
        
        for(int i = 0; i < max_COL_z.size(); i++){
            max_COL_z[i] = -pow(1,-10);            
            min_COL_z[i] = pow(1,-10);
        }
        
    };
   
    ~getCOL_z(){};
         
    int col_model;
    
    int col_widths;
    
    int col_evo;
    
    bool hoppet_col;

    vector<double> COL_z, max_COL_z, min_COL_z;
    
    double factor, PolyNorm;
    
    double Col_fav, Col_dis, N_fav, N_dis;
    
    void eval(const double & av_z, const double & Q2, const int & charge_in, const int & model_in, const int & evo_in, const std::vector<double> & FF_in, const std::vector<double> &par);
    
    void plot(const std::string outname_in, const double & Q2, const int & charge_in, const int & model_in, const int & evo_in, const std::vector<double> &par);
    
    void getPolyNorm(const double & A, const double & B);
    
    double maximum(const double & a, const double & b);
    
    double maximum(const double & a, const double & b, const double & c);
     
};


}
