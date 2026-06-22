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

#include "dss_wrapper.h"
//#include "FCN.h"
#include "FragFunct.h"

//HOPPET
#include "/st100-gr4/codes/hoppet-1.1.5-modified_cf-build/include/hoppet_v1_collins.h"
#include "/st100-gr4/codes/hoppet-1.1.5-modified_cff-build/include/hoppet_v1_collins2.h"

    
namespace COL{

    using namespace std;
    
class COLLINS{

 public:
     
    COLLINS(){
        // (const int & model_in, const int & widths_in, const int & evo_in){
           
        // col_widths = widths_in;
        
        // col_model = model_in;
    
        // col_evo = evo_in;

        vector<double> (13, 0).swap(COL_z);
        
        vector<double> (13, 0).swap(COL_first);
        
        vector<double> (13, 0).swap(COL_an_pow);
        
        max_COL_z.resize(13);
        
        min_COL_z.resize(13);
        
        Col_fav = 0;
        
        Col_dis = 0;
        
        N_fav = 0;
        
        N_dis = 0;
        
        factor = 0;
        
        PolyNorm = 0;
        
        // hoppet_col = false;
        
        for(int i = 0; i < max_COL_z.size(); i++){

            max_COL_z[i] = -pow(1,-10);            
            min_COL_z[i] = pow(1,-10);
        }
        
    };
   
    ~COLLINS(){};
         
    // int col_model;
    
    // int col_widths;
    
    // int col_evo;
    
    // bool hoppet_col;
    
    bool hasFF;

    std::string model, widths, evo, header;

    vector<std::string> parnames;

    vector<double> COL_z, max_COL_z, min_COL_z, COL_first, COL_an_pow;
    
    double factor, PolyNorm;
    
    double Col_fav, Col_dis, N_fav, N_dis;
    
    void eval(const double & av_z, const double & Q2, const int & charge_in, const std::vector<double> & FF_in, const std::vector<double> &par);
    
    void eval_first_moment(const double & av_z, const double & Q2, const int & charge_in, const std::vector<double> & FF_in, const std::vector<double> &par);
    
    void eval_analysing_power(const double & av_z, const double & Q2, const int & charge_in, const std::vector<double> & FF_in, const std::vector<double> &par);
    
    void plot(const std::string outname_in, const double & Q2, const int & charge_in, const std::vector<double> &par);

    void plot_first_moment(const std::string outname_in, const double & Q2, const int & charge_in, const std::vector<double> &par);
    
    void plot_analysing_power(const std::string outname_in, const double & Q2, const int & charge_in, const std::vector<double> &par);
        
    void has_FF(bool hasFF_in);

    void set_model(const string model_in);

    void set_header();

    void set_widths(const string widths_in);

    void set_evolution(const string evo_in);

    void getPolyNorm(const double & A, const double & B);
    
    double maximum(const double & a, const double & b);
    
    double maximum(const double & a, const double & b, const double & c);
     
};


}
