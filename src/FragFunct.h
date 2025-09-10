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

#include <LHAPDF/LHAPDF.h>
#include "dss_wrapper.h"

namespace FRAG{

class FF{
    
    public:
        
        FF(){
            
            theFF.resize(13);
            tmpFF.resize(13);
            tmpFF2.resize(13);
        }
        
        FF(const std::string FFset_in, const std::string FForder_in){
            
            FFset = FFset_in;
            
            FForder = FForder_in;
            
            if(FForder == "LO"){
                
                FForder_i = 0;
                iset_ff = 1;
            }
            else if(FForder == "NLO"){
                
                FForder_i = 1;
                iset_ff = 2;
            }
            
            if(FFset == "Kretzer"){
                
                FForder = "LO";
                FForder_i = 0;
                
            }
            
            theFF.resize(13);
            tmpFF.resize(13);
            tmpFF2.resize(13);
        };
   
        ~FF(){};

        std::vector<double> theFF, tmpFF, tmpFF2;
        
        std::string FFset, FForder;
        
        int FForder_i, iset_ff;

        bool useLHAPDF;

        LHAPDF::PDF* ff;
        LHAPDF::PDF* ff_pip;
        LHAPDF::PDF* ff_pim;
        LHAPDF::PDF* ff_kap;
        LHAPDF::PDF* ff_kam;
        LHAPDF::PDF* ff_prp;
        LHAPDF::PDF* ff_prm;
        LHAPDF::PDF* ff_pisum;
        LHAPDF::PDF* ff_kasum;
        LHAPDF::PDF* ff_prsum;

        void use_LHAPDF(bool use_LHAPDF_in);
        
        void set_FFset(const std::string FFset_in);

        void set_LHAPDF_FFset(map<string, string> FFsets);
        
        void set_FForder(const std::string FForder_in);
        
        void set_FF(const std::string FFset_in, const std::string FForder_in);
     
        void FF_eval(const int & hadron_in, const int & charge_in, const double & av_z_in, const double & Q2_in);

        void FF_plot(const std::string frag_name, const int & hadron_in, const int & charge_in, const double & Q2_in);
    
};





}
