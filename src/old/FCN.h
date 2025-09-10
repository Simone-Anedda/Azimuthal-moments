//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string>
#include <cstring>
#include <gsl/gsl_math.h>
#include <vector>
#include <mutex>          // std::mutex for locking threads. We use it for Minuit.

//Minuit libraries
#include "Minuit2/Minuit2/FCNBase.h"
#include "Minuit2/Minuit2/FunctionMinimum.h"
#include "Minuit2/Minuit2/MnMigrad.h"
#include "Minuit2/Minuit2/MnUserParameters.h"
#include "Minuit2/Minuit2/MnPrint.h"

//CUBA
#include <cuba.h>
#include "mycubasettings.h"
//HOPPET
#include <hoppet_v1.h>
#include <hoppet_v1_collins.h>
#include <hoppet_v1_collins2.h>
#include "hoppet_settings.h"
//Collinear PDFs
#include "CollPDF.h"
//Fragmentation Functions
#include "FragFunct.h"

#include "Data.h"

#include "TRANSV_x.h"
#include "COL_z.h"

#include "transv_wrapper.h"


//following from chi2 prob tables for std =2 , i.e. 95% confidence level
#define _2_STD_ERRDEF_0PAR      1.00             //standard for Minuit. If NumberOfParameters is not set to something diff	
#define _2_STD_ERRDEF_1PAR      4.00 	       //than initial value of 0, then this default value of 1 is to be taken.
#define _2_STD_ERRDEF_2PAR      6.18 	
#define _2_STD_ERRDEF_3PAR      8.02 	
#define _2_STD_ERRDEF_4PAR      9.72 	
#define _2_STD_ERRDEF_5PAR      11.31	
#define _2_STD_ERRDEF_6PAR      12.85	
#define _2_STD_ERRDEF_7PAR      14.34	
#define _2_STD_ERRDEF_8PAR      15.79	
#define _2_STD_ERRDEF_9PAR      17.21	
#define _2_STD_ERRDEF_10PAR     18.61	
#define _2_STD_ERRDEF_11PAR     19.99	
#define _2_STD_ERRDEF_12PAR     21.35	
#define _2_STD_ERRDEF_13PAR     22.69	
#define _2_STD_ERRDEF_14PAR     24.03	
#define _2_STD_ERRDEF_15PAR     25.34	
#define _2_STD_ERRDEF_16PAR     26.65	
#define _2_STD_ERRDEF_17PAR     27.95	
#define _2_STD_ERRDEF_18PAR     29.24	
#define _2_STD_ERRDEF_19PAR     30.52	
#define _2_STD_ERRDEF_20PAR     31.80	

///for reference only, extended table

// sigma 	1σ 	1.28 	1.64 	1.96 	2σ 	2.58 	3σ 	3.29 	4σ
// conf int % 	68.3% 	80% 	90% 	95% 	95.45% 	99% 	99.73% 	99.9% 	99.99%
// p-value 	0.317 	0.20 	0.10 	0.05 	0.0455 	0.01 	0.0027 	0.001 	0.00006
// chi2(k=1) 	1.00 	1.64 	2.71 	3.84 	4.00 	6.63 	9.00 	10.83 	16.00
// chi2(k=2) 	2.30 	3.22 	4.61 	5.99 	6.18 	9.21 	11.83 	13.82 	19.33
// chi2(k=3) 	3.53 	4.64 	6.25 	7.81 	8.02 	11.34 	14.16 	16.27 	22.06
// chi2(k=4) 	4.72 	5.99 	7.78 	9.49 	9.72 	13.28 	16.25 	18.47 	24.50
// chi2(k=5) 	5.89 	7.29 	9.24 	11.07 	11.31 	15.09 	18.21 	20.52 	26.77
// chi2(k=6) 	7.04 	8.56 	10.64 	12.59 	12.85 	16.81 	20.06 	22.46 	28.91
// chi2(k=7) 	8.18 	9.80 	12.02 	14.07 	14.34 	18.48 	21.85 	24.32 	30.96
// chi2(k=8) 	9.30 	11.03 	13.36 	15.51 	15.79 	20.09 	23.57 	26.12 	32.93
// chi2(k=9) 	10.42 	12.24 	14.68 	16.92 	17.21 	21.67 	25.26 	27.88 	34.85
// chi2(k=10) 	11.54 	13.44 	15.99 	18.31 	18.61 	23.21 	26.90 	29.59 	36.72
// chi2(k=11) 	12.64 	14.63 	17.28 	19.68 	19.99 	24.72 	28.51 	31.26 	38.54
// chi2(k=12) 	13.74 	15.81 	18.55 	21.03 	21.35 	26.22 	30.10 	32.91 	40.33
// chi2(k=13) 	14.84 	16.98 	19.81 	22.36 	22.69 	27.69 	31.66 	34.53 	42.09
// chi2(k=14) 	15.94 	18.15 	21.06 	23.68 	24.03 	29.14 	33.20 	36.12 	43.82
// chi2(k=15) 	17.03 	19.31 	22.31 	25.00 	25.34 	30.58 	34.71 	37.70 	45.52
// chi2(k=16) 	18.11 	20.47 	23.54 	26.30 	26.65 	32.00 	36.22 	39.25 	47.20
// chi2(k=17) 	19.20 	21.61 	24.77 	27.59 	27.95 	33.41 	37.70 	40.79 	48.86
// chi2(k=18) 	20.28 	22.76 	25.99 	28.87 	29.24 	34.81 	39.17 	42.31 	50.50
// chi2(k=19) 	21.36 	23.90 	27.20 	30.14 	30.52 	36.19 	40.63 	43.82 	52.13
// chi2(k=20) 	22.44 	25.04 	28.41 	31.41 	31.80 	37.57 	42.08 	45.31 	53.73



namespace ROOT {

namespace Minuit2{


using namespace std;
  
class FCN : public FCNBase{

 public:

     
     FCN(){
       
       theErrorDef.push_back(_2_STD_ERRDEF_0PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_1PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_2PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_3PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_4PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_5PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_6PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_7PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_8PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_9PAR );
       theErrorDef.push_back(_2_STD_ERRDEF_10PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_11PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_12PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_13PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_14PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_15PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_16PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_17PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_18PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_19PAR);
       theErrorDef.push_back(_2_STD_ERRDEF_20PAR);

       
       
       
       set_NumberOfParameters(0);

       
       
    };
    
    ~FCN(){};
    
    
    ///to set error
    void set_NumberOfParameters(int paranum){
      
     
      NumberOfParameters = paranum;
      
           
    }
    

    ///for Minuit
    double virtual  Up() const {
     
        double output = 0.0;
     
        output = theErrorDef[NumberOfParameters];
        
//         cout<<"Up output: "<<output<<endl;
     
        return output;
     
    }
//      virtual double operator ()(const std::vector<double>&) const{}; WE HAD WRITTEN {}, WHICH MEANS WE WERE DEFINING THE FUNCTION, WHICH IS WRONG.
    virtual double operator ()(const std::vector<double>&par) const;//MODIFIED: LINE HAS NO {}, AS IT SHOULD BE FOR A DECLARATION

    
    
 private:

   
   vector<double> theErrorDef;
   
   int		 NumberOfParameters;
 
   
  
 
  
  
  
   
};


    void plot_transversity(const string & outname, const double & Q2);




}

}
