#ifndef _getData_H_
#define _getData_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>

#include "rapidcsv.h"
    
namespace DATA{

    using namespace std;
    
class getData{

 public:
   
   getData(const char* fname_in, int hadron_in, int charge_in, int target_in, int process_in, int epem_sign_in){
     
     sprintf(fname,"%s",fname_in);
     
     getNoRows(fname_in);
     
     getNoCols(fname_in);
     
     readData(fname_in);
     
     string name(fname_in);
     
//      CSV = new rapidcsv::Document;
    
//      CSV->Load(name);
     
     vector<double> (Nrow, 0).swap(theory);
     
     vector<double> (Nrow, 0).swap(theory_0);
     
     vector<double> (Nrow, 0).swap(theory_sigma);
     
//      vector<double> (Nrow, 0).swap(theory_plus);
     
//      vector<double> (Nrow, 0).swap(theory_minus);
     
     vector<double> (Nrow, -1e10).swap(theory_max);
     
     vector<double> (Nrow, +1e10).swap(theory_min);
     
//      vector<double> (Nrow, 0).swap(theory_hesse_max);
     
//      vector<double> (Nrow, 0).swap(theory_hesse_min);     
     
     vector<double> (Nrow, 0).swap(chi2);
     
     vector<double> (Nrow, 0).swap(xvals);
     
     valid.resize(Nrow);
     
     getCol(5); xvals.swap(col_temp);
        
     getCol(11); meas.swap(col_temp);
     
     getCol(13); errs.swap(col_temp);
     
     check_pts();  //initialize valid and check npts 
     
     eval_nbins();
          
     tot_chi2 = 0;
     
     hadron = hadron_in;
     
     charge = charge_in;
          
     target = target_in;
     
     process = process_in;
     
     epem_sign = epem_sign_in;
     
  };
   
   ~getData(){};
     
    //initialize
   
     void readData(const char *name);
     
     void getNoCols(const char *name);
     
     void getNoRows(const char *name);
     
     void getCol(int ncol);
     
     void getRow(int nrow);
     
     //print
     
     void print_bin(int mode);
     
     void eval_npts();
     
     void check_pts();
     
     void set_xrange(double min, double max);
    
     void set_zrange(double min, double max); 
     
     void set_ptrange(double min, double max); 
     
     void set_Q2range(double min, double max); 
     
     void eval_chi2();
     
     void print();

     void print(const string filename, const string results, const int resdirLen, const int process);  //to print on a specific file

     void eval_nbins();
     
     void eval_nbins_z1(const int process);
     
     void split_in_bins(const string filename, const string results, const int resdirLen, const int process); //to obtain different bins to plot
     
     void get_xvals();
     
     int Nrow;
     
     int Ncol;
     
     int Nbins;
     
     char fname[500], filename[500];
     
     string results;
     
     vector<double> everything;
     
     vector<double> bin;
     
     vector<double> col_temp;
     
     vector<double> theory;
     
     vector<double> chi2;
    
     vector<double> errs;
    
     vector<double> meas; 
     
     vector<bool> valid;
     
     vector<double> xvals;
     
     vector<double> theory_max, theory_min, theory_0, theory_sigma;
     
     double tot_chi2;
     
     int npts;
     
     int target;
     
     int process;  //to distinguish between SIDIS and e+e-
     int epem_sign;  //to distinguish between UL (epem_sign=1) and UC (epem_sign=2) asymmetries
     
     //for fDSS    
     int hadron;
     int charge;
     
     
};


}

#endif
