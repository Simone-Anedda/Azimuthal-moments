#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>

#include "getData.h"


extern bool hessian_flag, MC_flag;

namespace DATA{

using namespace std;

void getData::readData(const char *name){
    
    vector<double> data;
    ifstream input (name);
    string lineData;
    
    //check in file opening
    if(!input) cerr<<"Error in opening your file "<<name<<" (check filename)\n";

    while(getline(input, lineData,'\n')){
        if (lineData[0] == '#') continue;
        char d[10];
        stringstream lineStream(lineData);
        
        while (lineStream >> d){
            double num = atof(d);
            data.push_back(num);           
        }        
    }
    
    everything.swap(data);
}


void getData::getNoCols(const char *name){
    
    int cols=0;
    char n[6];
    ifstream input(name);
    
    while (input.peek() != '\n' && input >> n){
            if (input.peek() == '#') continue;
            cols++;
    }
    
    Ncol=cols-1;
    
}


void getData::getNoRows(const char *name){
    
    int rows=0;
    ifstream input (name);
    string lineData;
    
    //check in file opening
    if(!input) cerr<<"Error in opening "<<fname<<" (check filename)\n";
    
    while(getline(input, lineData,'\n')){
    if (lineData[0] == '#') continue;
        rows++;
    }    
    
    Nrow=rows;
}


void getData::getCol(int ncol){

    //jogh: modify function to read from 'everything'

    int m = ncol-1;
    vector<double> col;
       if(ncol>Ncol){ 
        cerr<<"ERROR: File "<<fname<<" has only "<< Ncol <<" columns\n"; 
        cout<<"Choose a different column\n";
        cin >> ncol;
        m = ncol - 1;
    }
        
    col.push_back(everything[m]);
        for (int i=1; i<Nrow; i++){                  //LOOP STRUCTURE SHOULD BE OK
            col.push_back(everything[m+i*Ncol]);
        }
        
        col_temp.swap(col);

}        


void getData::getRow(int nrow){
    
  //jogh: modify function to read from 'everything'

    int n = nrow-1;
    vector<double> row;
    if(nrow>Nrow){
        cerr<<"ERROR: File "<<fname<<" has only "<< Nrow <<" rows\n"; 
        cout<<"Choose a different row\n";
        cin >> nrow;
        n = nrow - 1;
    }

    for (int i=(Ncol*n); i<(Ncol*n+Ncol); i++){                   //LOOP STRUCTURE SHOULD BE OK
            row.push_back(everything[i]);       
        }
        
        bin.swap(row);//jogh
        
}


void getData::print_bin(int mode){

    for (int i = 0; i < Ncol; i++) printf("%.4e\t", bin[i]);  
    if (mode != 1) printf("\n");
}


void getData::check_pts(){
    
    for(int i = 0; i < Nrow; i++){
        if(errs[i] != 0) valid[i] = true;
        else valid[i] = false;
    }
    
    eval_npts();
}

void getData::set_xrange(double min, double max){

    getCol(5);
    
    for(int i = 0; i < Nrow; i++){
        bool out_of_range = col_temp[i] <= min || col_temp[i] >= max;
        if(out_of_range) valid[i] = false;
//         else valid[i]=true;
    }

    eval_npts();
}


void getData::set_zrange(double min, double max){

    getCol(7);
    
    for(int i = 0; i < Nrow; i++){
        bool out_of_range = col_temp[i] <= min || col_temp[i] >= max;
        if(out_of_range) valid[i] = false;
//         else valid[i]=true;
    }

    eval_npts();
}

void getData::set_ptrange(double min, double max){

    getCol(8);
    
    for(int i=0; i<Nrow; i++){
        bool out_of_range = col_temp[i] <= min || col_temp[i] >= max;
        if(out_of_range) valid[i]=false;
//         else valid[i]=true;
    }

    eval_npts();
}

void getData::set_Q2range(double min, double max){

    getCol(10);
    
    for(int i=0; i<Nrow; i++){
        bool out_of_range = col_temp[i] <= min || col_temp[i] >= max;
        if(out_of_range) valid[i]=false;
//         else valid[i]=true;
    }

    eval_npts();
}


void getData::eval_npts(){

    npts=0;
    for(int i=0; i<valid.size();i++){
        if(valid[i]) npts++;
    }
}


void getData::get_xvals(){
    
    vector<double> xmin, xmax, av_x;
    getCol(3); xmin.swap(col_temp);
    getCol(4); xmax.swap(col_temp);
    
    for(int i=0; i<Nrow; i++){

        av_x[i] = (xmin[i]+xmax[i])/2;
    }

    xvals.swap(av_x);
}


void getData::eval_chi2(){
    
    tot_chi2=0;
    
    for(int i=0; i<Nrow; i++){
        
        if(valid[i]){
            chi2[i]=pow(meas[i]-theory[i],2)/pow(errs[i],2);
        }
        else chi2[i]=0;
    
        tot_chi2+=chi2[i];
    }
    
}


void getData::print(){
    
     for(int bin_i=1; bin_i<=Nrow; bin_i++){ 
         int i=bin_i-1;
        if(valid[i]){
        getRow(bin_i); print_bin(1);
//         printf("%.4e\t%.4e\n",theory[i],chi2[i]);
        printf("%.4e\t%.4e\t%.4e\t%.4e\n",theory[i],chi2[i],theory_min[i],theory_max[i]);
        }
    }
    
}


void getData::print(const string filename, const string results, const int resdirLen, const int process){
    
    ofstream file (filename);
    
    double chisq=0;
    
    //SIDIS
    if(process==0) file<<"#Q2min\t\tQ2max\t\tvarmin\t\tvarmax\t\t<x>\t\t<y>\t\t<z>\t\t<pT>\t\t<W>\t\t<Q^2>\t\tAN\t\terstat\t\tertot\t\ttheory\t\tchi2\t\tmin_th\t\tmax_th\n#\n";
    //e+e-(PT)
    if(process==1) file<<"#Q2min\t\tQ2max\t\tpt0min\t\tpt0max\t\t<z1>\t\t<z2>\t\tA0th\t\t<pT>\t\t<W>\t\t<Q^2>\t\tA0UL(C)\t\terstat\t\tertot\t\ttheory\t\tchi2\t\tmin_th\t\tmax_th\n#\n";
    //e+e-(z1z2)
    if(process==2) file<<"#Q2min\t\tQ2max\t\tz1min\t\tz1max\t\t<z1>\t\tz2min\t\tz2max\t\t<z2>\t\tA0th\t\t<Q^2>\t\tA0UL(C)\t\terstat\t\tertot\t\ttheory\t\tchi2\t\tmin_th\t\tmax_th\n#\n";
    
    cout.precision(4);
    
    for(int bin_i=1; bin_i<=Nrow; bin_i++){ 
         
         int i=bin_i-1;
        
         if(valid[i]){
            getRow(bin_i); //print_bin(1);
            for(int j=0;j<bin.size();j++){
            
                file<<scientific<<bin[j]<<"\t";
            }
            chisq+=chi2[i];
//             file<<scientific<<theory[i]<<"\t"<<chi2[i]<<"\t"<<theory_min[i]<<"\t"<<theory_max[i]<<endl;
            file<<scientific<<theory_0[i]<<"\t"<<chi2[i]<<"\t"<<theory_sigma[i]<<"\t"<<theory_sigma[i]<<endl;
        }
    }
    eval_npts();
    file<<scientific<<"\n#tot_chi2 = "<<chisq<<"\t tot_chi2/ndata = "<<chisq/npts<<endl;
    if(process==0 || process==1) eval_nbins();
    if(process==2) eval_nbins_z1(process);
    split_in_bins(filename,results,resdirLen,process);

    file.close();

}


void getData::eval_nbins(){

    getCol(1);
    Nbins=1;
    for(int j=1; j<col_temp.size()-1; j++){
        
        if(col_temp[j+1]!=col_temp[j]) Nbins++;
    }
    
}

void getData::eval_nbins_z1(const int process){

    if(process==0) cerr<<"z1 binning not used for PT-dependent SIDIS processes\n";
    if(process==1) cerr<<"z1 binning not used for PT-dependent e+e- processes\n";

    if(process==2){
        getCol(3);
        Nbins=1;
        for(int j=1; j<col_temp.size()-1; j++){
            
            if(col_temp[j+1]!=col_temp[j]) Nbins++;
        }
    }
    
}


void getData::split_in_bins(const string filename, const string results, const int resdirLen, const int process){

    int bin_i=1, k=0, pts_in_bin=0, npts=0;
    vector<double> check, tmp;
    double chi2_set=0, chisq=0;
    
    cout.precision(4);
    
    string fname = filename, resname = filename, index, outname;
    
    if(fname.size()>0) fname.resize(fname.size()-4); //to erase ".dat" from filename
    
    resname.erase(resname.begin(), resname.begin()+resdirLen+1);
    
    ofstream res(results, std::ios::app); 
    res<<endl;
    
    while(bin_i<=Nrow && k<=Nbins){
        stringstream a;
        a << k+1;
        index=a.str();
        outname=fname+"_"+index+".dat"; 
        
        ofstream out(outname.c_str());
        if(!out.is_open()) cerr<<"Error in opening your file"<<outname<<"\n";
        
        while(out.is_open()){
//             if(out.is_open()) cout<<"Your file "<<outname<<" is open\n";
            if(bin_i<Nrow){
                getRow(bin_i+1); 
                check.swap(bin); 
            }
            getRow(bin_i); 
            tmp.swap(bin);
            if(valid[bin_i-1]){
                
                for(int j=0;j<tmp.size();j++)  out<<scientific<<tmp[j]<<"\t";
                
//                 out<<theory[bin_i-1]<<"\t"<<chi2[bin_i-1]<<"\t"<<theory_min[bin_i-1]<<"\t"<<theory_max[bin_i-1]<<"\n";
                out<<theory_0[bin_i-1]<<"\t"<<chi2[bin_i-1]<<"\t"<<theory_sigma[bin_i-1]<<"\t"<<theory_sigma[bin_i-1]<<"\n";
                
                chi2_set+=chi2[bin_i-1]; //evaluating chi2 for the bin
                chisq+=chi2[bin_i-1]; npts++;
                pts_in_bin++;   //counting n_pts in bin
            
            }
            bin_i++;     
            
            if(process!=2){
                if(tmp[0]!=check[0] || bin_i==Nrow+1){ //Q^2 binning
                    res<<resname+"_"+index<<"\t"<<scientific<<chi2_set<<"\t"<<pts_in_bin<<"\t"<<scientific<<chi2_set/pts_in_bin<<endl; 
                    chi2_set=0;
                    pts_in_bin=0;
                    out.close();
//                  if(!out.is_open()) cout<<"Your file "<<outname<<" is closed\n";
                    out.clear();                
                }
            }
            if(process==2){
                if(tmp[3]!=check[3] || bin_i==Nrow+1){ //z1 binning
                    res<<resname+"_"+index<<"\t"<<scientific<<chi2_set<<"\t"<<pts_in_bin<<"\t"<<scientific<<chi2_set/pts_in_bin<<endl; 
                    chi2_set=0;
                    pts_in_bin=0;
                    out.close();
//                  if(!out.is_open()) cout<<"Your file "<<outname<<" is closed\n";
                    out.clear();                
                }
            }

        }        
        
        k++;
    }
        res<<"tot_chi2_set = "<<chisq<<"\t\tnpts = "<<npts<<"\ntot_chi2_set/npts = "<<chisq/npts<<endl;
        chisq=0; npts=0;
        res.close();
}

    
    
}
