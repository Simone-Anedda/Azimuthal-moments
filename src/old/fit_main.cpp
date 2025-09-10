#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>

//strings for data files
#include "datanames_Collins_csv.h"

//Collinear distributions
#include "CollPDF.h"
//Fragmentation Functions
#include "FragFunct.h"

//HOPPET
#include <hoppet_v1.h>
#include <hoppet_v1_collins.h>
#include <hoppet_v1_collins2.h>
#include "hoppet_settings.h"

//Minuit libraries
#include "Minuit2/Minuit2/FunctionMinimum.h"
#include "Minuit2/Minuit2/MnMigrad.h"
#include "Minuit2/Minuit2/MnMinos.h"
#include "Minuit2/Minuit2/MinosError.h"
#include "Minuit2/Minuit2/MnUserParameters.h"
#include "Minuit2/Minuit2/MnPrint.h"

#include "FCN.h"
#include "TRANSV_x.h"
#include "COL_z.h"
#include "CovMat.h"
#include "ParVec.h"

#include "Data.h"

using namespace ROOT::Minuit2;
    
    //format of DATA::Data data(fname, hadron, charge, target, process, epem_sign)
    // hadron: [1] pion, [2] kaon, [3] proton, [4] hadron; - DSS Parameter
    // charge: [0] neutral, [1] positive, [-1] negative - DSS Parameter
    // target: [0] proton, [1] neutron, [2] deuteron - to use the proper PDFs
    // process: "SIDIS", "e+e-(PT)", "e+e-(z1z2)"
    // epem_sign: "none" (for SIDIS), "A0UL" "A0UC" "A12UL" "A12UC" (for e+e-)
        
    DATA::Data\
    
    JLab11_Npp_x(jlab1, "pion", "pos", "neutron", "SIDIS", "none"), JLab11_Npm_x(jlab1, "pion", "neg", "neutron", "SIDIS", "none"),\
    
    Herm10_pip_x(herm1, "pion", "pos", "proton", "SIDIS", "none"), Herm10_pip_z(herm2, "pion", "pos", "proton", "SIDIS", "none"), Herm10_pip_pt(herm3, "pion", "pos", "proton", "SIDIS", "none"),\
    Herm10_pim_x(herm4, "pion", "neg", "proton", "SIDIS", "none"), Herm10_pim_z(herm5, "pion", "neg", "proton", "SIDIS", "none"), Herm10_pim_pt(herm6, "pion", "neg", "proton", "SIDIS", "none"),\
    
    Herm20_pip(herm203d_1, "pion", "pos", "proton", "SIDIS", "none"), Herm20_pim(herm203d_2, "pion", "neg", "proton", "SIDIS", "none"),\
    
    Comp12_Ppp_x(comp12_1, "pion", "pos", "proton", "SIDIS", "none"), Comp12_Ppp_z(comp12_2, "pion", "pos", "proton", "SIDIS", "none"), Comp12_Ppp_pt(comp12_3, "pion", "pos", "proton", "SIDIS", "none"),\
    Comp12_Ppm_x(comp12_4, "pion", "neg", "proton", "SIDIS", "none"), Comp12_Ppm_z(comp12_5, "pion", "neg", "proton","SIDIS", "none"), Comp12_Ppm_pt(comp12_6, "pion", "neg", "proton", "SIDIS", "none"),\
    Comp12_Dpp_x(comp12_7, "pion", "pos", "deuteron", "SIDIS", "none"), Comp12_Dpp_z(comp12_8, "pion", "pos", "deuteron", "SIDIS", "none"), Comp12_Dpp_pt(comp12_9, "pion", "pos", "deuteron", "SIDIS", "none"),\
    Comp12_Dpm_x(comp12_10, "pion", "neg", "deuteron", "SIDIS", "none"), Comp12_Dpm_z(comp12_11, "pion", "neg", "deuteron", "SIDIS", "none"), Comp12_Dpm_pt(comp12_12, "pion", "neg", "deuteron", "SIDIS", "none"),\
    
    Comp14_Ppp_x(comp14_1, "pion", "pos", "proton", "SIDIS", "none"), Comp14_Ppp_z(comp14_2, "pion", "pos", "proton", "SIDIS", "none"), Comp14_Ppp_pt(comp14_3, "pion", "pos", "proton", "SIDIS", "none"),\
    Comp14_Ppm_x(comp14_4, "pion", "neg", "proton", "SIDIS", "none"), Comp14_Ppm_z(comp14_5, "pion", "neg", "proton", "SIDIS", "none"), Comp14_Ppm_pt(comp14_6, "pion", "neg", "proton", "SIDIS", "none"),\
    
    Comp17_01hp_x(comp17_1, "hadron", "pos", "proton", "SIDIS", "none"), Comp17_01hp_z(comp17_2, "hadron", "pos", "proton", "SIDIS", "none"), Comp17_01hp_pt(comp17_3, "hadron", "pos", "proton", "SIDIS", "none"),\
    Comp17_02hp_x(comp17_4, "hadron", "pos", "proton", "SIDIS", "none"), Comp17_02hp_z(comp17_5, "hadron", "pos", "proton", "SIDIS", "none"), Comp17_02hp_pt(comp17_6, "hadron", "pos", "proton", "SIDIS", "none"),\
    Comp17_01hm_x(comp17_7, "hadron", "neg", "proton", "SIDIS", "none"), Comp17_01hm_z(comp17_8, "hadron", "neg", "proton", "SIDIS", "none"), Comp17_01hm_pt(comp17_9, "hadron", "neg", "proton", "SIDIS", "none"),\
    Comp17_02hm_x(comp17_10, "hadron", "neg", "proton", "SIDIS", "none"), Comp17_02hm_z(comp17_11, "hadron", "neg", "proton", "SIDIS", "none"), Comp17_02hm_pt(comp17_12, "hadron", "neg", "proton", "SIDIS", "none"),\
    
    Belle12_A0UL_z1z2(belle2012_1, "pion", "pos", "proton", "e+e-(z1z2)", "A0UL"), Belle12_A0UC_z1z2(belle2012_2, "pion", "pos", "proton", "e+e-(z1z2)", "A0UC"),\
    Belle12_A12UL_z1z2(belle2012_3, "pion", "pos", "proton", "e+e-(z1z2)", "A12UL"), Belle12_A12UC_z1z2(belle2012_4, "pion", "pos", "proton", "e+e-(z1z2)", "A12UC"),\
    
    Babar13_A0UL_z1z2(babar2013_1, "pion", "pos", "proton", "e+e-(z1z2)", "A0UL"), Babar13_A0UC_z1z2(babar2013_2, "pion", "pos", "proton", "e+e-(z1z2)", "A0UC"),\
    Babar13_A12UL_z1z2(babar2013_3, "pion", "pos", "proton", "e+e-(z1z2)", "A12UL"), Babar13_A12UC_z1z2(babar2013_4, "pion", "pos", "proton", "e+e-(z1z2)", "A12UC"),\
    
    Babar14_A0UL_PT(babar2014_1, "pion", "pos", "proton", "e+e-(PT)", "A0UL"), Babar14_A0UC_PT(babar2014_2, "pion", "pos", "proton", "e+e-(PT)", "A0UC"),\
    
    Babar15_A0UL_z1z2(babar2015_1, "pion", "pos", "proton", "e+e-(z1z2)", "A0UL"), Babar15_A0UC_z1z2(babar2015_2, "pion", "pos", "proton", "e+e-(z1z2)", "A0UC"),\
    Babar15_A12UL_z1z2(babar2015_3, "pion", "pos", "proton", "e+e-(z1z2)", "A12UL"), Babar15_A12UC_z1z2(babar2015_4, "pion", "pos", "proton", "e+e-(z1z2)", "A12UC"),\
    
    BESIII16_A0UL_PT(besIII2016_1, "pion", "pos", "proton", "e+e-(PT)", "A0UL"), BESIII16_A0UC_PT(besIII2016_2, "pion", "pos", "proton", "e+e-(PT)", "A0UC");

    //format of TRANSVX::getTRANSV_x myTransv(model,widths,evo) and of COLZ::getCOL_z myCol(model,widths,evo)
    // model: [0] 2015 fit [1] 2015 fit - Bernstein polynomials for Collins [2] 2018 fit - REWEIGHTING - unbiased
    // widhts (<k^2_\perp>, <p^2_\perp>): [0] (0.57,0.12) & (0.60, 0.20) (fit LO 2014 [1] (0.25,0.20) [2] (0.58, 0.12) & (0.52, 0.20) (fit NLO)
    // evo: [0] no evolution, [1] HOPPET evolution
    
    TRANSVX::getTRANSV_x myTransv(0, 0, 1); //note: remember to set model parameter equal for transversity and Collins functions
    COLZ::getCOL_z myCol(0, 0, 0);
    
    FRAG::FF myFF;

    COLLPDF::CollPDF myPDF;

    // int pdfID = pdf->lhapdfID();
    
    // int iset_ff = 2; //SET THIS INT TO 1 IF YOU WANT LO FFs SET, 2 FOR NLO FFs SET
    
    bool diffsame_flag = true; //SET THIS FLAG TRUE IF YOU WANT TO USE DIFFERENT WIDTHS FOR COMPASS AND HERMES
    
    bool fixed_Q2 = false, dssv_init = true, grv_flag = false;

    bool SB_flag; //to be set false if one does not want to use limits on NT parameters for transversity (see input file)
    
    bool hessian_flag = false, MC_flag = false;
    
    int Npoints = 0, Ndof, NfixedPars = 0;
    
    time_t now = time(0);
    
    tm *ltm = localtime(&now);
    
    string year = to_string(1900 + ltm->tm_year), month = to_string(1 + ltm->tm_mon), day = to_string(ltm->tm_mday);
    
    string today = day + "-" + month + "-" + year, Sivers_print_name;

int main(int argc,char *argv[]){

    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double zmin = atof(argv[5]);
    double zmax = atof(argv[6]);
    double Q2min = atof(argv[8]);
    double Q2max = atof(argv[9]);
    const int thefactscheme =  atoi(argv[11]);
    
    stringstream resfile, lhapdfstring, ffstring, fforderstring, modelstring, SBstring, fixedQ2string;
    char *res_directory;
    
    resfile << argv[13];
    res_directory = argv[15];
    lhapdfstring << argv[17];
    ffstring << argv[19];
    fforderstring << argv[21];
    SBstring << argv[23];
    fixedQ2string <<argv[24];

    
    string results = resfile.str();
    string res_errs = results;
    string res_pars = results;
    string resdir = res_directory;
    string theLHAPDFset = lhapdfstring.str(); 
    string theFFset = ffstring.str(); 
    string theFForder = fforderstring.str(); 
    string theSBstring = SBstring.str();
    string thefixedQ2string = fixedQ2string.str();

    if(theSBstring == "yes") SB_flag = true;
    if(theSBstring == "no") SB_flag = false;

    if(thefixedQ2string == "yes") fixed_Q2 = true;
    if(thefixedQ2string == "no") fixed_Q2 = false;


    if(SB_flag){
    
        if(results.size() > 0) results.resize(results.size()-4);
        results = results + "-usingSB.dat";
        if(res_errs.size() > 0) res_errs.resize(res_errs.size()-4);
        res_errs = res_errs + "-usingSB.dat";
        if(res_pars.size() > 0) res_pars.resize(res_pars.size()-4);
        res_pars = res_pars + "-usingSB.dat";

        strcat(res_directory,"-usingSB");
        resdir += "-usingSB";
    
    }

    else if(!SB_flag){

        if(results.size() > 0) results.resize(results.size()-4);
        results = results + "-noSB.dat";
        if(res_errs.size() > 0) res_errs.resize(res_errs.size()-4);
        res_errs = res_errs + "-noSB.dat";
        if(res_pars.size() > 0) res_pars.resize(res_pars.size()-4);
        res_pars = res_pars + "-noSB.dat";

        strcat(res_directory,"-noSB");
        resdir += "-noSB";

    }
    
    myPDF.set_PDFset(theLHAPDFset);
    myFF.set_FF(theFFset, theFForder);
    
    cout<<"Collinear PDF set: "<<myPDF.pdf->set().description()<<endl;
    cout<<"LHAPDF ID: "<<myPDF.pdf->lhapdfID()<<", member: "<<myPDF.pdf->memberID()<<endl;
    cout<<"Collinear FF set: "<<myFF.FFset<<", order: "<<myFF.FForder<<endl<<endl;

    // cout<<"Collinear PDF set: "<<pdf->set().description()<<endl;
    // cout<<"LHAPDF ID: "<<pdf->lhapdfID()<<", member: "<<pdf->memberID()<<endl;
//     cout<<"Collinear FF set: "<<myFF.FFset<<", order: "<<myFF.FForder<<endl<<endl; 

    
    
    if(myTransv.transv_evo == 1){ 
        
        if(results.size() > 0) results.resize(results.size()-4);
        results = results + "-h1-evo.dat";
        if(res_errs.size() > 0) res_errs.resize(res_errs.size()-4);
        res_errs = res_errs + "-h1-evo.dat";
        if(res_pars.size() > 0) res_pars.resize(res_pars.size()-4);
        res_pars = res_pars + "-h1-evo.dat";

        strcat(res_directory,"-h1-evo");
//         cout<<"At transv hoppet: res_directory = "<<res_directory<<endl;
        resdir += "-h1-evo";
        
        //starting HOPPET and choosing alpha_s scheme
//         hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, intorder, factscheme);
        hoppetStartExtended(ymaxT, dy, QminT, QmaxT, dlnlnQ, nloop, intorder, factscheme);
        hoppetSetFFN(fixed_nf);
        hoppetSetPoleMassVFN(mcharm, mbottom, mtop);
//         hoppetPreEvolve(asQ0, Q0alphas, nloop, muR_Q, Q0pdf);

    }
    
    if(myCol.col_evo == 1){ 
        
        if(results.size() > 0) results.resize(results.size()-4);
        results = results + "-H1p-evo.dat";
        if(res_errs.size() > 0) res_errs.resize(res_errs.size()-4);
        res_errs = res_errs + "-H1p-evo.dat";
        if(res_pars.size() > 0) res_pars.resize(res_pars.size()-4);
        res_pars = res_pars + "-H1p-evo.dat";

//         strcat(res_directory,"-H1p-evo");
//         cout<<"At collins hoppet: res_directory = "<<res_directory<<endl;
//         resdir += "-H1P-evo";
        
        //starting HOPPET and choosing alpha_s scheme for Collins pi+ and pi-
//         hoppetStartExtendedcf(ymaxCol, dy, QminC, QmaxC, dlnlnQ, nloop, intorder, factscheme);
        hoppetStartExtendedcf(ymaxCol, dy, QminC, QmaxC, dlnlnQ, nloop, intorder, factscheme);
        hoppetSetFFNcf(fixed_nf);
        hoppetSetPoleMassVFNcf(mcharm, mbottom, mtop);
        
//         hoppetStartExtendedcff(ymaxCol, dy, Qmin, Qmax, dlnlnQ, nloop, intorder, factscheme);
        hoppetStartExtendedcff(ymaxCol, dy, QminC, QmaxC, dlnlnQ, nloop, intorder, factscheme);
        hoppetSetFFNcff(fixed_nf);
        hoppetSetPoleMassVFNcff(mcharm, mbottom, mtop);

    }
    
//     strcat(res_directory,"/");
//     resdir += "/";
//     cout<<res_directory<<"\t"<<resdir<<endl;
    
    if(res_errs.size() > 0) res_errs.resize(res_errs.size()-4);
    res_errs = res_errs+"-errs.dat";
    if(res_pars.size() > 0) res_pars.resize(res_pars.size()-4);
    res_pars = res_pars + "-pars.dat";

    int resdirLen = resdir.length();
//     cout<<resdirLen<<endl;
    
    ofstream res(results), errs(res_errs), pars(res_pars);
    
    if(myCol.col_model == 0){
        res<<"Collins fit 2015"<<endl;
    }

    if(myCol.col_model == 1){
        res<<"Collins fit 2015 - Bernstein polynomials"<<endl;
    }

    if(myCol.col_model == 2){
        res<<"Collins fit 2018 - REWEIGHTING - Unbiased form"<<endl;
    }
    
    if(!res.is_open()) cerr<<"Error in opening results file\n";
    
    // res<<"\nLHAPDF ID: "<<pdf->lhapdfID()<<", member: "<<pdf->memberID()<<"\n\n";
    
    res<<"Collinear PDF set: "<<myPDF.pdf->set().description()<<endl;
    res<<"LHAPDF ID: "<<myPDF.pdf->lhapdfID()<<", member: "<<myPDF.pdf->memberID()<<endl;
    res<<"Collinear FF set: "<<myFF.FFset<<", order: "<<myFF.FForder<<endl<<endl;

    
    FCN function;
    MnUserParameters upar;
    
    if(myCol.col_model == 0 && myTransv.transv_model == 0){ //fit 2015 
//         upar.Add("NT_uval",0.4626649640857, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_dval",-0.9999994517746, 0.01);
//         upar.Add("alpha",0.6681874036618, 0.01);
//         upar.Add("beta", 2.039303368588, 0.01);
//         upar.Add("NC_fav", 0.8934879992544, 0.1);
//         upar.Add("NC_dis",-0.3687940726268, 0.1);
//         upar.Add("gamma", 2.014212741211, 0.001);
//         upar.Add("delta",0.0, 0.001);
//         upar.Add("sqMC",  0.2763349233157, 0.01);
   
//         upar.Add("NT_uval", 2.66, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_dval", -9.09, 0.01);
//         upar.Add("alpha", 0.744, 0.01);
//         upar.Add("beta", 1.77, 0.01);
//         upar.Add("NC_fav", 0.369, 0.1);
//         upar.Add("NC_dis", -0.326, 0.1);
//         upar.Add("gamma", 1.18, 0.001);
//         upar.Add("delta", -0.0492, 0.001);
//         upar.Add("sqMC", 0.273, 0.01);
        
        upar.Add("NT_uval",0.61, 0.01);   //(name,initial value,"error size")
        upar.Add("NT_dval",-1.0, 0.01);
        upar.Add("alpha",0.7, 0.01);
        upar.Add("beta",1.8, 0.01);
        upar.Add("NC_fav",0.9, 0.1);
        upar.Add("NC_dis",-0.37, 0.1);
        upar.Add("gamma",2.02, 0.001);
        upar.Add("delta",0.0, 0.001);
        upar.Add("sqMC",0.28, 0.01);
        
//         upar.Add("NT_uval",0.461781, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_dval",-1.48426, 0.01);
//         upar.Add("alpha",0.749479, 0.01);
//         upar.Add("beta",2.24748, 0.01);
//         upar.Add("NC_fav",0.867332, 0.1);
//         upar.Add("NC_dis",-0.377584, 0.1);
//         upar.Add("gamma",1.98177, 0.001);
//         upar.Add("delta",0.0, 0.001);
//         upar.Add("sqMC",0.27706, 0.01);
        
        if(SB_flag){
            upar.SetLimits("NT_uval",-1.0, 1.0);
            upar.SetLimits("NT_dval",-1.0, 1.0);
        }
        upar.SetLimits("NC_fav",-1.0, 1.0);
        upar.SetLimits("NC_dis",-1.0, 1.0);
//         upar.SetLimits("NC_fav",-10.0, 10.0);
//         upar.SetLimits("NC_dis",-10.0, 10.0);
        upar.SetLowerLimit("alpha",0.0);
        upar.SetLowerLimit("gamma", 0.0);
        upar.SetLowerLimit("beta",0.0);
//         upar.SetLimits("alpha", 0.0, 3.0);
//         upar.SetLimits("beta", 0.0, 3.0);
//         upar.SetLimits("gamma", 0.0, 3.0);
//         upar.SetLimits("delta", -3.0, 3.0);
        upar.SetLowerLimit("delta",0.0);
        upar.SetLowerLimit("sqMC", 0.0);
        
//         upar.Fix("NT_dval"); NfixedPars++;  
//         upar.Fix("NT_uval"); NfixedPars++;
//         upar.Fix("alpha"); NfixedPars++;
//         upar.Fix("beta"); NfixedPars++;
//         upar.Fix("NC_fav"); NfixedPars++;
//         upar.Fix("NC_dis"); NfixedPars++;
//         upar.Fix("gamma");  NfixedPars++;
        upar.Fix("delta");  NfixedPars++;
//         upar.Fix("sqMC"); NfixedPars++;
        
    }

    if(myCol.col_model == 100 && myTransv.transv_model == 100){ //fit 2015 beta_u != beta_d
        
//         cout<<"adding parameters for model == "<<myCol.col_model<<endl;
        upar.Add("NT_uval", 0.6, 0.01);   //(name,initial value,"error size")
        upar.Add("NT_dval", -1.0, 0.01);
        upar.Add("alpha", 0.82, 0.01);
        upar.Add("beta_uv", 2.1, 0.01);
        upar.Add("beta_dv", 3.0, 0.01);
        upar.Add("NC_fav", 0.87, 0.01);
        upar.Add("NC_dis", -0.38, 0.01);
        upar.Add("gamma", 2.0, 0.001);
        upar.Add("delta", 0.0, 0.001);
        upar.Add("sqMC",  0.28, 0.01);

        
        if(SB_flag){
            upar.SetLimits("NT_uval",-1.0, 1.0);
            upar.SetLimits("NT_dval",-1.0, 1.0);
        }
        upar.SetLowerLimit("alpha",0.0);
        upar.SetLowerLimit("beta_uv",0.0);
        upar.SetLowerLimit("beta_dv",0.0);
        upar.SetLimits("NC_fav",-1.0,1.0);
        upar.SetLimits("NC_dis",-1.0,1.0);
        upar.SetLowerLimit("gamma",0.0);
        upar.SetLowerLimit("delta",0.0);
        upar.SetLowerLimit("sqMC",0.0);
        
//          upar.Fix("NT_uval"); NfixedPars++;
//          upar.Fix("NT_dval"); NfixedPars++;  
//          upar.Fix("alpha"); NfixedPars++;
//          upar.Fix("beta_uv"); NfixedPars++;
//          upar.Fix("beta_dv"); NfixedPars++;
//          upar.Fix("NC_fav"); NfixedPars++;
//          upar.Fix("NC_dis"); NfixedPars++;
//          upar.Fix("gamma");  NfixedPars++;
         upar.Fix("delta");  NfixedPars++;
//          upar.Fix("sqMC"); NfixedPars++;
        
    }
    
    
    if(myCol.col_model == 1 && myTransv.transv_model == 1){ //fit 2015 - Bernstein polynomials for Collins
        upar.Add("NT_uval",0.58, 0.01);   //(name,initial value,"error size")
        upar.Add("NT_dval",-1.0, 0.01);
        upar.Add("alpha",0.79, 0.01);
        upar.Add("beta",1.44, 0.01);
        upar.Add("a_fav",-0.02, 0.01);
        upar.Add("a_dis",-1.0, 0.01);
        upar.Add("b_fav",0.66, 0.01);
        upar.Add("b_dis",0.12, 0.01);
        upar.Add("sqMC",0.27, 0.01);

        if(SB_flag){
            upar.SetLimits("NT_uval",-1.0,1.0);
            upar.SetLimits("NT_dval",-1.0,1.0);
        }
        upar.SetLowerLimit("alpha",0.0);
        upar.SetLowerLimit("beta",0.0);        
        upar.SetLimits("a_fav",-1.0,1.0);
        upar.SetLimits("a_dis",-1.0,1.0);
        upar.SetLimits("b_fav",-1.0,1.0);
        upar.SetLimits("b_dis",-1.0,1.0);
//         upar.SetLowerLimit("sqMC",0.0);
        
//         upar.Fix("NT_uval"); NfixedPars++;
//         upar.Fix("NT_dval"); NfixedPars++;  
//         upar.Fix("alpha"); NfixedPars++;
//         upar.Fix("beta"); NfixedPars++;
//         upar.Fix("a_fav"); NfixedPars++;
//         upar.Fix("a_dis"); NfixedPars++;
//         upar.Fix("b_fav");  NfixedPars++;
//         upar.Fix("b_dis");  NfixedPars++;
//         upar.Fix("sqMC"); NfixedPars++;
        
    }
    
    if(myCol.col_model == 2 && myTransv.transv_model == 2){ //fit 2018 - REWEIGHTING - unbiased
//         upar.Add("NT_u",2.7, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_d",-1.3, 0.01);
//         upar.Add("a_u",0.24, 0.01);
//         upar.Add("a_d",-0.1, 0.01);
//         upar.Add("b_u",3.8, 0.01);
//         upar.Add("b_d",5.1, 0.01);
//         upar.Add("NC_fav",0.6, 0.01);
//         upar.Add("NC_dis",-0.5, 0.01);
//         upar.Add("gamma",-1.0, 0.01);
//         upar.Add("delta",2.0, 0.01);
//         upar.Add("Col_p_width",0.3, 0.01);
        
        upar.Add("NT_u",0.2, 0.01);   //(name,initial value,"error size")
        upar.Add("NT_d",-0.9, 0.01);
        upar.Add("a_u",0.7, 0.01);
        upar.Add("a_d",1.0, 0.01);
        upar.Add("b_u",1.8, 0.01);
        upar.Add("b_d",3.0, 0.01);
        upar.Add("NC_fav",0.23, 0.01);
        upar.Add("NC_dis",-0.03, 0.01);
        upar.Add("gamma",-2.06, 0.01);
        upar.Add("delta",2.28, 0.01);
        upar.Add("Col_p_width",0.08, 0.01);

//         upar.SetLimits("NT_u",-1.0,1.0);
//         upar.SetLimits("NT_d",-1.0,1.0);
//         upar.SetLimits("a_u",-1.0,1.0);
//         upar.SetLimits("a_d",-1.0,1.0);     
//         upar.SetLimits("b_u",-1.0,1.0);
//         upar.SetLimits("b_d",-1.0,1.0);      
//         upar.SetLimits("NC_fav",-1.0,1.0);
//         upar.SetLimits("NC_dis",-1.0,1.0); 
//         upar.SetLimits("a_fav",-1.0,1.0);
//         upar.SetLimits("a_dis",-1.0,1.0);
//         upar.SetLimits("b_fav",-1.0,1.0);
//         upar.SetLimits("b_dis",-1.0,1.0);
//         upar.SetLowerLimit("a_u",0.0);
//         upar.SetLowerLimit("a_d",0.0);
//         upar.SetLowerLimit("b_u",0.0);
//         upar.SetLowerLimit("b_d",0.0);
//         upar.SetLowerLimit("gamma",0.0);
//         upar.SetLowerLimit("delta",0.0);
//         upar.SetLowerLimit("Col_p_width",0.0);
//         upar.SetLowerLimit("Col_p_width",0.0);
        upar.SetLowerLimit("Col_p_width",0.0);
        
//         upar.Fix("NT_u"); NfixedPars++;
//         upar.Fix("NT_d"); NfixedPars++;  
//         upar.Fix("a_u"); NfixedPars++;
//         upar.Fix("a_d"); NfixedPars++;
//         upar.Fix("b_u"); NfixedPars++;
//         upar.Fix("b_d"); NfixedPars++;
//         upar.Fix("NC_fav"); NfixedPars++;
//         upar.Fix("NC_dis"); NfixedPars++;
//         upar.Fix("gamma");  NfixedPars++;
//         upar.Fix("delta");  NfixedPars++;
//         upar.Fix("Col_p_width"); NfixedPars++;
        
    }
    
    if(myCol.col_model == 3 && myTransv.transv_model == 3 || myCol.col_model == 10 && myTransv.transv_model == 10){ //fit 2018 - REWEIGHTING - unbiased - all different
//         upar.Add("NT_u",2.7, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_d",-1.3, 0.01);
//         upar.Add("a_u",0.24, 0.01);
//         upar.Add("a_d",-0.1, 0.01);
//         upar.Add("b_u",3.8, 0.01);
//         upar.Add("b_d",5.1, 0.01);
//         upar.Add("NC_fav",0.6, 0.01);
//         upar.Add("NC_dis",-0.5, 0.01);
//         upar.Add("gamma",-1.0, 0.01);
//         upar.Add("delta",2.0, 0.01);
//         upar.Add("Col_p_width",0.3, 0.01);
        
        upar.Add("NT_u",1.0, 0.01);   //(name,initial value,"error size")
        upar.Add("NT_d",-3.0, 0.01);
        upar.Add("a_u",0.0, 0.01);
        upar.Add("a_d",0.0, 0.01);
        upar.Add("b_u",2.4, 0.01);
        upar.Add("b_d",3.6, 0.01);
        upar.Add("NC_fav",0.3, 0.01);
        upar.Add("NC_dis",-0.1, 0.01);
        upar.Add("g_fav",-2.0, 0.01);
        upar.Add("g_dis",-2.5, 0.01);
        upar.Add("d_fav",2.00, 0.01);
        upar.Add("d_dis",0.3, 0.01);
        upar.Add("Col_p_width",0.08, 0.01);

        
//         upar.Add("NT_u",0.7, 0.01);   //(name,initial value,"error size")
//         upar.Add("NT_d",-1.0, 0.01);
//         upar.Add("a_u",0.7, 0.01);
//         upar.Add("a_d",0.0, 0.01);
//         upar.Add("b_u",1.0, 0.01);
//         upar.Add("b_d",0.0, 0.01);
//         upar.Add("NC_fav",0.4, 0.01);
//         upar.Add("NC_dis",-0.5, 0.01);
//         upar.Add("g_fav",2.0, 0.01);
//         upar.Add("g_dis",2.0, 0.01);
//         upar.Add("d_fav",0.02, 0.01);
//         upar.Add("d_dis",0.02, 0.01);
//         upar.Add("Col_p_width",0.3, 0.01);

//         upar.SetLimits("NT_u",-1.0,1.0);
//         upar.SetLimits("NT_d",-1.0,1.0);
//         upar.SetLowerLimit("a_u",0.0);
//         upar.SetLowerLimit("a_d",0.0);  
//         upar.SetLowerLimit("b_u",0.0);
//         upar.SetLowerLimit("b_d",0.0);        
//         upar.SetLimits("NC_fav",-1.0,1.0);
//         upar.SetLimits("NC_dis",-1.0,1.0);
//         upar.SetLowerLimit("g_fav",0.0);
//         upar.SetLowerLimit("g_dis",0.0);
//         upar.SetLowerLimit("d_fav",0.0);
//         upar.SetLowerLimit("d_dis",0.0);
//         upar.SetLowerLimit("Col_p_width",0.0);
        
//         upar.Fix("NT_u"); NfixedPars++;
//         upar.Fix("NT_d"); NfixedPars++;  
        upar.Fix("a_u"); NfixedPars++;
        upar.Fix("a_d"); NfixedPars++;
//         upar.Fix("b_u"); NfixedPars++;
//         upar.Fix("b_d"); NfixedPars++;
//         upar.Fix("NC_fav"); NfixedPars++;
//         upar.Fix("NC_dis"); NfixedPars++;
//         upar.Fix("g_fav");  NfixedPars++;
//         upar.Fix("g_dis");  NfixedPars++;
//         upar.Fix("d_fav");  NfixedPars++;
//         upar.Fix("d_dis");  NfixedPars++;
//         upar.Fix("Col_p_width"); NfixedPars++;
        
    }
    
    
    //set x- and z-range
//     Comp17_01hp_x.set_xrange(xmin,xmax);
//     Comp17_01hp_z.set_xrange(xmin,xmax); 
//     Comp17_01hp_pt.set_xrange(xmin,xmax);
//     Comp17_02hp_x.set_xrange(xmin,xmax);
//     Comp17_02hp_z.set_xrange(xmin,xmax); 
//     Comp17_02hp_pt.set_xrange(xmin,xmax);
//     JLab_pim_x.set_xrange(xmin,xmax);
    
//     Comp17_01hm_x.set_zrange(zmin,zmax);
//     Comp17_01hm_z.set_zrange(zmin,zmax);    
//     Comp17_01hm_pt.set_zrange(zmin,zmax);   

//     Comp17_01hp_x.set_Q2range(Q2min,Q2max);
//     Comp17_01hp_z.set_Q2range(Q2min,Q2max);    
//     Comp17_01hp_pt.set_Q2range(Q2min,Q2max);   
//     Comp17_01hm_x.set_Q2range(Q2min,Q2max);
//     Comp17_01hm_z.set_Q2range(Q2min,Q2max);    
//     Comp17_01hm_pt.set_Q2range(Q2min,Q2max); 
    
//     Comp17_02hp_x.set_Q2range(Q2min,Q2max);
//     Comp17_02hp_z.set_Q2range(Q2min,Q2max);    
//     Comp17_02hp_pt.set_Q2range(Q2min,Q2max);   
//     Comp17_02hm_x.set_Q2range(Q2min,Q2max);
//     Comp17_02hm_z.set_Q2range(Q2min,Q2max);    
//     Comp17_02hm_pt.set_Q2range(Q2min,Q2max);  
    
//     myFF.FF_plot("./DSSLOPipQ2-112.dat", "DSS", 1, 1, 1, 0, 112.);
//     myFF.FF_plot("./DSSLOPimQ2-112.dat", "DSS", 1, 1, -1, 0, 112.);
    
    
    //minimization + calculation of No. dof
    vector<double> par = upar.Params(), fitted_pars(par.size());
    int Npar = par.size() - NfixedPars;
    
    function.set_NumberOfParameters(Npar);
    
    MnMigrad migrad(function, upar);
    FunctionMinimum min = migrad();
    
//     MnMinos minos(function, min);
//     MinosError beta = minos.Minos(3, 1000, 0.1);
    
    Ndof = Npoints - Npar;
//     cout<<Npoints<<"\t"<<Ndof<<endl;
    
    cout<<"minimum: "<<min<<endl<<endl;
    cout<<"No. of pts = "<<Npoints<<endl;
    cout<<"chi2/Npts="<<min.Fval()/Npoints<<"\nchi2/dof="<<min.Fval()/Ndof<<endl;
    
    res<<"minimum: "<<min<<"\n\nNo. of pts="<<Npoints<<"\nchi2/Npts="<<min.Fval()/Npoints<<"\nchi2/dof="<<min.Fval()/Ndof<<endl;
    
    if(!upar.Parameter(0).IsFixed()) pars<<min.UserParameters().Value(0);
    for(int i = 1; i < par.size(); i++){
        if(!upar.Parameter(i).IsFixed()) pars<<"\t"<<min.UserParameters().Value(i);
    }
    pars<<endl;
    
    errs<<min.UserCovariance();
//     cout<<min.UserCovariance();
    res<<"\n\nSet\t\t\t\t\tchi2_set\tnpts\tchi2/npts\n";
    
    //close files
    res.close();
    errs.close();
    pars.close();
    
    for(int i = 0; i < fitted_pars.size(); i++){
        fitted_pars[i] = min.UserParameters().Value(i);
//         if(fitted_pars[i] != par[i]) cerr<<"Fitted_pars["<<i<<"] is different from the initial one\n";
    }
        
    //creating directory for results
    const int dir_err = mkdir(res_directory, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (dir_err == -1) cout<<"Error in creating results directory (maybe folder already exists).\n";
    
    //print transversity and Collins
    if(!SB_flag){
        
        myTransv.plot(resdir+"/transv-" + today + "-CTEQ66-DSSnlo-noSB-h1-evo-2.4.dat", 2.4, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
        myTransv.plot(resdir+"/transv-" + today + "-CTEQ66-DSSnlo-noSB-h1-evo-4.0.dat", 4.0, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
    }
    
    else{
     
        myTransv.plot(resdir+"/transv-" + today + "-CTEQ66-DSSnlo-usingSB-h1-evo-2.4.dat", 2.4, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
        myTransv.plot(resdir+"/transv-" + today + "-CTEQ66-DSSnlo-usingSB-h1-evo-4.0.dat", 4.0, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
    }

    myTransv.plot_SB(resdir+"/SB-CTEQ66-DSSnlo-Q2-2.4.dat", 2.4, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
    myTransv.plot_SB(resdir+"/SB-CTEQ66-DSSnlo-Q2-4.0.dat", 4.0, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
//     myTransv.plot_SB("./SB-cteq66+dssv08.dat", 2.4, myTransv.transv_model, myTransv.transv_evo, fitted_pars);
//     myCol.plot(resdir+"/Collins-test-31-10-2019-fit2015-h1evot-H1pevot-1.0.dat", 1.0, 1, myCol.col_model, myCol.col_evo, fitted_pars);
//     myCol.plot(resdir+"/Collins-test-31-10-2019-fit2015-h1evot-H1pevot-2.4.dat", 2.4, 1, myCol.col_model, myCol.col_evo, fitted_pars);
//     myCol.plot(resdir+"/Collins-test-31-10-2019-fit2015-h1evot-H1pevot-112.dat", 112.0, 1, myCol.col_model, myCol.col_evo, fitted_pars);
    
    //print results for every experiment
    JLab11_Npp_x.print(resdir+"/res_jlab2011_col_aut_Npp_x.dat", results, resdirLen);
    JLab11_Npm_x.print(resdir+"/res_jlab2011_col_aut_Npm_x.dat", results, resdirLen);

    Herm10_pip_x.print(resdir+"/res_herm2010_col_aut_pip_x.dat", results, resdirLen);
    Herm10_pip_z.print(resdir+"/res_herm2010_col_aut_pip_z.dat", results, resdirLen);
    Herm10_pip_pt.print(resdir+"/res_herm2010_col_aut_pip_pt.dat",results,resdirLen);
    Herm10_pim_x.print(resdir+"/res_herm2010_col_aut_pim_x.dat",results,resdirLen);
    Herm10_pim_z.print(resdir+"/res_herm2010_col_aut_pim_z.dat",results,resdirLen);
    Herm10_pim_pt.print(resdir+"/res_herm2010_col_aut_pim_pt.dat",results,resdirLen);
    
    Herm20_pip.print(resdir+"/res_herm2020_col_aut_pip.dat",results,resdirLen);
    Herm20_pim.print(resdir+"/res_herm2020_col_aut_pim.dat",results,resdirLen);
    
//     Comp12_Ppp_x.print(resdir+"/res_comp2012_col_aut_Ppp_x.dat",results,resdirLen);
//     Comp12_Ppp_z.print(resdir+"/res_comp2012_col_aut_Ppp_z.dat",results,resdirLen);
//     Comp12_Ppp_pt.print(resdir+"/res_comp2012_col_aut_Ppp_pt.dat",results,resdirLen);
//     Comp12_Ppm_x.print(resdir+"/res_comp2012_col_aut_Ppm_x.dat",results,resdirLen);
//     Comp12_Ppm_z.print(resdir+"/res_comp2012_col_aut_Ppm_z.dat",results,resdirLen);
//     Comp12_Ppm_pt.print(resdir+"/res_comp2012_col_aut_Ppm_pt.dat",results,resdirLen);
    
    Comp12_Dpp_x.print(resdir+"/res_comp2012_col_aut_Dpp_x.dat",results,resdirLen);
    Comp12_Dpp_z.print(resdir+"/res_comp2012_col_aut_Dpp_z.dat",results,resdirLen);
    Comp12_Dpp_pt.print(resdir+"/res_comp2012_col_aut_Dpp_pt.dat",results,resdirLen);
    Comp12_Dpm_x.print(resdir+"/res_comp2012_col_aut_Dpm_x.dat",results,resdirLen);
    Comp12_Dpm_z.print(resdir+"/res_comp2012_col_aut_Dpm_z.dat",results,resdirLen);
    Comp12_Dpm_pt.print(resdir+"/res_comp2012_col_aut_Dpm_pt.dat",results,resdirLen);
    
    Comp14_Ppp_x.print(resdir+"/res_comp2014_col_aut_Ppp_x.dat",results,resdirLen);
    Comp14_Ppp_z.print(resdir+"/res_comp2014_col_aut_Ppp_z.dat",results,resdirLen);
    Comp14_Ppp_pt.print(resdir+"/res_comp2014_col_aut_Ppp_pt.dat",results,resdirLen);
    Comp14_Ppm_x.print(resdir+"/res_comp2014_col_aut_Ppm_x.dat",results,resdirLen);
    Comp14_Ppm_z.print(resdir+"/res_comp2014_col_aut_Ppm_z.dat",results,resdirLen);
    Comp14_Ppm_pt.print(resdir+"/res_comp2014_col_aut_Ppm_pt.dat",results,resdirLen);
    
//     Comp17_01hp_x.print(resdir+"/res_comp2017_Zgt01_p_x_Col.dat",results,resdirLen);
//     Comp17_01hp_z.print(resdir+"/res_comp2017_Zgt01_p_z_Col.dat",results,resdirLen);
//     Comp17_01hp_pt.print(resdir+"/res_comp2017_Zgt01_p_pt_Col.dat",results,resdirLen);
//     Comp17_02hp_x.print(resdir+"/res_comp2017_Zgt02_p_x_Col.dat",results,resdirLen);
//     Comp17_02hp_z.print(resdir+"/res_comp2017_Zgt02_p_z_Col.dat",results,resdirLen);
//     Comp17_02hp_pt.print(resdir+"/res_comp2017_Zgt02_p_pt_Col.dat",results,resdirLen);
    
//     Comp17_01hm_x.print(resdir+"/res_comp2017_Zgt01_n_x_Col.dat",results,resdirLen);
//     Comp17_01hm_z.print(resdir+"/res_comp2017_Zgt01_n_z_Col.dat",results,resdirLen);
//     Comp17_01hm_pt.print(resdir+"/res_comp2017_Zgt01_n_pt_Col.dat",results,resdirLen);
//     Comp17_02hm_x.print(resdir+"/res_comp2017_Zgt02_n_x_Col.dat",results,resdirLen);
//     Comp17_02hm_z.print(resdir+"/res_comp2017_Zgt02_n_z_Col.dat",results,resdirLen);
//     Comp17_02hm_pt.print(resdir+"/res_comp2017_Zgt02_n_pt_Col.dat",results,resdirLen);
    
    Belle12_A0UL_z1z2.print(resdir+"/res_belle2012_col_A0UL_z1z2_pp.dat",results,resdirLen);
    Belle12_A0UC_z1z2.print(resdir+"/res_belle2012_col_A0UC_z1z2_pp.dat",results,resdirLen);
//     Belle12_A12UL_z1z2.print(resdir+"/res_belle2012_col_A12UL_z1z2_pp.dat",results,resdirLen);
//     Belle12_A12UC_z1z2.print(resdir+"/res_belle2012_col_A12UC_z1z2_pp.dat",results,resdirLen);
    
    Babar13_A0UL_z1z2.print(resdir+"/res_babar2013_col_A0UL_z1z2_pp.dat",results,resdirLen);
    Babar13_A0UC_z1z2.print(resdir+"/res_babar2013_col_A0UC_z1z2_pp.dat",results,resdirLen);
//     Babar13_A12UL_z1z2.print(resdir+"/res_babar2013_col_A12UL_z1z2_pp.dat",results,resdirLen);
//     Babar13_A12UC_z1z2.print(resdir+"/res_babar2013_col_A12UC_z1z2_pp.dat",results,resdirLen);
    
    Babar14_A0UL_PT.print(resdir+"/res_babar2014_col_A0UL_pt_pp.dat",results,resdirLen);
    Babar14_A0UC_PT.print(resdir+"/res_babar2014_col_A0UC_pt_pp.dat",results,resdirLen);
    
    Babar15_A0UL_z1z2.print(resdir+"/res_babar2015_col_A0UL_z1z2_pp.dat",results,resdirLen);
    Babar15_A0UC_z1z2.print(resdir+"/res_babar2015_col_A0UC_z1z2_pp.dat",results,resdirLen);
//     Babar15_A12UL_z1z2.print(resdir+"/res_babar2015_col_A12UL_z1z2_pp.dat",results,resdirLen);
//     Babar15_A12UC_z1z2.print(resdir+"/res_babar2015_col_A12UC_z1z2_pp.dat",results,resdirLen);
    
    BESIII16_A0UL_PT.print(resdir+"/res_besIII2016_col_A0UL_pt_pp.dat",results,resdirLen);
    BESIII16_A0UC_PT.print(resdir+"/res_besIII2016_col_A0UC_pt_pp.dat",results,resdirLen);

    
    return 0;
  
}

