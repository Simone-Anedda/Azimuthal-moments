//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#include "FCN.h"

    using namespace std;

    //add more extern objects to add more data    
    extern DATA::Data\
    
    JLab11_Npp_x, JLab11_Npm_x,\
    
    Herm10_pip_x, Herm10_pip_z, Herm10_pip_pt,\
    Herm10_pim_x, Herm10_pim_z, Herm10_pim_pt,\
    
    Herm20_pip_x_1, Herm20_pip_x_2, Herm20_pip_x_3, Herm20_pip_x_4,\
    Herm20_pim_x_1, Herm20_pim_x_2, Herm20_pim_x_3, Herm20_pim_x_4,\

    Comp12_Ppp_x, Comp12_Ppp_z, Comp12_Ppp_pt,\
    Comp12_Ppm_x, Comp12_Ppm_z, Comp12_Ppm_pt,\
    Comp12_Dpp_x, Comp12_Dpp_z, Comp12_Dpp_pt,\
    Comp12_Dpm_x, Comp12_Dpm_z, Comp12_Dpm_pt,\
    
    Comp14_Ppp_x, Comp14_Ppp_z, Comp14_Ppp_pt,\
    Comp14_Ppm_x, Comp14_Ppm_z, Comp14_Ppm_pt,\
    
    Comp17_01hp_x, Comp17_01hp_z, Comp17_01hp_pt,\
    Comp17_02hp_x, Comp17_02hp_z, Comp17_02hp_pt,\
    Comp17_01hm_x, Comp17_01hm_z, Comp17_01hm_pt,\
    Comp17_02hm_x, Comp17_02hm_z, Comp17_02hm_pt,\
    
    Belle12_A0UL_z1z2, Belle12_A0UC_z1z2,\
    Belle12_A12UL_z1z2, Belle12_A12UC_z1z2,\
    
    Babar13_A0UL_z1z2, Babar13_A0UC_z1z2,\
    Babar13_A12UL_z1z2, Babar13_A12UC_z1z2,\
    
    Babar14_A0UL_PT, Babar14_A0UC_PT,\
    
    Babar15_A0UL_z1z2, Babar15_A0UC_z1z2,\
    Babar15_A12UL_z1z2, Babar15_A12UC_z1z2,\
    
    BESIII16_A0UL_PT, BESIII16_A0UC_PT,\
    
    Belle19_A12UL_z1z2, Belle19_A12UC_z1z2,\
    
    Comp_14_Ppp_x_SS, Comp_14_Ppp_z_SS, Comp_14_Ppp_pt_SS,\
    Comp_14_Ppm_x_SS, Comp_14_Ppm_z_SS, Comp_14_Ppm_pt_SS,\
    Belle12_A0UL_z1z2_SS, Belle12_A0UC_z1z2_SS,\
    
    Comp24_hp_x, Comp24_hp_z, Comp24_hp_pt,\
    Comp24_hm_x, Comp24_hm_z, Comp24_hm_pt;
    
    extern COLLPDF::CollPDF myPDF;

    extern TRANSV::TRANSVERSITY myTransv;
    extern COL::COLLINS myCol;
    
    extern FRAG::FF myFF;
        
    std::vector<double> evo_par(20), int_par(20);
    std::vector<double> PDF(13);
    std::vector<double> COL_z(13), COL_ppz1(13), COL_ppz2(13), COL_pmz1(13), COL_pmz2(13);
    std::vector<double> FF(13), FF_ppz1(13), FF_ppz2(13), FF_pmz1(13), FF_pmz2(13),\
        FF_ppz1_fixedQ2(13), FF_ppz2_fixedQ2(13), FF_pmz1_fixedQ2(13), FF_pmz2_fixedQ2(13),\
        FF_evo(13);
        
    double f[13], h1[13]; //for HOPPET and transversity
    
    const double charges[13] = {-2./3., 1./3., -2./3., 1./3., -2./3., 1./3., 0, -1./3., 2./3., -1./3., 2./3., -1./3., 2./3.};
    
    double k_width, p_width, M1;
    double PT_width, pperpC_width, PTC_width;
    double FC_ee, rho_C;
    
    double output = 0, F_UU = 0, F_UT = 0,\
        numU = 0, denU = 0,\
        numL = 0, denL = 0,\
        numC = 0, denC = 0,\
        numFact = 0, denFact = 0;

    double tmp = 0, tmp1 = 0, tmp2 = 0;
    int icount = 0, hessian_count = 0;
    int ivalid;
    
    int comp, nregions, neval, fail;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    int int_charge, int_order, int_hadron;//, int_epem_sign, int_process;
    double int_av_pt, int_av_Q2;
    string int_fname, int_epem_sign, int_process;

    int evo_charge, evo_hadron, evo_process;
    
    extern int Npoints, NfixedPars, hessian_sign, i_hesse, i_rep;

    extern string mode;
    
    extern double target_chi2;
    
    bool evo_flag, poly_flag = true;
    extern bool diffsame_flag, fixed_Q2, dssv_init, grv_flag; 
    
    extern bool hessian_flag, MC_flag, hesse_fit;  

    bool check_fit;
    // ofstream chi2_iter("./chi2_iter.dat");

double maximum(const double & a, const double & b){ return a < b ? b : a; }

double maximum(const double & a, const double & b, const double & c){ return maximum(maximum(a, b), c); }

    
void transv_subroutine(const double & x, const double & Q, double * res)
{
    double Q02 = pow(Q,2);
     
//     cout<<"in transv_subroutine: \t";
//     for(int i = 0; i < evo_par.size(); i++) cout<<i<<": "<<evo_par[i]<<"\t";
//     cout<<endl;
        
    myTransv.eval(x, Q02, evo_par);
    
    for(int i = 0; i < 6; i++) res[i] = x * myTransv.TRANSV_x[i];
    res[6] = 0;
    for(int i = 7; i < 13; i++) res[i] = x * myTransv.TRANSV_x[i];

}


void Collins_subroutine_pip(const double & z, const double & Q, double * res)
{
    double Q02 = pow(Q,2);
    
    int hadron = 1, charge = 1;
    
    vector<double> FF(13);
    
//     if(z == 1) z -= EPSILON;
    
//     cout<<"in Collins_subroutine_pip: z = "<<z<<endl;
//     for(int i = 0; i < evo_par.size(); i++) cout<<i<<": "<<evo_par[i]<<"\t";
//     cout<<endl;
    
//     cout<<"Q02 pip = "<<Q02<<endl;
    
    // myFF.FF_eval("DSS", iset_ff, hadron, charge, int_order, z, Q02);
    myFF.FF_eval(hadron, charge, z, Q02);
    for(int i = 0; i < FF.size(); i++) FF[i] = myFF.theFF[i] / z; //to scale FFs: routines give z*D(z)  

    myCol.eval(z, Q02, charge, FF, evo_par);
        
    for(int i = 0; i < 6; i++) res[i] = z * myCol.COL_z[i];
    res[6] = 0;
    for(int i = 7; i < 13; i++) res[i] = z * myCol.COL_z[i];

}

void Collins_subroutine_pim(const double & z, const double & Q, double * res)
{
    double Q02 = pow(Q,2);
    
    int hadron = 1, charge = -1;
    
    vector<double> FF(13);

//     if(z == 1) z -= EPSILON;
//     cout<<"in Collins_subroutine_pim: z = "<<z<<endl;
//     cout<<"in Collins_subroutine_pim: \t";
//     for(int i = 0; i < evo_par.size(); i++) cout<<i<<": "<<evo_par[i]<<"\t";
//     cout<<endl;
//     cout<<"Q02 pim = "<<Q02<<endl;

    // myFF.FF_eval("DSS", iset_ff, hadron, charge, int_order, z, Q02);
    myFF.FF_eval(hadron, charge, z, Q02);
    for(int i = 0; i < FF.size(); i++) FF[i] = myFF.theFF[i] / z; //to scale FFs: routines give z*D(z)  

    myCol.eval(z, Q02, charge, FF, evo_par);
    
    for(int i = 0; i < 6; i++) res[i] = z * myCol.COL_z[i];
    res[6] = 0;
    for(int i = 7; i < 13; i++) res[i] = z * myCol.COL_z[i];

}


void Collins_PT(const double & av_z, const double & av_Q2, const std::vector<double> & par, const string & dataname){

    // int model = myCol.model; 
    // int widths = myCol.widths; //choice of (<k^2_\perp>,<p^2_\perp>)
    std::size_t found = dataname.find("comp");        
    
    double sqMC;
    
    //COLLINS 2015 fit
    if(myTransv.model == "fit-2015" || myTransv.model == "fit-2015-bernstein" || myCol.model == "fit-2013-polynomial") sqMC = par[8];
    if(myTransv.model == "fit-2015-bu-bd") sqMC = par[9];

     //fit 2018 - REWEIGHTING - unbiased
    if(myTransv.model == "ref-fit-unbiased-1"){ 
        pperpC_width = par[10];
        sqMC= -(pperpC_width * p_width) / (pperpC_width - p_width);
    }
        
    if(myTransv.model == "ref-fit-unbiased-2" || myTransv.model == "fit-2015-full"){ 
        pperpC_width = par[12];
        sqMC= -(pperpC_width * p_width) / (pperpC_width - p_width);
    }

    if(myTransv.widths == "latest"){
        
        k_width = 0.57;
        p_width = 0.12;
    
        if (found != std::string::npos && diffsame_flag){
            k_width = 0.60;
            p_width = 0.20;
        } 
    }
    
    if(myTransv.widths == "old"){
        
        k_width = 0.25;
        p_width = 0.20;
    }

    if(myTransv.widths == "latest-bis"){
        
        k_width = 0.58;
        p_width = 0.12;
    
        if (found != std::string::npos && diffsame_flag){
            k_width = 0.52;
            p_width = 0.21;
        } 
    }
    
    
    //set widths
    PT_width = p_width + pow(av_z,2) * k_width; //cout<<"PTwidth="<<PT_width<<endl;
    
    if(myTransv.model == "fit-2015" || myTransv.model == "fit-2015-bernstein" || myTransv.model == "fit-2015-bu-bd" || myTransv.model == "fit-2015-full" || myCol.model == "fit-2013-polynomial"){

        pperpC_width = (sqMC * p_width) / (sqMC + p_width); //cout<<"pperpC="<<pperpC_width<<endl; fit 2015 - MC^2 is the parameter
    } 
        
    PTC_width = pperpC_width + pow(av_z, 2) * k_width; //cout<<"PTCwidth="<<PTC_width<<endl;
    
    //factors for fit 2015 and for z1-z2 integrated asymmetries
    rho_C = sqMC / (sqMC + p_width); //cout<<"rhoC="<<rho_C<<endl;
    
    FC_ee = 2.0 * EulerConst * pow(rho_C, 2) * (1. - rho_C); //cout<<"FCee="<<FC_ee<<endl;
    
}

// vector<double> CollinsHoppetEval(const double z, int charge, double f[]){
// 
//     vector<double> COL(13);
//     
// //     cout<<"in CollinsHoppetEval:  ";
// //     for(int i = 5; i <= 9; i++){
// //         if (i == 5 || i == 7 || i == 8 || i == 9)  cout<<i<<": "<<f[i]<<"\t";
// //     }
// //     cout<<"fav="<<z*myCol.Col_fav<<"\tunf="<<z*myCol.Col_dis;
// //     cout<<endl;
//     
//     COL[0] = 0;
//     COL[1] = 0;
//     COL[2] = 0;
//             
//     for(int i = 3; i <= 9; i++){
//         
//         if(i == 3 || i == 9) COL[i] = f[i] / z; //s or sbar are always disfavoured for pions 
//         
//         if(charge == 1){ //pi+
//         
//             if(i == 5 || i == 8) COL[i] = f[8] / z; //favoured for u and dbar
//             if(i == 4 || i == 7) COL[i] = f[7] / z; //disfavoured for d and ubar
//         }
//         
//         if(charge == -1){ //pi-
// 
//             if(i == 5 || i == 8) COL[i] = f[7] / z; //disfavoured for u and dbar
//             if(i == 4 || i == 7) COL[i] = f[8] / z; //favoured for d and ubar
//         }
//     }
//     
//     COL[6] = 0; //no gluon
//     COL[10] = 0;
//     COL[11] = 0;
//     COL[12] = 0; 
//     
//     return COL;
//     
// }


// void Collins_FF(const int & iset_ff_in, const int & hadron_in, const int & charge_in,const int & order_in, double & av_z1, double & av_z2, double & av_Q2){
void Collins_FF(const int & hadron_in, const int & charge_in, double & av_z1, double & av_z2, double & av_Q2){

    // int iset_ff = iset_ff_in, hadron = hadron_in, charge = charge_in, order = order_in;
    int hadron = hadron_in, charge = charge_in;
    
    double u, ub, d, db, s, sb, c, b, gl;
    
    //call to the Fragmentation Functions
    if(!fixed_Q2){
    
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z1, av_Q2);
        myFF.FF_eval(hadron, charge, av_z1, av_Q2);
        for(int i = 0; i < FF_ppz1.size(); i++) FF_ppz1[i] = myFF.theFF[i] / av_z1;
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z2, av_Q2);
        myFF.FF_eval(hadron, charge, av_z2, av_Q2);
        for(int i = 0; i < FF_ppz2.size(); i++) FF_ppz2[i] = myFF.theFF[i] / av_z2;
        
        charge*=-1; //to call pi- FFs
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z1, av_Q2);
        myFF.FF_eval(hadron, charge, av_z1, av_Q2);
        for(int i = 0; i < FF_pmz1.size(); i++) FF_pmz1[i] = myFF.theFF[i] / av_z1;
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z2, av_Q2);
        myFF.FF_eval(hadron, charge, av_z2, av_Q2);
        for(int i = 0; i < FF_pmz2.size(); i++) FF_pmz2[i] = myFF.theFF[i] / av_z2;
    
    }
    
    //call to the Fragmentation Functions at fixed Q^2
    if(fixed_Q2){
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z1, FIXEDQ2);
        myFF.FF_eval(hadron, charge, av_z1, FIXEDQ2);
        for(int i = 0; i < FF_ppz1.size(); i++) FF_ppz1_fixedQ2[i] = myFF.theFF[i] / av_z1;
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z2, FIXEDQ2);
        myFF.FF_eval(hadron, charge, av_z2, FIXEDQ2);
        for(int i = 0; i < FF_ppz2.size(); i++) FF_ppz2_fixedQ2[i] = myFF.theFF[i] / av_z2;
        
        charge*=-1; //to call pi- FFs
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z1, FIXEDQ2);
        myFF.FF_eval(hadron, charge, av_z1, FIXEDQ2);
        for(int i = 0; i < FF_pmz1.size(); i++) FF_pmz1_fixedQ2[i] = myFF.theFF[i] / av_z1;
        
        // myFF.FF_eval("DSS", iset_ff, hadron, charge, order, av_z2, FIXEDQ2);
        myFF.FF_eval(hadron, charge, av_z2, FIXEDQ2);
        for(int i = 0; i < FF_pmz2.size(); i++) FF_pmz2_fixedQ2[i] = myFF.theFF[i] / av_z2;
    
    }
}


void Collins_epem_loop(const std::vector<double> &COL_ppz1_in, const std::vector<double> &COL_ppz2_in, const std::vector<double> &COL_pmz1_in, const std::vector<double> &COL_pmz2_in, const std::vector<double> &FF_ppz1_in, const std::vector<double> &FF_ppz2_in, const std::vector<double> &FF_pmz1_in, const std::vector<double> &FF_pmz2_in){
    
    for(int i = 3; i <= 9; i++){
            
        if(i == 3){ //sb
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i+6] + COL_pmz1_in[i] * COL_ppz2_in[i+6]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i+6] + FF_pmz1_in[i] * FF_ppz2_in[i+6]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i+6] + COL_pmz1_in[i] * COL_pmz2_in[i+6]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i+6] + FF_pmz1_in[i] * FF_pmz2_in[i+6]);
        }                                
                                         
        if(i == 4){ //ub                 
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i+4] + COL_pmz1_in[i] * COL_ppz2_in[i+4]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i+4] + FF_pmz1_in[i] * FF_ppz2_in[i+4]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i+4] + COL_pmz1_in[i] * COL_pmz2_in[i+4]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i+4] + FF_pmz1_in[i] * FF_pmz2_in[i+4]);
        }                                
                                         
        if(i == 5){ //db                 
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i+2] + COL_pmz1_in[i] * COL_ppz2_in[i+2]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i+2] + FF_pmz1_in[i] * FF_ppz2_in[i+2]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i+2] + COL_pmz1_in[i] * COL_pmz2_in[i+2]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i+2] + FF_pmz1_in[i] * FF_pmz2_in[i+2]);
        }                                
                                         
        if(i == 6){ //g                  
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i] + COL_pmz1_in[i] * COL_ppz2_in[i]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i] + FF_pmz1_in[i] * FF_ppz2_in[i]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i] + COL_pmz1_in[i] * COL_pmz2_in[i]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i] + FF_pmz1_in[i] * FF_pmz2_in[i]);
        }                                
                                         
        if(i == 7){ //d                  
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i-2] + COL_pmz1_in[i] * COL_ppz2_in[i-2]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i-2] + FF_pmz1_in[i] * FF_ppz2_in[i-2]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i-2] + COL_pmz1_in[i] * COL_pmz2_in[i-2]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i-2] + FF_pmz1_in[i] * FF_pmz2_in[i-2]);            
        }                                
                                         
        if(i == 8){ //u                  
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i-4] + COL_pmz1_in[i] * COL_ppz2_in[i-4]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i-4] + FF_pmz1_in[i] * FF_ppz2_in[i-4]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i-4] + COL_pmz1_in[i] * COL_pmz2_in[i-4]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i-4] + FF_pmz1_in[i] * FF_pmz2_in[i-4]);
        }                                
                                         
        if(i == 9){ //s                  
            numU += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_pmz2_in[i-6] + COL_pmz1_in[i] * COL_ppz2_in[i-6]);
            denU += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_pmz2_in[i-6] + FF_pmz1_in[i] * FF_ppz2_in[i-6]);
            numL += pow(charges[i], 2) * (COL_ppz1_in[i] * COL_ppz2_in[i-6] + COL_pmz1_in[i] * COL_pmz2_in[i-6]);
            denL += pow(charges[i], 2) * (FF_ppz1_in[i] * FF_ppz2_in[i-6] + FF_pmz1_in[i] * FF_pmz2_in[i-6]);
        }  
    }
    
}


int Integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata){

#define f0 ff[0]
#define f1 ff[1]
#define f2 ff[2]
#define f3 ff[3]
    
    double z1 = 0.15 + 0.75 * xx[0];
    double z2 = 0.15 + 0.75 * xx[1];
    
//     cout<<"In integrand: z1 = "<<z1<<"\tz2 = "<<z2<<endl;
    
    double f[13];
    
    numU = 0, denU = 0, numL = 0, denL = 0, numC = 0, denC = 0, numFact = 0, denFact = 0;
    
    Collins_PT(z1, int_av_Q2, int_par, int_fname); //note: av_z1 is not used here, just needed to call the function to fix FC_ee factor
    
    // Collins_FF(iset_ff, int_hadron, int_charge, int_order, z1, z2, int_av_Q2);
    Collins_FF(int_hadron, int_charge, z1, z2, int_av_Q2);
    
//     cout<<"In integrand: int_had = "<<int_hadron<<"\t int_ch = "<<int_charge<<endl;

    if(!fixed_Q2){
       
        //call for the Collins functions
        if(myCol.evo == "DGLAP" || myCol.evo == "none"){ 
            
            myCol.eval(z1, int_av_Q2, int_charge, FF_ppz1, int_par);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, int_av_Q2, int_charge, FF_ppz2, int_par);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = myCol.COL_z[i];
            
            int_charge *= -1; //to call pi- Collins
                
            myCol.eval(z1, int_av_Q2, int_charge, FF_pmz1, int_par);
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, int_av_Q2, int_charge, FF_pmz2, int_par);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = myCol.COL_z[i];
        }
        
        if(myCol.evo == "CT3"){
        
            //call for the Collins functions
            hoppetEvalcf(z1, sqrt(int_av_Q2), f);
//             COL_ppz1 = CollinsHoppetEval(z1, int_charge, f);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = f[i] / z1;

            hoppetEvalcf(z2, sqrt(int_av_Q2), f);
//             COL_ppz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = f[i] / z2;
            
            int_charge *= -1; //to call pi- Collins
                
            hoppetEvalcff(z1, sqrt(int_av_Q2), f);      //modified to use proper hoppet calls for pip and pim
//             COL_pmz1 = CollinsHoppetEval(z1, int_charge, f);    
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = f[i] / z1;

            hoppetEvalcff(z2, sqrt(int_av_Q2), f);
//             COL_pmz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = f[i] / z2;

        }
        
        //calling the loop to calculate numerator and denominator of A0
        Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1, FF_ppz2, FF_pmz1, FF_pmz2);
    }
    
    if(fixed_Q2){
        
        if(myCol.evo == "DGLAP" || myCol.evo == "none"){
            //call for the Collins functions
            myCol.eval(z1, FIXEDQ2, int_charge, FF_ppz1_fixedQ2, int_par);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, FIXEDQ2, int_charge, FF_ppz2_fixedQ2, int_par);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = myCol.COL_z[i];
            
            int_charge *= -1; //to call pi- Collins
                
            myCol.eval(z1, FIXEDQ2, int_charge, FF_pmz1_fixedQ2, int_par);
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = myCol.COL_z[i];
            
            myCol.eval(z2, FIXEDQ2, int_charge, FF_pmz2_fixedQ2, int_par);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = myCol.COL_z[i];
        
        }
        
        if(myCol.evo == "CT3"){
        
            //call for the Collins functions
            hoppetEvalcf(z1, sqrt(FIXEDQ2), f);
//             COL_ppz1 = CollinsHoppetEval(z1, int_charge, f);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = f[i] / z1;

            hoppetEvalcf(z2, sqrt(FIXEDQ2), f);
//             COL_ppz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = f[i] / z2;
            
            int_charge *= -1; //to call pi- Collins
                
            hoppetEvalcff(z1, sqrt(FIXEDQ2), f);      //modified to use proper hoppet calls for pip and pim
//             COL_pmz1 = CollinsHoppetEval(z1, int_charge, f);    
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = f[i] / z1;

            hoppetEvalcff(z2, sqrt(FIXEDQ2), f);
//             COL_pmz2 = CollinsHoppetEval(z2, int_charge, f);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = f[i] / z2;

        }
        
        //calling the loop to calculate numerator and denominator of A0
        Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1_fixedQ2, FF_ppz2_fixedQ2, FF_pmz1_fixedQ2, FF_pmz2_fixedQ2);
        
    }
    
    //construct the four functions to be integrated; see 'notes-epem-A0-z1z2integrated'            
    //note: no explicit constant factor (zmax-zmin)^2 in the expressions, it simplyfies in the ratio (but it's important for the tensor charge!!)
    denFact = (pow(z2, 2) / (pow(z1, 2) + pow(z2, 2))) * exp((-pow(int_av_pt, 2) / p_width) * (pow(z2, 2) / (pow(z1, 2) + pow(z2, 2))));
    numFact = ((pow(z2, 5) * z1) / pow(pow(z1, 2) + pow(z2, 2), 3)) * exp(-(pow(int_av_pt, 2) / (rho_C * p_width)) * (pow(z2, 2) / (pow(z1, 2) + pow(z2, 2))));
    if(myCol.model == "ref-fit-unbiased-1" || myCol.model == "ref-fit-unbiased-2") numFact *= z1 * z2;
    
    
    f0 = denU * denFact;
    f1 = numU * numFact;

//     if(int_epem_sign == 1){ //A0_UL
    if(int_epem_sign == "A0UL"){ //A0_UL
        
        f2 = denL * denFact;
        f3 = numL * numFact;
    }
        
//     if(int_epem_sign == 2){ //A0_UC
    if(int_epem_sign == "A0UC"){ //A0_UC
        
        numC = numU + numL;
        denC = denU + denL;
        f2 = denC * denFact;
        f3 = numC * numFact;
    }
    
    return 0;
}

void SetTarget(const int target_in){
    
    int target = target_in;

    
    if(target == 1){         //neutron

        if(grv_flag){
        
            tmp = PDF[8];
            PDF[8] = PDF[7];
            PDF[7] = tmp;
            tmp = PDF[4];
            PDF[4] = PDF[5];
            PDF[5] = tmp;

        }
        
        if(!grv_flag){

            tmp = myPDF.thePDF[8];
            myPDF.thePDF[8] = myPDF.thePDF[7];
            myPDF.thePDF[7] = tmp;
            tmp = myPDF.thePDF[4];
            myPDF.thePDF[4] = myPDF.thePDF[5];
            myPDF.thePDF[5] = tmp;
        }

        tmp = h1[8];
        h1[8] = h1[7];
        h1[7] = tmp;
        tmp = h1[4];
        h1[4] = h1[5];
        h1[5] = tmp;
    }
        
    if(target == 2){ //deuteron
        
        if(grv_flag){
            
            tmp1 = PDF[5]; 
            tmp2 = PDF[4];
            PDF[4] += tmp1;
            PDF[5] += tmp2;
            tmp1 = PDF[8]; 
            tmp2 = PDF[7];
            PDF[7] += tmp1;
            PDF[8] += tmp2;
            PDF[3] *= 2;
            PDF[9] *= 2;

        }

        if(!grv_flag){
        
            tmp1 = myPDF.thePDF[5]; 
            tmp2 = myPDF.thePDF[4];
            myPDF.thePDF[4] += tmp1;
            myPDF.thePDF[5] += tmp2;
            tmp1 = myPDF.thePDF[8]; 
            tmp2 = myPDF.thePDF[7];
            myPDF.thePDF[7] += tmp1;
            myPDF.thePDF[8] += tmp2;
            tmp1 = myPDF.thePDF[3];
            myPDF.thePDF[3] += myPDF.thePDF[9];
            myPDF.thePDF[9] += tmp1;
        
        }

        tmp1 = h1[5]; 
        tmp2 = h1[4];
        h1[4] += tmp1;
        h1[5] += tmp2;
        tmp1 = h1[8]; 
        tmp2 = h1[7];
        h1[7] += tmp1;
        h1[8] += tmp2;
    }
    
}


// double Collins_Asy(const std::vector<double> &bin, const int hadron_in, const int charge_in, const int target_in, const string process_in, const std::vector<double> &par, const string dataname_in, const string epem_sign_in){
double Collins_Asy(const std::map<string, double> &bin, const int hadron_in, const int charge_in, const int target_in, const string process_in, const std::vector<double> &par, const string dataname_in, const string epem_sign_in){

    int hadron = hadron_in, charge = charge_in, target = target_in;
    std::string process(process_in), epem_sign(epem_sign_in);
    
    string dataname = dataname_in; 
    std::size_t found = dataname.find("comp");
    std::size_t found_SS = dataname.find("StringSpinner");
    
    // cout << dataname <<endl;
    
    double sqMC;
    
    //COLLINS 2015 fit
    if(myCol.model == "fit-2015" || myCol.model == "fit-2015-bernstein"|| myCol.model == "fit-2013-polynomial") sqMC = par[8];
    if(myCol.model == "fit-2015-bu-bd" || myCol.model == "fit-2015-full") sqMC = par[9];
    
    double D_NN = 0;
    
    tmp = 0, tmp1 = 0, tmp2 = 0;
    
//     double f[13], h1[13]; //for HOPPET and transversity
    
    output = 0, F_UU = 0, F_UT = 0, numU = 0, denU = 0, numL = 0, denL = 0, numC = 0, denC = 0;
    double Q2min, Q2max, varmin, varmax, av_x, av_y, av_z, av_pt, av_W, av_Q2;
    double pt0min, pt0max, z1min, z1max, av_z1, z2min, z2max, av_z2, ssth; //kinematics
    
    //different choices for kinematics according to data type: [0] SIDIS, [1] e+e-(PT), [2] e+e- (z1z2)

    if(process == "SIDIS"){ 
        //SIDIS data
        // Q2min = bin[0], Q2max = bin[1], varmin = bin[2], varmax = bin[3], av_x = bin[4], av_y = bin[5], av_z = bin[6], av_pt = bin[7], av_W = bin[8], av_Q2 = bin[9];
        // Q2min = bin.at("Q2min"), Q2max = bin.at("Q2max"), 
        av_x = bin.at("x"), av_y = bin.at("y"), av_z = bin.at("z"), av_pt = bin.at("pT"), av_Q2 = bin.at("Q2");
        
//         cout<<av_x<<"\t"<<av_y<<"\t"<<av_z<<"\t"<<av_pt<<endl;
    }
    
    if(process == "e+e-(PT)"){  
        //e+e- (PT) data
        // Q2min = bin[0], Q2max = bin[1], pt0min = bin[2], pt0max = bin[3], av_z1 = bin[4], av_z2 = bin[5], ssth = bin[6], av_pt = bin[7], av_W = bin[8], av_Q2 = bin[9];
        // Q2min = bin.at("Q2min"), Q2max = bin.at("Q2max"), 
        pt0min = bin.at("pt0min"), pt0max = bin.at("pt0max");
        av_z1 = bin.at("z1"), av_z2 = bin.at("z2"), ssth = bin.at("A0th"), av_pt = bin.at("pT"), av_Q2 = bin.at("Q2");
    }
    
    if(process == "e+e-(z1z2)"){  
        //e+e- (z1z2) data
        // Q2min = bin[0], Q2max = bin[1], z1min = bin[2], z1max = bin[3], av_z1 = bin[4], z2min = bin[5], z2max = bin[6], av_z2 = bin[7], ssth = bin[8], av_Q2 = bin[9];
        // Q2min = bin.at("Q2min"), Q2max = bin.at("Q2max"), 
        z1min = bin.at("z1min"), z1max = bin.at("z1max"), av_z1 = bin.at("z1");
        z2min = bin.at("z2min"), z2max = bin.at("z2max"), av_z2 = bin.at("z2") , av_Q2 = bin.at("Q2");
        
        if(epem_sign == "A0UL" || epem_sign == "A0UC") ssth = bin.at("A0th");
        if(epem_sign == "A12UL" || epem_sign == "A12UC") ssth = bin.at("A12th");

    }
    
    double av_Q = sqrt(av_Q2);
    int_av_pt = av_pt;
    int_av_Q2 = av_Q2;
    
    double uv = 0, dv = 0, ub = 0, db = 0, st = 0, gl = 0; //for grv98
    int lo = 1;
    
    if(process == "SIDIS"){ //SIDIS
        
        //call and normalize PDFs
        if(!grv_flag) myPDF.PDF_eval(av_x, av_Q2);
        
        else if(grv_flag){
            //call GRV98 from its fortran code
            grv98pa_(&lo, &av_x, &av_Q2, &uv, &dv, &ub, &db, &st, &gl); 
            PDF[0] = 0;
            PDF[1] = 0;
            PDF[2] = 0;
            PDF[3] = st;
            PDF[4] = ub;
            PDF[5] = db;
            PDF[6] = gl;
            PDF[7] = dv + db;
            PDF[8] = uv + ub;
            PDF[9] = st;
            PDF[10] = 0;
            PDF[11] = 0;
            PDF[12] = 0;

            for(int i = 0; i < PDF.size(); i++) PDF[i] /= av_x; //to scale PDFs: LHAPDF gives x*f(x)
        }
        

        //call and normalize FFs
        myFF.FF_eval(hadron, charge, av_z, av_Q2);
        for(int i = 0; i < FF.size(); i++) FF[i] = myFF.theFF[i] / av_z; //to scale FFs: routines give z*D(z)
         
        if(fixed_Q2){
//             cout<<"fixing Q^2 = "<<FIXEDQ2<<endl;
            av_Q2 = FIXEDQ2;
            av_Q = sqrt(FIXEDQ2);
//             int_av_Q2 = FIXEDQ2;
        }
        
        //Transversity
        if(myTransv.evo == "DGLAP" || myTransv.evo == "none"){
            
            myTransv.eval(av_x, av_Q2, par);
            for(int i = 0; i <= 6; i++) h1[i] = 0; //modified to set ubar and dbar transversity to zero
            for(int i = 7; i < 9; i++) h1[i] = myTransv.TRANSV_x[i]; //we already divide by x in TRANSV_x.cpp
            for(int i = 9; i < 13; i++) h1[i] = 0;
        
        }
        else if(myTransv.evo == "CT3"){ 
            
//             hoppetEval(av_x, av_Q, h1);
            hoppetEval(av_x, av_Q, f);
            for(int i = 0; i <= 6; i++) h1[i] = 0; //modified to set ubar and dbar transversity to zero
//             for(int i = 7; i < 9; i++) h1[i] /= av_x; //HOPPET gives x*f
            h1[7] = f[7] / av_x;
            h1[8] = f[8] / av_x;
            for(int i = 9; i < 13; i++) h1[i] = 0;
//             cout<<"in FCN.cpp: "<<h1[7]<<"\t"<<h1[8]<<endl;
        }
        
        //Collins z and PT dependence
        if(myCol.evo == "DGLAP" || myCol.evo == "none"){
            
            myCol.eval(av_z, av_Q2, charge, FF, par);
            for(int i = 0; i < COL_z.size(); i++) COL_z[i] = myCol.COL_z[i];
//             for(int i=0; i<COL_z.size(); i++) cout<<"asy: COL_z["<<i<<"]="<<COL_z[i]<<endl;
        
        }
        
        else if(myCol.evo == "CT3"){ //TO DO: CHECK HOPPET EVOLUTION FOR COLLINS

            if(charge == 1) hoppetEvalcf(av_z, av_Q, f);
            if(charge == -1) hoppetEvalcff(av_z, av_Q, f);
            
            for(int i = 0; i < COL_z.size(); i++) COL_z[i] = f[i] / av_z;
//             COL_z = CollinsHoppetEval(av_z, charge, f);
//             for(int i=0; i<COL_z.size(); i++) cout<<"asy: COL_z["<<i<<"]="<<COL_z[i]<<endl;

        }
        
        Collins_PT(av_z, av_Q2, par, dataname); 
        
        SetTarget(target);

        //loop on quarks
        for(int i = 3; i <= 9; i++){
            
            F_UT += pow(charges[i], 2) * h1[i] * COL_z[i];
            // F_UU += pow(charges[i], 2) * PDF[i] * FF[i];
            if(!grv_flag) F_UU += pow(charges[i], 2) * myPDF.thePDF[i] * FF[i];
            else if(grv_flag) F_UU += pow(charges[i], 2) * PDF[i] * FF[i];
        }
                
//         cout<<"before F_UT = "<<F_UT<<endl;

                
        if(myCol.model == "fit-2015" || myCol.model == "fit-2015-bernstein" || myCol.model == "fit-2015-bu-bd" || myCol.model == "fit-2015-full" || myCol.model == "fit-2013-polynomial"){ //fit 2015
            
            F_UT *= (sqrt(2.0 * EulerConst)) * (av_pt / sqrt(sqMC)) * (1.0 - av_y) * (pow(pperpC_width, 2) / p_width) * (exp(-pow(av_pt, 2) / PTC_width) / (pow(PTC_width, 2))); 
//             cout<<"after F_UT = "<<F_UT<<endl;
        }
        if(myCol.model == "ref-fit-unbiased-1" || myCol.model == "ref-fit-unbiased-2"){ //fit 2018 - REWEIGHTING - unbiased
            
            F_UT *= (4.0 * av_z * av_pt * Mpi) * (1.0 - av_y) * (exp(-pow(av_pt, 2) / PTC_width) / (pow(PTC_width, 2))); 
        }
            
        F_UU *= (1.0 + pow((1.0 - av_y), 2)) * (exp(-pow(av_pt, 2) / PT_width) / (PT_width));
        
        D_NN = 2.0 * (1.0 - av_y) / (1.0 + pow(1.0 - av_y, 2));
        
        output = F_UT / F_UU;
        
        if(found != std::string::npos && found_SS == std::string::npos) output *= - 1 / D_NN; //COMPASS \pi mismatch (they measure sin(\phi_h+\phi_S +- \pi) that is -sin(\phi_h+\phi_S)) and divide out the depolarisation factor        
        
//         if(found_SS != std::string::npos) output *= - 1; //COMPASS \pi mismatch (they measure sin(\phi_h+\phi_S +- \pi) that is -sin(\phi_h+\phi_S)) 

    }
    
    
    if(process == "e+e-(PT)"){ //e+e-(PT)
        
        Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
            EPSREL, EPSABS, VERBOSE | LAST,
            MINEVAL, MAXEVAL, KEY,
            STATEFILE, SPIN,
            &nregions, &neval, &fail, integral, error, prob);
        
        if(fail != 0){
            
            cout<<"CUHRE RESULT:\tnregions "<<nregions<<"\tneval "<<neval<<"\tfail "<<fail<<endl;
//         for( comp = 0; comp < NCOMP; comp++ ) cout<<"CUHRE RESULT:\t"<<integral[comp]<<" +- "<<error[comp]<<"\tp = "<<prob[comp]<<endl;
            cout<<endl;
            for (int i = 0; i < par.size(); i++) cout<<par[i]<<"\t";
            cout<<endl;
        }
        
        output = (integral[1] / integral[0]) - (integral[3] / integral[2]);
        
        //multiply for the overall factor; see 'notes-epem-A0-z1z2integrated'        
        if(myCol.model == "fit-2015" || myCol.model == "fit-2015-bernstein" || myCol.model == "fit-2015-bu-bd" || myCol.model == "fit-2015-full" || myCol.model == "fit-2013-polynomial"){ 
            
             output *= (EulerConst / 2) * ssth * (pow(av_pt, 2) / (sqMC + p_width)); //z1-z2-integrated asymmetries - fit2015
        }
        
        if(myCol.model == "ref-fit-unbiased-1" || myCol.model == "ref-fit-unbiased-2"){
            
            output *= 4.0 * pow(Mpi, 2) * pow(av_pt, 2) * ssth * (p_width / pow(pperpC_width, 3)); //z1-z2-integrated asymmetries - fit2018 - REWEIGHTING - unbiased
        }
//         cout<<"A0_PT = "<<output<<endl;
    }
    
    
    if(process == "e+e-(z1z2)"){ //e+e-(z1z2)
         //Asymmetry expressions taken from  arXiv:1809.09500v1 [hep-ph]
        Collins_PT(av_z1, av_Q2, par, dataname); //note: av_z1 is not used here, just needed to call the function to fix FC_ee factor

        //call for the Fragmentation Functions
        // Collins_FF(iset_ff, hadron, charge, order, av_z1, av_z2, av_Q2);
        Collins_FF(hadron, charge, av_z1, av_z2, av_Q2);
                
        //call for the Collins functions
        if(myCol.evo == "DGLAP" || myCol.evo == "none"){ 
            
            myCol.eval(av_z1, av_Q2, charge, FF_ppz1, par);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = myCol.COL_z[i];
            myCol.eval(av_z2, av_Q2, charge, FF_ppz2, par);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = myCol.COL_z[i];
            
            charge *= -1; //to call pi- Collins
            
            myCol.eval(av_z1, av_Q2, charge, FF_pmz1, par);
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = myCol.COL_z[i];
            myCol.eval(av_z2, av_Q2, charge,  FF_pmz2, par);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = myCol.COL_z[i];
        }
        
        if(myCol.evo == "CT3"){
                    
            //call for the Collins functions
            hoppetEvalcf(av_z1, av_Q, f);
//             COL_ppz1 = CollinsHoppetEval(av_z1, charge, f);
            for(int i = 0; i < COL_ppz1.size(); i++) COL_ppz1[i] = f[i] / av_z1;

            hoppetEvalcf(av_z2, av_Q, f);
//             COL_ppz2 = CollinsHoppetEval(av_z2, charge, f);
            for(int i = 0; i < COL_ppz2.size(); i++) COL_ppz2[i] = f[i] / av_z2;
            
            charge *= -1; //to call pi- Collins
//             evo_charge = charge;
                
            hoppetEvalcff(av_z1, av_Q, f);
//             COL_pmz1 = CollinsHoppetEval(av_z1, charge, f);
            for(int i = 0; i < COL_pmz1.size(); i++) COL_pmz1[i] = f[i] / av_z1;

            hoppetEvalcff(av_z2, av_Q, f);
//             COL_pmz2 = CollinsHoppetEval(av_z2, charge, f);
            for(int i = 0; i < COL_pmz2.size(); i++) COL_pmz2[i] = f[i] / av_z2;
        
        }
        
        //calling the loop to calculate numerator and denominator of A0
        Collins_epem_loop(COL_ppz1, COL_ppz2, COL_pmz1, COL_pmz2, FF_ppz1, FF_ppz2, FF_pmz1, FF_pmz2);
        
//         if(epem_sign == 1){ //A0_UL
        if(epem_sign == "A0UL" || epem_sign == "A12UL"){ //A0_UL
            
            output = (numU / denU) - (numL / denL);
        }
        
//         if(epem_sign == 2){ //A0_UC
        if(epem_sign == "A0UC" || epem_sign == "A12UC"){ //A0_UC
            
            numC = numU + numL;
            denC = denU + denL;
            output = (numU / denU) - (numC / denC);
        }
        
        //multiply for the overall factor; see eqq. (47)-(50) of Phys.Rev. D98 (2018) no.9, 094023 (arXiv:1809.09500v1 [hep-ph])
        //Asymmetry expressions taken from  arXiv:1809.09500v1 [hep-ph]
        if(myCol.model == "fit-2015" || myCol.model == "fit-2015-bernstein" || myCol.model == "fit-2015-bu-bd" || myCol.model == "fit-2015-full" || myCol.model == "fit-2013-polynomial"){

//             Collins_PT(av_z1, av_Q2, par, dataname); //note: av_z1 is not used here, just needed to call the function to fix FC_ee factor
         
            if(epem_sign == "A0UL" || epem_sign == "A0UC") output *= (1./4.) * ssth * ((av_z1 * av_z2) / (pow(av_z1, 2) + pow(av_z2, 2))) * FC_ee; //PT1-integrated A0 asymmetries 
            if(epem_sign == "A12UL" || epem_sign == "A12UC"){
                
                output *= (1./4.) * ssth * pi * EulerConst / 2. * pow(pperpC_width, 3) / (sqMC * pow(p_width, 2)); //PT1-integrated A12 asymmetries  //update 19/03/24: A12 expressions should be ok now // update 15/07/24 factor 1/4 was missing (def. of A12 asymmetry)
//                 if(dataname.find("belle_2019") != std::string::npos) output *= 4;
                
            }
        }
        
        if(myCol.model == "ref-fit-unbiased-1" || myCol.model == "ref-fit-unbiased-2"){
            
            output *= (4.0 * pow(Mpi, 2) / pperpC_width) * ssth * (pow(av_z1 * av_z2, 2) / (pow(av_z1, 2) + pow(av_z2, 2))); //PT1-integrated asymmetries
        }
//         cout<<"A0_z1z2 = "<<output<<"\tdata = "<<bin[10]<<"\t spread = "<<abs(output - bin[10])<<endl;
    }
    
    
    
    return output;
    
}

double Collins_Asy_loop(DATA::Data &data, const std::vector<double> &par){

    double tot_chi2;
    int nrow = data.Nrow; 
    int_charge = data.charge;
    int_hadron = data.hadron;
    int_fname = data.fname;
    int_process = data.the_process;
    int_epem_sign = data.the_epem_sign;

    Npoints += data.npts;
    data.is_fitted = true;
    
//     cout << data.name << "\t" << data.has_normalization << "\t N = " <<data.normalization<< "\tdelta_N = " <<data.delta_N  << endl;

    // cout << "i_rep = "<< i_rep <<endl;
//     cout<<"In asy_loop: evo_proc = "<<evo_process<<"\t evo_ch = "<<evo_charge<<"\tevo_had = "<<evo_hadron<<"\n";

    if(int_process == "e+e-(PT)"){ //e+e- PT (to get parameters for the numerical integration)

        int_par.resize(par.size());
        int_par.shrink_to_fit();
    
        for(int i = 0; i < int_par.size(); i++){ 
            int_par[i] = par[i];
            if(fabs(int_par[i] - par[i]) > EPSILON) cout<<"WARNING in copying parameters:\nint_par["<<i<<"]="<<int_par[i]<<"\tpar["<<i<<"]="<<par[i]<<endl;    
        }
    }

    if(!hessian_flag){
        
        for(int i = 0; i < data.Nrow; i++){
            

            if(data.valid[i]){

                // data.getRow(i);
                data.getRowMap(i);
                
//                 data.theory[i] = Collins_Asy(data.bin, data.hadron, data.charge, data.target, data.process, par, data.fname, data.epem_sign);
                // data.theory[i] = Collins_Asy(data.bin, data.hadron, data.charge, data.target, int_process, par, data.fname, int_epem_sign);
                data.theory[i] = Collins_Asy(data.bin_map, data.hadron, data.charge, data.target, int_process, par, data.fname, int_epem_sign);
                if(mode == "replicas_bands" || mode == "MC_bands") data.theory_rep_vals[i][i_rep] = data.theory[i];
                if(data.theory[i] > data.theory_max[i]) data.theory_max[i] = data.theory[i];
                if(data.theory[i] < data.theory_min[i]) data.theory_min[i] = data.theory[i];
            }
        }
        
        for(int i = 0; i < data.theory.size(); i++) data.theory_0[i] = data.theory[i];
        data.eval_chi2();
    }

    if(hesse_fit)  check_fit = hessian_flag;
    if(!hesse_fit) check_fit = hessian_flag && i_hesse <= (par.size() - NfixedPars);

    if(check_fit){ //Hessian method from Eur.Phys.J.C 75 (2015) 132 (LHAPDF6)
    
        for(int i = 0; i < data.Nrow; i++){
            
            if(data.valid[i]){
                
                // data.getRow(i);
                data.getRowMap(i);
                
                if(hessian_sign == 1){
                    
                    data.theory_plus[i] = Collins_Asy(data.bin_map, data.hadron, data.charge, data.target,  int_process, par, data.fname, int_epem_sign);
                    data.delta_i0_plus[i] = data.theory_plus[i] - data.theory_0[i]; 
                    
                    data.eval_chi2(data.theory_plus);
                }
                
                if(hessian_sign == -1){
                    
                    data.theory_minus[i] = Collins_Asy(data.bin_map, data.hadron, data.charge, data.target,  int_process, par, data.fname, int_epem_sign);
                    data.delta_i0_minus[i] = data.theory_minus[i] - data.theory_0[i];

                    data.eval_chi2(data.theory_minus);
                }
                
                if(hessian_count % 2 == 0 && hessian_count != 0){
                    
                    //asymmetric errors
                    data.theory_hesse_max[i] += pow(maximum(data.delta_i0_plus[i], data.delta_i0_minus[i], 0.0), 2); 
                    data.theory_hesse_min[i] += pow(maximum(-data.delta_i0_plus[i], -data.delta_i0_minus[i], 0.0), 2);  //delta_0i = -delta_i0
                    cout<<data.theory_hesse_min[i]<<"\t"<<data.theory_hesse_max[i]<<endl;
                    //symmetric errors
//                     data.theory_hesse_max[i] += pow(data.theory_plus[i] - data.theory_minus[i], 2); 
//                     data.theory_hesse_min[i] += pow(data.theory_plus[i] - data.theory_minus[i], 2);
                }
                                
            }
        }
            
        if(hessian_count == 2 * (par.size() - NfixedPars)){

            for(int i = 0; i < data.theory.size(); i++){
//                                 
                //asymmetric errors
                data.theory_hesse_max[i] = data.theory_0[i] + sqrt(data.theory_hesse_max[i]);
                data.theory_hesse_min[i] = data.theory_0[i] - sqrt(data.theory_hesse_min[i]);
                
                //symmetric errors
//                 data.theory_hesse_max[i] = data.theory_0[i] + 0.5 * sqrt(data.theory_hesse_max[i]);
//                 data.theory_hesse_min[i] = data.theory_0[i] - 0.5 * sqrt(data.theory_hesse_min[i]);
            }
        }
        
    }

//     if (data.has_normalization) cout<<"in main loop in FCN: chi2 = "<<data.tot_chi2<<endl;

    return data.tot_chi2;
}

namespace ROOT{

namespace Minuit2{

double FCN::operator()(const std::vector<double>& par) const{
    
    icount++;
    ivalid = 0, Npoints = 0;
    double output = 0.0;

    // for(int i = 0; i < par.size(); i++) cout<<par[i]<<"\t";
    // cout<<endl;
    
    if(hessian_flag){ 
        
        hessian_count++; 
        if (hessian_count % 2 == 0 && i_hesse <= par.size() - NfixedPars) i_hesse++;
    }

    if(myTransv.evo ==  "CT3" || myCol.evo == "CT3")
    {
        evo_flag = true;
        evo_par.resize(par.size());
        evo_par.shrink_to_fit();
//         hoppetAssign(transv_subroutine);
    
        for(int i = 0; i < evo_par.size(); i++){ 
            
            evo_par[i] = par[i];
            if(fabs(evo_par[i] - par[i]) > EPSILON) cout<<"WARNING in copying parameters:\nevo_par["<<i<<"]="<<evo_par[i]<<"\tpar["<<i<<"]="<<par[i]<<endl;
//             cout<<evo_par[i]<<"\t";
        }
//         cout<<endl;
    }
    
    //evolve functions and create grids for HOPPET
    if(evo_flag){
        
        if(myTransv.evo == "CT3") hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, transv_subroutine, Q0pdf);
//         if(myTransv.transv_evo == 1) hoppetCachedEvolve(transv_subroutine);
        // && evo_flag){ 
        if(myCol.evo == "CT3"){ 
            
            hoppetEvolvecf(asQ0, Q0alphas, nloop, muR_Q, Collins_subroutine_pip, Q0pdf);
            hoppetEvolvecff(asQ0, Q0alphas, nloop, muR_Q, Collins_subroutine_pim, Q0pdf);
        
        }
        
        evo_flag = false;
    }
    
//     if(poly_flag){
    //add new definition for different data files  

// //JLab 2011 neutron     
    output += Collins_Asy_loop(JLab11_Npp_x, par);
    output += Collins_Asy_loop(JLab11_Npm_x, par);
    
//  //HERMES 2010 proton     
//     output += Collins_Asy_loop(Herm10_pip_x, par); 
//     output += Collins_Asy_loop(Herm10_pip_z, par); 
//     output += Collins_Asy_loop(Herm10_pip_pt, par);
//     output += Collins_Asy_loop(Herm10_pim_x, par); 
//     output += Collins_Asy_loop(Herm10_pim_z, par); 
//     output += Collins_Asy_loop(Herm10_pim_pt, par);

// // HERMES 2020 proton
    output += Collins_Asy_loop(Herm20_pip_x_1, par);
    output += Collins_Asy_loop(Herm20_pip_x_2, par);
    output += Collins_Asy_loop(Herm20_pip_x_3, par);
    output += Collins_Asy_loop(Herm20_pip_x_4, par);
    
    output += Collins_Asy_loop(Herm20_pim_x_1, par);
    output += Collins_Asy_loop(Herm20_pim_x_2, par);
    output += Collins_Asy_loop(Herm20_pim_x_3, par);
    output += Collins_Asy_loop(Herm20_pim_x_4, par);

//     //COMPASS 2012 proton    
//     output += Collins_Asy_loop(Comp12_Ppp_x, par); 
//     output += Collins_Asy_loop(Comp12_Ppp_z, par); 
//     output += Collins_Asy_loop(Comp12_Ppp_pt, par);
//     output += Collins_Asy_loop(Comp12_Ppm_x, par); 
//     output += Collins_Asy_loop(Comp12_Ppm_z, par); 
//     output += Collins_Asy_loop(Comp12_Ppm_pt, par);
    
// //COMPASS 2012 Deuteron     
    output += Collins_Asy_loop(Comp12_Dpp_x, par); 
    output += Collins_Asy_loop(Comp12_Dpp_z, par); 
    output += Collins_Asy_loop(Comp12_Dpp_pt, par);
    output += Collins_Asy_loop(Comp12_Dpm_x, par); 
    output += Collins_Asy_loop(Comp12_Dpm_z, par); 
    output += Collins_Asy_loop(Comp12_Dpm_pt, par);

// //COMPASS 2014 proton     
    output += Collins_Asy_loop(Comp14_Ppp_x, par); 
    output += Collins_Asy_loop(Comp14_Ppp_z, par); 
    output += Collins_Asy_loop(Comp14_Ppp_pt, par);
    output += Collins_Asy_loop(Comp14_Ppm_x, par); 
    output += Collins_Asy_loop(Comp14_Ppm_z, par); 
    output += Collins_Asy_loop(Comp14_Ppm_pt, par);
    
//     // COMPASS 2017 proton 
//     output += Collins_Asy_loop(Comp17_02hp_x, par); 
//     output += Collins_Asy_loop(Comp17_02hp_z, par); 
//     output += Collins_Asy_loop(Comp17_02hp_pt, par);
//     output += Collins_Asy_loop(Comp17_02hm_x, par); 
//     output += Collins_Asy_loop(Comp17_02hm_z, par); 
//     output += Collins_Asy_loop(Comp17_02hm_pt, par);

//  //BELLE 2012
    output += Collins_Asy_loop(Belle12_A0UL_z1z2, par);
    output += Collins_Asy_loop(Belle12_A0UC_z1z2, par);
//     output += Collins_Asy_loop(Belle12_A12UL_z1z2, par);
//     output += Collins_Asy_loop(Belle12_A12UC_z1z2, par);

//     //BABAR 2013 
//     output += Collins_Asy_loop(Babar13_A0UL_z1z2, par); 
//     output += Collins_Asy_loop(Babar13_A0UC_z1z2, par); 
    // output += Collins_Asy_loop(Babar13_A12UL_z1z2, par);
    // output += Collins_Asy_loop(Babar13_A12UC_z1z2, par);

// //BABAR 2014
    output += Collins_Asy_loop(Babar14_A0UL_PT, par); 
    output += Collins_Asy_loop(Babar14_A0UC_PT, par);
    
// //BABAR 2015     
    output += Collins_Asy_loop(Babar15_A0UL_z1z2, par); 
    output += Collins_Asy_loop(Babar15_A0UC_z1z2, par); 
//     output += Collins_Asy_loop(Babar15_A12UC_z1z2, par);
//     output += Collins_Asy_loop(Babar15_A12UL_z1z2, par); 

//     //BESIII 2016 A0 PT-dependent
    output += Collins_Asy_loop(BESIII16_A0UL_PT, par);
    output += Collins_Asy_loop(BESIII16_A0UC_PT, par);
    
// //  BELLE 2019 A12
//     output += Collins_Asy_loop(Belle19_A12UL_z1z2, par);
//     output += Collins_Asy_loop(Belle19_A12UC_z1z2, par);
    
//     //STRING-SPINNER pseudo-data
//     output += Collins_Asy_loop(Comp_14_Ppm_x_SS, par);
//     output += Collins_Asy_loop(Comp_14_Ppm_z_SS, par);
//     output += Collins_Asy_loop(Comp_14_Ppm_pt_SS, par);
//     output += Collins_Asy_loop(Comp_14_Ppm_x_SS, par);
//     output += Collins_Asy_loop(Comp_14_Ppm_z_SS, par);
//     output += Collins_Asy_loop(Comp_14_Ppm_pt_SS, par);
//     
//     output += Collins_Asy_loop(Belle12_A0UL_z1z2_SS, par);
//     output += Collins_Asy_loop(Belle12_A0UC_z1z2_SS, par);
    
//     //COMPASS 2024 deuteron
    // output += Collins_Asy_loop(Comp24_hp_x, par); 
    // output += Collins_Asy_loop(Comp24_hp_z, par); 
    // output += Collins_Asy_loop(Comp24_hp_pt, par);
    // output += Collins_Asy_loop(Comp24_hm_x, par); 
    // output += Collins_Asy_loop(Comp24_hm_z, par); 
    // output += Collins_Asy_loop(Comp24_hm_pt, par);

    
// //     }

    // chi2_iter<<icount<<"\t"<<output<<endl;
//     transv_pars<<par[0]<<"\t"<<par[1]<<"\t"<<par[2]<<"\t"<<par[3]<<"\t"<<output<<endl;
//     poly_flag=true;
    // cout<<"beta_u = "<<par[3]<<"  beta_d = "<<par[4];
    cout<<"\n#chi2="<<output;
    
    if(!hessian_flag){
        
        if(myCol.evo != "CT3" && myTransv.evo != "CT3") cout<<"\tcount="<<icount<<endl;
        if(myCol.evo == "CT3" || myTransv.evo == "CT3") cout<<"\tevo_count="<<icount<<endl;
        else cout<<endl;
    }
    else cout<<"\thessian_count = "<<hessian_count<<endl<<endl;
//     iweights++;
//     cout<<"\tiweights = "<<iweights;

    return output;
  

}



}


}
