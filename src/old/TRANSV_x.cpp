//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#include "TRANSV_x.h"

    extern bool dssv_init, grv_flag;
    extern COLLPDF::CollPDF myPDF;
    // extern LHAPDF::PDF* pdf;

namespace TRANSVX{

using namespace std;

vector<double> getTRANSV_x::eval_SB(const double &av_x, const double &Q2){
    
    int iset_grsv = 2; //LO grids for Soffer bound GRSV
//     int iset_grsv = 1; //NLO grids for Soffer bound GRSV
    
    vector<double> PDF(13), SB(13);
    double u, d, ub, db, st; //for parpolt_
    double uv, dv, ubar, dbar, strange, gl; //for DSSV
    double x = av_x, sqQ = Q2;
    
//     if(dssv_init) grv_flag = false;
//     if(grv_flag) dssv_init = false;
    
    if(dssv_init && !grv_flag){ 
//         cout<<"initialising DSSV\n";
        dssvini_();
        dssv_init = false;
    }
//     
    if(!grv_flag){
//         cout<<"creating SB with DSSV\n";
        myPDF.PDF_eval(x, sqQ);
        for (int i = 0; i < PDF.size(); i++) PDF[i] = x * myPDF.thePDF[i]; //NOTE: here we need x*f(x) to build the SB
        dssvfit_(&x, &sqQ, &uv, &dv, &ub, &db, &st, &gl);
    }
//     cout<<x<<"\t"<<sqQ<<"\n";
//     <<uv<<"\t"<<dv<<"\t"<<ubar<<"\t"<<dbar<<"\t"<<strange<<"\t"<<gl<<endl;

    if(grv_flag) parpolt_(&iset_grsv, &x, &sqQ, &u, &d, &ub, &db, &st);
    
//     cout<<"CTEQ6L1 + DSSV: x = "<<x<<"\tQ2 = "<<sqQ<<"\t uv = "<<(uv+(PDF[8]-PDF[4]))/2<<"\t dv = "<<(dv+(PDF[7]-PDF[5]))/2<<endl;
//     cout<<"PARPOLT:        x = "<<x<<"\tQ2 = "<<sqQ<<"\t uv = "<<u-ub<<"\t dv = "<<d-db<<endl;
//     cout<<"RATIO (CTQ6L1+DSSV)/PARPOLT:\t\t uv = "<<(uv+(PDF[8]-PDF[4]))/(2*(u-ub))<<"\t dv = "<<(dv+(PDF[7]-PDF[5]))/(2*(d-db))<<"\n\n";
    
    SB[0] = 0;
    SB[1] = 0;
    SB[2] = 0;
    SB[3] = st;
    SB[4] = ub;
    SB[5] = db;
    SB[6] = 0;
    
    if(grv_flag){
//         cout<<"SB with GRV and GRSV\n";
        SB[7] = d; //SB with GRV and GRSV
        SB[8] = u; //SB with GRV and GRSV
    }
    else if(!grv_flag){
//         cout<<"SB with DSSV\n";
        SB[7] = (dv + (PDF[7] - PDF[5])) / 2.0; //dval SB DSSV
        SB[8] = (uv + (PDF[8] - PDF[4])) / 2.0; //uval SB DSSV
    }
    
    SB[9] = st;
    SB[10] = 0;
    SB[11] = 0;
    SB[12] = 0;
    
    for(int i = 0; i < SB.size(); i++) SB[i] /= x; //parpolt_ gives the Soffer bound multiplied by x!!! (checked with Alexei and Osvaldo)
//     for(int i=0; i<h1.size(); i++) cout<<"h1["<<i<<"]="<<h1[i]<<"\t";
//     cout<<endl;

    return SB;
}

void getTRANSV_x::eval(const double & av_x, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par){

    int model = model_in, evo = evo_in;
    double x = av_x, sqQ = Q2;
    
    SB_x = eval_SB(x, sqQ);
    
    if(model == 0 || model == 1){ //fit 2015 ([1] for Bernstein polynomials in Collins function)
        
        double NT_uval = par[0], NT_dval = par[1], alpha = par[2], beta = par[3];
        
        TRANSV_x[0] = 0;
        TRANSV_x[1] = 0;
        TRANSV_x[2] = 0;
        TRANSV_x[3] = 0;
        TRANSV_x[4] = 0;
        TRANSV_x[5] = 0;
        TRANSV_x[6] = 0;
        if(grv_flag){
//             cout<<"transversity with GRV and GRSV\n";
//             TRANSV_x[7] = (SB_x[7] - SB_x[5]) * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for grv
//             TRANSV_x[8] = (SB_x[8] - SB_x[4]) * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for grv
            TRANSV_x[7] = (SB_x[7] - SB_x[5]) * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for grv
            TRANSV_x[8] = (SB_x[8] - SB_x[4]) * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for grv
        }
        else if(!grv_flag){          
//             cout<<"transversity with DSSV\n";
//             TRANSV_x[7] = SB_x[7] * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for DSSV
//             TRANSV_x[8] = SB_x[8] * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta);// * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for DSSV
            TRANSV_x[7] = SB_x[7] * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta) * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for DSSV
            TRANSV_x[8] = SB_x[8] * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta) * (pow(alpha + beta, alpha + beta) / (pow(alpha, alpha) * pow(beta, beta))); //for DSSV
        }
        
        TRANSV_x[9] = 0;
        TRANSV_x[10] = 0;
        TRANSV_x[11] = 0;
        TRANSV_x[12] = 0;
//         for(int i=0; i<TRANSV_x.size(); i++) cout<<"transv["<<i<<"]="<<TRANSV_x[i]<<endl;

    }
    
    
    if(model == 100){ //fit 2015 beta_u != beta_d
        
        double NT_uval = par[0], NT_dval = par[1], alpha = par[2], beta_uv = par[3], beta_dv = par[4];
        
//         beta_dv = beta_uv;
        
        TRANSV_x[0] = 0;
        TRANSV_x[1] = 0;
        TRANSV_x[2] = 0;
        TRANSV_x[3] = 0;
        TRANSV_x[4] = 0;
        TRANSV_x[5] = 0;
        TRANSV_x[6] = 0;
        
        if(grv_flag){
//             cout<<"transversity with GRV and GRSV\n";
            TRANSV_x[7] = (SB_x[7] - SB_x[5]) * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta_dv); //for grv
            TRANSV_x[8] = (SB_x[8] - SB_x[4]) * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta_uv); //for grv
        }
        else if(!grv_flag){          
//             cout<<"transversity with DSSV\n";
            TRANSV_x[7] = SB_x[7] * NT_dval * pow(av_x, alpha) * pow(1. - av_x, beta_dv); //for DSSV
            TRANSV_x[8] = SB_x[8] * NT_uval * pow(av_x, alpha) * pow(1. - av_x, beta_uv); //for DSSV
        }
        
        TRANSV_x[9] = 0;
        TRANSV_x[10] = 0;
        TRANSV_x[11] = 0;
        TRANSV_x[12] = 0;

//         cout<<"in TRANSV_x.cpp: ";
//         for(int i=7; i<9; i++) cout<<"transv["<<i<<"]="<<TRANSV_x[i]<<"\t";
//         cout<<endl;

    }
    
     if(model == 10){ 
         
        double NT_u = par[0], NT_d = par[1], a_u = par[2], a_d = par[3], b_u = par[4], b_d = par[5];
        
//         a_d = a_u, b_d = b_u;
        
        TRANSV_x[0] = 0;
        TRANSV_x[1] = 0;
        TRANSV_x[2] = 0;
        TRANSV_x[3] = 0;
        TRANSV_x[4] = 0;
        TRANSV_x[5] = 0;
        TRANSV_x[6] = 0;
        if(grv_flag){
//             cout<<"transversity with GRV and GRSV\n";
            TRANSV_x[7] = (SB_x[7] - SB_x[5]) * NT_d * pow(av_x, a_d) * pow(1. - av_x, b_d) * (pow(a_d + b_d, a_d + b_d) / (pow(a_d, a_d) * pow(b_d, b_d))); //for grv
            TRANSV_x[8] = (SB_x[8] - SB_x[4]) * NT_u * pow(av_x, a_u) * pow(1. - av_x, b_u) * (pow(a_u + b_u, a_u + b_u) / (pow(a_u, a_u) * pow(b_u, b_u))); //for grv
        }
        else if(!grv_flag){          
//             cout<<"transversity with DSSV\n";
            TRANSV_x[7] = SB_x[7] * NT_d * pow(av_x, a_d) * pow(1. - av_x, b_d) * (pow(a_d + b_d, a_d + b_d) / (pow(a_d, a_d) * pow(b_d, b_d))); //for DSSV
            TRANSV_x[8] = SB_x[8] * NT_u * pow(av_x, a_u) * pow(1. - av_x, b_u) * (pow(a_u + b_u, a_u + b_u) / (pow(a_u, a_u) * pow(b_u, b_u))); //for DSSV
        }
        
        TRANSV_x[9] = 0;
        TRANSV_x[10] = 0;
        TRANSV_x[11] = 0;
        TRANSV_x[12] = 0;
//         for(int i=0; i<TRANSV_x.size(); i++) cout<<"transv["<<i<<"]="<<TRANSV_x[i]<<endl;

    }
    
    if(model == 2 || model == 3){ //fit 2018 - REWEIGHTING
        
        double NT_u = par[0], NT_d = par[1], a_u = par[2], a_d = par[3], b_u = par[4], b_d = par[5];
        
//         a_d = a_u, b_d = b_u;
        
        TRANSV_x[0] = 0;
        TRANSV_x[1] = 0;
        TRANSV_x[2] = 0;
        TRANSV_x[3] = 0;
        TRANSV_x[4] = 0;
        TRANSV_x[5] = 0;
        TRANSV_x[6] = 0;
//         TRANSV_x[7] = NT_d*pow(av_x,a_d)*pow(1.-av_x,b_d)*(pow(a_d+b_d,a_d+b_d)/(pow(a_d,a_d)*pow(b_d,b_d)));
//         TRANSV_x[8] = NT_u*pow(av_x,a_u)*pow(1.-av_x,b_u)*(pow(a_u+b_u,a_u+b_u)/(pow(a_u,a_u)*pow(b_u,b_u)));
        TRANSV_x[7] = NT_d * pow(av_x, a_d) * pow(1. - av_x, b_d);
        TRANSV_x[8] = NT_u * pow(av_x, a_u) * pow(1. - av_x, b_u);
//         TRANSV_x[7] = NT_d*av_x*((1-a_d-b_d) +a_d*av_x+b_d*pow(av_x,2));
//         TRANSV_x[8] = NT_u*av_x*((1-a_u-b_u) +a_u*av_x+b_u*pow(av_x,2));
        TRANSV_x[9] = 0;
        TRANSV_x[10] = 0;
        TRANSV_x[11] = 0;
        TRANSV_x[12] = 0;
//         for(int i=0; i<TRANSV_x.size(); i++) cout<<"transv["<<i<<"]="<<TRANSV_x[i]<<endl;

    }
    
    for(int i = 0; i < TRANSV_x.size(); i++){
        
        if(TRANSV_x[i] > max_TRANSV_x[i]) max_TRANSV_x[i] = TRANSV_x[i];
        if(TRANSV_x[i] < min_TRANSV_x[i]) min_TRANSV_x[i] = TRANSV_x[i];
    }
    
}

void getTRANSV_x::getPolyNorm(const double & A, const double & B){
        
    double a = A, b = B;
    double xp = 0.0, xm = 0.0, Fxp = 0.0, Fxm = 0.0;
    
    if(pow(a, 2) - 3. * b * (1. - a - b) < 0){
        
        PolyNorm = 1;
        cout<<"Warning: imaginary root\n";
    }
    
    if(pow(a, 2) - 3. * b * (1. - a - b) >= 0){
        
        xp = (1. / (3. * b)) * (- a + sqrt(pow(a, 2) - 3. * b * (1. - a - b)));
        xm = (1. / (3. * b)) * (- a - sqrt(pow(a, 2) - 3. * b * (1. - a - b)));
        
        Fxp = xp * (1. - a - b + a * xp + b * pow(xp, 2));
        Fxm = xm * (1. - a - b + a * xm + b * pow(xm, 2));
    
        PolyNorm = maximum(1, abs(Fxp), abs(Fxm));
    }
}

double getTRANSV_x::maximum(const double & a, const double & b){
    
    return a < b ? b : a;
}

double getTRANSV_x::maximum(const double & a, const double & b, const double & c){
    
    return maximum(maximum(a, b), c);
    
}

void getTRANSV_x::plot(const std::string outname_in, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par){

//     cout.precision(4);
    
    ofstream transv(outname_in);
    
    double f[13];
    
    double Q = sqrt(Q2);
    
    vector<double> x_vals = {.001, .002, .003, .004, .005, .006, .007, .008, .009, .01, .015, .02, .025, .03, .035, .04, .045, .05, .055, .06, .065, .07, .075, .08, .085, .09, .095, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, .999};    
    
    cout<<"printing transversity on "<<outname_in<<endl;
    if(!transv.is_open()) cerr<<"Error in opening "<<outname_in<<endl;
    
    transv<<"#Q^2="<<Q2<<" GeV^2"<<endl;
    transv<<"#x\tDeltau\t\tDeltad\t\tDeltaub\t\tDeltadb\t\tDeltast\n";
    
    for(int i = 0; i < x_vals.size(); i++){
        
        if(evo_in == 0){
            eval(x_vals[i], Q2, model_in, evo_in, par);
            transv<<x_vals[i]<<"\t"<<TRANSV_x[8]<<"\t"<<TRANSV_x[7]<<"\t"<<TRANSV_x[4]<<"\t"<<TRANSV_x[5]<<"\t"<<TRANSV_x[9]<<endl;
        }
        
        if(evo_in == 1){
            hoppetEval(x_vals[i], Q, f);
            transv<<x_vals[i]<<"\t"<<f[8] / x_vals[i]<<"\t"<<f[7] / x_vals[i]<<"\t"<<f[4] / x_vals[i]<<"\t"<<f[5] / x_vals[i]<<"\t"<<f[9] / x_vals[i]<<endl;
        }
    
    }
    
    transv.close();
    
}
    
void getTRANSV_x::plot_SB(const std::string outname_in, const double & Q2, const int & model_in, const int & evo_in, const std::vector<double> &par){

//     cout.precision(4);
    
    ofstream sb(outname_in);
    
    double f[13];
    
    double Q = sqrt(Q2);
    
    vector<double> x_vals = {.001, .002, .003, .004, .005, .006, .007, .008, .009, .01, .015, .02, .025, .03, .035, .04, .045, .05, .055, .06, .065, .07, .075, .08, .085, .09, .095, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, .999};    
    
    cout<<"printing Soffer Bound on "<<outname_in<<endl;
    if(!sb.is_open()) cerr<<"Error in opening "<<outname_in<<endl;
    
    sb<<"Collinear PDF set: "<<myPDF.pdf->set().description()<<endl;
    sb<<"LHAPDF ID: "<<myPDF.pdf->lhapdfID()<<", member: "<<myPDF.pdf->memberID()<<endl;
    sb<<"#Soffer Bound: lhapdfID:"<<myPDF.pdf->lhapdfID()<<" + DSSV\n";
    sb<<"#Q^2 = "<<Q2<<" GeV^2"<<endl;
    sb<<"#x\tSBu\t\tSBd\t\tSBub\t\tSBdb\t\tSBst\n";
    
    for(int i = 0; i < x_vals.size(); i++){
        
        SB_x = eval_SB(x_vals[i], Q2);
        sb<<x_vals[i]<<"\t"<<SB_x[8]<<"\t"<<SB_x[7]<<"\t"<<SB_x[4]<<"\t"<<SB_x[5]<<"\t"<<SB_x[9]<<endl;

    }
    
    sb.close();
    
}    
    
    
}
