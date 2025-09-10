//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#include "COL_z.h"

    extern int iset_ff;

namespace COLZ{

using namespace std;

void getCOL_z::eval(const double & av_z, const double & Q2, const int & charge_in, const int & model_in, const int & evo_in, const std::vector<double> & FF_in, const std::vector<double> &par){

    int model = model_in, evo = evo_in, charge = charge_in;    
    N_fav = 0, N_dis = 0;
     
    //get the fragmentation function
    std::vector<double> FF(13);
    for(int i = 0; i < FF.size(); i++) FF[i] = FF_in[i];
    
    if(av_z != 1.){
    
        if(model == 0){ //fit 2015 [0] same beta 
            
            double NC_fav = par[4], NC_dis = par[5], gamma = par[6], delta = par[7]; 
            
//             N_fav = NC_fav * pow(av_z, gamma) * pow(1. - av_z, delta);// * (pow(gamma + delta, gamma + delta) / (pow(gamma , gamma) * pow(delta, delta)));
            N_fav = NC_fav * pow(av_z, gamma) * pow(1. - av_z, delta) * (pow(gamma + delta, gamma + delta) / (pow(gamma , gamma) * pow(delta, delta)));
    //         N_dis=NC_dis*pow(av_z,gamma)*pow(1.-av_z,delta)*(pow(gamma+delta,gamma+delta)/(pow(gamma,gamma)*pow(delta,delta)));
            N_dis = NC_dis;
            
            if(isinf(pow(1. - av_z, delta))){
                cout<<"in COL_z: N_fav = "<<N_fav<<"\tz = "<<av_z<<"\tdelta = "<<delta<<"\t(1-z)^delta = "<<pow(1. - av_z, delta)<<endl;
            }   
                
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = 2.0 * N_dis * FF[i];//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    
                    Col_fav = 2.0 * N_fav * FF[8];
                    Col_dis = 2.0 * N_dis * FF[7];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for d and ubar
                    
                }
                if(charge == -1){ //pi-
                    
                    Col_fav = 2.0 * N_fav * FF[7];
                    Col_dis = 2.0 * N_dis * FF[8];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
        }
        
        if(model == 100){ //fit 2015 [100] beta_u != beta_d
            
            double NC_fav = par[5], NC_dis = par[6], gamma = par[7], delta = par[8]; 
            
            N_fav = NC_fav * pow(av_z, gamma) * pow(1. - av_z, delta) * (pow(gamma + delta, gamma + delta) / (pow(gamma, gamma) * pow(delta, delta)));
    //         N_dis = NC_dis*pow(av_z,gamma)*pow(1.-av_z,delta)*(pow(gamma+delta,gamma+delta)/(pow(gamma,gamma)*pow(delta,delta)));
            N_dis = NC_dis;

            
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = 2.0 * N_dis * FF[i];//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    
                    Col_fav = 2.0 * N_fav * FF[8];
                    Col_dis = 2.0 * N_dis * FF[7];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for d and ubar
                }                                                 
                                                                
                if(charge == -1){ //pi-                           
                                                                
                    Col_fav = 2.0 * N_fav * FF[7];                    
                    Col_dis = 2.0 * N_dis * FF[8];                    
                                                                
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
        }
        
        
        if(model == 1){//fit 2015 - Bernstein polynomials for Collins
            
            double a_fav = par[4], a_dis = par[5], b_fav = par[6], b_dis = par[7]; 
            
    //         cout<<"model="<<model<<"\tpars:\t"<<a_fav<<"\t"<<a_dis<<"\t"<<b_fav<<"\t"<<b_dis<<endl;
            
            N_fav = a_fav * (1.0 - av_z) + b_fav * av_z;
            N_dis = a_dis * (1.0 - av_z) + b_dis * av_z;
            
    //         cout<<"a_fav="<<a_fav<<"\ta_dis="<<a_dis<<"\tb_fav="<<b_fav<<"\tb_dis="<<b_dis<<endl;
    //         cout<<"N_fav="<<N_fav<<"\tN_dis="<<N_dis<<"\tav_z="<<av_z<<endl;
                    
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = 2.0 * N_dis * FF[i];//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    
                    Col_fav = 2.0 * N_fav * FF[8];
                    Col_dis = 2.0 * N_dis * FF[7];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for d and ubar
                }
                if(charge == -1){ //pi-
                    
                    Col_fav = 2.0 * N_fav * FF[7];
                    Col_dis = 2.0 * N_dis * FF[8];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
        }
        
        if(model == 10){
            
            double NC_fav = par[6], NC_dis = par[7], ga_fav = par[8], ga_dis = par[9], de_fav = par[10], de_dis = par[11]; 
            
    //         cout<<"model="<<model<<"\tpars:\t"<<a_fav<<"\t"<<a_dis<<"\t"<<b_fav<<"\t"<<b_dis<<endl;
            
            N_fav = NC_fav * pow(av_z, ga_fav) * pow(1. - av_z, de_fav) * (pow(ga_fav + de_fav, ga_fav + de_fav) / (pow(ga_fav, ga_fav) * pow(de_fav, de_fav)));
            N_dis = NC_dis * pow(av_z, ga_dis) * pow(1. - av_z, de_dis) * (pow(ga_dis + de_dis, ga_dis + de_dis) / (pow(ga_dis, ga_dis) * pow(de_dis, de_dis)));
    //         N_dis=NC_dis;
            
    //         cout<<"NC_fav="<<NC_fav<<"\tNC_dis="<<NC_dis<<"\tga_fav="<<ga_fav<<"\tga_dis="<<ga_dis<<"\tde_fav="<<de_fav<<"\tde_dis="<<de_dis<<endl;
    //         cout<<"N_fav="<<N_fav<<"\tN_dis="<<N_dis<<"\tav_z="<<av_z<<endl;
                    
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = 2.0 * N_dis * FF[i];//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    
                    Col_fav = 2.0 * N_fav * FF[8];
                    Col_dis = 2.0 * N_dis * FF[7];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for d and ubar
                }
                
                if(charge == -1){ //pi-
                    
                    Col_fav = 2.0 * N_fav * FF[7];
                    Col_dis = 2.0 * N_dis * FF[8];
                    
                    if(i == 5 || i == 8) COL_z[i] = 2.0 * N_dis * FF[i]; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = 2.0 * N_fav * FF[i]; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
        }
        
        if(model == 2){ //fit 2018 - REWEIGHTING
            
            double NC_fav = par[6], NC_dis = par[7], gamma = par[8], delta = par[9]; 
            
            N_fav = NC_fav * pow(av_z, gamma) * pow(1.-av_z, delta);
    //         N_dis=NC_dis*pow(av_z,gamma)*pow(1.-av_z,delta);
            N_dis = NC_dis;
                    
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = N_dis;//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    if(i == 5 || i == 8) COL_z[i] = N_fav; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = N_dis; //disfavoured for d and ubar
                }
                
                if(charge == -1){ //pi-
                    if(i == 5 || i == 8) COL_z[i] = N_dis; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = N_fav; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
            Col_fav = N_fav;
            Col_dis = N_dis;
        }
        

        if(model == 3){ //fit 2018 - REWEIGHTING
            
            double NC_fav = par[6], NC_dis = par[7], ga_fav = par[8], ga_dis = par[9], de_fav = par[10], de_dis = par[11]; 
            
            N_fav = NC_fav * pow(av_z, ga_fav) * pow(1.-av_z, de_fav);
            N_dis = NC_dis * pow(av_z, ga_dis) * pow(1.-av_z, de_dis);
    //         N_dis=NC_dis;
                    
            COL_z[0] = 0;
            COL_z[1] = 0;
            COL_z[2] = 0;
            for(int i = 3; i <= 9; i++){
                
                if(i == 3 || i == 9) COL_z[i] = N_dis;//s or sbar are always disfavoured for pions 
                
                if(charge == 1){ //pi+
                    
                    if(i == 5 || i == 8) COL_z[i] = N_fav; //favoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = N_dis; //disfavoured for d and ubar
                }
                
                if(charge == -1){ //pi-
                    
                    if(i == 5 || i == 8) COL_z[i] = N_dis; //disfavoured for u and dbar
                    if(i == 4 || i == 7) COL_z[i] = N_fav; //favoured for d and ubar
                }
            }
            
            COL_z[6] = 0; //no gluon
            COL_z[10] = 0;
            COL_z[11] = 0;
            COL_z[12] = 0;  
            
            Col_fav = N_fav;
            Col_dis = N_dis;
        }
    }
    
    else if(av_z == 1.){
//         cout<<"z = "<<av_z<<"\n Setting Collins to zero\n";
        for(int i = 0; i < COL_z.size(); i++) COL_z[i] = 0.;
    }
    
    for(int i = 0; i < COL_z.size(); i++){
        
        if(COL_z[i] > max_COL_z[i]) max_COL_z[i] = COL_z[i];
        if(COL_z[i] < min_COL_z[i]) min_COL_z[i] = COL_z[i];
    }
}

// void getCOL_z::getPolyNorm(const double & A, const double & B){
//         
//     double a=A, b=B;
//     double xp=0.0, xm=0.0, Fxp=0.0, Fxm=0.0;
//     
//     if(pow(a,2)-3.*b*(1.-a-b)<0){
//         PolyNorm=1;
//         cout<<"Warning: imaginary root\n";
//     }
//     if(pow(a,2)-3.*b*(1.-a-b)>=0){
//         xp=(1./(3.*b))*(-a+sqrt(pow(a,2)-3.*b*(1.-a-b)));
//         xm=(1./(3.*b))*(-a-sqrt(pow(a,2)-3.*b*(1.-a-b)));
//         
//         Fxp=xp*(1.-a-b +a*xp + b*pow(xp,2));
//         Fxm=xm*(1.-a-b +a*xm + b*pow(xm,2));
//     
//         PolyNorm=maximum(1,abs(Fxp),abs(Fxm));
//     }
// }

double getCOL_z::maximum(const double & a, const double & b){
    
    return a < b ? b : a;
}

double getCOL_z::maximum(const double & a, const double & b, const double & c){
    
    return maximum(maximum(a,b),c);
    
}
    
void getCOL_z::plot(const std::string outname_in, const double & Q2, const int & charge_in, const int & model_in, const int & evo_in, const std::vector<double> &par){
    
    FRAG::FF myFF;
    
    std::vector<double> FFpp(13), FFpm(13);
    
    double f[13];
    
    cout<<"printing Collins on "<<outname_in<<endl;
    
    vector<double> z_vals = {.05, .055, .06, .065, .07, .075, .08, .085, .09, .095, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, 1.0};
    
    ofstream col(outname_in);
    if(!col.is_open()) cerr<<"Error in opening "<<outname_in<<endl;

    
//     col<<"#"<<FFset_in<<", iset="<<iset_ff_in<<", order="<<order_in<<endl;
    col<<"#Q^2="<<Q2<<" GeV^2"<<endl;
    col<<"#z\tfav(z)\tunf(z)\n";
    
    for(int i = 0; i < z_vals.size(); i++){
        col<<z_vals[i]<<"\t";
        
        if(evo_in == 0){
            
            // myFF.FF_eval("DSS", iset_ff, 1, 1, 0, z_vals[i], Q2);
            myFF.FF_eval(1, 1, z_vals[i], Q2);
            for(int j = 0; j < FFpp.size(); j++) FFpp[j] = myFF.theFF[j] / z_vals[i]; //to scale FFs: routines give z*D(z)
            eval(z_vals[i], Q2, 1,  model_in, evo_in, FFpp, par);
            col<<Col_fav<<"\t";
            
            // myFF.FF_eval("DSS", iset_ff, 1, -1, 0, z_vals[i], Q2);
            myFF.FF_eval(1, -1, z_vals[i], Q2);
            for(int j = 0; j < FFpm.size(); j++) FFpm[j] = myFF.theFF[j] / z_vals[i]; //to scale FFs: routines give z*D(z)
            eval(z_vals[i], Q2, -1,  model_in, evo_in, FFpm, par);
            col<<Col_dis<<endl;
        }
        
        if(evo_in == 1){
            
            hoppetEvalcf(z_vals[i], sqrt(Q2), f);
            col<<f[8] / z_vals[i]<<"\t"<<f[7] / z_vals[i]<<"\t"<<f[9] / z_vals[i]<<"\t";
            hoppetEvalcff(z_vals[i], sqrt(Q2), f);
//             hoppetEvalcff(z_vals[i], sqrt(Q2), f);
            col<<f[8] / z_vals[i]<<"\t"<<f[7] / z_vals[i]<<"\t"<<f[9] / z_vals[i]<<"\n";
        }
    }

    col.close();
    
}    
    
    
}
