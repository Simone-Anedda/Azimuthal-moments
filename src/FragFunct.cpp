//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#include "FragFunct.h"

namespace FRAG{

using namespace std;

void FF::use_LHAPDF(bool use_LHAPDF_in){
    
    useLHAPDF = use_LHAPDF_in;
}
    

void FF::set_FFset(const std::string FFset_in){
    
    FFset = FFset_in;
}

void FF::set_FForder(const std::string FForder_in){
    
    FForder = FForder_in;
}

void FF::set_FF(const std::string FFset_in, const std::string FForder_in){
    
    set_FFset(FFset_in);
    set_FForder(FForder_in);
    
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
}

void FF::set_LHAPDF_FFset(map<string, string> FFsets){
    
    useLHAPDF = true;

    ff_pip = LHAPDF::mkPDF(FFsets["pip"] + "/0");
    ff_pim = LHAPDF::mkPDF(FFsets["pim"] + "/0");
    ff_kap = LHAPDF::mkPDF(FFsets["kap"] + "/0");
    ff_kam = LHAPDF::mkPDF(FFsets["kam"] + "/0");
    ff_pisum = LHAPDF::mkPDF(FFsets["pisum"] + "/0");
    if(!FFsets["kasum"].empty()) ff_kasum = LHAPDF::mkPDF(FFsets["kasum"] + "/0");
    if(!FFsets["prp"].empty()) ff_prp = LHAPDF::mkPDF(FFsets["prp"] + "/0");
    if(!FFsets["prm"].empty()) ff_prm = LHAPDF::mkPDF(FFsets["prm"] + "/0");
    if(!FFsets["prsum"].empty()) ff_prsum = LHAPDF::mkPDF(FFsets["prsum"] + "/0");
}



void FF::FF_eval(const int & hadron_in, const int & charge_in, const double & av_z_in, const double & Q2_in){
    
    int hadron = hadron_in, charge = charge_in;
    
    int khffset, khffcharge; //for Kretzer
    double uff[2], dff[2], sff[2], cff[2], bff[2], gff; 

    // cout << "calling FF_eval\n";
    
    int pion = 1, kaon = 2;
    double up, ubp, dp, dbp, sp, sbp, cp, bp, glp;  
    double uk, ubk, dk, dbk, sk, sbk, ck, bk, glk; 
    
    double u, ub, d, db, s, sb, c, b, gl;
    
    double av_z = av_z_in, av_Q2 = Q2_in;
    
    if(!useLHAPDF){
        
        if(hadron == 1) khffset = 1;
        if(hadron == 2) khffset = 3;
        if(hadron == 4) khffset = 5;
        if(charge == 0) khffcharge = 3;
        if(charge == -1) khffcharge = 2;
        
        if(FFset == "DSS"){
        
            if(hadron != 4){
                if(iset_ff == 1) fdss_(&hadron, &charge, &FForder_i, &av_z, &av_Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
                if(iset_ff == 2) fdss2017nlo_(&iset_ff, &hadron, &charge, &FForder_i, &av_z, &av_Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
                theFF[0] = 0; 
                theFF[1] = b;
                theFF[2] = c;
                theFF[3] = sb;
                theFF[4] = ub;
                theFF[5] = db;
                theFF[6] = gl;
                theFF[7] = d;
                theFF[8] = u;
                theFF[9] = s;
                theFF[10] = c;
                theFF[11] = b;
                theFF[12] = 0;
                cout<<iset_ff<<"\t"<<hadron<<": d="<<theFF[7]<<"\tu="<<theFF[8]<<endl;
                
            }
        
            if(hadron == 4){
                if(iset_ff == 1){
                    fdss_(&hadron, &charge, &FForder_i, &av_z, &av_Q2, &u, &ub, &d, &db, &s, &sb, &c, &b, &gl);
                    theFF[0] = 0; 
                    theFF[1] = b;
                    theFF[2] = c;
                    theFF[3] = sb;
                    theFF[4] = ub;
                    theFF[5] = db;
                    theFF[6] = gl;
                    theFF[7] = d;
                    theFF[8] = u;
                    theFF[9] = s;
                    theFF[10] = c;
                    theFF[11] = b;
                    theFF[12] = 0;
                    // cout<<hadron<<": d="<<theFF[7]<<"\tu="<<theFF[8]<<endl;
                }
            
                if(iset_ff == 2){ 
                    fdss2017nlo_(&iset_ff, &pion, &charge, &FForder_i, &av_z, &av_Q2, &up, &ubp, &dp, &dbp, &sp, &sbp, &cp, &bp, &glp);
                    fdss2017nlo_(&iset_ff, &kaon, &charge, &FForder_i, &av_z, &av_Q2, &uk, &ubk, &dk, &dbk, &sk, &sbk, &ck, &bk, &glk);
                    theFF[0] = 0; 
                    theFF[1] = bp + bk;
                    theFF[2] = cp + ck;
                    theFF[3] = sbp + sbk;
                    theFF[4] = ubp + ubk;
                    theFF[5] = dbp + dbk;
                    theFF[6] = glp + glk;
                    theFF[7] = dp + dk;
                    theFF[8] = up + uk;
                    theFF[9] = sp + sk;
                    theFF[10] = cp + ck;
                    theFF[11] = bp + bk;
                    theFF[12] = 0;    
                    // cout<<"dp = "<<dp<<"    up="<<up<<"   dk = "<<dk<<"    uk="<<uk<<endl;
                    // cout<<hadron<<": d="<<theFF[7]/dp<<"\tu="<<theFF[8]/up<<endl;
                }
            }
            
        }
        
        if(FFset == "Kretzer"){
        
            pkhff_(&hadron, &khffset, &khffcharge, &av_z, &av_Q2, uff, dff, sff, cff, bff, &gff);
            theFF[0] = 0; 
            theFF[1] = bff[1];
            theFF[2] = cff[1];
            theFF[3] = sff[1];
            theFF[4] = uff[1];
            theFF[5] = dff[1];
            theFF[6] = gff;
            theFF[7] = dff[0];
            theFF[8] = uff[0];
            theFF[9] = sff[0];
            theFF[10] = cff[0];
            theFF[11] = bff[0];
            theFF[12] = 0;
            
            if(khffcharge == 3) for(int i = 0; i < theFF.size(); i++) theFF[i] /= 2.0;
        }
    }
    
    else if (useLHAPDF){

        if(charge == 1){  

            if(hadron == 1) ff_pip->xfxQ2(av_z, av_Q2, theFF);
            if(hadron == 2) ff_kap->xfxQ2(av_z, av_Q2, theFF);            
            if(hadron == 4){

                ff_pip->xfxQ2(av_z, av_Q2, theFF);
                ff_kap->xfxQ2(av_z, av_Q2, tmpFF);
                if(FFset == "NNFF10"){
                    
                     ff_prp->xfxQ2(av_z, av_Q2, tmpFF2);
                     for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i] + tmpFF2[i];
                }
                else if(FFset == "DEHSS" || FFset == "MAPFF10") for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i];

                }
            }
 

        if(charge == -1){

            if(hadron == 1) ff_pim->xfxQ2(av_z, av_Q2, theFF);
            if(hadron == 2) ff_kam->xfxQ2(av_z, av_Q2, theFF);
            if(hadron == 4){

                ff_pim->xfxQ2(av_z, av_Q2, theFF);
                ff_kam->xfxQ2(av_z, av_Q2, tmpFF);
                
                if(FFset == "NNFF10"){
                    
                     ff_prm->xfxQ2(av_z, av_Q2, tmpFF2);
                     for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i] + tmpFF2[i];
                }
                else if(FFset == "DEHSS" || FFset == "MAPFF10") for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i];
                // cout<<hadron<<": d="<<tmpFF[7]/theFF[7]<<"\tu="<<tmpFF[8]/theFF[8]<<endl;

            }

        }

        if(charge == 0){
            
            if(hadron == 1) ff_pisum->xfxQ2(av_z, av_Q2, theFF);

            if(hadron == 2){

                if(FFset == "DEHSS"){

                    ff_kap->xfxQ2(av_z, av_Q2, theFF);
                    ff_kam->xfxQ2(av_z, av_Q2, tmpFF);
                    for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i];
                }
                else ff_kasum->xfxQ2(av_z, av_Q2, theFF);
            
            } 
            
            if(hadron == 4){

                ff_pisum->xfxQ2(av_z, av_Q2, theFF);
                if(FFset != "DEHSS"){
                    
                    ff_kasum->xfxQ2(av_z, av_Q2, tmpFF);
                    if(FFset == "NNFF10") ff_prsum->xfxQ2(av_z, av_Q2, tmpFF2);
                }
                else if(FFset == "DEHSS"){
                    
                    ff_kap->xfxQ2(av_z, av_Q2, tmpFF);
                    ff_kam->xfxQ2(av_z, av_Q2, tmpFF2);
                }

                if(FFset != "MAPFF10") for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i] + tmpFF2[i];
                if(FFset == "MAPFF10") for(int i = 0; i < theFF.size(); i++) theFF[i] += tmpFF[i];
                
            }

            for(int i = 0; i < theFF.size(); i++) theFF[i] / 2.;

        }
    
    }
    
    
}


void FF::FF_plot(const std::string frag_name, const int & hadron_in, const int & charge_in, const double & Q2_in){
    
    vector<double> z_vals = {.01, .015, .02, .025, .03, .035, .04, .045, .05, .055, .06, .065, .07, .075, .08, .085, .09, .095, .1, .11, .12, .13, .14, .15, .16, .17, .18, .19, .2, .21, .22, .23, .24, .25, .26, .27, .28, .29, .3, .31, .32, .33, .34, .35, .36, .37, .38, .39, .4, .41, .42, .43, .44, .45, .46, .47, .48, .49, .5, .51, .52, .53, .54, .55, .56, .57, .58, .59, .6, .61, .62, .63, .64, .65, .66, .67, .68, .69, .7, .71, .72, .73, .74, .75, .76, .77, .78, .79, .8, .81, .82, .83, .84, .85, .86, .87, .88, .89, .9, .91, .92, .93, .94, .95, .96, .97, .98, .99, 1.0};
    
    ofstream frag(frag_name);
    
    frag<<"#"<<FFset<<", iset="<<iset_ff<<", order="<<FForder_i<<endl;
//     frag<<"#"<<FFset_in<<", iset="<<iset_ff_in<<", order="<<order_in<<endl;
    frag<<"#hadron="<<hadron_in<<", charge="<<charge_in<<", Q^2="<<Q2_in<<" GeV^2"<<endl;
    frag<<"#z\tz*u(z)\tz*d(z)\tz*ub(z)\tz*db(z)\tz*s(z)\tz*g(z)\n";
    
    for(int i = 0; i < z_vals.size(); i++){
        
        FF_eval(hadron_in, charge_in, z_vals[i], Q2_in);
//         FF_eval(FFset_in, iset_ff_in, hadron_in, charge_in, order_in, z_vals[i], Q2_in);
        frag<<z_vals[i]<<"\t"<<theFF[8]<<"\t"<<theFF[7]<<"\t"<<theFF[4]<<"\t"<<theFF[5]<<"\t"<<theFF[9]<<"\t"<<theFF[6]<<endl;
    }

    frag.close();
    
}

}
