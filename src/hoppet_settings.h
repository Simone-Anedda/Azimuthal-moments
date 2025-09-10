#ifndef HOPPETSETTINGS_H
#define HOPPETSETTINGS_H

//HOPPET parameters
#define  ymax        log(1.0/0.001)     //< highest value of ln1/x user wants to access  
#define  ymaxT       log(1.0/0.001)     //< highest value of ln1/x user wants to access  
#define  ymaxCol     log(1.0/0.05) 
#define  dy                0.2          //< internal ln1/x grid spacing: 0.1-0.25 is a sensible range
#define  Qmin              1.0       //< lower limit of Q range
#define  Qmax              sqrt(40.0)         //< upper limit of Q range
#define  dlnlnQ            dy/16.0       //< internal table spacing in lnlnQ (e.g. dy/4)
#define  nloop             1            //< the maximum number of loops we'll want (<=3)
#define  intorder         -5            //< order of numerical interpolation (e.g. -6)
#define  factscheme        5        //Transversity kernel

#define  QminT             0.9/*1.0*/     //< lower limit of Q range
#define  QmaxT             10000.0      //< upper limit of Q range
#define  QminC             1.0     //< lower limit of Q range
#define  QmaxC             1000.0     //< upper limit of Q range

#define fixed_nf 4 //cteq66
#define asQ0 0.117982 //cteq66 alhpa_s(M_Z)
#define Q0alphas 91.1876 //cteq66 M_Z
#define Q0pdf 1.3 //cteq66

//GRV98 scheme -- see 'notes-unpol-xsecs' 
//#define fixed_nf 3
//#define asQ0 0.44
//#define Q0alphas 1.0
//#define Q0pdf 1.0


//#define fixed_nf 4
//#define asQ0 0.12979
//#define Q0alphas 91.2
//#define Q0pdf sqrt(1.2)
    
#define muR_Q 1.0
    
//hoppet default (pole masses)
#define mcharm 1.414213563
#define mbottom 4.5
#define mtop 175.0
//MSbar masses PDG    
#define mcPDG 1.275
#define mbPDG 4.18
#define mtPDG 160.0

#endif 
