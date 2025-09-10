//
// Author: Carlo Flore <carlo.flore@ijclab.in2p3.fr>
//

#pragma once

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>



using namespace std;



#ifdef __cplusplus
extern "C" {
  #endif
    
  void pkhff_(int *ih, int *iset, int *icharge, double *z, double *Q2, 
              double uff[2], double dff[2], double sff[2], 
              double cff[2], double bff[2], double *gff);
  ///old dss
  void fdss_(int *hadron, int *charge, int *order, double *z, double *Q2,
	     double *u, double *ub, double *d, double *db, double *s, double *sb,
	     double *c, double *b, double *gl );
  
  ///new dss
  void fdssh_(int *is, int *il, int *charge, int *order, double *z, double *Q2,
	     double *u, double *ub, double *d, double *db, double *s, double *sb,
	     double *c, double *b, double *gl );
  ///new dss
  void fdss17_(int *is, int *hadron, int *charge, int *order, double *z, double *Q2,
	     double *u, double *ub, double *d, double *db, double *s, double *sb,
	     double *c, double *b, double *gl );
  ///new dss
  void fdss2017lo_(int *is, int *hadron, int *charge, int *order, double *z, double *Q2,
	     double *u, double *ub, double *d, double *db, double *s, double *sb,
	     double *c, double *b, double *gl );
  ///new dss
  void fdss2017nlo_(int *is, int *hadron, int *charge, int *order, double *z, double *Q2,
	     double *u, double *ub, double *d, double *db, double *s, double *sb,
	     double *c, double *b, double *gl );
  
  // Access to FORTRAN FRAGINI COMMON block (see file fDSS.f or pkhff.f).
  extern struct {
//     int ifini;
    int ifinip;
    int ifinik;
    int ifinih;
  } fraginipkhff_;
 
  extern struct {
    int finip; 
    int finik;
    int finipr;
    int finih;
  } fragini_;

  extern struct {
    int finip; 
    int finik;
    int finipr;
    int finih;
  } fraginilo2017_;
  
  extern struct {
    int finip; 
    int finik;
    int finipr;
    int finih;
  } fragininlo2017_;
 
  
  #ifdef __cplusplus
}
#endif


/*
 * *******************************************************************
 *                                                                  *
 *        fdss  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
 *  D.de Florian, R.Sassot, M.Stratmann   hep-ph/0703242)           *
 *    Phys.Rev.D.75:114010,2007
 *                                                                  *
 *     CALL fdss (IH,IC,IO, X, Q2, U, UB, D, DB, S, SB, C, B, GL)   *
 *                                                                  *	
 *  INPUT:                                                          *
 *  IH = hadron type    1: PION                                     *
 *                      2: KAON                                     *
 *                      3: PROTON                                   *
 *                      4: CHARGED HADRONS                          *
 *                                                                  *
 *  IC = Hadron Charge  0: 0 (as average of + and -)                *
 *                      1: +                                        *
 *                     -1: -                                        *
 *                                                                  *
 *  IO= Order           0: LO                                       *
 *                      1: NLO                                      *
 *                                                                  *
 *            X                    (between  0.05   and  1.0)       *
 *            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
 *             (for values outside the allowed range the program    *
 *              writes a warning and extrapolates to the x and      *
 *              Q2 values requested)                                *
 *                                                                  *
 *   OUTPUT: U, UB, D, DB, S, SB,   C,           B,       GL        *
 *           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon      *
 *           Always X times the distribution is returned            *
 *                                                                  *
 *                                                                  *
 *   COMMON:  The main program or the calling routine has to have   *
 *            a common block  COMMON / FRAGINI / FINI , and  FINI   *
 *            has always to be zero when DSS is called for the      *
 *            first time or when the SET has been changed.          *
 *                                                                  *
 ********************************************************************
 */