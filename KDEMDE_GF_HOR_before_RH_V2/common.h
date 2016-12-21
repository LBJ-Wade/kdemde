//
//  common.h
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/13/16.
//  Copyright (c) 2016 c. All rights reserved.
//

#ifndef KinDecEMDE_common_h
#define KinDecEMDE_common_h


#endif


//==================================================================================================
// Standards
//==================================================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

#include <sys/types.h>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>// GSL_SUCCESS ...
#include <gsl/gsl_math.h> // cos(x), sin(x) ...
#include <gsl/gsl_odeiv2.h>// ODE solver
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>
#include <exception>

#define Lmax 15
#define Nmax 1100

//set these to Lmax=1, Nmax=0 when decoupled

#define _SUCCESS_ 0
#define _FAILURE_ 1

#define GAMMA gsl_sf_gamma_inc

using namespace std;


//==================================================================================================
// Flag Switches
//==================================================================================================

const bool numjac=0;//1 For numerical evaluation of Jacobian.
const bool decoupled=0;//1;// if decoupled 1 then gammakds=0. When working on the Dcpl limit also set Lmax=1 and Nmax=0.


//==================================================================================================
// DIMENSION of SYSTEM for DM perturbations ONLY
//==================================================================================================


const extern long DIM;

const int Shift=0;//If its 9 it corresponds to placing DM perturbations at the END of the column rather than at the start.

//==================================================================================================
// Constants. Units GeV
//==================================================================================================

const double InvGeVTocm=0.1973300000e-13;
const double pcTocm=3.068e18;
const double cmTopc=1.0/pcTocm;
const double InvGeVTopc=InvGeVTocm*cmTopc;
const double GeVToInvpc=1.0/InvGeVTopc;
const double InvPcToGev=1.0/GeVToInvpc;

const double mPL=1.22e19;
const double mPL2=1.4884e38;

//=== Cosmological Parameters Planck 15 =======/

const float  OmegaM=0.308;
const float  h=0.678;
const double  rho_C=8.0992e-47*pow(h, 2);
const double T_today= 2.348e-13;//GeV

//=== Model and derrived Params ===============/

const float alpha=3.0/8.0;//See V&G paper
const int nu=4.0;//See V&G paper
const int n=2;// p wave scattering
const int beta=4+2-nu;
const double lambda=(2.0-alpha)/(alpha*beta);
const float fdm=0.0;

const double a_HOR=100;//100;//*pow(20.0/6.0,8.0/3.0/n);//pow(20.0/6.0,8.0/3.0/n);//100;//*pow(20/2.81165,8.0/3.0/n);
const float a_I=1.0;
const float xini=log10(a_I);
const double astar =1.0e-10;// arbitrary, chosen to represent the a when T_L = T_DM in the fully decoupled limit. Should be zero.

const double kTilde=pow(1.0*a_HOR/a_I, -0.5)*a_I;


struct Params
{
    double a_fin;
    double c; //c is defined as sqrt(2T_kds/m_DM)
    double kTokRH1;//Related to GammaTilde
    double T_REH;
    double T_kdS;
    double m_DM;
    double a_REH;
    double a_HOR;
    //add latter derrived parameters and make a function that INITIALIZES THE PARAMETERS
    
};

typedef struct Params* pParams;



struct Params_Full
{
    double a_fin;
    double c; //c is defined as sqrt(2T_kds/m_DM)
    double kTokRH1;//Related to GammaTilde
    double T_REH;
    double T_kdS;
    double a_RH;
    //add latter derrived parameters and make a function that INITIALIZES THE PARAMETERS
    
};


struct Background
{
    double rhophi;
    double rhom;
    double rhor;
    double TDM;
    double TL;
};

typedef struct Background* pBackground;

//==================================================================================================
// index map
//==================================================================================================
struct var_nl
{
    int n;
    int l;
};


struct ExtraIndices
{
    int BkgP;
    int BkgR;
    int BkgM;
    int BkgT;
    int deltaP;
    int thetaP;
    int deltaR;
    int thetaR;
    int Phi;
};

typedef struct ExtraIndices* pExtraIndices;




//extern ExtraIndices ExtInd;



extern int column[Nmax+1][Lmax+1];// This will store a precomputed map from (n,l) to Column index. 


extern vector<var_nl> indices; // This will store a precomputed map of Column to (n,l) via indices[col].n and indices[col].l


extern Params *glob_Par;//=NULL;

extern ExtraIndices *pExtInd;//Those are the 9 extra Column indices for the Background (4) and for the Radiation (2), Scalar (2), and Grav Potential (1) perturbations.

struct Output_Files
{
    ofstream FileBkg;
    ofstream FilePert;
    ofstream FilePertDM;
    ofstream FileFullPertDM;
};

typedef struct Output_Files* pOutFiles;

struct Growth_Function_File
{
    ofstream GFFile;
};

typedef struct Growth_Function_File* pGFFile;

/*** Debug variables below *********/

extern int test_func;
extern int test_jac;
