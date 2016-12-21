//
//  main.cpp
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/13/16.
//  Copyright (c) 2016 c. All rights reserved.
//

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


//==================================================================================================
// ode solver & routines
//==================================================================================================
#include "ODE_solver_Rec.h"
#include "Definitions.h"
#include "routines.h"


//==================================================================================================
// Functions for computing the terms in Eq.24 in Bertchinger 06 adapted to EMDE coded in functions.cpp.
// EMDEPerturbations.cpp is where the ODEs are coded for both Background and Perturbation Evoution. 
//==================================================================================================

#include "common.h"
#include "functions.h"
#include "EMDEPerturbations.h"



const long DIM = (Nmax+1) * (Lmax +1); //Nmax and Lmax are defined in common.h Defines size of system for DM pert only.


using namespace std;
using namespace ODE_solver_Rec;
//==================================================================================================
// Extern Parameters and structures
//==================================================================================================


extern vector<var_nl> indices;
extern Params *glob_Par;//=NULL;



//==================================================================================================
// Main Code
//==================================================================================================

int main(int argc, const char * argv[]) {
    
    //==========Parameters===================
    double a_fin;
    double c=atof(argv[2]);//0.5;//0.24;//1.5e-1;//0.01; // c is defined as sqrt(2T_d/m_DM)
    double kTokRH1;
    double kTokRH1_min=atof(argv[1]), kTokRH1_max=5.1;//sqrt(10);//sqrt(10)/100;//sqrt(10)/100;//100;
    int nk_grid=1;
    double T_REH=5;
    double T_kdS=10;//19.94079648;//10;//5.460626859;// corresponds to a_KD=10 19.94079648; // a_KD=100 if TRH=5 for  T_kdS~20
    double m_DM = 1.e9;//This parameter is only "relevant" in the decoupled limit. Needs to be chosen large enough to supress the kcT term in Eq.18, at the largest value of a considered.
  
    vector<double> kvals(nk_grid);
    
    
    //==============================================================
    // initialize kappa grid
    //==============================================================
    init_xarr(kTokRH1_min, kTokRH1_max, &kvals[0], nk_grid, 0, 1); // 5th argument ZERO FOR LINEAR grid

    time_t  t0, t1; /* time_t is defined on <time.h> and <sys/types.h> as long */
    clock_t c0, c1; /* clock_t is defined on <time.h> and <sys/types.h> as int */
    
    
    t0 = time(NULL);
    c0 = clock();
    
    printf ("\tbegin (wall):            %ld\n", (long) t0);
    printf ("\tbegin (CPU):             %d\n\n", (int) c0);


    
    ExtraIndices EI={DIM,DIM+1,DIM+2,DIM+3,DIM+4,DIM+5,DIM+6,DIM+7,DIM+8};//This amounts to placing the DM pert at the begining of the unknown functions column.
    pExtInd=&EI;
    
    create_index_map();//computes and stores an index map from column to n,l. The way it is written DO NOT call more than once, as it appends copies of itself. 

    Params P ={a_fin, c, kTokRH1, T_REH, T_kdS,m_DM};
    glob_Par=&P;
    
    Output_Files Files;
    Growth_Function_File GF_File;
    
    create_output_file_growth_function(&GF_File);
    
    for(int i=0; i < kvals.size(); i++)
    {
        kTokRH1= kvals[i];
        a_fin=100*a_HOR*pow(kTokRH1, 2); // set a_fin to be 100*a_RH for all k modes!
//        Params P ={a_fin, c, kTokRH1, T_REH, T_kdS,m_DM};
//        glob_Par=&P;
        
        P.a_fin=a_fin;
        P.kTokRH1=kTokRH1;

        create_output_files(&Files);
        cout<<"kTokRH1 From Glob="<<glob_Par->kTokRH1<<endl;
        Background InitialBkgDens;
        ICBkgRho_TDM(&InitialBkgDens, &P);//Initializes the Bkg Energy Densities and Temp.
    
        cout<<"GammaT="<<GammaTilde(&P)<<" UpsiIni="<<UpsiPar(&InitialBkgDens, &P)<<" kTilde="<<kTilde<<endl;
        printf("a_KD=%f a_HOR= %f and a_REH=%f\n", a_KD(&InitialBkgDens, &P), a_HOR, a_REH(&P) );
    
        setup_kevol(&InitialBkgDens, &P, &Files, &GF_File);
        close_output_files(&Files);

        
    }
    
//    close_output_files(&Files);
    close_GF_file(&GF_File);
    
    t1 = time(NULL);
    c1 = clock();
    
    printf ("\xfin (wall):              %ld\n", (long) t1);
    printf ("\xfin (CPU);               %d\n", (int) c1);
    printf ("\telapsed wall clock time: %ld\n", (long) (t1 - t0));
    printf ("\telapsed CPU time:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
    
    
    
    //=========== TEST AREA =======================
    
    
    return 0;
}
