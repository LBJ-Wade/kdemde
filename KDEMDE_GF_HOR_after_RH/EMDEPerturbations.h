//
//  EMDEPerturbations.h
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/18/16.
//  Copyright (c) 2016 c. All rights reserved.
//

#ifndef __KinDecEMDE__EMDEPerturbations__
#define __KinDecEMDE__EMDEPerturbations__

#include <stdio.h>

#endif /* defined(__KinDecEMDE__EMDEPerturbations__) */

//==================================================================================================
// ode solver & routines
//==================================================================================================
#include "ODE_solver_Rec.h"
#include "Definitions.h"
#include "routines.h"


using namespace std;
using namespace ODE_solver_Rec;

void test_extern_params();

void fcn_k_eval(int *neq, double *t, double *y, double *f, int col);
void jac_k_eval(int *neq, double *t, double *y, double *f, int col);
void fcn_k(int *neq, double *t, double *y, double *f, int col);

void output_current_pert_k(ofstream &ofile, const ODE_solver_Solution &Sz);


void output_current_BkgSol(ofstream &ofile, const ODE_solver_Solution &Sz);


void output_current_DMpert_k(ofstream &ofile, const ODE_solver_Solution &Sz);



double setup_kevol(struct Background *bkgdens, void * parameters, struct Output_Files *pFiles, struct Growth_Function_File *pGFFile);
