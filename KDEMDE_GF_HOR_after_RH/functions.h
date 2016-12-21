//
//  functions.h
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/13/16.
//  Copyright (c) 2016 c. All rights reserved.
//

#ifndef __KinDecEMDE__functions__
#define __KinDecEMDE__functions__

#include <stdio.h>

#endif /* defined(__KinDecEMDE__functions__) */




long double GAMMA_split(double a, double x);

int kronekerDetlta(int x ,int  y);

double deriv(double x, double (* func)(double x, void * params), void * params);

double gstarS(double x, void *params);
double gstar(double x, void *params);
double dlngstar_dlnT(double x, void *params);
double Fgstar(double a, void *params);

double T_kd(void *parameters);

double GammaTilde(void *parameters);

//double a_REH(void *parameters);

double gammakdS(void *parameters);

double a_Today(void *parameters);

double Gamma(void *parameters);

double H1(void *parameters);

double rho_critI(void *parameters);

double k_RH(void *parameters);

double cm(void *parameters);

double a_KD(struct Background *bkg_ivals, void *parameters);

double TFromRhoT(double r, void *parameters);
double dT_drhoT(double r, void *parameters);

double H(double a, double b, double c, void *parameters);

double E(double a, double b, double c);

double Upsi(double p, double m, double r, void *parameters);


double UpsiPar(void *bkgdens, void *parameters);

double gammaa(double r, double x, void *parameters);

double Rcr(double p, double m, double r,double x, void *parameters);

double R(double p, double m, double r,double a, void *parameters);

int ICBkgRho_TDM(struct Background *bkgdens, void *parameters);

void create_index_map();

void create_nl_tocol_map();

void create_output_files();

void create_Extra_indices_map();

int Initialize_Parameters(struct Params *parameters);

void create_output_files(struct Output_Files *pFiles);
void create_output_file_growth_function(struct Growth_Function_File *ptrFile);

void close_output_files(struct Output_Files *pFiles);
void close_GF_file(struct Growth_Function_File *ptrFile);


