//==================================================================================================
// We code a simplified version of Eq. 24 in Bertchiner '06: arXiv:astro-ph/0607319
// Namely we ignore the fourth term on the lhs (the one starting with k) and the delta_l0 term on
// the rhs. Those terms are codded as well but commented out.
// The integration variable is chosen as x=tau/tau_dec
// When using the Solver with numjac=1, i.e. numerical, compared with the GSL solver where I have
// coded the Jacobian analytically, for the system with Nmax=4 and Lmax=4 the performance is:
// 31s CPU time (This version)vs 0.7s CPU time (gsl_odeiv2_step_msbdf). Therefore we need to implement
// the sparse linear algebra methods to improve speed.
//==================================================================================================




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
// Definitions used for computing the terms in Eq.24 in Bertchinger 06
//==================================================================================================

#include <sys/types.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_errno.h>// GSL_SUCCESS ...
#include <gsl/gsl_math.h> // cos(x), sin(x) ...
#include <gsl/gsl_odeiv2.h>// ODE solver
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_gamma.h>

#define Lmax 15
#define Nmax 400



bool numjac=0;

#define GAMMA gsl_sf_gamma_inc

const extern size_t DIM;
const size_t DIM = (Nmax+1) * (Lmax +1); //Nmax and Lmax are defined in common.h Defines size of system.


using namespace std;
using namespace ODE_solver_Rec;

//==================================================================================================
// parameters. I have folded kappa (which is k in units of tau_decoupling inside the Params struct
//==================================================================================================


struct Params
{
    double xini, xfin;
    double c; //c is defined as sqrt(2T_d/m_DM)
    double kappa;//kappa = k*tau_dec
};

typedef struct Params* pParams;

struct Output_File
{
    ofstream FileGF;
};

typedef struct Output_File* pOutFile;




Params *glob_P=NULL;// could not use const anymore, because of the way I coded the functions

double glob_kappa;//kept it in case we want to return to glob_kappa & Params w.o. it
double glob_c;

//==================================================================================================
// create index map
//==================================================================================================
struct var_nl
{
    int n;
    int l;
};

vector<var_nl> indices;

void create_index_map()
{
    var_nl dum;
    
    for(int l=0; l<Lmax+1; l++)
    {
        dum.l=l;
        
        for(int n=0; n<Nmax+1; n++)
        {
            dum.n=n;
            indices.push_back(dum);
        }
    }
    
    return;
}

//==================================================================================================
// necessary functions
//==================================================================================================
long double GAMMA_split(double a, double x){
    
    /*
     
     GAMMA function that has included expansion at large x. gsp_sf_gamma_inc(a,x) returns only double.
     
     For a=3/4 if x> 0.195^-4~ 705 then GAMMA(3/4,x) can no longer be respresented as adouble
     
     since it becomes smaller than 2.3E-308, which is the absolute value of the smallest double.
     
     
     */
    
    int i;
    double uk;
    long double sum;
    sum =0.0;
    uk=1.0;
    
    if (x < pow(0.2,-4)) // 0.2 is valid only for GAMMA(3/4,x)!! It corresponds, roughly, to the value
        // where Gamma becomes smaller than 2.3xE-308, which is the smallest absolute
        // value of a double on my system  as checked with DBL_MIN
        
        return gsl_sf_gamma_inc(a,x);
    else{
        for (i =1;i < 5; i++) {
            sum += uk*pow(x, (-i+1));
            uk  *= (a-i);
            
        }
        return pow(x, (a-1.0))*expl(-x)*sum;
    }
}
//==================================================================================================

double deriv(double x, double (* func)(double x, void * params), void * params){
    
    /****************************************************************/
    /*    Takes x and a pointer to a function and returns deriv     */
    /****************************************************************/
    
    
    gsl_function F;
    double result, abserr;
    F.function =  func;
    F.params = params;
    
    gsl_deriv_central(&F, x, 1e-8, &result, &abserr);
    
    return result;
}
//==================================================================================================

double gamma_a_taudec(double x, void * params){
    
    return 2.0*pow(x,-5); //rad domination used here!
}
//==================================================================================================

double R_taudec(double x, void * parameters){
    
    //Params par = *(Params *)parameters;
    
    return 1/2.0/x; //rad domination used here!
}
//==================================================================================================

double Phi (double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    double t1 = par.kappa * x;
    double t2 = sqrt(0.3e1);
    double t4 = t2 * t1 / 0.3e1;
    double t5 = sin(t4);
    double t7 = cos(t4);
    double t13 = par.kappa * par.kappa;
    double t16 = x * x;
    
    return -(0.9e1 * t5 - 0.3e1 * t1 * t2 * t7) * t2 / t13 / par.kappa / t16 / x;
}
//==================================================================================================

double Psi (double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    double t1 = par.kappa * x;
    double t2 = sqrt(0.3e1);
    double t4 = t2 * t1 / 0.3e1;
    double t5 = sin(t4);
    double t7 = cos(t4);
    double t13 = par.kappa * par.kappa;
    double t16 = x * x;
    
    return -(0.9e1 * t5 - 0.3e1 * t1 * t2 * t7) * t2 / t13 / par.kappa / t16 / x;
}
//==================================================================================================

double dPsi_dx(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /*                                                              */
    /*   When expanding tho the general case we will have an indep  */
    /*     set of ODEs taht give the sols for Psi, and Phi. We can  */
    /*    then use it to find the derrivative (its the function that*/
    /*     defines the  system.)                                    */
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    gsl_function F;
    double result, abserr;
    F.function = & Psi;
    F.params = & par;
    
    gsl_deriv_central(&F, x, 1e-8, &result, &abserr);
    
    return result;
}
//==================================================================================================

double dPsi_dx_th(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    //    double kappa = par->kappa;
    
    
    double t1 = par.kappa * x;
    double t2 = sqrt(0.3e1);
    double t4 = t2 * t1 / 0.3e1;
    double t5 = sin(t4);
    double t8 = x * x;
    double t14 = cos(t4);
    double t20 = par.kappa * par.kappa;
    double t23 = t8 * t8;
    
    return -0.3e1 * t5 * t2 / par.kappa / t8 + (0.27e2 * t5 - 0.9e1 * t1 * t2 * t14) * t2 / t20 / par.kappa / t23;
}
//==================================================================================================

double dPhi_dx_th(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    //    double kappa = par->kappa;
    
    
    double t1 = par.kappa * x;
    double t2 = sqrt(0.3e1);
    double t4 = t2 * t1 / 0.3e1;
    double t5 = sin(t4);
    double t8 = x * x;
    double t14 = cos(t4);
    double t20 = par.kappa * par.kappa;
    double t23 = t8 * t8;
    
    return -0.3e1 * t5 * t2 / par.kappa / t8 + (0.27e2 * t5 - 0.9e1 * t1 * t2 * t14) * t2 / t20 / par.kappa / t23;
}
//==================================================================================================

double dPhi_dx(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06   */
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    gsl_function F;
    double result, abserr;
    F.function = & Phi;
    F.params = & par;
    
    gsl_deriv_central(&F, x, 1e-8, &result, &abserr);
    
    
    
    return result;
}
//==================================================================================================

double d2Psi_dx2_th(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06   */
    /*   WILL NEED TO THINK HOW TO GET THIS NUMERICALLY for Psi     */
    /*     and Phi given numerically by a call to the solver!        */
    /****************************************************************/
    
    //    params_struct par = *(params_struct *)params;
    
    
    Params * par = (Params *)parameters;
    
    
    double kappa = par->kappa;
    
    double t1 = kappa * x;
    double t2 = sqrt(0.3e1);
    double t4 = t1 * t2 / 0.3e1;
    double t5 = sin(t4);
    double t6 = kappa * kappa;
    double t9 = t6 * kappa;
    double t11 = cos(t4);
    double t12 = t11 * t2;
    double t16 = 0.1e1 / t9;
    double t17 = x * x;
    double t19 = 0.1e1 / t17 / x;
    double t32 = t17 * t17;
    
    return  -(0.3e1 * t5 * t6 + t9 * x * t12) * t2 * t16 * t19 + 0.18e2 * t5 * t2 / kappa * t19
    - (0.108e3 * t5 - 0.36e2 * t1 * t12) * t2 * t16 / t32 / x;
}
//==================================================================================================

double Tdec_to_Tl(double x, void *parameters){
    
    /****************************************************************/
    /*                                                              */
    /*                                                              */
    /*  Assumes tau does not cross g* jumps due to phase transtions */
    /*     and a~tau, i.e. RD!!                                     */
    /****************************************************************/
    
    return x;
}
//==================================================================================================

double Tdm_to_Tl(double x){
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*         This is Eq 12 in Bert 06. See if it needs changed    */
    /*         It assumes gamma ~T^6 and RD!!                       */
    /****************************************************************/
    
    double t1 = x * x;
    double t2 = t1 * t1;
    double t3 = 0.1e1 / t2;                // 0.1e1 --> 1.0???
    double t4 = pow(t3, 0.1e1 / 0.4e1);    // 0.1e1 / 0.4e1--> 0.25???
    double t5 = exp(t3);
    double t7 = GAMMA(0.3e1 / 0.4e1, t3);  // 0.3e1 / 0.4e1 --> 0.75???
    
    if (x > 0.2) // here 0.2 is the value where the calculation breaks down and we need to use asympt
        return  t4 * t5 * t7;
    else
        return (1.0-pow(x, 4)/4.0+5.0*pow(x, 8)/16.0 - 45.0*pow(x,12)/64.0);// using Gama(3/4,s)=1-1/(4s)+5/(16s^2)+... expansion
}
//==================================================================================================

double Tdm_to_Tl_orig(double x){
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*         This is Eq 12 in Bert 06. See if it needs changed    */
    /*         It assumes gamma ~T^6 and RD!!                       */
    /****************************************************************/
    
    double t1 = x * x;
    double t2 = t1 * t1;
    double t3 = 0.1e1 / t2;
    double t4 = pow(t3, 0.1e1 / 0.4e1);
    long double t5 = expl(t3);
    long double t7 = GAMMA_split(0.3e1 / 0.4e1, t3);
    
    return  t4 * t5 * t7;
}
//==================================================================================================

double A_n(double x, int n){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /*         This is Eq 25 in Bert 06. See if it needs changed    */
    /****************************************************************/
    
    return pow((1-Tdm_to_Tl(x)),n);
}
//==================================================================================================

double A_np(double x, void * params){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /*         This is Eq 25 in Bert 06. See if it needs changed    */
    /****************************************************************/
    
    int n = *(int *)params;
    
    return pow((1-Tdm_to_Tl(x)),n);
}
//==================================================================================================

double B_n(double x, int n){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /*         This is Eq 25 in Bert 06. See if it needs changed    */
    /****************************************************************/
    
    return n * Tdm_to_Tl(x) * pow((1-Tdm_to_Tl(x)),n-1);
}
//==================================================================================================

double B_np(double x, void * parameters){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /*         This is Eq 25 in Bert 06. See if it needs changed    */
    /****************************************************************/
    
    int n = *(int *)parameters;
    
    return n * Tdm_to_Tl(x) * pow((1-Tdm_to_Tl(x)),n-1);
}
//==================================================================================================

double dA_n_dx_zero(double x, int i){
    
    /*Check with Maple */
    
    double t1 = x * x;
    double t2 = t1 * t1;
    double t4 = t2 * t2;
    double t10 = pow(t2 / 0.4e1 - 0.5e1 / 0.6e1 * t4 + 0.15e2 / 0.8e1 * t4 * t2, (double) (i - 1));
    
    return(0.6e1 - 0.40e2 * t2 + 0.135e3 * t4) * t1 * x * (double) i * t10 / 0.6e1;
}
//==================================================================================================

double dA_n_dx(double x, int i){
    
    
    double t1 = x * x;
    double t2 = t1 * t1;
    double t3 = 0.1e1 / t2;
    double t4 = pow(t3, 0.1e1 / 0.4e1);
    double t5 = exp(t3);
    double t7 = GAMMA(0.3e1 / 0.4e1, t3);
    double t11 = pow(0.1e1 - t4 * t5 * t7, (double) (i - 1));
    double t16 = exp(-t3);
    double t18 = t4 * t4;
    double t19 = t18 * t4;
    double t23 = t2 * t2;
    
    if (x>0.2) {
        return - t11 * (double) i * t5 * (-t7 * t2 - 0.4e1 * t7 + 0.4e1 * t16 * t2 * t19) / t23 / x / t19;
    }
    else
        return dA_n_dx_zero(x, i) ;
    
}
//==================================================================================================

double dB_n_dx(double x, int i){
    
    /* Need to include the x>0.2 case as well as for dA_n_dx-->worth thinking about ~Sheridan */
    
    
    double t1 = x * x;
    double t2 = t1 * t1;
    double t3 = 0.1e1 / t2;
    double t4 = GAMMA(0.3e1 / 0.4e1, t3);
    double t7 = exp(-t3);
    double t9 = pow(t3, 0.1e1 / 0.4e1);
    double t10 = t9 * t9;
    double t11 = t9 * t10;
    double t14 = exp(t3);
    double t15 = t14 * t9;
    double t16 = t4 * t4;
    double t17 = (double) i * t16;
    double t30 = pow(0.1e1 - t4 * t15, (double) (i - 2));
    double t31 = t2 * t2;
    
    return 0.1e1 / t11 / x / t31 * t30 * t14 * (double) i * (-t2 * t4 - 0.4e1 * t4 + 0.4e1 * t11 * t2 * t7 + t2 * t17 * t15 + 0.4e1 * t17 * t15 - 0.4e1 * (double) i * t4);
}
//==================================================================================================

double ul_to_taudec(double x, void *parameters){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /* This is Eq 26b in Bert 06.  needs changed for general case   */
    /****************************************************************/
    
    Params * par = (Params *)parameters;
    
    double kappa = par->kappa;
    
    
    double t1 =  (kappa * kappa);
    double t3 =  (x * x);
    double t9 = sqrt(0.3e1);
    double t11 = kappa * x * t9 / 0.3e1;
    double t12 = sin(t11);
    double t14 = cos(t11);
    double t18 = 0.1e1 / kappa / x;
    double calH_ul = -t18 * ((double) (3 - 18 / t1 / t3) * t12 + 0.6e1 * t14 * t9 * t18) * t9 / 0.2e1;//Eq 26b Bert06
    
    double calH_taudec = pow(x,-1); //Assumes RD and calH ~ 1/tau
    
    return calH_ul/calH_taudec;
}
//==================================================================================================

double deltaL(double x , void *parameters){
    /****************************************************************/
    /*			Chked with Maple                                    */
    /* This is Eq 26c in Bert 06.  needs changed for general case   */
    /****************************************************************/
    Params * par = (Params *)parameters;
    
    double kappa = par->kappa;
    
    return 0.3e1 / 0.2e1 / kappa / x * (0.3e1 * sin(kappa * x * sqrt(0.3e1) / 0.3e1) - kappa * x * sqrt(0.3e1) * cos(kappa * x * sqrt(0.3e1) / 0.3e1)) * sqrt(0.3e1) + 0.3e1 / 0.2e1 * (0.3e1 * (0.1e1 - 0.6e1 * pow(kappa, -0.2e1) * pow(x, -0.2e1)) * sin(kappa * x * sqrt(0.3e1) / 0.3e1) + 0.6e1 * cos(kappa * x * sqrt(0.3e1) / 0.3e1) / kappa / x * sqrt(0.3e1)) / kappa / x * sqrt(0.3e1);
}
//==================================================================================================

double deltaTl_to_Tdm(double x, void * parameters){
    
    /****************************************************************/
    /*         Code checked with Maple                              */
    /*                                                              */
    /*  RD assumed -> perf fluid trans functions as Eq.26 Bert 06 	*/
    /****************************************************************/
    
    Params par = *(Params *)parameters;
    
    double deltaTL_over_TL;
    
    deltaTL_over_TL = deltaL(x, &par)/3.0; //Eqn 26.c in Bertschinger 06
    return deltaTL_over_TL/(Tdm_to_Tl(x));
}





//==================================================================================================

//==================================================================================================
// derivative function wrapper
//==================================================================================================
void fcn_k_eval(int *neq, double *t, double *y, double *f, int col)
{
    double x=*t, c=glob_P->c, kappa=glob_P->kappa;
    
    //    =================================================================
    //     Below we code a simplified version of Eq 24 in Bert'06
    //    =================================================================
    
    //    *****************************************************************
    //
    //     This func codes the ODEs to solve as following:
    //     f  functions that define the system for y(x) and params
    //     dydx[0]=f[0], etc. Not to be confused with fnl we solve for
    //
    //     we treat f_nl(x) as a row column ordered 1D array
    //     of dim (1+Nmax)*(1+Lmax): y[l*(Nmax+1)+n]=f[n][l]
    //     Since Nmax>>Lmax its more efficient to have N in the inner loop
    //
    //    *****************************************************************
    
    int l, n; // indices labeling n and l of Eq. 26 Bert '06
    int inl, inpl, inml, inlp, inlm, inmlp, inplm; // 1D indices corresponding to (n,l), (n+1,l), (n-1,l), etc. touples.
    
    
    double C_nl, C_nml;    // The coefficients of the terms on the RHS of Eq. 24 Bert'06. E.g. C_nml= 2*n*R (moved to rhs)
    double C_k_lp, C_k_lm; // The term prop to k on Eq. 24
    double D_l0, D_l1;     // The driving terms prop to delta_l0 and delta_l1 respectively
    
    //=====================================================================
    // evaluate all time-dependent variables outside loops. For efficiency
    // this is important and it also avoids double evaluations...
    //=====================================================================
    double R=R_taudec(x, glob_P);
    double gamma_a=gamma_a_taudec(x, glob_P);
    double dPsi=dPsi_dx_th(x, glob_P);
    double Phiv=Phi(x, glob_P);
    double ul=ul_to_taudec(x, glob_P);
    double Tdec_Tl=Tdec_to_Tl(x, glob_P);
    double deltaTl_Tdm=deltaTl_to_Tdm(x, glob_P);
    double kcT=-kappa*c*pow(Tdec_Tl, -1.0/2.0);
    
    for(l=0; l<Lmax+1; l++){
        for ( n=0; n<Nmax+1; n++){
            
            inl   = l*(Nmax+1)+n;
            inpl  = l*(Nmax+1)+n+1;
            inml  = l*(Nmax+1)+n-1;
            inlp  = (l+1)*(Nmax+1)+n;
            inlm  = (l-1)*(Nmax+1)+n;
            inmlp = inlp -1;
            inplm = inlm +1;
            
            C_nl  = -(2.0*n+l)*(gamma_a + R); //1st term of Eq.24 Bert'06
            C_nml =   2.0*n*R;                // 2nd term of Eq.24
            C_k_lp=C_k_lm=0.0;
            
            // if statements cost close to nothing but make things so much easier to read...
            if(l<Lmax) {
                
                C_k_lp+=(1.5+n+l)*y[inlp];
                if(n>0) C_k_lp+=-n*y[inmlp];
                
                C_k_lp*=kcT*(1.0+l)/(2.0*l+1.0);
            }
            
            if(l>0) {
                
                C_k_lm+=-y[inlm];
                if(n<Nmax) C_k_lm+=y[inplm];
                
                C_k_lm*=(kcT*l)/(2.0*l+1.0);
            }
            
            // driving terms ('if' here is also much better than evaluating things but then multiplyin by zero...)
            D_l0 = (l==0 ?  3.0*dPsi*A_n(x, n) - 2.0*(dPsi+gamma_a*deltaTl_Tdm)*B_n(x,n) : 0.0); //Driving term prop to delta_l0
            D_l1 = (l==1 ? (kappa/3.0)*(2.0/c)*sqrt(Tdec_Tl)*(Phiv+gamma_a*ul)*A_n(x, n) : 0.0); //Driving term prop to delta_l1
            
            // ode terms
            f[inl]=C_nl*y[inl] + C_nml*y[inml] + C_k_lp + C_k_lm + D_l0 + D_l1; //The full system
        }
    }
    
    return;
}

//==================================================================================================
// jacobian function wrapper
//==================================================================================================
void jac_k_eval(int *neq, double *t, double *y, double *f, int col)
{
    double x=*t, c=glob_P->c, kappa=glob_P->kappa;
    
    //=====================================================================
    // get indices from index map
    //=====================================================================
    int n=indices[col].n, l=indices[col].l;
    int inl, inpl, inml, inlp, inlm, inmlp, inplm; // 1D indices corresponding to (n,l), (n+1,l), (n-1,l), etc. touples.
    
    //=====================================================================
    // comment: one could make sure that the time-update is only done once
    // since the jacobian routine will be called many times with the same x
    // This could be fully consolidated with fcn call too...
    //=====================================================================
    double R=R_taudec(x, glob_P);
    double gamma_a=gamma_a_taudec(x, glob_P);
    double Tdec_Tl=Tdec_to_Tl(x, glob_P);
    double kcT=-kappa*c*pow(Tdec_Tl, -1.0/2.0);
    
    //=====================================================================
    // comment: one could precompute a map from (n, l) to col and avoid all
    // these operations...
    //=====================================================================
    inl   = l*(Nmax+1)+n;
    inpl  = l*(Nmax+1)+n+1;
    inml  = l*(Nmax+1)+n-1;
    inlp  = (l+1)*(Nmax+1)+n;
    inlm  = (l-1)*(Nmax+1)+n;
    inmlp = inlp -1;
    inplm = inlm +1;
    
    // jacobian terms
    f[inl]= -(2.0*n+l)*(gamma_a + R);  //1st term of Eq.24 Bert'06
    if(n<Nmax) f[inpl]= 2.0*(1.0+n)*R;  // 2nd term of Eq.24 but shifted
    
    if(l>0)
    {
        f[inlm] =(kcT*l)/(2.0*l-1.0)*(0.5+n+l);
        if(n<Nmax) f[inplm]=-(kcT*l)/(2.0*l-1.0)*(1.0+n);
    }
    
    if(l<Lmax)
    {
        if(n>0) f[inmlp]=kcT*(1.0+l)/(2.0*l+3.0);
        f[inlp]=-kcT*(1.0+l)/(2.0*l+3.0);
    }
    
    return;
}

//==================================================================================================
// jacobian and fcn evaluation depending on col
//==================================================================================================
void fcn_k(int *neq, double *t, double *y, double *f, int col)
{
    //==============================================================================================
    // reset everything
    //==============================================================================================
    for(int i=0; i<*neq; i++) f[i]=0.0;
    
    //==============================================================================================
    // when choosing numerical jacobian setup, just return directly afterwards...
    //==============================================================================================
    if(numjac){ fcn_k_eval(neq, t, y, f, col); return; }
    
    //==============================================================================================
    // normal fcn evaluation
    //==============================================================================================
    if(col<0) fcn_k_eval(neq, t, y, f, col);
    
    //==============================================================================================
    // jacobian evaluation
    //==============================================================================================
    else jac_k_eval(neq, t, y, f, col);
    
    return;
}

//==================================================================================================
// output of solution
//==================================================================================================

void output_current_solution_k(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    ofile << Sz.z << " ";
    
    for(int k=0; k<(int)Sz.y.size(); k++) ofile << Sz.y[k] << " "; 
    for(int k=0; k<(int)Sz.y.size(); k++) ofile << Sz.dy[k] << " ";
    
    ofile << endl;
    
    return;
}


//==================================================================================================
// output of selected solution
//==================================================================================================

void output_selected_solution_k(ofstream &ofile, const ODE_solver_Solution &Sz)
{
	//important: selected solution file outputs x, kappa, delta DM, which is what we need!
    double x = Sz.z, Three_CalH_u= 4.5*glob_P->c/glob_P->kappa*pow(x,-1.5)*Sz.y[(Nmax+1)], ThetaT=1.5*glob_P->c*glob_P->kappa*pow(x,-0.5)*Sz.y[(Nmax+1)];
    ofile << Sz.z << " "<<glob_P->kappa<<" "<<Sz.y[0]<<" "<<Three_CalH_u<<" "<< Sz.y[0]+Three_CalH_u<<" "<<deltaL(x, glob_P)<<" "<<ThetaT; //for each x, this outputs all of the density perturbations
    ofile << endl;
    
    return;
}


//==================================================================================================
// output Growth Function
//==================================================================================================

void output_growth_function(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    double x = Sz.z, Three_CalH_u= 4.5*glob_P->c/glob_P->kappa*pow(x,-1.5)*Sz.y[(Nmax+1)];
    ofile <<glob_P->kappa<<" "<< Sz.y[0]+Three_CalH_u<<" "<<Sz.y[0]<<" "<< x;
    ofile << endl;
    
    return;
}

void create_output_file_growthfunction(struct Output_File *ptrFile){
    
    cout << " Creating Output File for Growth function"<<endl;
    char pGFFile[100];
    double c=glob_P->c;
    
   
    sprintf(pGFFile, "./DATAFILES/DM_GrFunc_KDRD_Nmax%dLmax%dc_%0.E.dat",Nmax, Lmax,glob_c);
    
    ptrFile->FileGF.open(pGFFile);
    
    ptrFile->FileGF.precision(6);
    ptrFile->FileGF <<"# kappa"<<" nu "<<" delta_DM "<<" x "<<endl;
    
}

void close_growth_function_file(struct Output_File *ptrFile){
    ptrFile->FileGF.close();
}


//==================================================================================================
// setup for background solution
//==================================================================================================
double setup_kmode(void * parameters, struct Output_File *pFile)
{
    //==============================================================
    
    Params  P = *(Params *)parameters;
    
    //    glob_P=&P;
    double kappa=P.kappa;
    
    int nsteps=200;
    
    int neq=DIM;
    
    vector<double> xvals(nsteps);
    
    double xs=P.xini, xend;
    
    //==============================================================
    // initialize time axis
    //==============================================================
    init_xarr(P.xini, P.xfin, &xvals[0], nsteps, 1, 1); // 5th argument ZERO FOR LINEAR grid
    
    //==============================================================
    // solution setup
    //==============================================================
    ODE_solver_Solution Sx;
    ODE_solver_accuracies tols;
    set_verbosity(-1);
    
    Sx.y .resize(neq, 0.0);        // this sets all to zero initially...
    Sx.dy.resize(neq, 0.0);
    
    tols.rel.resize(neq, 1.0e-4);
    tols.abs.resize(neq, 1.0e-8);
    
    //==============================================================
    // initial solution: Eq 27 in Bert 06
    //==============================================================
    
    Sx.z=P.xini;  // Because of the way ODE_solver_Solution class is written cannot use Sx.x which would be more natural. Fix later
    Sx.y[0]= deltaL(P.xini, &P);
    Sx.y[(Nmax+1)]=(2.0*kappa/3.0)*ul_to_taudec(P.xini, &P)*(1.0/P.c)*sqrt(Tdec_to_Tl(P.xini,&P));
    Sx.y[1]=-deltaL(P.xini,&P)/3.0;
    
    //==============================================================
    ODE_Solver_data ODE_Solver_info;
    ODE_Solver_info.Snewptr=NULL;
    ODE_Solver_set_up_solution_and_memory(Sx, tols, ODE_Solver_info,
                                          fcn_k, numjac);
    
    cout << " computing evolution for kappa= " << kappa << endl;
    
    char buffer_DataFile[100], buffer_FullDataFile[100];
    
    sprintf(buffer_DataFile, "./DATAFILES/DM_Pert_KDRD_Nmax%dLmax%dkappa_%.1Ec_%0.E.dat",Nmax, Lmax, kappa,P.c); //not sure if we need every single data file but overwriting is a waste
    sprintf(buffer_FullDataFile, "./DATAFILES/DM_FullPert_KDRD_Nmax%dLmax%dkappa_%.0Ec_%0.E.dat",Nmax, Lmax, kappa,P.c); //this one overwrites and takes a lot of space.... comment out while not being used?
    
    ofstream ofile(buffer_DataFile);
    ofile.precision(6);
    ofile<<"# x"<< " kappa " <<"deltaDM"<<"  3Cal_Hu"<<" nu"<<" deltaL"<<" ThetaDM_T"<<endl;
    
    ofstream ofileFull(buffer_FullDataFile);
    ofileFull.precision(6);
    ofileFull <<"# x"<<" f_nl with n first then l"<<"df_nl/dx "<<endl;
    
    
    
    create_index_map();
    
    for(int i=1; i<nsteps; i++)
    {
        

        
        xend=xvals[i];
        
        ODE_Solver_Solve_history(xs, xend,
                                 1.0e-5*fabs(xs-xend), // minimal time-step
                                 5.0e-1*fabs(xs-xend), // maximal time-step
                                 Sx, ODE_Solver_info);
        
        output_current_solution_k(ofileFull,Sx);
        output_selected_solution_k(ofile,Sx);
        
        cout << "Moments computed for t = " << xend << endl;

        xs=xend;
    }
    
    output_growth_function(pFile->FileGF,Sx);
    ofile.close();
    ofileFull.close();
    
    clear_memory(Sx,tols);
    
    //==============================================================
    return 0;
}

//==================================================================================================
// call code
//==================================================================================================
int main(int narg, char *args[])
{
    double xini, xfin= pow(10.0,6.0); //72;
    double kappa_i=0.2, kappa_f=pow(10.0,2); //pow(10.0,-1.5), kappa_f=pow(10.0,2);
    double nk_grid=100;//200;
    double kappa;
    Params P;
    glob_P=&P;
    
    
    if (narg<2){
        //Expecting at least 2 arguments: the program name, and c
        cerr<<"Usage: "<<args[0]<<" c "<<args[1]<<endl;
    }
    
    double c = atof(args[1]); // if c=0.24 in EMDE, then c=0.2617 here
    glob_c=c;
    
    vector<double> kvals(nk_grid);
    
    
    
    //==============================================================
    // initialize kappa grid
    //==============================================================
    init_xarr(kappa_i, kappa_f, &kvals[0], nk_grid, 1, 1); // 5th argument ZERO FOR LINEAR grid
    
    
    time_t  t0, t1; /* time_t is defined on <time.h> and <sys/types.h> as long */
    clock_t c0, c1; /* clock_t is defined on <time.h> and <sys/types.h> as int */
    
    
    t0 = time(NULL);
    c0 = clock();
    
    printf ("\tbegin (wall):            %ld\n", (long) t0);
    printf ("\tbegin (CPU):             %d\n\n", (int) c0);
    

    Output_File File;
    create_output_file_growthfunction(&File);
    
    for(int i=0; i < kvals.size(); i++){
        kappa=kvals[i];
        xini=min(0.1, 1.0/kappa);
        P.xini=xini;
        P.xfin=xfin;
        P.c=c;
        P.kappa=kappa;
        
        setup_kmode(glob_P,&File);
    }
    
    close_growth_function_file(&File);


   
    
    t1 = time(NULL);
    c1 = clock();
    
    printf ("\xfin (wall):              %ld\n", (long) t1);
    printf ("\xfin (CPU);               %d\n", (int) c1);
    printf ("\telapsed wall clock time: %ld\n", (long) (t1 - t0));
    printf ("\telapsed CPU time:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
    
    
    //==================================================================================================
    // Testing area below
    //==================================================================================================
    
    
    
    
    
    
    return 0;
}

//==================================================================================================
//==================================================================================================

//==================================================================================================
//==================================================================================================


