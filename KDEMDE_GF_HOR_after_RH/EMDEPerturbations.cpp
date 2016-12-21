//
//  EMDEPerturbations.cpp
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/18/16.
//  Copyright (c) 2016 c. All rights reserved.
//

#include "EMDEPerturbations.h"
#include<vector>

#include "functions.h"
#include "common.h"

//==================================================================================================
//  Glob Parameters and indices
//==================================================================================================

Params *glob_Par;
ExtraIndices *pExtInd;


vector<var_nl> indices;

int column[Nmax+1][Lmax+1];

int countRD=0;





//==================================================================================================
// derivative function wrapper for ODE System coded using a as integration variable.
//==================================================================================================
void fcn_k_eval(int *neq, double *t, double *y, double *f, int col)
{
    
    
    double a = *t;
    
    //==================================================================================
    // Setup Indices and parameters and common combinations of functions that enter the ODEs
    //==================================================================================
    
    int iBkgP   = pExtInd->BkgP, iBkgR= pExtInd->BkgR, iBkgM=pExtInd->BkgM, iTDM=pExtInd->BkgT;
    int ideltaP = pExtInd->deltaP, ideltaR= pExtInd->deltaR, ithetaP=pExtInd->thetaP, ithetaR=pExtInd->thetaR;
    int iPhi =pExtInd->Phi;
    int ideltaM= column[0][0];//column[n][l] is the generic column index for f_nl.
    int ithetaM= column[0][1];//column[n][l] is the generic column index for f_nl.
    
    double c=glob_Par->c, kT=glob_Par->kTilde,TkdS=glob_Par->T_kdS, m_DM=glob_Par->m_DM;
    double GT = GammaTilde(glob_Par);
    double a2 = pow(a,2);
    
    double rhoTp = pow(a, -3)*y[iBkgP], rhoTr=pow(a, -4)*y[iBkgR], rhoTm=pow(a, -3)*y[iBkgM];
    //converting back to the tilde Energy densities. All functions in functions.cpp are coded in terms of the 'tilde' Energy densities.
    double T_DM = y[iTDM]/a2;// We absorb the a^-2 scaling for the DM temperature
    
    double Upsilon = Upsi(rhoTp, rhoTm, rhoTr, glob_Par);
    double T_L = TFromRhoT(rhoTr, glob_Par);
    //cout << "T_L: " << T_L << " rhoTr: " << rhoTr << " rho_critI: " <<rho_critI(&glob_Par) <<  endl;
    double TDM_to_TL=T_DM/T_L, TL_to_TDM=1.0/TDM_to_TL;
    double Fgstr = Fgstar(rhoTr, glob_Par);
    
    double E = sqrt(rhoTp+rhoTr+rhoTm);
    double Ea2 = E*a2;
    double Over_Ea2 = 1.0/Ea2;
    double kT2=pow(kT, 2);
    double R_To_H1Ea2 = (1.0+0.125*Fgstr*(GT*rhoTp/rhoTr/E - 4.0))/a;
    double G = (rhoTp/rhoTr)*GT/E/a;
    
    double dPhi_da = -1.0/a*((1+pow(kT/E/a, 2)/3)*y[iPhi] -0.5*pow(E,-2)*(rhoTp*y[ideltaP]+rhoTr*y[ideltaR]+rhoTm*y[ideltaM]));
    
    double thetaTm = 1.5*kT*c*sqrt(T_L/TkdS)*y[ithetaM];
    //cout << "CURRENT VALUE OF C: " <<  c << endl;
    
    double kcT  = (decoupled ? -kT*Over_Ea2*sqrt(2.0*T_L/m_DM) : -kT*c*Over_Ea2*sqrt(T_L/TkdS));// full coeff of kc*sqrt{TL} term in Eq 18 of notes
    double kcT1 = (decoupled? (kT/3.0)*2.0*Over_Ea2*sqrt(m_DM/T_L) : (kT/3.0)*(2.0/c)*Over_Ea2*sqrt(TkdS/T_L));// full coef of deltal1 term In Eq18.
    
    
    double C_nl, C_nml;    // The coefficients of the terms on the RHS of Eq. 24 Bert'06. E.g. C_nml= 2*n*R (moved to rhs)
    double C_k_lp, C_k_lm; // The term prop to k on Eq. 24
    double D_l0, D_l1;     // The driving terms prop to delta_l0 and delta_l1 respectively
    
   
    
    vector<double> A(Nmax+1);
    vector<double> B(Nmax+1);
    
    for (int i=0; i<Nmax+1; i++)
    {
        A[i]=(decoupled? 1:pow((1-TDM_to_TL),i));
        B[i]=(decoupled? 0: (i==0? 0: i*TDM_to_TL*pow((1-TDM_to_TL),i-1)));
    }
    
    
    //============================================================
    //  Background Boltzmann ODEs. Eqns 8-10 & 12 in pdf notes "EMDE_ODEs..."
    //============================================================
    
    f[iBkgP] = -GT/a*y[iBkgP]/E;//phi
    
    f[iBkgR] =  GT*y[iBkgP]/E;//r
    
    f[iBkgM] = 0.0;//m
    
    f[iTDM] = -2.0*a*Upsilon*(T_DM-T_L);//DM Temp scaled with a^2.
    
    //============================================================
    //  Petrurbations Lepton, Scalar and the Grav Potential. Eqns 13-17 pdf notes
    //============================================================
    f[ideltaP]= -Over_Ea2*y[ithetaP] -3.0*dPhi_da + GT/E/a*y[iPhi];
    
    f[ithetaP]= -y[ithetaP]/a-kT2*Over_Ea2*y[iPhi];
    
    f[ideltaR]= -4.0/3.0*Over_Ea2*y[ithetaR] - 4.0*dPhi_da + G*(y[ideltaP]-y[ideltaR]-y[iPhi]);
    
    f[ithetaR]= -kT2*Over_Ea2*(y[iPhi]-0.25*y[ideltaR]) +G*(0.75*y[ithetaP]-y[ithetaR])
                 +(decoupled? 0: 0.75*rhoTm/rhoTr*Upsilon*(thetaTm - y[ithetaR])/a);
    
    f[iPhi]= dPhi_da;
    
    //============================================================
    // switch radiation evolution off when it becomes small
    //============================================================
//    if(abs(y[ithetaR])<1.0e-5 && abs(y[iPhi])<1.0e-10) f[ideltaR]=f[ithetaR]=0.0;
    
    if(abs(y[ithetaR])<1.0e-5 && abs(y[iPhi])<1.0e-10)
    {
        f[ideltaR]=f[ithetaR]=0.0;
      
    }
    
    //============================================================
    // switch rhoPhi evolution when it becomes small
    //============================================================
    if(abs(y[iBkgP])<1.0e-20) // 1.e-20  is 10 orders of magnitude below the constant y[iBkgM] so at that point the scalar field is essentially off!
    {
        f[iBkgP]=0;
    }

    
    //============================================================
    //  Petrurbations DM Eqn 26 Bert 06. See Eq. 18 of pdf notes
    //============================================================
    
    for(int l=0; l<Lmax+1; l++){
        for (int  n=0; n<Nmax+1; n++){
            
            
            C_nl  = -(2.0*n+l)*(Upsilon/a + R_To_H1Ea2); //1st term of Eq.24 Bert'06
            C_nml =   (n==0 ? 0.0: 2.0*n*R_To_H1Ea2*y[column[n-1][l]]);                // 2nd term of Eq.24
            C_k_lp=C_k_lm=0.0;
            
            //TRUNCATIONS
            if(l<Lmax) {
                
                C_k_lp+=(1.5+n+l)*y[column[n][l+1]];
                if(n>0) C_k_lp+=-n*y[column[n-1][l+1]];
                
                C_k_lp*=kcT*(1.0+l)/(2.0*l+1.0);
            }
            
            if(l>0) {
                
                C_k_lm+=-y[column[n][l-1]];
                if(n<Nmax){
			C_k_lm+=y[column[n+1][l-1]];                
                } else { //when n=Nmax
                        C_k_lm+=( (2.0*y[column[n][l-1]]) - y[column[n-1][l-1]]);
                        }
                //new truncation scheme above
                C_k_lm*=(kcT*l)/(2.0*l+1.0);
            }
            // driving terms
            D_l0 = (l==0 ?  -3.0*dPhi_da*A[n] - 2.0*(-dPhi_da + Upsilon/a*0.25*Fgstr*y[ideltaR]*TL_to_TDM)*B[n] : 0.0); //Driving term prop to delta_l0.
            D_l1 = (l==1 ? kcT1*(-y[iPhi]+Upsilon*E*a*y[ithetaR]/kT2)*A[n] : 0.0); //Driving term prop to delta_l1
            
            // ode terms
            f[column[n][l]]=C_nl*y[column[n][l]] + C_nml + C_k_lp + C_k_lm + D_l0 + D_l1; //The DM Pert ODEs in Eq 18 of pdf notes
            
            
        }
    }
    
    return;
}


//==================================================================================================
// jacobian function wrapper
//==================================================================================================
double loc_df_dx_2point(const double &fp1, const double &fm1, const double &h, double eps)
{
    if(fm1==0.0) return (fp1-fm1)/(2.0*h);
    double dum=fp1/fm1-1.0;
    return fabs(dum)<=eps ? 0.0 : dum*fm1/(2.0*h);
}

void jac_k_eval(int *neq, double *t, double *y, double *f, int col)
{
    //==============================================================================================
    // Setup Indices
    //==============================================================================================
    int iBkgP = pExtInd->BkgP, iBkgR= pExtInd->BkgR, iBkgM=pExtInd->BkgM, iTDM=pExtInd->BkgT;
    int idelP = pExtInd->deltaP, ithP=pExtInd->thetaP, idelR=pExtInd->deltaR, ithR=pExtInd->thetaR, iPhi=pExtInd->Phi;
    
    //==============================================================================================
    // for all the messy variables (which are a few), just do numerical derivatives
    //==============================================================================================
    if (col==iBkgP || col==iBkgR|| col==iBkgM || col==iTDM ||
        //col==idelP || col==ithP || col==idelR || col==ithR || col==iPhi)
        col==idelP || col==idelR || col==ithR || col==iPhi)
    {
        //===============================================================================
        // r[i] is reset here
        // y[i] should contain the solution at z
        //===============================================================================
        double y0=y[col], Dyj=y[col]*1.0e-8, eps=1.0e-12;
        if(y0==0.0) Dyj=1.0e-8;
        
        vector<double> fp1(*neq), fm1(*neq);
        
        //===============================================================================
        // derivatives with respect to Xj.
        // A two-point formula is used. It has accuracy O(h^2)
        //===============================================================================
        // get f(yj+Dyj)
        y[col]=y0+Dyj;
        fcn_k_eval(neq, t, &y[0], &fp1[0], col);
        // get f(yj-Dyj)
        y[col]=y0-Dyj;
        fcn_k_eval(neq, t, &y[0], &fm1[0], col);
        // restore y again
        y[col]=y0;
        
        //===============================================================================
        // define numerical derivative
        //===============================================================================
        for(int k=0; k<*neq; k++) f[k]=loc_df_dx_2point(fp1[k], fm1[k], Dyj, eps);
        
        //===============================================================================
        // redo some of the derivatives to get better values
        //===============================================================================
        double a=*t, kT=glob_Par->kTilde;
        double rhoTp = pow(a, -3)*y[iBkgP], rhoTr=pow(a, -4)*y[iBkgR], rhoTm=pow(a, -3)*y[iBkgM];
        double E = sqrt(rhoTp+rhoTr+rhoTm);
        double Gp = GammaTilde(glob_Par)/E/a;
        
        if(col==idelP)
        {
            double deldPhidx_delP=0.5*pow(E,-2)*rhoTp/a;
            
            f[idelP]=-3.0*deldPhidx_delP;
            f[ithP] = 0.0;
            
            double G = (rhoTp/rhoTr)*Gp;
            f[idelR]=G-4.0*deldPhidx_delP;
            f[ithR] = 0.0;
            
            f[iPhi] = deldPhidx_delP;
        }

        if(col==iPhi)
        {
            double deldPhidx_dPhi=-(1.0+pow(kT/E/a, 2)/3.0)/a;
            double xi=kT*kT/E/a/a;
            
            f[idelP]=Gp-3.0*deldPhidx_dPhi;
            f[ithP] =-xi;
            
            double G = (rhoTp/rhoTr)*Gp;
            f[idelR]=-G-4.0*deldPhidx_dPhi;
            f[ithR]=-xi;
            
            f[iPhi] = deldPhidx_dPhi;
        }

        return;
    }
    
    // this is bulk of the matrix and it is done analytically
    else{

        //=====================================================================
        // Setup Parameters and common combinations
        //=====================================================================
        
        double a=*t;
        double c=glob_Par->c, kT=glob_Par->kTilde, TkdS=glob_Par->T_kdS, m_DM=glob_Par->m_DM;
        double GT = GammaTilde(glob_Par);
        
        double rhoTp = pow(a, -3)*y[iBkgP], rhoTr=pow(a, -4)*y[iBkgR], rhoTm=pow(a, -3)*y[iBkgM];
        double Upsilon = Upsi(rhoTp, rhoTm, rhoTr, glob_Par);
        double T_L = TFromRhoT(rhoTr, glob_Par);
        double Fgstr = Fgstar(rhoTr, glob_Par);
        
        double E = sqrt(rhoTp+rhoTr+rhoTm);
        double R_To_H1Ea2 = (1.0+0.125*Fgstr*(GT*rhoTp/rhoTr/E - 4.0))/a;
        double kcT  = (decoupled ? -kT*sqrt(2.0*T_L/m_DM) : -kT*c*sqrt(T_L/TkdS))/E/a/a;// full coeff of kc*sqrt{TL} term in Eq 18 of notes
        
        //=====================================================================
        // Rad, Scalar and Phi Perturbations
        //=====================================================================
        int idelM= column[0][0];//column[n][l] is the generic column index for f_nl.
        int itheM= column[0][1];//column[n][l] is the generic column index for f_nl.
        
        if(col==idelM)// && abs(y[iPhi])>1.0e-10)
        {
            double deldPhidx_delM=0.5*pow(E,-2)*rhoTm/a;
            f[idelP] = -3.0*deldPhidx_delM;
            f[idelR] = -4.0*deldPhidx_delM;
            f[iPhi]  = deldPhidx_delM;
        }
        
        if(col==itheM && !decoupled)
        {
            double thetaTm_dtheM = 1.5*kT*c*sqrt(T_L/TkdS);
            f[ithR]=0.75*rhoTm/rhoTr*Upsilon*thetaTm_dtheM/a;
        }

        if(col==ithP)
        {
            f[ithP]=-1.0/a;
            f[idelP]=-1.0/E/a/a;
            
            double G = (rhoTp/rhoTr)*GT/E/a;
            f[ithR]=0.75*G;
            return;
        }

        //============================================================
        //  Perturbations DM Eqn 26 Bert 06
        //============================================================
        // jacobian terms
        
        int n=indices[col].n, l=indices[col].l;
        
        f[column[n][l]]= -(2.0*n+l)*(Upsilon/a + R_To_H1Ea2);  // 1st term of Eq.24 Bert'06
        
        if(n<Nmax) f[column[n+1][l]]= 2.0*(1.0+n)*R_To_H1Ea2;      // 2nd term of Eq.24 but shifted
        
        if(l>0)
        {
            f[column[n][l-1]]=(kcT*l)/(2.0*l-1.0)*(0.5+n+l);
            if(n<Nmax) f[column[n+1][l-1]]=-(kcT*l)/(2.0*l-1.0)*(1.0+n);
        }
        
        if(l<Lmax)
        {
            if(n>0) f[column[n-1][l+1]]=kcT*(1.0+l)/(2.0*l+3.0);
            f[column[n][l+1]]=-kcT*(1.0+l)/(2.0*l+3.0);
        }
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
// Output current Lepton Scalar and Phi perturbations file.
//==================================================================================================

void output_current_pert_k(ofstream &ofile, const ODE_solver_Solution &Sz)
{

    
    ofile << Sz.z << " "<< Sz.y[pExtInd->deltaP] << " " << Sz.y[pExtInd->thetaP] << " ";
    ofile << Sz.y[pExtInd->deltaR] << " " << Sz.y[pExtInd->thetaR]<< " ";
    ofile << Sz.y[pExtInd->Phi];
    
    ofile << endl;
    
    return;
}


//==================================================================================================
// Output current DM perturbations to file
//==================================================================================================


void output_growth_function(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    double kTilde=glob_Par->kTilde;
    double a = Sz.z;
    double rhoTR=pow(a, -4)*Sz.y[pExtInd->BkgR], rhoTP=pow(a, -3)*Sz.y[pExtInd->BkgP], rhoTM=rhoTP=pow(a, -3)*Sz.y[pExtInd->BkgM];
    double TL=TFromRhoT(rhoTR, glob_Par);
    double E=sqrt(rhoTM+rhoTP+rhoTR);
    double thetaTDM=1.5*kTilde*Sz.y[column[0][1]]*(decoupled ? sqrt(2.0*TL/glob_Par->m_DM): glob_Par->c*sqrt(TL/glob_Par->T_kdS));
    double Three_CalH_u= 3*thetaTDM*E*a/pow(kTilde,2);
    
//    cout<<"TCH="<<Three_CalH_u<<" TDM="<<thetaTDM<<" E="<<E<<endl;
//    cout<<"a="<<a<<" nu="<<Sz.y[column[0][0]]+Three_CalH_u<<endl;
    
    ofile <<glob_Par->kTokRH1<<" "<< Sz.y[column[0][0]]+Three_CalH_u<<" "<<Sz.y[column[0][0]]<<" "<< a;
    ofile << endl;
    
    return;
}

void output_current_DMpert_k(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    
    double rhoTR=pow(Sz.z, -4)*Sz.y[pExtInd->BkgR];
    double TL=TFromRhoT(rhoTR, glob_Par);
    
    ofile << Sz.z << " "<< Sz.y[column[0][0]] << " "; //THIS IS IT! OUTPUTS X AND PERTURBATION
    ofile << 1.5*glob_Par->kTilde*Sz.y[column[0][1]]*(decoupled ? sqrt(2.0*TL/glob_Par->m_DM): glob_Par->c*sqrt(TL/glob_Par->T_kdS));

//    for (int l=0; l<Lmax+1; l++) {
//        for (int n=0; n<Nmax+1; n++) {
//            ofile <<Sz.y[column[n][l]]<< " ";
//        }
//    }
    

    ofile << endl;
    
    return;
}



void output_current_DMFullPert_k(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    ofile << Sz.z << " ";
    
    
        for (int l=0; l<Lmax+1; l++) {
            for (int n=0; n<Nmax+1; n++) {
                ofile <<Sz.y[column[n][l]]<<" ";
            }
        }
    
    ofile << endl;
    
    return;
}



//==================================================================================================
// Output current DM perturbations to screen
//==================================================================================================


void cout_current_DMpert_k(const ODE_solver_Solution &Sz)
{
    

    cout<<"a="<< Sz.z;
    for (int l =0; l<min(Lmax+1, 5); l++) {
        for (int n=0; n<min(Nmax+1, 5); n++) {
            cout<<"  f_"<<n<<"_"<<l<<"="<<Sz.y[column[n][l]]<<" ";
        }
    }
    
    cout<<endl;
    return;
}



//==================================================================================================
// Output current Bakcground solution to file
//==================================================================================================

void output_current_BkgSol(ofstream &ofile, const ODE_solver_Solution &Sz)
{
    double rhoTP, rhoTR, rhoTM, T_DM;
    
//    rhoTP  = pow(Sz.z, -3)*Sz.y[pExtInd->BkgP];
//    rhoTR  = pow(Sz.z, -4)*Sz.y[pExtInd->BkgR];
//    rhoTM = pow(Sz.z, -3)*Sz.y[pExtInd->BkgM];
//    T_DM =  pow(Sz.z, -2)*Sz.y[pExtInd->BkgT];
    
    rhoTP  = Sz.y[pExtInd->BkgP];
    rhoTR  = Sz.y[pExtInd->BkgR];
    rhoTM = Sz.y[pExtInd->BkgM];
    T_DM =  Sz.y[pExtInd->BkgT];
    
    ofile << Sz.z << " "<<rhoTP<<" "<<rhoTR<<" "<<rhoTM<<" "<<T_DM<<" ";
    
    ofile << TFromRhoT(pow(Sz.z, -4)*Sz.y[pExtInd->BkgR], glob_Par) << " ";
    
    ofile << endl;
    
    return;

    
}



//==================================================================================================
// setup for k solution
//==================================================================================================
double setup_kevol(struct Background *bkgdens, void * parameters, struct Output_Files *pFiles, struct Growth_Function_File *pGFFile)
{
    //==============================================================
    
    if (decoupled and ((Nmax!=0 or Lmax!=1))){
        
        cout<< "You are attempting to run decoupled mode with Nmax>0 and Lmax>1. Please set them to 0 and 1 respectivelly in common.h file"<<endl;
        
        throw exception();
        
        return _FAILURE_;
    }
    
    else {
        
        Params  P = *(Params *)parameters;
        
        int nsteps=3000;//500;
        long neq=DIM+9;// +9 from the Bckg (4), Lepton and Scalar Perturbations (4) and Grav Potential (1).

        vector<double> xvals(nsteps);
        
        double xs=a_I, xend;
        
        //==============================================================
        // initialize a axis
        //==============================================================
        // we will  need to evolve system to a larger than 1e8 for this configuration, but with numjac=1 it becoms slow at large a.
        init_xarr(a_I, P.a_fin, &xvals[0], nsteps, 1, 1);
        
        //==============================================================
        // solution setup
        //==============================================================
        ODE_solver_Solution Sx;
        ODE_solver_accuracies tols;
        set_verbosity(-1);
        
        Sx.y .resize(neq, 0.0);
        Sx.dy.resize(neq, 0.0);
        
        //==============================================================
        // relative error setup (one should play with the settings)
        //==============================================================
        double rel_eps=5.0e-4;
        tols.rel.resize(neq, rel_eps);
        //tols.rel.resize(neq, 1.0e-2);
        
        tols.rel[0]=tols.rel[1]=tols.rel[2]=rel_eps;

        tols.rel[pExtInd->BkgP]=tols.rel[pExtInd->BkgM]=rel_eps;
        tols.rel[pExtInd->BkgR]=tols.rel[pExtInd->BkgT]=rel_eps;
        
        tols.rel[pExtInd->deltaP]=tols.rel[pExtInd->thetaP]=rel_eps;
        tols.rel[pExtInd->thetaR]=tols.rel[pExtInd->deltaR]=rel_eps;
        tols.rel[pExtInd->Phi]=rel_eps;
        
        //==============================================================
        // absolute error setup (one should play with the settings)
        //==============================================================
        tols.abs.resize(neq, 1.0e-6);
        
        tols.abs[pExtInd->BkgP] = 1.0e-30;//1.0e-25;// the code is very sensitive to this value. Need to figure a way to 'automate' adjustment based on other params
        tols.abs[pExtInd->BkgR] = 1.0e-10;//1.0e-3;
        tols.abs[pExtInd->BkgM] = 1.0e-15;
        tols.abs[pExtInd->BkgT] = 1.0e-3;
        
        tols.abs[pExtInd->deltaP] = 1.0e-3;
        tols.abs[pExtInd->thetaP] = 1.0e-3;
        tols.abs[pExtInd->deltaR] = 1.0e-3;
        tols.abs[pExtInd->thetaR] = 1.0e-3;
        tols.abs[pExtInd->Phi]    = 1.0e-20;

        //==============================================================
        // initial Conditions Background
        //==============================================================
        
        Sx.z=a_I;
        
        Sx.y[pExtInd->BkgP]= pow(a_I, 3)*(bkgdens->rhophi);
        Sx.y[pExtInd->BkgR]= pow(a_I, 4)*(bkgdens->rhor);
        Sx.y[pExtInd->BkgM]= pow(a_I, 3)*(bkgdens->rhom);
        Sx.y[pExtInd->BkgT]= pow(a_I, 2)*bkgdens->TDM;//This handles both cpld and decoupled cases
        
        //========================================================================
        // initial Conditions Lepton Scalar and Phi perturbations.Eqns 12 E&S 2011
        //========================================================================
        
        double PhiIni=1.0;
        
        Sx.y[pExtInd->Phi]= PhiIni;
        Sx.y[pExtInd->deltaR] = PhiIni+ 46.0/63.0*pow(glob_Par->kTilde,2)*PhiIni*a_I;
        Sx.y[pExtInd->thetaR] = -2.0/3.0*pow(glob_Par->kTilde, 2)*PhiIni*sqrt(a_I);
        Sx.y[pExtInd->deltaP] = 2*PhiIni+2.0/3.0*pow(glob_Par->kTilde, 2)*PhiIni*a_I;
        Sx.y[pExtInd->thetaP] = -2.0/3.0*pow(glob_Par->kTilde, 2)*PhiIni*sqrt(a_I);
        
        
        //=======================================================================================
        // initial Conditions DM. Adiabatic to set delta_m and Tight Coupling for the other fnls
        //=======================================================================================
        
        double Inv_Ti=1.0/TFromRhoT(bkgdens->rhor, &P);
        
        create_nl_tocol_map();

        Sx.y[column[0][0]] = 2*Sx.y[pExtInd->deltaR];
        
        Sx.y[column[0][1]] = 2.0*Sx.y[pExtInd->thetaR]/(3.0*glob_Par->kTilde)*(decoupled? sqrt(P.m_DM*Inv_Ti/2.0): sqrt(P.T_kdS*Inv_Ti)/P.c);
        
        if (Nmax>0)
        {
            Sx.y[column[1][0]] = -1.0/8.0*Sx.y[column[0][0]];
        }
        
        
        
        //==============================================================
        
        
        ODE_Solver_data ODE_Solver_info;
        ODE_Solver_info.Snewptr=NULL;
        ODE_Solver_set_up_solution_and_memory(Sx, tols, ODE_Solver_info,
                                              fcn_k, numjac);


        
        for(int i=1; i<nsteps; i++)
        {
	      cout << i << "\n";            
            
            //output_current_solution
            
//            time_t  t0, t1; /* time_t is defined on <time.h> and <sys/types.h> as long */
//            clock_t c0, c1; /* clock_t is defined on <time.h> and <sys/types.h> as int */
//            
//            
//            t0 = time(NULL);
//            c0 = clock();
            

            output_current_BkgSol(pFiles->FileBkg, Sx);
            
            output_current_pert_k(pFiles->FilePert, Sx);
            
            output_current_DMpert_k(pFiles->FilePertDM, Sx);
//            output_current_DMFullPert_k(pFiles->FileFullPertDM, Sx);
            
//            cout_current_DMpert_k(Sx);
            

            xend=xvals[i];
         
            
            // advance solution
            
            ODE_Solver_Solve_history(xs, xend,
                                     1.0e-5*fabs(xs-xend), // minimal time-step
                                     5.0e-1*fabs(xs-xend), // maximal time-step. ADJUST Make smaller. was 0.5*()
                                     Sx, ODE_Solver_info);

//            t1 = time(NULL);
//            c1 = clock();
//            float WT=(float)(t1-t0);
//            float CPUT=(float)(c1 - c0)/CLOCKS_PER_SEC;
//            cout<<""<<endl;
//            cout<<"Moments computed for a="<<xs<<" in CPU Time:"<<CPUT<<" Wall Time:"<<WT<<endl;
//            printf ("\telapsed wall clock time: %ld\n", (long) (t1 - t0));
//            //            printf ("\telapsed CPU time:        %f\n", (float) (c1 - c0)/CLOCKS_PER_SEC);
            
            xs=xend;
        }
        
        
        output_growth_function(pGFFile->GFFile, Sx);
        
        clear_memory(Sx,tols);

       
//        pFiles->FileBkg.close();
//        pFiles->FilePert.close();
//        pFiles->FilePertDM.close();
        
        //==============================================================
        return 0;
    }
}










