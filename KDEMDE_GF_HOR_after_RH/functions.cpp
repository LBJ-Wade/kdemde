//
//  functions.cpp
//  KinDecEMDE
//
//  Created by Cosmin Ilie on 4/13/16.
//  Copyright (c) 2016 c. All rights reserved.
//

#include "functions.h"
#include "common.h"
#include "Definitions.h"

using namespace std;


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




double gstarS(double x, void *params){
    (void) (params);
    return 10.75; //assmed constant
}

double gstar(double x, void *params){
    
    (void) (params);
    return 10.75; //assumed constant
}


double dlngstar_dlnT(double x, void *params){
//Placeholder for now. To code as a function of rhoTr for changing g* case
    
    (void)(params);
    
    return 0;

}

double Fgstar(double a, void *params){
    //This combination appears often: (1+1/4dlng*/dlnT)^-1. Coded here
    
    return (pow(1+0.25*dlngstar_dlnT(a, params),-1));
    
}

double T_kd(void *parameters){
    Params par = *(Params *)parameters;

    return sqrt(5.0/2.0)*(pow(par.T_kdS,2)/par.T_REH);//Eqn 23 E15 assuming gstar is a constant. Later put root finder. ALSO ASSUMES n=2. Pwave
}


double GammaTilde(void *parameters){
    
    Params par = *(Params *)parameters;
    double a_REH=par.a_REH;
    return pow(a_REH/a_I,-3.0/2.0);
    
    //return pow(a_HOR/a_I, -3.0/2.0)*pow(par.kTokRH1, -3);//THIS FORM ASSUMES HOR entry prior to REHEATING//

    
}


//double a_REH(void *parameters){
//    
//    Params par = *(Params *)parameters;
//    
//    return a_HOR*pow(par.kTokRH1, 2);//THIS FORM ASSUMES HOR entry prior to REHEATING//
//}

double gammakdS(void *parameters){
    
    Params par = *(Params *)parameters;
    
    if (decoupled)
        return 0.0;
    else
    return sqrt(8*PI3/90.0/mPL2*gstar(par.T_kdS,&par)*pow(par.T_kdS, 4));

    
}


double a_Today(void *parameters){
    
    Params par = *(Params *)parameters;
    
    return par.a_REH*pow(1.087,0.25)*pow(gstarS(par.T_REH,&par)/gstarS(T_today,&par), 1.0/3.0)*(par.T_REH/T_today);
//    return a_REH(&par)*pow(1.087,0.25)*pow(gstarS(par.T_REH,&par)/gstarS(T_today,&par), 1.0/3.0)*(par.T_REH/T_today);
}


double Gamma(void *parameters){
    Params par = *(Params *)parameters;
    
    return sqrt(8*PI3*gstar(par.T_REH,&par)/90)*pow(par.T_REH, 2)/mPL;
}

double H1(void *parameters){
    Params par = *(Params *)parameters;
    
    return Gamma(&par)/GammaTilde(&par);//THIS FORM ASSUMES HOR entry prior to REHEATING THROUGH GammaTilde//


}

double rho_critI(void *parameters){
    
    Params par = *(Params *)parameters;
    return 3*pow(H1(&par)*mPL, 2)/2/FOURPI;
    
}

double k_RH(void *parameters){
    
    Params par = *(Params *)parameters;
    
    return pow(gstarS(T_today,&par)/gstarS(par.T_REH,&par), 1.0/3.0)*T_today*par.T_REH/mPL*sqrt(8.0*PI3/90.0*gstar(par.T_REH,&par));
}


double cm(void *parameters){
    
    //rho_m=cm/a^3
    
    Params par = *(Params *)parameters;
    
    return OmegaM*rho_C*pow(a_Today(&par),3);
}


double TFromRhoT(double x, void *parameters){
    // first argument x is the radiation density in units of critical density initially, i.e. rhoTilde
    //For now coded assuming a constant gstar!
    
    Params par = *(Params *)parameters;

    return pow(x*rho_critI(&par)*30/PI2/gstar(1.0,&par), 0.25);
    
    
}


double dT_drhoT(double r, void *parameters){
    Params par = *(Params *) parameters;
    
    return 0.25*TFromRhoT(r, &par)/r*Fgstar(r, &par);
}


double H(double a, double b, double c, void *parameters){
//Hubble rate in terms of rho_tilde quantities a, b, c for three component fluid
    
    Params par = *(Params *)parameters;

    return H1(&par)*sqrt(a+b+c);

}

double E(double a, double b, double c){
    //Hubble rate in terms of rho_tilde quantities a, b, c for three component fluid
    
    
    return sqrt(a+b+c);
    
}

double Upsi(double p, double m, double r, void *parameters){
    
    Params par = *(Params *)parameters;
    
    return (decoupled? 0.0:(gammakdS(&par)/H(p, m, r, &par))*pow(TFromRhoT(r, &par)/par.T_kdS,4+n));
    
}

double UpsiPar(void *bkgdens, void *parameters){
    
    Background bkg = *(Background *)bkgdens;
    Params par = *(Params *)parameters;
    
    double rhophi= bkg.rhophi;
    double rhom=bkg.rhom;
    double rhor=bkg.rhor;
    
    return (gammakdS(&par)/H(rhophi, rhom, rhor, &par))*pow(TFromRhoT(rhor, &par)/par.T_kdS,4+n);
    
}


double gammaa(double r, double x, void *parameters){
    
     Params par = *(Params *)parameters;
    
    return gammakdS(&par)*pow(TFromRhoT(r, &par)/par.T_kdS, n+4)*pow(10, x);
    
}


double Rcr(double p, double m, double r,double x, void *parameters){
    
    //This is a more general form, coded withoug assuming g* is constant
    Params par = *(Params *)parameters;
    
    double F=(1+1.0/2.0*dT_drhoT(r, &par)/TFromRhoT(r, &par)*
              (-4*r+GammaTilde(&par)*p/sqrt(p+m+r)));
    
    return H(p, m, r, &par)*pow(10, x)*F;
}


double R(double p, double m, double r,double a, void *parameters){
    
    //This is coded  for a general g* //
    
    Params par = *(Params *)parameters;
    

    double F=(1.0 +Fgstar(a, &par)*(-0.5+GammaTilde(&par)/(8*sqrt(p+m+r))*(p/r)));

    
    return H(p, m, r, &par)*a*F;
}




int ICBkgRho_TDM(struct Background *bkg_ivals, void *parameters){

    //This will generate the Initial values for the background densities in units of rho_crit,I and store them
    //Can be adapted to include all cases: Rescaled, lna, a scaling
    
    Params par = *(Params *)parameters;
    
    double rhophiI;
    double rhomI;
    double rhorI;
    double RatioIni;//The initial ratio of Dm to Lepton Temperature. See V&G paper Eq 44.
    double UpsiIni;
    double si;
    
    RatioIni=1.0; //placeholder
    
    rhomI=cm(&par)/rho_critI(&par);
    rhophiI=(1-rhomI)/(2.0/5.0*GammaTilde(&par)+1);
    rhorI=2.0/5.0*GammaTilde(&par)*rhophiI;
    
    UpsiIni=Upsi(rhophiI, rhomI, rhorI, &par);
    si=2*UpsiIni/(alpha*beta);
    cout<<"si="<<si<<endl;
    RatioIni=(UpsiIni>100? 1 : exp(si)*pow(si, lambda)*GAMMA(1-lambda,si));
    
    cout<<"RatioIni="<<RatioIni<<endl;
    
    bkg_ivals->rhom=rhomI;
    bkg_ivals->rhophi=rhophiI;
    bkg_ivals->rhor=rhorI;
    bkg_ivals->TL =TFromRhoT(rhomI, &par);
    bkg_ivals->TDM= (decoupled? TFromRhoT(rhorI, &par)*pow(a_I/astar, alpha-2.0):TFromRhoT(rhorI, &par)*RatioIni);//THIS INCLUDES CASE FOR PARTIAL COUPLED a_KD<1!!!! Test for decoupled to set it to 0 and also remove the ODE for TDM, as in that limit TDM is zero.

    
    
    return _SUCCESS_;
    
}


int Initialize_Parameters(struct Params *parameters){
//  Placeholder. Add all Derrived parameters //

    
//    double xfin=10.0;
//    double c=0.01; // c is defined as sqrt(2T_d/m_DM)
//    double kTokRH1=100;
//    double T_REH=5;
//    double T_kdS=19.94079648; // a_KD=100 if TRH=5 for this T_kdS
//    
//
//    parameters->xfin=xfin;
//    parameters->c=c;
//    parameters->kTokRH1=kTokRH1;
//    parameters->T_REH=T_REH;
//    parameters->T_kdS=T_kdS;
//    parameters->a_REH=a_REH(parameters);
//    
//    printf("a_reh=%f",a_REH(parameters));
    
    
    return _SUCCESS_;

}


double a_KD(struct Background *bkg_ivals, void *parameters){
    //======TO CHK if g* const assumption was used!======//
   Params par = *(Params *)parameters;
    cout<<"gammakds="<<gammakdS(&par)<<" n="<<n<<" H1="<<H1(&par)<<" Tini="<<TFromRhoT(bkg_ivals->rhor, &par)<<" Tkds="<<par.T_kdS<<endl;
    return a_I*pow(gammakdS(&par)/H1(&par),8.0/3.0/n)*pow(TFromRhoT(bkg_ivals->rhor, &par)/par.T_kdS,8.0*(4+n)/3.0/n);
    
}


void create_index_map()
{
    var_nl dum;
    indices.clear();
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


void create_nl_tocol_map(){
    for(int l=0; l<Lmax+1; l++)
    {
                for(int n=0; n<Nmax+1; n++)
        {
            column[n][l]=Shift+n+l*(Nmax+1);//
        }
    }
    
}

void create_output_files(struct Output_Files *ptrFiles){
    
    cout << " Creating Output Files for kTokRH1= " << glob_Par->kTokRH1 <<"c="<<glob_Par->c<<endl;
    double aRH=glob_Par->a_REH, c=glob_Par->c, T_RH=glob_Par->T_REH, Tkds=glob_Par->T_kdS, kTokRH1=glob_Par->kTokRH1, mDM=glob_Par->m_DM;

    
    char buffer_kTokRH1[10] , pFT_Bkg[100], pFT_DM_Pert[100], pFT_Pert[100], pFT_FULL_DM_Pert[100];
    sprintf(buffer_kTokRH1, "%.0E",glob_Par->kTokRH1);
    
    sprintf(pFT_Bkg, "./DATAFILES/B_Bkg_CPLD_Nmax%dLmax%daRH%0.ETRH%.0ETkds%0.EkTokRH1_%0.E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1);

    if (decoupled)
    {
        
        sprintf(pFT_Pert, "./DATAFILES/B_DcpldPert_Nmax%dLmax%daRH%0.ETRH%0.ETkds%0.EkTokRH1_%0.EmDM%.0E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,mDM);
        sprintf(pFT_DM_Pert,"./DATAFILES/B_DcpldPertDM_Nmax%dLmax%daRH%0.ETRH%0.ETkds%0.EkTokRH1_%0.EmDM%.0E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,mDM);
        sprintf(pFT_FULL_DM_Pert,"./DATAFILES/B_DcpldPertDM_FULL_Nmax%dLmax%daRH%0.ETRH%0.ETkds%0.EkTokRH1_%0.EmDM%.0E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,mDM);


    }
    
    else
    {
        
        sprintf(pFT_Pert, "./DATAFILES/B_Pert_Nmax%dLmax%daRH%0.ETRH%0.ETkds%0.EkTokRH1_%0.Ec%.0E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,c);
        sprintf(pFT_DM_Pert,"./DATAFILES/B_PertDM_Nmax%dLmax%daRH%2.2ETRH%1.ETkds%0.EkTokRH1_%2.2Ec%2.2E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,c);
        sprintf(pFT_FULL_DM_Pert,"./DATAFILES/B_PertDM_FULL_Nmax%dLmax%daRH%0.ETRH%.0ETkds%0.EkTokRH1_%0.Ec%.0E.dat",Nmax,Lmax,aRH,T_RH,Tkds,kTokRH1,c);

    }
    

    ptrFiles->FilePert.open(pFT_Pert);
    ptrFiles->FilePert.precision(6);
    ptrFiles->FilePert << " # c=" << c << "  kTokRH1="<<kTokRH1<<" Tkds="<<Tkds<<" T_RH="<<T_RH<<" aRH="<<aRH<<endl;
    ptrFiles->FilePert << " # x" << " delta_p " << " theta_p " << " delta_r " << " theta_r " << " Phi " << endl;
    
    
    ptrFiles->FilePertDM.open(pFT_DM_Pert);
    ptrFiles->FilePertDM.precision(6);
    ptrFiles->FilePertDM << " # c=" << c << "  kTokRH1="<<kTokRH1<<" Tkds="<<Tkds<<" T_RH="<<T_RH<<" aRH="<<aRH<<endl;
    ptrFiles->FilePertDM << " # x" <<" delta_m "<< "theta_m "<< endl;
    
    ptrFiles->FileFullPertDM.open(pFT_FULL_DM_Pert);
    ptrFiles->FileFullPertDM.precision(6);
    ptrFiles->FileFullPertDM << " # c=" << c << "  kTokRH1="<<kTokRH1<<" Tkds="<<Tkds<<" T_RH="<<T_RH<<endl;

    ptrFiles->FileFullPertDM << " # x" <<" f_nl with n running first "<< endl;
    ptrFiles->FileFullPertDM << " # x";
    for (int l=0; l<Lmax+1; l++) {
        for (int n=0; n<Nmax+1; n++) {
            ptrFiles->FileFullPertDM <<"f_"<<n<<"_"<<l<<" ";
        }
    }
    ptrFiles->FileFullPertDM << endl;
    

    ptrFiles->FileBkg.open(pFT_Bkg);
    ptrFiles->FileBkg.precision(6);
    ptrFiles->FileBkg << " # c=" << c << "  kTokRH1="<<kTokRH1<<" Tkds="<<Tkds<<" T_RH="<<T_RH<<endl;
    ptrFiles->FileBkg << " # x" << " rhoT_p " << " rhoT_r " << " rhoT_m " << " T_DM " << " T_L " << endl;
 
    

}

void create_output_file_growth_function(struct Growth_Function_File *ptrFile)
{
    
    cout << " Creating  Growth Function File for = " <<"c="<<glob_Par->c<<endl;
    double c=glob_Par->c, T_RH=glob_Par->T_REH, Tkds=glob_Par->T_kdS, kTokRH1=glob_Par->kTokRH1, mDM=glob_Par->m_DM;
    
    
    char buffer_kTokRH1[10]  ,pGF_DM[100];
    sprintf(buffer_kTokRH1, "%.0E",glob_Par->kTokRH1);
    
    
    if (decoupled)
    {
        sprintf(pGF_DM, "./DATAFILES/DcpldPert_GF_2_DM_Nmax%dLmax%dTRH%0.ETkds%0.EkTokRH1_%0.EmDM%.0E.dat",Nmax,Lmax,T_RH,Tkds,kTokRH1,mDM);
    }
    
    else
    {
        sprintf(pGF_DM,"./DATAFILES/PertDM_GF_2_Nmax%dLmax%dTRH%.0ETkds%0.Ec%.0E.dat",Nmax,Lmax,T_RH,Tkds,c);
    }
    
    
    
    
    ptrFile->GFFile.open(pGF_DM);
    ptrFile->GFFile.precision(6);
    ptrFile->GFFile << " # c=" << c <<" Tkds="<<Tkds<<" T_RH="<<T_RH<<endl;
    ptrFile->GFFile <<"# kTokRH1"<<" nu "<<" delta_DM "<<" a "<<endl;
    
    
    
}


void close_output_files(struct Output_Files *ptrFiles)
{
    ptrFiles->FilePert.close();
    ptrFiles->FilePertDM.close();
    ptrFiles->FileBkg.close();
    ptrFiles->FileFullPertDM.close();
}

void close_GF_file(struct Growth_Function_File *ptrFile)
{
    ptrFile->GFFile.close();
}

//void create_Extra_indices_map(){
//
//    ExtInd.BkgP=0;
//    ExtInd.Bkgr=1;
//    ExtInd.Bkgm=2;
//    ExtInd.BkgT=3;
//    
//}

