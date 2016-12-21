Currently the code handles the kTokRH1 <1 and kToKRH1 >1 scenarios with two separate versions, named: KDEMDE_GF_HOR_after_RH and KDEMDE_GF_HOR_before_RH_V2 respectively. 

In common.h one sets Lmax Nmax and a number of flag switches. Additionally for the kTokRH1>1 case a_HOR is set in there as well. 

The code computes the perturbations for a k grid with nk_grid points and outputs those in the folder DATAFILES. Each data file has a header describing it and the name is descriptive, including the values of parameters used in the file name. The files that start with PertDM_GF_ in the title store the values of the DM density perturbation and the gauge invariant quantity nu (defined in Bertschinger 06) for  

For each mode the perturbations are evolved from a_I to a_FIN. 

The parameters for the solver such as accuracies and number of steps for the time evolution are set in the setup_kevol function that is found in the EMDEPerturbations.cpp file. That file contains all the perturbation equations, the initial conditions, as well as the functions that code how to output the results. 

The 'main' file is called KDEMDE.cpp and it is the place where T_kdS, T_RH, c and the kToKRH1 values are defined. 