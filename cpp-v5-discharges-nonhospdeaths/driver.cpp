#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "prms.h"
#include "essentials.h"
#include "parseargs.h"

#include <ctime>
#include <sys/time.h>

using namespace std;


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double b,d,s,c0,c2,c3;	// this are defined, i believe, in rkf.cpp
double* yic;

FILE* OutFile = NULL;

prms* ppc;  // need to allocate space for this in main

double G_CLO_BETA = 0.3;
double G_CLO_INTRODUCTION_TIME = -1;
int G_CLO_INTRODUCTION_COUNT;
string G_CLO_LOCATION = "RI";

double G_CLO_TF = 365.0;
int G_CLO_STEPS_PER_DAY=100;
double G_CLO_P_HOSP_TO_ICU = 0.30;

double G_CLO_SYMP_FRAC = 0.25;
double G_CLO_HOSPFRAC_YOUNG_DEV = 1.0;      //DEPRECATED
double G_CLO_HOSPFRAC_MID_DEV = 0.5;        //DEPRECATED
double G_CLO_HOSPFRAC_OLD_DEV = 0.5;        //DEPRECATED

double G_CLO_HOSPFRAC_10 = 0.021; // numbers below are the raw RI data (June 8)
double G_CLO_HOSPFRAC_20 = 0.017;
double G_CLO_HOSPFRAC_30 = 0.031;
double G_CLO_HOSPFRAC_40 = 0.05;
double G_CLO_HOSPFRAC_50 = 0.1;         //0.125;
double G_CLO_HOSPFRAC_60 = 0.15;         //0.222;
double G_CLO_HOSPFRAC_70 = 0.15;         //0.307;
double G_CLO_HOSPFRAC_80 = 0.15;         //0.198;

double G_CLO_DEATHPROB_HOME_60 = 0.00;
double G_CLO_DEATHPROB_HOME_70 = 0.05;
double G_CLO_DEATHPROB_HOME_80 = 0.22;

double G_CLO_RELSUSC_0 = 0.6;
double G_CLO_RELSUSC_10 = 0.6;
double G_CLO_RELSUSC_20 = 1.0;
double G_CLO_RELSUSC_30 = 1.0;
double G_CLO_RELSUSC_40 = 1.0;
double G_CLO_RELSUSC_50 = 1.0;
double G_CLO_RELSUSC_60 = 1.0;
double G_CLO_RELSUSC_70 = 1.0;
double G_CLO_RELSUSC_80 = 1.0;

double G_CLO_TIME_SYMPTOHOSP = 7.0;
double G_CLO_SELFISOLATION_FACTOR = 1.0;

double G_CLO_EARLYMARCH_HOSPRATE = 1.0;     // scaling factor showing how much more likely hospitalization was in early March than later in the epidemic
double G_CLO_EARLYMARCH_ENDDAY = 84.5;      // the day that the high-hosp rate ended (probably because the patient load became too high)


double G_CLO_DEV_LEN_HOSPSTAY = 1.0;

double G_CLO_ICUFRAC_DEV = 1.0;
double G_CLO_PROB_ICU_TO_VENT = 0.75;

double G_CLO_VENTDEATH_MID_DEV = 0.7;
double G_CLO_MEANTIME_ON_VENT_SURV = 10.8;   //  mean time on a ventilator for survivors

double G_CLO_RELATIVE_BETA_HOSP = 0.2;

bool G_B_DIAGNOSTIC_MODE = false;
bool G_B_CHECKPOP_MODE = false;
bool G_B_USESOCIALCONTACTMATRIX = true;
bool G_B_REMOVE_99COLUMNS = false;

double G_C_COMM[NUMAC][NUMAC]; // this is the contact rate in the community between susc (first index) and infected (second index)                              - this matrix is symmetric
double G_C_HOSP[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected hospitalized patient (second index)     - this matrix is not symmetric
double G_C_ICU[NUMAC][NUMAC];  // this is the contact rate btw the community (first index, susceptible) and an infected ICU patient (second index)              - this matrix is not symmetric
double G_C_VENT[NUMAC][NUMAC]; // this is the contact rate btw the community (first index, susceptible) and an infected Ventilated patient (second index)       - this matrix is not symmetric
 

// END ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void InitializeContactMatrices( void );


//
//
int main(int argc, char* argv[])
{
    

    //
    // ###  1.  ALLOCATE SPACE FOR A PARAMETERS CLASS AND FOR THE MAIN STATE VARIABLES
    //
    ppc = new prms; assert( ppc );
    yic = new double[STARTK+NUMAC];
    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out


    // assign some default values to the v_beta and v_betatimes arrays
    ppc->v_betas.push_back( 1.0 );
    ppc->v_betas.push_back( 0.9 );
    ppc->v_betas.push_back( 0.8 );
    ppc->v_betatimes.push_back( 0.0 );
    ppc->v_betatimes.push_back( 60.0 );
    ppc->v_betatimes.push_back( 95.0 );
    
//     string s1(argv[0]); string s2(argv[1]); string s3("3.8.44"); 
//     
//     printf("\n\t %s float outcome %d", s1.c_str() , isFloat(s1)?1:0 );
//     printf("\n\t %s float outcome %d", s2.c_str() , isFloat(s2)?1:0 );
//     printf("\n\t %s float outcome %d", s3.c_str() , isFloat(s3)?1:0 );
//     printf("\n\n");
    
    
    // this just sets all the matrix values to one
    InitializeContactMatrices();
    
    ParseArgs( argc, argv );
    
    

    
    
    /*for(int i=0; i < ppc->v_betas.size(); i++)
    {
        printf("\n\t time-start: %1.3f   \t beta-val: %1.3f", ppc->v_betatimes[i], ppc->v_betas[i] );
    }
    printf("\n\t time final is %1.1f\n", G_CLO_TF);
    
   return 1;*/
    
    
    
    //
    // ###  2.  INITIALIZE PARAMETERS - these are the default/starting values
    //
    ppc->v[ i_len_incub_period ]                            = 6.0;  // this is an average of the Lauer et al estimate (5.5d) and Backer et al (6.5d)
    ppc->v[ i_len_symptomatic_infectious_period_phase_1 ]   = G_CLO_TIME_SYMPTOHOSP;
    ppc->v[ i_len_symptomatic_infectious_period_phase_2 ]   = 7.0;
    ppc->v[ i_len_medicalfloor_hospital_stay ]              = 10.7 * G_CLO_DEV_LEN_HOSPSTAY; //   take from Lewnard et al, survivors only
    ppc->v[ i_mean_time_vent ]                              = G_CLO_MEANTIME_ON_VENT_SURV;   //   mean time on ventilator for survivors
    
    
    // params below are for relative infectiousness of certain individuals; I1 and I2 individuals have infectiousness = 1.0
    ppc->v[ i_phi_asymp ]           = 0.5;  // set to same value as Imperial college models, and JosephWu NatMed paper assumption (sensitivity analysis performed)
    ppc->v[ i_phi_symp_phase2 ]     = 0.75; // the second half of the symptomatic period is somewhat less infectious (TODO get evidence here)
    ppc->v[ i_phi_incub ]           = 0.5;  // currently a complete unknown
    ppc->v[ i_phi_hosp ]            = 1.0;
    ppc->v[ i_phi_hosp_recovering ] = 0.2;  // likely to be at late stage of infection with much less culturable virus
    ppc->v[ i_phi_icu ]             = 1.0;
    ppc->v[ i_phi_vent ]            = 0.1;  // to discuss
    
    
    ppc->v[ i_selfisolation_factor ] = G_CLO_SELFISOLATION_FACTOR;  // set to 1.0 by default; this means that infected and symptomatic individuals in the community 
                                                                    // do not self-isolate any more than anyone else; set this to something smaller than one
                                                                    // to make these individuals isolate more than non-symp people
    
    
    // params below are for relative contact levels of hospitalized, ICUed, and vented individuals
    ppc->v[ i_beta_hosp ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // these relbeta's do not change with social distancing measures
    ppc->v[ i_beta_icu ]  = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // because the contact rate of a hospitalized patient is unaffected by the 
    ppc->v[ i_beta_vent ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // social distancing policies set for healthy individuals
    
    
    
    // set the relative susceptibilities of the different age groups
    for(int ac=0;ac<NUMAC;ac++) 
    {
        ppc->v_rel_susc[0] = G_CLO_RELSUSC_0;
        ppc->v_rel_susc[1] = G_CLO_RELSUSC_10;
        ppc->v_rel_susc[2] = G_CLO_RELSUSC_20;
        ppc->v_rel_susc[3] = G_CLO_RELSUSC_30;
        ppc->v_rel_susc[4] = G_CLO_RELSUSC_40;
        ppc->v_rel_susc[5] = G_CLO_RELSUSC_50;
        ppc->v_rel_susc[6] = G_CLO_RELSUSC_60;
        ppc->v_rel_susc[7] = G_CLO_RELSUSC_70;
        ppc->v_rel_susc[8] = G_CLO_RELSUSC_80;
    }

    // set the fraction of individuals who progress to asymptomatic infection
    double a=G_CLO_SYMP_FRAC;  //NOTE this number has to be somewhere between 0.1 and 0.3, closer to 0.3 probably; default is 0.25;
    //
    ppc->v_prob_E_A[0] = 1.0 - a * 0.05; // 98.5% to 99.5% asymp       // these relative params are all taken from Joseph Wu et al, Nat Med, 2020
    ppc->v_prob_E_A[1] = 1.0 - a * 0.08; // 97.6% to 99.2% asymp
    ppc->v_prob_E_A[2] = 1.0 - a * 0.41; // 87.7% to 95.9% asymp
    ppc->v_prob_E_A[3] = 1.0 - a * 1.00; // 70% to 90% asymp           // 30-39, the reference group
    ppc->v_prob_E_A[4] = 1.0 - a * 1.35; // 59.5% to 86.5% asymp
    ppc->v_prob_E_A[5] = 1.0 - a * 1.99; // 40.3% to 80.1% asymp 
    ppc->v_prob_E_A[6] = 1.0 - a * 2.85; // 14.5% to 71.5% asymp 
    ppc->v_prob_E_A[7] = 1.0 - a * 3.05; //  8.5% to 69.5% asymp 
    ppc->v_prob_E_A[8] = 1.0 - a * 2.49; // 25.3% to 75.1% asymp 

    ppc->v_prob_E_A[0] = 0.83;      // these values are taken from a range of studies, and approximated 
    ppc->v_prob_E_A[1] = 0.83;      // their current state (2020-06-09) is for a RI manual fit exercise
    ppc->v_prob_E_A[2] = 0.5; 
    ppc->v_prob_E_A[3] = 0.5; 
    ppc->v_prob_E_A[4] = 0.45; 
    ppc->v_prob_E_A[5] = 0.35; 
    ppc->v_prob_E_A[6] = 0.35; 
    ppc->v_prob_E_A[7] = 0.25; 
    ppc->v_prob_E_A[8] = 0.2; 

    // set the fraction of individuals who are hospitalized immediately after I_2
    //
//     double b1 = G_CLO_HOSPFRAC_YOUNG_DEV;   // default is 1.8 :: REASON is that we want these rates to match the hosp-age-dist in the Lewnard paper
//     double b2 = G_CLO_HOSPFRAC_MID_DEV;     // default is 
//     double b3 = G_CLO_HOSPFRAC_OLD_DEV;     // default is 
    ppc->v_prob_I2_H[0] = G_CLO_HOSPFRAC_10;    // NOTE because there is not likely to be data on hosp rates for 0-9 y.o., we adopt the hosp rate for 10-19 year olds here
    ppc->v_prob_I2_H[1] = G_CLO_HOSPFRAC_10;
    ppc->v_prob_I2_H[2] = G_CLO_HOSPFRAC_20;
    ppc->v_prob_I2_H[3] = G_CLO_HOSPFRAC_30;    // this is the hospitalization probability for 30-39 year-olds
    ppc->v_prob_I2_H[4] = G_CLO_HOSPFRAC_40;
    ppc->v_prob_I2_H[5] = G_CLO_HOSPFRAC_50;
    ppc->v_prob_I2_H[6] = G_CLO_HOSPFRAC_60;
    ppc->v_prob_I2_H[7] = G_CLO_HOSPFRAC_70;
    ppc->v_prob_I2_H[8] = G_CLO_HOSPFRAC_80;
    //
    // IMPORTANT-NOTE:  now that you have assigned the hospitalization probabilities, adjust them upwards for the 
    //                  early March period when hosp rates were higher than normal
    ppc->earlymarch_highhosp_factor = G_CLO_EARLYMARCH_HOSPRATE;
    ppc->earlymarch_highhosp_endday = G_CLO_EARLYMARCH_ENDDAY;
    //ppc->apply_earlymarch_hosprates();



    // set the probability of death for I4 for all 9 age classes
    ppc->v_prob_I4_D[0] = 0.0; // COMPLETELY UNKNOWN
    ppc->v_prob_I4_D[1] = 0.0;
    ppc->v_prob_I4_D[2] = 0.0;
    ppc->v_prob_I4_D[3] = 0.0;
    ppc->v_prob_I4_D[4] = 0.0;
    ppc->v_prob_I4_D[5] = 0.0;
    ppc->v_prob_I4_D[6] = G_CLO_DEATHPROB_HOME_60;
    ppc->v_prob_I4_D[7] = G_CLO_DEATHPROB_HOME_70;
    ppc->v_prob_I4_D[8] = G_CLO_DEATHPROB_HOME_80;

    // set the probability of progression from the "HA_1" -state to the CA state; this is the only way to go from hospitalizatio to ICU
    ppc->v_prob_HA_CA[0] = 0.304; 
    ppc->v_prob_HA_CA[1] = 0.293;     // NOTE the estimated Lewnard et al (medRxiv, Apr 16) probabilities are used here
    ppc->v_prob_HA_CA[2] = 0.2825;    //      averaged over male/female equally
    ppc->v_prob_HA_CA[3] = 0.301;
    ppc->v_prob_HA_CA[4] = 0.463;
    ppc->v_prob_HA_CA[5] = 0.4245;
    ppc->v_prob_HA_CA[6] = 0.460;
    ppc->v_prob_HA_CA[7] = 0.4835;
    ppc->v_prob_HA_CA[8] = 0.416;
    
    // the probabilities above range from: 
    double c1 = G_CLO_ICUFRAC_DEV;
    for(int ac=0; ac<NUMAC; ac++)
    {
        ppc->v_prob_HA_CA[ac] *= c1;
        //ppc->v_prob_HA_CA[ac] = 0.0;
    }
    
    // set the probability of death for HA4 for all 9 age classes - 
    ppc->v_prob_HA4_D[0] = 0.0;
    ppc->v_prob_HA4_D[1] = 0.0;
    ppc->v_prob_HA4_D[2] = 0.0;
    ppc->v_prob_HA4_D[3] = 0.0;
    ppc->v_prob_HA4_D[4] = 0.0;
    ppc->v_prob_HA4_D[5] = 0.0;   // NOTE there are no data right now for these numbers
    ppc->v_prob_HA4_D[6] = 0.0;
    ppc->v_prob_HA4_D[7] = 0.03;
    ppc->v_prob_HA4_D[8] = 0.05;

    
    // set the probability of ventilation for CA-individuals for all 9 age classes - 
    ppc->v_prob_CA_V[0] = G_CLO_PROB_ICU_TO_VENT; // these can be set to about .75 for the higher age classes; from the Seattle ICU paper on 24 patients
    ppc->v_prob_CA_V[1] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[2] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[3] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[4] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[5] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[6] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[7] = G_CLO_PROB_ICU_TO_VENT;
    ppc->v_prob_CA_V[8] = G_CLO_PROB_ICU_TO_VENT;

    // from the Seattle ICU data on 24 patients, this probability is 60% -- obviously, it's a small sample size of older patients
    // these are being set to the ICU-to-Death probabilities since it's very difficult to get good data on death when on and not on a ventilator (for ICU patients)
    //
    double vd = G_CLO_VENTDEATH_MID_DEV; // default set to 1.0
    ppc->v_prob_V_D[0] = 0.03125;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[1] = 0.05119;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[2] = 0.15;          // this range should be between 14% (Lewnard) and 16.7%  (Graselli)
    ppc->v_prob_V_D[3] = 0.15;          // this range should be between 13% (Lewnard) and 17%  (Graselli), but Yang LRM observed 0%
    ppc->v_prob_V_D[4] = 0.400*vd;          // this range should be between 31% and 50%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[5] = 0.460*vd;          // this range should be between 26% and 70%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[6] = 0.585*vd;         // this range should be between 39% and 72%  (Lewnard, Graselli, Yang LRM, Bhatraju) 
    ppc->v_prob_V_D[7] = 0.70;          // this range should be between 60% and 88%  (Lewnard, Graselli, Yang LRM, Bhatraju)
    ppc->v_prob_V_D[8] = 0.90;          // this range should be between 60% and 100% (Lewnard, Graselli, Yang LRM, Bhatraju)
    

    // set the probability of death for CA-individuals for all 9 age classes 
    //
    // if you are in the ICU, and you don't progress to death, and your don't progress to ventilation, that means you are back on the medical-floor level of care
    //
    // very little data here; Seattle study on 24 ICU patients says this should be about 0.125
    //
    ppc->v_prob_CA_D[0] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];         
    ppc->v_prob_CA_D[1] = (1.0-ppc->v_prob_CA_V[1]) * ppc->v_prob_V_D[1];                                     
    ppc->v_prob_CA_D[2] = (1.0-ppc->v_prob_CA_V[2]) * ppc->v_prob_V_D[2];
    ppc->v_prob_CA_D[3] = (1.0-ppc->v_prob_CA_V[3]) * ppc->v_prob_V_D[3];
    ppc->v_prob_CA_D[4] = (1.0-ppc->v_prob_CA_V[4]) * ppc->v_prob_V_D[4];
    ppc->v_prob_CA_D[5] = (1.0-ppc->v_prob_CA_V[5]) * ppc->v_prob_V_D[5];
    ppc->v_prob_CA_D[6] = (1.0-ppc->v_prob_CA_V[6]) * ppc->v_prob_V_D[6];  // 0.146 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[7] = (1.0-ppc->v_prob_CA_V[7]) * ppc->v_prob_V_D[7];  // 0.175 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[8] = (1.0-ppc->v_prob_CA_V[8]) * ppc->v_prob_V_D[8];  // 0.225 --- so it's not close to the average of 0.125 obersved in Seattle 24-patent study; but these patients are >80
    
    
    
    
    
    
    //
    // ###  3.  RUN THE MODEL
    //
    generate_trajectories( 0.01, 3.0, 0.0, G_CLO_TF, 0.0 );
    if( OutFile != NULL ) fclose(OutFile);




    
    
    //
    // ###  4.  OUTPUT DIAGNOSTICS
    //
    if( G_B_DIAGNOSTIC_MODE )
    {
        int ac;
        
        printf("\n\t\t\t\t 0-9 \t\t 10-19  \t 20-29  \t 30-39  \t 40-49  \t 50-59  \t 60-69  \t 70-79  \t 80+ ");
        printf("\n\tTotal Symp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTJ + ac]+0.5) );
        }
        printf("\n\tTotal Deaths    ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTD + ac]+0.5) + (int)(yic[STARTDHOSP + ac]+0.5));
        }
        printf("\n\tCFR             ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t %1.3f%%   ", 100.0 * (yic[STARTD + ac]+yic[STARTDHOSP + ac]) / yic[STARTJ + ac]  );
        }

        printf("\n\n\tTotal Hosp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTK + ac]+0.5) );
        }
        printf("\n\tHosp FR         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTDHOSP + ac] / yic[STARTK + ac]  );
        }

        double total_num_hospitalized=0.0;
        for(ac=0; ac<NUMAC; ac++) total_num_hospitalized += yic[STARTK + ac];
        
        printf("\n\n\t%% Hosp by Age         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTK + ac] / total_num_hospitalized  );
        }
        printf("\n\tLewnard et al  \t\t  0.1%% \t\t  0.2%% \t\t  3.6%% \t\t  9.8%% \t\t  14.6%% \t  20.5%% \t  22.2%% \t  16.6%% \t  12.3%%");
        
    
        printf("\n\n");

        /*for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
        {
            for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
            {
                printf("\t %1.2f ", G_C_COMM[acs][aci]);
            }
            printf("\n");
        }
        
        printf("\n\n");*/
        
    }
    

    delete[] yic;
    delete ppc;
    return 0;
}



void InitializeContactMatrices( void )
{
    // ### GENERAL POPULATION CONTACT MATRIX
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            G_C_COMM[acs][aci] = 1.0;
        }
    }

    
    // ### CONTACT MATRIX BTW INFECTED HOSPITALIZED INDIVIDUALS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs==0 )        // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.0; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; HEALTHY CHILDREN DO NOT VISIT THE HOSPITAL
            }
            else if( acs==1 )   // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.0; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; HEALTHY TEENS VISIT THE HOSPITAL RARELY
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_HOSP[acs][aci] = 0.2; // NO HOSPITAL VISITS OF ANY KIND EXCEPT EOL CIRCUMSTANCES; SOME MAY WORK IN THE HOSPITAL
            }
            else
            {
                G_C_HOSP[acs][aci] = 1.0;
            }
        }
    }

    
    // ### CONTACT MATRIX BTW INFECTED ICU PATIENTS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs<=1 )        // age-class of the susceptible contact
            {
                G_C_ICU[acs][aci] = 0.0; // HEALTHY CHILDREN & TEENS DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_ICU[acs][aci] = 0.0; // UNINFECTED ELDERLY DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else
            {
                G_C_ICU[acs][aci] = 1.0;
            }
        }
    }
    
    
    // ### CONTACT MATRIX BTW INFECTED-AND-VENTILATED ICU PATIENTS AND THE REST OF THE POPULATION
    for(int acs=0; acs<NUMAC; acs++) // age-class of the susceptible individual
    {
        for(int aci=0; aci<NUMAC; aci++) // age-class of the infected individual   
        {
            if( acs<=1 )        // age-class of the susceptible contact
            {
                G_C_VENT[acs][aci] = 0.0; // HEALTHY CHILDREN & TEENS DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else if( acs>=7 )   // age-class of the susceptible contact
            {
                G_C_VENT[acs][aci] = 0.0; // UNINFECTED ELDERLY DO NOT VISIT THE ICU, AND DO NOT WORK IN THE ICU
            }
            else
            {
                G_C_VENT[acs][aci] = 1.0;
            }
        }
    }
    
    
    
}


