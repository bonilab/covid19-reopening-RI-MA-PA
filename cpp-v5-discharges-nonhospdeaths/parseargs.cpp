#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "essentials.h"
#include "prms.h"
#include "assert.h"



// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double* yic;  
extern prms* ppc;
extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;
extern FILE* OutFile;
extern double G_CLO_TF;
extern int G_CLO_STEPS_PER_DAY;
extern string G_CLO_LOCATION;
extern double G_CLO_P_HOSP_TO_ICU;
extern bool G_B_DIAGNOSTIC_MODE;
extern bool G_B_CHECKPOP_MODE;
extern bool G_B_USESOCIALCONTACTMATRIX;
extern bool G_B_REMOVE_99COLUMNS;

extern double G_CLO_SYMP_FRAC;
extern double G_CLO_HOSPFRAC_YOUNG_DEV; //DEPRECATED
extern double G_CLO_HOSPFRAC_OLD_DEV;   //DEPRECATED
extern double G_CLO_HOSPFRAC_MID_DEV;   //DEPRECATED

extern double G_CLO_HOSPFRAC_10;
extern double G_CLO_HOSPFRAC_20;
extern double G_CLO_HOSPFRAC_30;
extern double G_CLO_HOSPFRAC_40;
extern double G_CLO_HOSPFRAC_50;
extern double G_CLO_HOSPFRAC_60;
extern double G_CLO_HOSPFRAC_70;
extern double G_CLO_HOSPFRAC_80;

extern double G_CLO_DEATHPROB_HOME_60;
extern double G_CLO_DEATHPROB_HOME_70;
extern double G_CLO_DEATHPROB_HOME_80;

extern double G_CLO_RELSUSC_0;
extern double G_CLO_RELSUSC_10;
extern double G_CLO_RELSUSC_20;
extern double G_CLO_RELSUSC_30;
extern double G_CLO_RELSUSC_40;
extern double G_CLO_RELSUSC_50;
extern double G_CLO_RELSUSC_60;
extern double G_CLO_RELSUSC_70;
extern double G_CLO_RELSUSC_80;

extern double G_CLO_TIME_SYMPTOHOSP;
extern double G_CLO_SELFISOLATION_FACTOR;

extern double G_CLO_EARLYMARCH_HOSPRATE;     // scaling factor showing how much more likely hospitalization was in early March than later in the epidemic
extern double G_CLO_EARLYMARCH_ENDDAY;       // the day that the high-hosp rate ended (probably because the patient load became too high)




extern double G_CLO_ICUFRAC_DEV;
extern double G_CLO_PROB_ICU_TO_VENT;
extern double G_CLO_VENTDEATH_MID_DEV;
extern double G_CLO_MEANTIME_ON_VENT_SURV;
extern double G_CLO_RELATIVE_BETA_HOSP;

extern double G_CLO_DEV_LEN_HOSPSTAY;

extern double G_C_COMM[NUMAC][NUMAC];

//  END  ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void PrintUsageModes();
void SetLocationData( string loc );

bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

void PrintUsageModes()
{
    printf("\n\tUSAGE: ./odesim  outfilename   [-options] \n\n");
}


//
// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=2;

    if( argc<2 )
    { 
        PrintUsageModes(); 
	exit(-1);
    }
	
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to outfile\n\n");
    }
    else 
    {
        if( argv[1][0]=='-'){
            fprintf(stderr, "You probably meant to use '%s' as an option, not the output file.\n", argv[1]);
            exit(-1);
        }
        OutFile = fopen( argv[1], "w" );
    }

    string str;

    i=start;
    
    // read in options
    while(i<argc)
    {
	str = argv[i];
        
        // ### 1 ### IF BLOCK FOR BETA
        if( str == "-beta" )
        {
            ppc->v_betas.clear();
            ppc->v_betatimes.clear();
            i++;
            bool bFirstBetaPushedBack = false;
            
            //BEGIN LOOPING THROUGH THE LIST OF BETAS
            while(i<argc)
            {
                string s( argv[i] );
                if( isFloat(s) ) // if the current string is a floating point number, write it into the v_betas array
                {
                    // if the command line argument is <0, just set it back to zero
                    double d = atof( argv[i] );
                    if( d < 0.0 ) d = 0.0;
                    
                    ppc->v_betas.push_back( d );
                    
                    // if it's the first beta value being read in, then push it back twice, as it will be used as 
                    // the Jan 1 to March 1 beta (when there were almost no cases) as well as the beta for the 
                    // first time period we are invesigating
                    if(!bFirstBetaPushedBack) {ppc->v_betas.push_back( d );bFirstBetaPushedBack=true;}
                    
                    // increment and move on in this sub-loop
                    i++;
                }
                else
                {
                    // if the current string is NOT a float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }
                 
            } 
            //END OF LOOPING THROUGH THE LIST OF BETAS
             
            // MAKE SURE AT LEAST ONE BETA VALUE WAS READ IN 
            if( ppc->v_betas.size() == 0 )
            {
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-beta' option must be followed by at least one floating point number.\n\n");
            }
            else
            {
                double total_time_period = G_CLO_TF - 60.0;
                double time_step = total_time_period / ( (double) ppc->v_betas.size() - 1 );
                
                ppc->v_betatimes.push_back(0.0);
                
                for(int jj=1; jj<ppc->v_betas.size(); jj++)
                {
                    ppc->v_betatimes.push_back( 60.0 + time_step*((double)(jj-1)) );
                }
                
            }
                         
        }
        else if( str == "-checkpop" )
        {
            G_B_CHECKPOP_MODE = true;
        }
        /*else if( str == "-dev-hosp-young" )
        {
            G_CLO_HOSPFRAC_YOUNG_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_YOUNG_DEV > 3.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 3.0;
            if( G_CLO_HOSPFRAC_YOUNG_DEV < 0.0 ) G_CLO_HOSPFRAC_YOUNG_DEV = 0.0;
        }
        else if( str == "-dev-hosp-mid" )
        {
            G_CLO_HOSPFRAC_MID_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_MID_DEV > 3.0 ) G_CLO_HOSPFRAC_MID_DEV = 3.0;
            if( G_CLO_HOSPFRAC_MID_DEV < 0.0 ) G_CLO_HOSPFRAC_MID_DEV = 0.0;
        }
        else if( str == "-dev-hosp-old" )
        {
            G_CLO_HOSPFRAC_OLD_DEV = atof( argv[++i] );
            if( G_CLO_HOSPFRAC_OLD_DEV > 1.3 ) G_CLO_HOSPFRAC_OLD_DEV = 1.3;
            if( G_CLO_HOSPFRAC_OLD_DEV < 0.0 ) G_CLO_HOSPFRAC_OLD_DEV = 0.0;
        }*/
        else if( str == "-dev-icu-frac" )
        {
            G_CLO_ICUFRAC_DEV = atof( argv[++i] );
            if( G_CLO_ICUFRAC_DEV > 5.0 ) G_CLO_ICUFRAC_DEV = 5.0;
            if( G_CLO_ICUFRAC_DEV < 0.0 ) G_CLO_ICUFRAC_DEV = 0.0;
        }
        else if( str == "-dev-len-hospstay" )
        {
            G_CLO_DEV_LEN_HOSPSTAY = atof( argv[++i] );
            if( G_CLO_DEV_LEN_HOSPSTAY > 2.5 ) G_CLO_DEV_LEN_HOSPSTAY = 2.5;
            if( G_CLO_DEV_LEN_HOSPSTAY < 0.0 ) G_CLO_DEV_LEN_HOSPSTAY = 0.0;
        }
        else if( str == "-dev-ventdeath-mid" )
        {
            G_CLO_VENTDEATH_MID_DEV = atof( argv[++i] );
            if( G_CLO_VENTDEATH_MID_DEV > 1.7 ) G_CLO_VENTDEATH_MID_DEV = 1.7;
            if( G_CLO_VENTDEATH_MID_DEV < 0.0 ) G_CLO_VENTDEATH_MID_DEV = 0.0;
        }
        
        else if( str == "-death-prob-home-60" )   {     G_CLO_DEATHPROB_HOME_60 = atof( argv[++i] );    }
        else if( str == "-death-prob-home-70" )   {     G_CLO_DEATHPROB_HOME_70 = atof( argv[++i] );    }
        else if( str == "-death-prob-home-80" )   {     G_CLO_DEATHPROB_HOME_80 = atof( argv[++i] );    }
        else if( str == "-diag" )           {           G_B_DIAGNOSTIC_MODE = true;                     }
        else if( str == "-earlymarch-hosp-factor" ) {   G_CLO_EARLYMARCH_HOSPRATE = atof( argv[++i] );  }
        else if( str == "-introday" )       {           G_CLO_INTRODUCTION_TIME = atof( argv[++i] );    }
        else if( str == "-loc" )            {           G_CLO_LOCATION = argv[++i];                     }
        else if( str == "-hosp-frac-10" )   {           G_CLO_HOSPFRAC_10 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-20" )   {           G_CLO_HOSPFRAC_20 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-30" )   {           G_CLO_HOSPFRAC_30 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-40" )   {           G_CLO_HOSPFRAC_40 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-50" )   {           G_CLO_HOSPFRAC_50 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-60" )   {           G_CLO_HOSPFRAC_60 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-70" )   {           G_CLO_HOSPFRAC_70 = atof( argv[++i] );          }
        else if( str == "-hosp-frac-80" )   {           G_CLO_HOSPFRAC_80 = atof( argv[++i] );          }
        else if( str == "-mean-time-vent" ) {           G_CLO_MEANTIME_ON_VENT_SURV=atof( argv[++i] );  }        
        else if( str == "-printIndices" )
        {
            printf("NUMAC %d\n", NUMAC );       // number age groups
            printf("NUME %d\n", NUME );         // number of Exposed-class
            printf("STARTE %d\n", STARTE );     // start index of E
            printf("NUMA %d\n", NUMA );         // Asymptomatic-class
            printf("STARTA %d\n", STARTA );
            printf("NUMI %d\n", NUMI );         // I (symptomatic)-class
            printf("STARTI %d\n", STARTI );
            printf("NUMHA %d\n", NUMHA );       // Hospotalized Acute-class
            printf("STARTHA %d\n", STARTHA );
            printf("STARTCA %d\n", STARTCA );   // Critical care (ICU) Acute-class
            printf("NUMV %d\n", NUMV );         // mechanical Ventilator-class
            printf("STARTV %d\n", STARTV );
            printf("STARTCR %d\n", STARTCR );   // Critical care (ICU) Recovering-class
            printf("STARTHR %d\n", STARTHR );   // Hospitalized Recovering-class
            printf("STARTD %d\n", STARTD );             // Died-at-home class (age-stratified)
            printf("STARTDHOSP %d\n", STARTDHOSP );     // Died-in-hospital class (age-stratified)
            printf("STARTR %d\n", STARTR );             // Recovered-at-home-class (age-stratified)
            printf("STARTRHOSP %d\n", STARTRHOSP );         // Recovered-in-hospital-and-discharged class (age-stratified)
            printf("STARTJ %d\n", STARTJ );   // Cumulative I-class (age-stratified)
            printf("STARTK %d\n", STARTK );   // Cumulative hospitalization incidence-class (age-stratified)

            exit(0);
        }
        else if( str == "-prob-icu-vent" )  {           G_CLO_PROB_ICU_TO_VENT = atof( argv[++i] ); }
        else if( str == "-remove99" )       {           G_B_REMOVE_99COLUMNS = true;                }   
        else if( str == "-version" ) 
        {
            printf("%d\n", 5);
            exit(0);
        }
        else if( str == "-rel-beta-hosp" )
        {
            G_CLO_RELATIVE_BETA_HOSP = atof( argv[++i] );
            if( G_CLO_RELATIVE_BETA_HOSP > 1.0 ) G_CLO_RELATIVE_BETA_HOSP = 1.0;
            if( G_CLO_RELATIVE_BETA_HOSP < 0.0 ) G_CLO_RELATIVE_BETA_HOSP = 0.0;
        }
        else if( str == "-scm" )
        {
            G_B_USESOCIALCONTACTMATRIX = true;
            //G_C_COMM[0][3] = 0.5;
        }
        else if( str == "-self-isolation-factor" )  {            G_CLO_SELFISOLATION_FACTOR = atof( argv[++i] );    }
        else if( str == "-steps-per-day" )          {            G_CLO_STEPS_PER_DAY = atof( argv[++i] );           }
        else if( str == "-susc-0-20" )
        {
            G_CLO_RELSUSC_0  = atof( argv[++i] );
            G_CLO_RELSUSC_10 = G_CLO_RELSUSC_0;
        }
        else if( str == "-susc-60-100" )
        {
            G_CLO_RELSUSC_60 = atof( argv[++i] );
            G_CLO_RELSUSC_70 = G_CLO_RELSUSC_60;
            G_CLO_RELSUSC_80 = G_CLO_RELSUSC_60;            
        }
        else if( str == "-symp-frac" ) // what you're really setting here is the symp fraction for 30-39 year-olds
        {
            G_CLO_SYMP_FRAC = atof( argv[++i] );
            if( G_CLO_SYMP_FRAC > 0.325 ) G_CLO_SYMP_FRAC = 0.325;
            if( G_CLO_SYMP_FRAC < 0.0 ) G_CLO_SYMP_FRAC = 0.0;
        }
        else if( str == "-tf" )                     {            G_CLO_TF = atof( argv[++i] );                      }
        else if( str == "-time-symp-to-hosp" )      {            G_CLO_TIME_SYMPTOHOSP = atof( argv[++i] );         }
        // ### FINAL ### IF BLOCK FOR AN UNKNOWN COMMAND-LINE OPTIONS            
        else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
        
            
        
        
        
        
        //END OF MAIN WHILE-LOOP BLOCK; INCREMENT AND MOVE ON
        i++;
    }
    
    // after parsing all of the arguments, do the following
    
    // 1. set the specific age structure, population, and introduction time for this location
    SetLocationData( G_CLO_LOCATION );


    return;
}

void SetLocationData( string loc )
{
    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out
    
    if( loc=="RI" )
    {
        // 2019 est pop is 1,059,361
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 1059361.0;
        yic[0]  = ppc->v[i_N]   *   .105;   //  0-9
        yic[1]  = ppc->v[i_N]   *   .123;   //  10-19
        yic[2]  = ppc->v[i_N]   *   .140;   //  20-29
        yic[3]  = ppc->v[i_N]   *   .127;   //  30-39
        yic[4]  = ppc->v[i_N]   *   .124;   //  40-49
        yic[5]  = ppc->v[i_N]   *   .135;   //  50-59
        yic[6]  = ppc->v[i_N]   *   .120;   //  60-69
        yic[7]  = ppc->v[i_N]   *   .074;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 55.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;
        //G_CLO_EARLYMARCH_ENDDAY = 84.5; 
        
    }
    else if( loc=="MA" )
    {
        // 2018 est pop is 6,897,212
        

        // data from https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-total.html
        ppc->v[i_N] = 6897212.0;
        yic[0]  = ppc->v[i_N]   *   0.10586466;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12243686;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.14498786;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.13384234;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.12230812;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.14064248;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.11801015;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.06958116;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 61.0;
        
        G_CLO_INTRODUCTION_COUNT = 1;        
        
    }
    else if( loc=="PA" )
    {
        // 2019 est pop is 12,800,721
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 12800721.0;
        yic[0]  = ppc->v[i_N]   *   0.11160395;   //  0-9
        yic[1]  = ppc->v[i_N]   *   0.12229803;   //  10-19
        yic[2]  = ppc->v[i_N]   *   0.13156525;   //  20-29
        yic[3]  = ppc->v[i_N]   *   0.12581869;   //  30-39
        yic[4]  = ppc->v[i_N]   *   0.11809624;   //  40-49
        yic[5]  = ppc->v[i_N]   *   0.13878546;   //  50-59
        yic[6]  = ppc->v[i_N]   *   0.1270166;   //  60-69
        yic[7]  = ppc->v[i_N]   *   0.07657303;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );

        // if the introduction time was never set on the command-line, set it to this value
        if( G_CLO_INTRODUCTION_TIME < 0.0 ) G_CLO_INTRODUCTION_TIME = 62.0;
        
        G_CLO_INTRODUCTION_COUNT = 2;        
        
    }
    else
    {
        fprintf(stderr, "\n\tUnknown location [%s] entered after -loc on command line.\n\n", loc.c_str() );
        exit(-1);
    }
    
}





