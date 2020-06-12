//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "essentials.h"
#include "assert.h"
#include "prms.h"


// constructor
prms::prms()
{
    v.insert( v.begin(), num_params, 0.0 );
    assert( v.size()==num_params );
    
    v_rel_susc.insert( v_rel_susc.begin(), NUMAC, 0.0 );
    v_prob_E_A.insert( v_prob_E_A.begin(), NUMAC, 0.0 );
    v_prob_I2_H.insert(  v_prob_I2_H.begin(),  NUMAC, 0.0 );
    v_prob_I4_D.insert(  v_prob_I4_D.begin(),  NUMAC, 0.0 );
    v_prob_HA4_D.insert( v_prob_HA4_D.begin(), NUMAC, 0.0 );
    v_prob_HA_CA.insert( v_prob_HA_CA.begin(), NUMAC, 0.0 );
    v_prob_V_D.insert( v_prob_V_D.begin(), NUMAC, 0.0 );
    v_prob_CA_D.insert( v_prob_CA_D.begin(), NUMAC, 0.0 );
    v_prob_CA_V.insert( v_prob_CA_V.begin(), NUMAC, 0.0 );
    
    assert( v_rel_susc.size()==NUMAC );    
    assert( v_prob_E_A.size()==NUMAC );
    assert( v_prob_I2_H.size() ==NUMAC );
    assert( v_prob_I4_D.size() ==NUMAC );
    assert( v_prob_HA4_D.size()==NUMAC );
    assert( v_prob_HA_CA.size()==NUMAC );
    assert( v_prob_V_D.size()==NUMAC );
    assert( v_prob_CA_D.size()==NUMAC );
    assert( v_prob_CA_V.size()==NUMAC );
    
    v_betas.clear();
    v_betatimes.clear();
    assert( v_betas.size() == 0 );
    assert( v_betatimes.size() == 0 );
    
    earlymarch_highhosp_period = false;
    earlymarch_highhosp_factor = 1.0;
    earlymarch_highhosp_endday = -1.0;
    
    index_current_beta=-1;
    
}

// destructor
prms::~prms()
{
}


void prms::assign_new_beta( void )
{
    assert( v_betas.size() == v_betatimes.size() );
    assert( v_betas.size() > 0 );
    
    // if it's invalid, just set it to zero
    // this is what happens at initialization as well
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        index_current_beta = 0;
    }
    else // if it's valid, increment it if you can
    {
        if( index_current_beta != v_betas.size() - 1 ) // if it's NOT the last element, then increment it
        {
            index_current_beta++;
        }
    }
    
    // assign the possibly new beta value to the main parameter vector "v"
    v[i_beta] = v_betas[ index_current_beta ];
}


double prms::get_new_update_time( void )
{
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        return 1000000.0;
    }
    else // if the index is within range, return the next update time
    {
        if( index_current_beta != v_betas.size() - 1 ) // if the index does NOT point to the last element, then return the next element
        {
            return v_betatimes[ index_current_beta+1 ];
        }
        else // if it is the last element in the array, just return 1000000.0, as in the next update time is a million days from now
        {
            return 1000000.0;
        }
    }
    
    
}

void prms::apply_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) / earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = true;
}

void prms::end_earlymarch_hosprates( void )
{
    for(int ac=0; ac<NUMAC; ac++)
    {
        v_prob_I2_H[ac] = 1.0 - ( (1.0-v_prob_I2_H[ac]) * earlymarch_highhosp_factor );
    }
    earlymarch_highhosp_period = false;
}

