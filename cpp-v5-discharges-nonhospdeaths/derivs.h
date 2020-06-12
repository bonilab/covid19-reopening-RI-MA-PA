
#ifndef _DERIVS_
#define _DERIVS_


//void read_init(void);
void derivs( double t, double *y, double *dydt);

enum host_state {exposed, infected, severe, quarantined, hospitalized, endstate}; 
typedef enum host_state HOST_STATE;

#endif



