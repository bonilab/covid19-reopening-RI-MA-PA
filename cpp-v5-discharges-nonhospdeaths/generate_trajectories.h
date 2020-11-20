#include "essentials.h"
#include "prms.h"

#ifndef _GENERATE_TRAJECTORIES_
#define _GENERATE_TRAJECTORIES_

void generate_trajectories( double inital_number_of_cases, double param_beta, double t0, double tf, double h );

void write(double tt, double * yic, size_t dim);

void write_bin(double tt, double * yic, size_t dim);

#endif // _GENERATE_TRAJECTORIES_



