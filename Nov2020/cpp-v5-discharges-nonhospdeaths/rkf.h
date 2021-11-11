
#ifndef _RKF_
#define _RKF_

void rkck(double y[],double dydx[],int n,double x,double h, double yout[],
          double yerr[],void (*derivs)(double, double*, double*));
		  
int rkqs(double y[],double dydx[],int n, double *x,double htry,double eps,
          double yscal[],double *hdid,double *hnext,void (*derivs)(double, double*, double*));
		  
int odeint(double ystart[],int nvar,double x1,double x2,double eps,double *h1,
            double hmin,int *nok,int *nbad, void (*derivs)(double, double*, double*),
			int (*rkqs)(double y[],double dydx[],int n, double *x,double htry,double eps,
			double yscal[],double *hdid,double *hnext,void (*derivs)(double, double*, double*))); 
			
#endif

