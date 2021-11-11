#include <math.h>
#include <stdio.h>

#include "rkf.h"

#define TINY 1.0e-30
#define SAFTY 0.9
#define PGROW -0.2
#define PSHRINK -0.25
#define ERRCON 1.89e-4

// NOTE this needs to match the "dim" variable that is calculated at the beginning of the generate_trajectories function
#define Q 306


void rkck(double y[],double dydx[],int n,double x,double h, double yout[],
          double yerr[],void (*derivs)(double, double*, double*))
/* 5th order Runge-Kutta-Fehlberg method Num Res (2ed) p 719 */
{ int i;
  static double a2= 0.2, a3= 0.3, a4= 0.6, a5= 1.0, a6= 0.875,
                b21= 0.2, b31= 3.0/40.0, b32= 9.0/40.0,
                b41= 0.3, b42= -0.9, b43= 1.2,
                b51= -11.0/54.0, b52= 2.5, b53= -70.0/27.0, b54= 35.0/27.0,
                b61= 1631.0/55296.0, b62= 175.0/512.0, b63= 575.0/13824.0,
                b64= 44275.0/110592.0, b65= 253.0/4096.0,
                c1= 37.0/378.0, c3= 250.0/621.0, c4= 125.0/594.0,
                c6= 512.0/1771.0, dc5= -277.0/14336.0;
  double dc1= c1-2825.0/27648.0, dc3= c3-18575.0/48384.0,
         dc4= c4-13525.0/55296.0, dc6= c6-0.25;
  double ak2[Q],ak3[Q],ak4[Q],ak5[Q],ak6[Q],ytemp[Q];

  for(i=0;i<n;i++) ytemp[i] = y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2);

  for(i=0;i<n;i++) ytemp[i] = y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);

  for(i=0;i<n;i++) ytemp[i] = y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);

  for(i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);

  for(i=0;i<n;i++)
    ytemp[i] = y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);

  for(i=0;i<n;i++)
  { yout[i] = y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    yerr[i] = h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  }
}


int rkqs(double y[],double dydx[],int n, double *x,double htry,double eps,
          double yscal[],double *hdid,double *hnext,void (*derivs)(double, double*, double*))
{ int i;
  double errmax,h,htemp,xnew,yerr[Q],ytemp[Q];

  h=htry;
  for(;;)
  { rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    errmax= 0.0;
    for(i=0;i<n;i++)
       if(errmax<fabs(yerr[i]/yscal[i])) errmax=fabs(yerr[i]/yscal[i]);
    errmax /= eps;
    if(errmax>1.0)
     { htemp=SAFTY*h*pow(errmax,PSHRINK);
       if(h>=0.0) if(htemp>0.1*h) h=htemp; else h *= 0.1;
       if(h<0.0)  if(htemp<0.1*h) h=htemp; else h *= 0.1;
       xnew = (*x)+h;
       if(xnew==(*x))
       { fprintf(stderr,"stepsize underflow in rkqs\n");
         fprintf(stderr,"x= %f errmax = %f\n",*x,errmax); return(-1);}
     }
    else
     { if(errmax>ERRCON) *hnext=SAFTY*h*pow(errmax,PGROW);
        else *hnext=5.0*h;
       *x += (*hdid=h);
       for(i=0;i<n;i++) y[i] = ytemp[i];
       return(1);
     }
   }
}


int odeint(double ystart[],int nvar,double x1,double x2,double eps,double *h1,
            double hmin,int *nok,int *nbad, void (*derivs)(double, double*, double*),
			int (*rkqs)(double y[],double dydx[],int n, double *x,double htry,double eps,
          double yscal[],double *hdid,double *hnext,void (*derivs)(double, double*, double*)))
{     int i;
      double xsav,x,hnext,hdid,h;
      double yscal[Q],y[Q],dydx[Q];

      x=x1;
      if(x2>x1) h=(*h1); else h= -(*h1);
      *nok = (*nbad) = 0;
      for (i=0;i<nvar;i++) y[i]=ystart[i];
      for (;;)
      { (*derivs)(x,y,dydx);
        for (i=0;i<nvar;i++)
          yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
        if( ((*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs)) !=1)
          { fprintf(stderr,"leaving ODEINT\n"); return(-1);}
        if (hdid == h) ++(*nok); else ++(*nbad);
        if ((x-x2)*(x2-x1) >= 0.0)
        { for (i=0;i<nvar;i++) ystart[i]=y[i];
          (*h1)= hnext;
          return(1);
        }
        if (fabs(hnext) <= hmin)
        { fprintf(stderr,"Step size too small in ODEINT");
          return(-1);
        }
        h=hnext;
      }
}


