#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "error2.h"
#include "proba_utils.h"
#include "integration_nr.h"

#define NRANSI
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define NRANSI
#define MAXSTP 10000
#define TINY 1.0e-30
#ifndef MAX
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#endif

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


int rkck(   double y[],
            double dydx[],
            int n,
            double x,
            double h,
            double yout[],
	        double yerr[],
            int (*derivs)(  double,
                            double [],
                            double [],
                            void*),
            void *param)
{
    static char routine[] = "rkck";
    int status = FAILURE;
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.0/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double  *ak2 = NULL,
            *ak3 = NULL,
            *ak4 = NULL,
            *ak5 = NULL,
            *ak6 = NULL,
            *ytemp= NULL;

	ak2=malloc(n*sizeof(double));
    if(ak2==NULL) goto RETURN;
	ak3=malloc(n*sizeof(double));
    if(ak3==NULL) goto RETURN;
	ak4=malloc(n*sizeof(double));
    if(ak4==NULL) goto RETURN;
	ak5=malloc(n*sizeof(double));
    if(ak5==NULL) goto RETURN;
	ak6=malloc(n*sizeof(double));
    if(ak6==NULL) goto RETURN;
	ytemp=malloc(n*sizeof(double));
    if(ytemp==NULL) goto RETURN;

	for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+b21*h*dydx[i];
    }
	(*derivs)(x+a2*h,ytemp,ak2, param);
	for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    }
	(*derivs)(x+a3*h,ytemp,ak3, param);
	for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    }
	(*derivs)(x+a4*h,ytemp,ak4, param);
	for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    }
	(*derivs)(x+a5*h,ytemp,ak5, param);
	for (i=0;i<n;i++)
    {
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    }
	(*derivs)(x+a6*h,ytemp,ak6, param);
	for (i=0;i<n;i++)
    {
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    }
	for (i=0;i<n;i++)
    {
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    }

    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(ytemp) free(ytemp);
    if(ak6) free(ak6);
    if(ak5) free(ak5);
    if(ak4) free(ak4);
    if(ak3) free(ak3);
    if(ak2) free(ak2);
    return status;
}


int rkqs(   double y[],
            double dydx[],
            int n,
            double *x,
            double htry,
            double eps,
	        double yscal[],
            double *hdid,
            double *hnext,
	        int (*derivs)( double,
                            double [],
                            double [],
                            void *),
            void *param)
{
    static char routine[] = "rkqs";
    int status = FAILURE;
	int i;
	double  errmax,
            h,
            xnew,
            *yerr = NULL,
            *ytemp = NULL;

	yerr=malloc(n*sizeof(double));
    if(yerr==NULL) goto RETURN;
	ytemp=malloc(n*sizeof(double));
    if(ytemp==NULL) goto RETURN;
	h=htry;
	
    for (;;)
    {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs, param);
		errmax=0.0;
		for (i=0;i<n;i++)
        {
            errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
        }
		errmax /= eps;
		if (errmax > 1.0) {
			h=SAFETY*h*pow(errmax,PSHRNK);
			if (h < 0.1*h)
            {
                h *= 0.1;
            }
			xnew=(*x)+h;
			if (xnew == *x)
            {   
                DR_Error("stepsize underflow in rkqs");
            }
			continue;
		}
        else
        {
			if (errmax > ERRCON)
            {
                *hnext=SAFETY*h*pow(errmax,PGROW);
            }
			else
            {
                *hnext=5.0*h;
            }
			*x += (*hdid=h);
            for (i=0;i<n;i++)
            {
                y[i]=ytemp[i];
            }
			break;
		}
	}
    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);   
    }
	if(ytemp) free(ytemp);
	if(yerr) free(yerr);
    return status;
}


int odeint( double ystart[],
            int nvar,
            double x1,
            double x2,
            double eps,
            double h1,
	        double hmin,
            int maxNbPoints,
            int *nok,
            int *nbad,
	        int (*derivs)(  double,
                            double [],
                            double [],
                            void *),
	        int (*rkqs)(    double [],
                            double [],
                            int,
                            double *,
                            double,
                            double,
                            double [],
	                        double *,
                            double *,
                            int  (*)(   double,
                                        double [],
                                        double [],
                                        void *),
                            void *),
            void *param,
            DEBUGINFO *debugInfo)
{
    static char routine[] = "odeint";
    int status = FAILURE;
	int nstp,i;
	double xsav=0.0,x,hnext,hdid,h;
	double  *yscal = NULL,
            *y = NULL,
            *dydx = NULL;
    double *xp = NULL;
    double **yp = NULL;
    double dxsav =0.0;
    int kmax = 0;
    int kount;
    
    /* variables to store intermediate values */
    if(debugInfo)
    {
        kmax = debugInfo->kmax,kount;
        xp = debugInfo->xp,
        yp = debugInfo->yp,
        dxsav = debugInfo->dxsav;
    }
            
	yscal=malloc(nvar*sizeof(double));
    if(yscal==NULL) goto RETURN;
	y=malloc(nvar*sizeof(double));
    if(y==NULL) goto RETURN;
	dydx=malloc(nvar*sizeof(double));
    if(dydx==NULL) goto RETURN;

	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=0;nstp<maxNbPoints;nstp++)
    {
		(*derivs)(x,y,dydx, param);
		for (i=0;i<nvar;i++)
			yscal[i]=/*fabs(y[i])+*/fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav))
        {
			xp[kount]=x;
			for (i=0;i<nvar;i++) yp[i][kount]=y[i];
			xsav=x;
            kount++;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0)
        {
            h=x2-x;
        }
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs, param);
		if (hdid == h)
        {   
            ++(*nok);
        }
        else
        {
            ++(*nbad);
        }
		if ((x-x2)*(x2-x1) >= 0.0)
        {
			for (i=0;i<nvar;i++)
            {
                ystart[i]=y[i];
            }
			if (kmax)
            {
				xp[kount]=x;
				for (i=0;i<nvar;i++)
                {
                    yp[i][kount]=y[i];
                }
			}
            if(kmax)
            {
                debugInfo->nbPoints = kount;
            }
            status = SUCCESS;
			goto RETURN;
		}
		if (fabs(hnext) <= hmin)
        {
            DR_Error("Step size too small in odeint");
        }
		h=hnext;
	}
	DR_Error("Too many steps in routine odeint");
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);        
    }
    if(dydx) free(dydx);
	if(y) free(y);
    if(yscal) free(yscal);
    return status;
}

#undef MAXSTP
#undef TINY
#undef NRANSI
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software .@7M4]2+N3. */
