// S.Galluccio: 21 November 2000
// Added function for 1D Guassian quadrature with pointers only

/* ======================================================
   FILENAME:  num_f_gausslegendre.c
   
   PURPOSE:   integration by Gauss Legendre (NUMC p152)
   ====================================================== */
#include "utallhdr.h"
#include "num_h_allhdr.h"
#include "math.h"
#define LEGENDRE_PRECISION 3.0e-11


void gauleg(double x1,double x2,double x[], double w[],int n)

{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x1+x2);
	xl=0.5*(x2-x1);
	for(i=1;i<=m;i++){

		z=cos(SRT_PI*(i-0.25)/(n+0.5));
		do 
		{
			p1=1.0;
			p2=0.0;
			for(j=1;j<=n;j++) 
			{
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1)>LEGENDRE_PRECISION);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}


/*-------------------------------------------------------------------------------------------------------
---------- Same routine as before with pointers as actual arguments only -------------------------------*/


void GaussLeg(	double x1, 
				double x2, 
				double *x, 
				double *w, 
				int n)
{
	int m,j,i;
	double q=3.141592654/(n+0.5);
	double z1,z,xm,xl,pp,p3,p2,p1;
	static int s=0;
	static double*r=NULL;

	if (n>s) {
		r=realloc(r,n*sizeof(double));
		for (j=s;j<n;j++)
			r[j]=1.0/(j+1);
		s=n;
	} else if (n==0) {
		free(r);
		r=NULL;
		s=0;
	}
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(q*(i-0.25));
		do {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j+1)*z*p2-j*p3)*r[j];
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > LEGENDRE_PRECISION);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}





