/* ------------------------------------------------------------------
   FILENAME:   num_f_hermite.c
   
   PURPOSE:    provide the roots of the Hermite polynomials for a 
               quick integration of a function f on the gaussian 
			   density as a discrete sum on the w[]
               (cf Numerical Recipes in C p. 154)
   ------------------------------------------------------------------ */

#include "num_h_allhdr.h"
#include "num_h_hermite.h"
#include "utconst.h"
#include "math.h"

#define   HERMITE_EPS        1.0e-14
#define   PIM4               0.7511255444649425
#define   MAXIT              10

Err store_hermite_in_file(int n)
{
	FILE     *f;
	long     j;
	long     k;
	double   *x;
	double   *w;
	Err      err;

	f = fopen("c:\\temp\\hermite.dat", "w");
	for (j = 0; j <= n; j++)
	{
		x = dvector(1,j);
		w = dvector(1,j);
		err = gauss_hermite(x, w, j);
		if (err)
		{
			free_dvector(x, 1, j);
			free_dvector(w, 1, j);
			return err;
		}
		fprintf(f, " N = %d", j);
		fprintf(f, "\n");
		for ( k = 1; k <=j; k++)
		{
			fprintf(f,"%.16f\n",x[k]);
		}
		fprintf(f,"\n");
		for ( k = 1; k <=j; k++)
		{
			fprintf(f,"%.16f\n",w[k]);
		}
		fprintf(f,"\n");
		fprintf(f,"\n");
		
		free_dvector(x, 1, j);
		free_dvector(w, 1, j);
		
	}
	fclose(f);
	
	return NULL;
}

/* --------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------
   Given n, this routine returns arrays x[1..n] and w[1..n] containing the abscissas
   and weights of the n-point Gauss-Hermite quadrature formula. The largest abscissa
   is returned in x[1], the most negative in x[n]
   Mind you, the integration function is exp(-x*x) and not the real Gaussian density
   --------------------------------------------------------------------------------- */
Err gauss_hermite(double x[], double w[], int n)
{
	int     i,its,j,m;
	double  p1,p2,p3,pp,z,z1;

	m = (n+1)/2;
	for (i=1;i<=m;i++) 
	{
		if (i == 1) 
		{
			z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
		} 
		else 
		if (i == 2) 
		{
			z -= 1.14*pow((double)n,0.426)/z;
		}
		else 
		if (i == 3) 
		{
			z=1.86*z-0.86*x[1];
		} 
		else 
		if (i == 4) 
		{
			z=1.91*z-0.91*x[2];
		}
		else 
		{
			z=2.0*z-x[i-2];
		}
		
		for (its=1;its<=MAXIT;its++) 
		{
			p1 = PIM4;
			p2 = 0.0;
			for (j=1;j<=n;j++) 
			{
				p3 = p2;
				p2 = p1;
				p1 = z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
			}
			pp = sqrt((double)2*n)*p2;
			z1 = z;
			z = z1-p1/pp;
			if (fabs(z-z1) <= HERMITE_EPS) 
				break;
		}
		if (its > MAXIT) 
			return serror("Too many iterations in Gauss Hermite");
		x[i] = z;
		x[n+1-i] = -z;
		w[i] = 2.0/(pp*pp);
		w[n+1-i] = w[i];
                if (i == 0) {
                  z = sqrt((double)(2 * n + 1)) -
                      1.85575 * pow((double)(2 * n + 1), -0.16667);
                } else if (i == 1) {
                  z -= 1.14 * pow((double)n, 0.426) / z;
                } else if (i == 2) {
                  z = 1.86 * z - 0.86 * x[1];
                } else if (i == 3) {
                  z = 1.91 * z - 0.91 * x[2];
                } else {
                  z = 2.0 * z - x[i - 2];
                }

                for (its = 1; its <= MAXIT; its++) {
                  p1 = PIM4;
                  p2 = 0.0;
                  for (j = 1; j <= n; j++) {
                    p3 = p2;
                    p2 = p1;
                    p1 = z * sqrt(2.0 / j) * p2 -
                         sqrt(((double)(j - 1)) / j) * p3;
                  }
                  pp = sqrt((double)2 * n) * p2;
                  z1 = z;
                  z = z1 - p1 / pp;
                  if (fabs(z - z1) <= HERMITE_EPS) break;
                }
                if (its > MAXIT)
                  return serror("Too many iterations in Gauss Hermite");
                x[i] = z;
                x[n + 1 - i] = -z;
                w[i] = 2.0 / (pp * pp);
                w[n + 1 - i] = w[i];
	}

	return NULL;
}

#undef HERMITE_EPS
#undef PIM4
#undef MAXIT
/* (C) Copr. 1986-92 Numerical Recipes Software +135[)6=. */

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------
   Given n, this routine returns arrays x[1..n] and w[1..n] containing the abscissas
   and weights of the n-point Gauss-Hermite quadrature formula. The largest abscissa
   is returned in x[1], the most negative in x[n].
   If the polynomials roots have already been computed, it does not recalculate them
   Mind you, the integration function is exp(-x*x) and not the real Gaussian density
   --------------------------------------------------------------------------------- */

/* Static variables to store the coefficients one computed */
static long      _hermite_prev_n  = -33333;
static double   *_herm_x          = NULL ;
static double   *_herm_w          = NULL;

Err hermite_gauss_quick(int n, double x[], double w[])
{
Err             err       = NULL;
long            i;

	if (n < 0 )
	{
		return serror("Hermite polynomials of negative degree: U must be KIDDING !!!");
	}
	if (_hermite_prev_n == n)
	{
		for (i = 1; i <= n; i++)
		{
			x[i] = _herm_x[i];
			w[i] = _herm_w[i];
		}
	}
	else
	{
		err = gauss_hermite(x, w, n);
		if (err)
		{
			return err;
		}
		
		if (_herm_x) free_dvector(_herm_x, 1, _hermite_prev_n);
		if (_herm_w) free_dvector(_herm_w, 1, _hermite_prev_n);
		
		_hermite_prev_n = n;

		_herm_x = dvector(1, n);
		_herm_w = dvector(1, n);

		memcpy(&_herm_x[1], &x[1], n * sizeof(double));
		memcpy(&_herm_w[1], &w[1], n * sizeof(double));
	}

	return NULL;
}


Err HermiteStandard(double *x, double *w, int n)
{
int i;
	gauss_hermite(x, w, n); 
	x[0] = w[0] = 0;	
	for (i = 1; i <= n; i++)
	{
		x[i] *= SQRT_TWO;
		w[i] *= INV_SQRT_PI;
	}
	return NULL;
}

		

