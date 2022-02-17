/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file nonparametric_spline.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date	fev 2007
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/removenagwarning.h"

#include <cmath>
#include <complex>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"


#include "gpclosedforms/long_double.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/nonparametric_spline.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000


/// The precompute routine fill the vector y2 : 

/*
void cubicspline_precompute(const ARM_GP_Vector *x, const ARM_GP_Vector *y, const double yp1, const double ypn,
	ARM_GP_Vector &y2)
{
	int i,k;
	double p,qn,sig,un;

	int n=y2.size();
	ARM_GP_Vector u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/((*x)[1]-(*x)[0]))*(((*y)[1]-(*y)[0])/((*x)[1]-(*x)[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=((*x)[i]-(*x)[i-1])/((*x)[i+1]-(*x)[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=((*y)[i+1]-(*y)[i])/((*x)[i+1]-(*x)[i]) - ((*y)[i]-(*y)[i-1])/((*x)[i]-(*x)[i-1]);
		u[i]=(6.0*u[i]/((*x)[i+1]-(*x)[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/((*x)[n-1]-(*x)[n-2]))*(ypn-((*y)[n-1]-(*y)[n-2])/((*x)[n-1]-(*x)[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


/// assume that the vector y2a is precomputed;the second derivatives

double cubicspline_interpolate(const ARM_GP_Vector &xa, const ARM_GP_Vector &ya,const ARM_GP_Vector &y2a, const double x)
{
	int k;
	double h,b,a;

	int n=xa.size();
	if ((x<xa[0]) || (x>xa[n-1])) Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"cubicspline_interpolate: x ouside bounds " );
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"cubicspline_interpolate: bad xa " );
			
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double cubicspline_inter(const ARM_GP_Vector& xa, const ARM_GP_Vector& ya, const ARM_GP_Vector& y2a, double x)
{
	int n = xa.size();

	if(x >= xa[0] && x <= xa[n-1])
		return cubicspline_interpolate(xa,ya,y2a,x);

	if(x < xa[0])
	{
		double leftder = (ya[1]-ya[0])/(xa[1]-xa[0])-(xa[1]-xa[0])/6.*y2a[1];
		double leftval = ya[0];
		if (leftder>0){
			return leftval*(1.+tanh(leftder/leftval*(x-xa[0])));
		}else{
			return leftval+leftder*(x-xa[0]);
		}
	}
	else
	{
		double rightder = (ya[n-1]-ya[n-2])/(xa[n-1]-xa[n-2])+(xa[n-1]-xa[n-2])/6.*y2a[n-2];
		double rightval = ya[n-1];
		if (rightder<0){
			return rightval*(1+tanh(rightder/rightval*(x-xa[n-1])));
		}else{
			return rightval+rightder*(x-xa[n-1]);
		}
	}
}

double cubicspline_interder(const ARM_GP_Vector& xa, const ARM_GP_Vector& ya, const ARM_GP_Vector& y2a, double x)
{
	int n = xa.size();

	if(x >= xa[0] && x <= xa[n-1])
	{
		int klo,khi,k;
		double h,b,a;
		int n=xa.size();
		klo=0; 
		khi=n-1;
		while (khi-klo > 1) {
			k=(khi+klo) >> 1;
			if (xa[k] > x) khi=k;
			else klo=k;
		}
		h=xa[khi]-xa[klo];
		if (h == 0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SplineDensityFunctor::ComputeDVol: bad input");
		a=(xa[khi]-x)/h;
		b=(x-xa[klo])/h;
		return 1./h*(-ya[klo]+ya[khi]+(-(3*a*a-1)*y2a[klo]+(3*b*b-1)*y2a[khi])*(h*h)/6.0);
	}
	else if(x < xa[0])
	{
		double leftder = (ya[1]-ya[0])/(xa[1]-xa[0])-(xa[1]-xa[0])/6.*y2a[1];
		double leftval = ya[0];
		if (leftder>0){
			double res=cosh(leftder/leftval*(x-xa[0]));
			return leftder/res/res;
		}else{
			return leftder;
		}
	}
	else
	{
		double rightder = (ya[n-1]-ya[n-2])/(xa[n-1]-xa[n-2])+(xa[n-1]-xa[n-2])/6.*y2a[n-2];
		double rightval = ya[n-1];
		if (rightder<0){
			double res=cosh(rightder/rightval*(x-xa[n-1]));
			return rightder/res/res;
		}else{
			return rightder;
		}
	}
}
*/

void cubicspline_precompute(const double *x, const double *y, int n, const double yp1, const double ypn,
	double * y2)
{
	if(x == NULL || y == NULL || y2 == NULL)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a are NULL");

	int i,k;
	double p,qn,sig,un;

	if(n < 3) return;

	ARM_GP_Vector u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


/// assume that the vector y2a is precomputed;the second derivatives

double cubicspline_interpolate(const double * xa, const double * ya,const double * y2a, int n, const double x)
{
	if(xa == NULL || ya == NULL || y2a == NULL)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a are NULL");

	if(n < 3)
	{
		return n == 1 ? ya[0] : ya[0] + (x-xa[0])*(ya[1]-ya[0])/(xa[1]-xa[0]);
	}

	int k;
	double h,b,a;

	if ((x<xa[0]) || (x>xa[n-1])) Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"cubicspline_interpolate: x ouside bounds " );
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"cubicspline_interpolate: bad xa " );
			
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double cubicspline_inter(const double * xa, const double * ya, const double * y2a, int n, double x)
{
	if(xa == NULL || ya == NULL || y2a == NULL)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a are NULL");

	if(x >= xa[0] && x <= xa[n-1])
		return cubicspline_interpolate(xa,ya,y2a,n,x);

	if(x < xa[0])
	{
		if(n < 3) return ya[0];

		double leftder = (ya[1]-ya[0])/(xa[1]-xa[0])-(xa[1]-xa[0])/6.*y2a[1];
		double leftval = ya[0];
		if (leftder>0){
			return leftval*(1.+tanh(leftder/leftval*(x-xa[0])));
		}else{
			return leftval+leftder*(x-xa[0]);
		}
	}
	else
	{
		if(n < 3) return ya[n-1];

		double rightder = (ya[n-1]-ya[n-2])/(xa[n-1]-xa[n-2])+(xa[n-1]-xa[n-2])/6.*y2a[n-2];
		double rightval = ya[n-1];
		if (rightder<0){
			return rightval*(1+tanh(rightder/rightval*(x-xa[n-1])));
		}else{
			return rightval+rightder*(x-xa[n-1]);
		}
	}
}

double cubicspline_interder(const double * xa, const double * ya, const double * y2a, int n, double x)
{
	if(xa == NULL || ya == NULL || y2a == NULL)
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a are NULL");

	if(x >= xa[0] && x <= xa[n-1])
	{
		int klo,khi,k;
		double h,b,a;
		klo=0; 
		khi=n-1;
		while (khi-klo > 1) {
			k=(khi+klo) >> 1;
			if (xa[k] > x) khi=k;
			else klo=k;
		}
		h=xa[khi]-xa[klo];
		if (h == 0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SplineDensityFunctor::ComputeDVol: bad input");
		a=(xa[khi]-x)/h;
		b=(x-xa[klo])/h;
		return 1./h*(-ya[klo]+ya[khi]+(-(3*a*a-1)*y2a[klo]+(3*b*b-1)*y2a[khi])*(h*h)/6.0);
	}
	else if(x < xa[0])
	{
		double leftder = (ya[1]-ya[0])/(xa[1]-xa[0])-(xa[1]-xa[0])/6.*y2a[1];
		double leftval = ya[0];
		if (leftder>0){
			double res=cosh(leftder/leftval*(x-xa[0]));
			return leftder/res/res;
		}else{
			return leftder;
		}
	}
	else
	{
		double rightder = (ya[n-1]-ya[n-2])/(xa[n-1]-xa[n-2])+(xa[n-1]-xa[n-2])/6.*y2a[n-2];
		double rightval = ya[n-1];
		if (rightder<0){
			double res=cosh(rightder/rightval*(x-xa[n-1]));
			return rightder/res/res;
		}else{
			return rightder;
		}
	}
}

CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
