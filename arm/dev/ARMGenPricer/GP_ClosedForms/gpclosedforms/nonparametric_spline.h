/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file nonparametric_spline.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date fev 2007
 */
 
#ifndef _GP_CF_NONPARAMETRIC_SPLINE_H
#define _GP_CF_NONPARAMETRIC_SPLINE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"

#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

void cubicspline_precompute(const double *x, const double *y, int n,
							const double yp1, const double ypn,double * y2);

double cubicspline_interpolate(const double * xa, const double * ya,
							 const double * y2a, int n, const double x);

double cubicspline_inter(const double * xa, const double * ya, const double * y2a, int n, double x);

double cubicspline_interder(const double * xa, const double * ya, const double * y2a, int n, double x);


inline void cubicspline_precompute(const ARM_GP_Vector *x, const ARM_GP_Vector *y,
							const double yp1, const double ypn,ARM_GP_Vector &y2)
{
	if(x->size() != y->size() || x->size() != y2.size())
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a should have same size");

	cubicspline_precompute(x->begin(),y->begin(),x->size(),yp1,ypn,y2.begin());
}

inline double cubicspline_interpolate(const ARM_GP_Vector &xa, const ARM_GP_Vector &ya,
							 const ARM_GP_Vector &y2a, const double x)
{
	if(xa.size() != ya.size() || xa.size() != y2a.size())
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a should have same size");

	return cubicspline_interpolate(xa.begin(),ya.begin(),y2a.begin(),y2a.size(),x);
}


inline double cubicspline_inter(const ARM_GP_Vector& xa, const ARM_GP_Vector& ya, const ARM_GP_Vector& y2a, double x)
{
	if(xa.size() != ya.size() || xa.size() != y2a.size())
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a should have same size");

	return cubicspline_inter(xa.begin(), ya.begin(), y2a.begin(), y2a.size(), x);
}


inline double cubicspline_interder(const ARM_GP_Vector& xa, const ARM_GP_Vector& ya, const ARM_GP_Vector& y2a, double x)
{
	if(xa.size() != ya.size() || xa.size() != y2a.size())
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" cubicspline_precompute : x, y and y2a should have same size");

	return cubicspline_interder(xa.begin(), ya.begin(), y2a.begin(), y2a.size(), x);
}

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

