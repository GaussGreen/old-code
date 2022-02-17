/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: gplinalginterpol.h,v $
 * Revision 1.1  2004/12/08 16:47:59  ebenhamou
 * Initial revision
 */

/*! \file gplinalginterpol.h
 *
 *  \brief files to do all the interpolation for a given gplinalg element 
 *			that can be either a vector, a matrix or a tensor.
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPINFRA_GPLINALGINTERPOL_H
#define _INGPINFRA_GPLINALGINTERPOL_H

#include "port.h"
#include "env.h"


CC_BEGIN_NAMESPACE( ARM )

template <typename T> T LinearInterpolation( T* xElt, T* yElt, size_t xSize, T xVal)
{
	size_t i = CC_NS(std,lower_bound)(xElt,xElt+xSize,xVal)-xElt;

#if defined(__GP_STRICT_VALIDATION )
	CheckPointorVectorIncreasing<T>(xElt,xSize,"xElt", "LinearInterpolation", __LINE__, __FILE__  );
#endif

	if ( i==xSize-1 )
		return yElt[i];
	else
		if( xElt[i+1]-xElt[i] < K_NEW_DOUBLE_TOL )
			return yElt[i+1]
		else
			return yElt[i]+(yElt[i+1]-yElt[i])/(xElt[i+1]-xElt[i])*(xVal-xElt[i]);
}

template <typename T> T LinearInterpolation(T x, T x1, T y1, T x2, T y2)
{
    if (x1 == x2) return y1;
    return ((x2 - x)*y1 + (x - x1)*y2)/(x2 - x1);
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
