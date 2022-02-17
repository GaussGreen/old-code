/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: solver.cpp,v $
 * Revision 1.1  2004/09/22 10:15:09  ebenhamou
 * Initial revision
 *
 *
 *
 */


/*! \file nagsolver.cpp
 *
 *  \brief
 *
 *	\author  A. TRIKI
 *	\version 1.0
 *	\date February 2005
 */

#include "gpnumlib/nagsolver.h"

CC_BEGIN_NAMESPACE( ARM )

void  NAG_CALL WeightedSquared(double x, 
		double* fx,
	    Nag_Comm *comm)
{ 
	UnaryFuncWithDerivative<double,double>* solver = static_cast<UnaryFuncWithDerivative<double,double>* > ((comm->p));
	double result= (*solver)(x);
	result/=0.0000001;
    *fx = result*result;
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

