/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gaussian_integrals.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_INVERSE_H
#define _GP_CF_INVERSE_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include <functional>

#include "expt.h"   // for the exceptions


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000


///#define MAXIT 60
#define UNUSED (-1.11e30)


//////////////////////////////////////////////////////////////////////////
///
///  Inverse Functor
///
//////////////////////////////////////////////////////////////////////////

// C++ Functor
struct DoubleToDoubleFunc : public CC_NS(std,unary_function)<double,double> 
{
	virtual double operator()(double ) const = 0;
};

// C Functor wrapper
struct ptr_funDoubleToDouble : DoubleToDoubleFunc 
{
	typedef double (*funcDoubleToDouble)( double );
	ptr_funDoubleToDouble( funcDoubleToDouble ptrF ) : itsPFunc(ptrF ) {}
	virtual double operator()(double x) const {return (*itsPFunc)(x); }

private:
	funcDoubleToDouble itsPFunc;
};

// Inverse Class using Brent algorithm
class Inverse
{
public: 

	// Keywords to choose the bracketing strategy
	enum {REAL, ALWAYSPOSITIVE, CORRELATION, BOUNDEDBY0AND1};
	
	// Contructor with C++ functior
	Inverse( DoubleToDoubleFunc& f , int YConstraints=REAL)
		: itspFunc( &f ), hasCreatedObject(false),imageConstraints(YConstraints){}

	// Constructor with C functor
	typedef double (*funcDoubleToDouble)( double );
	Inverse(funcDoubleToDouble const ptrF, int YConstraints=REAL)
			: itspFunc( NULL ), hasCreatedObject(false),imageConstraints(YConstraints)
	{
		// Create C functor wrapper
		itspFunc = new ptr_funDoubleToDouble(ptrF);
		hasCreatedObject = true;
	}

	// Destructor
	virtual ~Inverse()
	{
		// Destroy the C functor wrapper
		if(hasCreatedObject)
			delete itspFunc;
	}

	// Eval Operator
	double operator() (double arg,double typical_x=0.0, double delta_x=0.1, double acc=3.0e-13);

private:
	DoubleToDoubleFunc* itspFunc;
	bool hasCreatedObject;
	int imageConstraints;
};

double returnImage(double x, int ImageConst);

double brentSolve(const DoubleToDoubleFunc& func, double ygoal, double lowerBound, double upperBound, double tol, int MAXITER = ARM_CF_MAXIT, int ImageConst = 0, double * best = NULL);


#undef ARM_CF_EPS
CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

