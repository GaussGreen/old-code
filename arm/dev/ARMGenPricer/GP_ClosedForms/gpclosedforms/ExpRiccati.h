/*!
 *
 * Copyright (c) IXIS CIB June 2007 Paris
 *
 *  function solution of the Riccati Equation:
 *
 *	dx/dt = d*x(t)² + a *(x0 - exp(-l*t) ) * x(t) + b * ( x+ - exp(-l*t) ) *( x- - exp(-l*t) )
 *
 *	if we put s = a * t, ys=xs/a
 *
 *	dy/dt = y(t)² + 2*epsilon*lambda * (1 - delta*eta/ ( epsilon*(2*epsilon-1) ) *exp(-lambda*t) ) * x(t) + lambda² * delta²(eta²/((2*epsilon-1)²)- 1/4 ) *exp(-2*lambda*t)
 *
 *   f(t0)= f(0)
 *  
 *  a,b,c,d are constants
 *
 *	\file ExpRiccati.h
 *
 *	\author  MBernardo
 */
 
#ifndef _GP_CF_EXPRICCATI_H
#define _GP_CF_EXPRICCATI_H


#include <complex>

#include "gpnumlib/odefunctions.h"
#include "gpbase/rootobject.h"

typedef std::complex<double> Complexe;

CC_BEGIN_NAMESPACE(ARM)


class ARM_ExpRiccati:public ARM_ODEFunc
{

public :

	ARM_ExpRiccati (		const Complexe&			alpha,
							const Complexe&			beta,
							const Complexe&			delta,
							const Complexe&			lambda,
							const Complexe&			x0,
							const Complexe&			x1,
							const Complexe&			x2,
							const Complexe&			y0,
							const Complexe&			t0);

	ARM_ExpRiccati(const ARM_ExpRiccati& rhs):	ARM_ODEFunc(rhs),
												itsLambda	(rhs.itsLambda	),
												itsAlpha	(rhs.itsAlpha	),
												itsBeta		(rhs.itsBeta	),
												itsGamma	(rhs.itsGamma	),
												itsKappa	(rhs.itsKappa	),
												itsEta		(rhs.itsEta		),
												itsSin		(rhs.itsSin		),
												itsCos		(rhs.itsCos		),
												itsX		(rhs.itsX		),
												itsY0		(rhs.itsY0		),
												itsT0		(rhs.itsT0		){}
	virtual ~ARM_ExpRiccati(){}
	
	ASSIGN_OPERATOR	( ARM_ExpRiccati );
	virtual ARM_Object* Clone() const{ return new ARM_ExpRiccati( *this); }
	virtual string toString(const string& indent="",const string& nextIndent="") const{ return ""; }
	virtual void derivs(double x, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;	
	virtual string	ExportShortName() const { return "LRICC";}

	/// value of the primitive
    virtual Complexe operator () (double t, double T) const;

    /// value of the solution
    virtual Complexe operator [] (double t) const;

public : 
		
	Complexe	itsLambda;
	Complexe	itsAlpha;
	Complexe	itsDelta;
	Complexe	itsX;
	Complexe	itsY0;
	Complexe	itsT0;

	Complexe	itsBeta;
	Complexe	itsGamma;
	Complexe	itsKappa;
	Complexe	itsOmega;
	Complexe	itsTheta;
	
	Complexe	itsEta;
	Complexe	itsSin;
	Complexe	itsCos;
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/