/*!
 *
 * Copyright (c) IXIS CIB June 2007 Paris
 *
 *  function solution of the Riccati Equation:
 *
 *	dx/dt = alpha*x(t)² + (beta0 + beta1*exp(lambda*t) )*x(t) + gamma0+gamma1*exp(lambda*t)+gamma2**exp(2*lambda*t)
 *
 *  x(t0)=x0
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
							const Complexe&			beta0,
							const Complexe&			beta1,
							const Complexe&			gamma0,
							const Complexe&			gamma1,
							const Complexe&			gamma2,
							const Complexe&			lambda,
							const Complexe&			x0,
							const Complexe&			t0);

	ARM_ExpRiccati(const ARM_ExpRiccati& rhs):	ARM_ODEFunc (rhs),
												itsLambda	(rhs.itsLambda	),
												itsAlpha	(rhs.itsAlpha	),
												itsBeta0	(rhs.itsBeta0	),
												itsBeta1	(rhs.itsBeta1	),
												itsMu		(rhs.itsMu		),
												itsNu		(rhs.itsNu		),
												itsEta		(rhs.itsEta		),
												itsConstant	(rhs.itsConstant){}
	virtual ~ARM_ExpRiccati(){};
	
	//ASSIGN_OPERATOR	( ARM_ExpRiccati );
	virtual ARM_Object* Clone() const{ return new ARM_ExpRiccati( *this); }
	virtual string toString(const string& indent="",const string& nextIndent="") const{ return ""; }
	virtual void derivs(double x, std::vector<double>* yt, std::vector<double>* dyt) const;	
	virtual string	ExportShortName() const { return "LRICC";}

	/// value of the primitive
    virtual Complexe operator () (double t, double T) const;

    /// value of the solution
    virtual Complexe operator [] (double t) const;

public : 
		
	Complexe	itsLambda;
	Complexe	itsAlpha;
	Complexe	itsBeta0;
	Complexe	itsBeta1;
	Complexe	itsMu;
	Complexe	itsNu;
	Complexe	itsEta;

	Complexe	itsConstant;
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/