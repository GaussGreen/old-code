
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Payoff Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFPAYOFFVISITOR_H
#define _INGPINFLATION_INFPAYOFFVISITOR_H

#pragma warning(disable : 4786) 


#include "gpinflation\infvisitor.h"

#include <gpbase\typedef.h>
#include <gpbase\assignop.h>
#include <map>
#include <string>


CC_USING_NS( std, map )
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )

typedef enum { 	CAP, DIG }	VolType;

typedef map<string, ARM_GP_CurvePtr >			ARM_MAP_Curve;
typedef ARM_MAP_Curve::const_iterator			Iter_Curve;
typedef map<string, VolType>					ARM_MAP_VolType;
typedef map<string, double>						ARM_MAP_Double;
typedef ARM_CountedPtr<ARM_MAP_Double> 			ARM_MAP_DbPtr;



class ARM_InfPayOffValue{
public:
	ARM_InfPayOffValue() {}
	virtual ~ARM_InfPayOffValue(){}
	virtual double operator() (	const  double			& , 
								const  ARM_MAP_Double	& , 
								const  ARM_MAP_Double	&   )=0;

	virtual VolType operator[] (	const  string & )=0;
};

typedef ARM_CountedPtr<ARM_InfPayOffValue>		ARM_InfPayOffValuePtr;



class ARM_InfHybridPayOff;
class ARM_InfHybridCap;
class ARM_InfHybridDigit;
class ARM_InfDoubleDigit;
class ARM_InfCorridor;
	
class	ARM_InfPayOffVisitor:	public ARM_AcyclicVisitor,
								public ARM_Visitor<ARM_InfHybridPayOff>,
								public ARM_Visitor<ARM_InfHybridCap>,
								public ARM_Visitor<ARM_InfHybridDigit>,
								public ARM_Visitor<ARM_InfDoubleDigit>,
								public ARM_Visitor<ARM_InfCorridor>{

public:
	ARM_InfPayOffVisitor(){}
	virtual ~ARM_InfPayOffVisitor(){}

	ARM_InfPayOffValuePtr GetInfPayOffValue() const{ return itsInfPayOffValue; }

	void Visit( ARM_InfHybridPayOff	& );
	void Visit( ARM_InfHybridCap	& );
	void Visit( ARM_InfHybridDigit	& );
	void Visit( ARM_InfDoubleDigit	& );
	void Visit( ARM_InfCorridor		& );

protected:
	ARM_InfPayOffValuePtr itsInfPayOffValue;
};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















