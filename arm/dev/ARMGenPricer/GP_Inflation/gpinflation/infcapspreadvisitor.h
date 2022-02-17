
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


#include "infvisitor.h"
#include <gpbase\typedef.h>
#include <map>
#include <string>


CC_USING_NS( std, map )
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )

typedef map<string, ARM_GP_CurvePtr >			ARM_MAP_Curve;
typedef ARM_MAP_Curve::const_iterator			Iter_Curve;
typedef map<string, double>						ARM_MAP_Double;
typedef ARM_CountedPtr<ARM_MAP_Double> 			ARM_MAP_DbPtr;

class ARM_InfPayOffValue{
public:
	ARM_InfPayOffValue() {}
	virtual ~ARM_InfPayOffValue(){}
	virtual double operator() (	const  double			& , 
								const  ARM_MAP_Double	& , 
								const  ARM_MAP_Double	&   )=0;
};

typedef ARM_CountedPtr<ARM_InfPayOffValue>		ARM_InfPayOffValuePtr;


class ARM_InfHybridPayOff:  public ARM_InfPayOff{

public:
	ARM_InfHybridPayOff	():ARM_InfPayOff(){};
	ARM_InfHybridPayOff	(	const ARM_MAP_Curve &,		const ARM_MAP_Curve &);
	ARM_InfHybridPayOff	(	const ARM_InfHybridPayOff & );
	virtual ~ARM_InfHybridPayOff(){};
	
	ASSIGN_OPERATOR		( ARM_InfHybridPayOff			);

	virtual ARM_Object* Clone	() const { return new ARM_InfHybridPayOff(*this); }
 	virtual string		toString(	const string& indent="", const string& nextIndent="") const;
	virtual void		Accept ( ARM_AcyclicVisitor &  );

	ARM_MAP_Curve		GetOptCurve () const{ return itsOptCurve; }
	ARM_MAP_Curve		GetCpnCurve () const{ return itsCpnCurve; }

protected:

	ARM_MAP_Curve	itsOptCurve;
	ARM_MAP_Curve	itsCpnCurve;
};

inline void ARM_InfHybridPayOff::Accept ( ARM_AcyclicVisitor & visitor ){   
		ARM_Visitor<ARM_InfHybridPayOff>* tmp =dynamic_cast< ARM_Visitor<ARM_InfHybridPayOff>* > (&visitor);
		tmp->Visit( *this ) ;
}

class ARM_InfHybridCap:  public ARM_InfHybridPayOff{

public:
	ARM_InfHybridCap	( ):ARM_InfHybridPayOff(){};
	ARM_InfHybridCap	(	const ARM_MAP_Curve &);
	ARM_InfHybridCap	(	const ARM_InfHybridCap & rhs):ARM_InfHybridPayOff(rhs){}
	virtual ~ARM_InfHybridCap(){};
	
	ASSIGN_OPERATOR		( ARM_InfHybridCap		);

	virtual ARM_Object* Clone	() const { return new ARM_InfHybridCap(*this); }
 	virtual void		Accept ( ARM_AcyclicVisitor &  );
};


inline void ARM_InfHybridCap::Accept ( ARM_AcyclicVisitor & visitor ){   
		ARM_Visitor<ARM_InfHybridCap>* tmp =dynamic_cast< ARM_Visitor<ARM_InfHybridCap>* > (&visitor);
		tmp->Visit( *this ) ;
}



class ARM_InfHybridPayOffValue:public ARM_InfPayOffValue{

public:
	ARM_InfHybridPayOffValue():ARM_InfPayOffValue(){}
	ARM_InfHybridPayOffValue( const ARM_InfHybridPayOff & payOff){	itsPayOff = payOff;	}
	virtual ~ARM_InfHybridPayOffValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double & 	);
	virtual double CptOpt		(	const double &, const ARM_MAP_Double &  );

protected:
	ARM_InfHybridPayOff itsPayOff;
};

class ARM_InfHybridCapValue:public ARM_InfHybridPayOffValue{

public:
	ARM_InfHybridCapValue(): ARM_InfHybridPayOffValue(){}
	ARM_InfHybridCapValue( const ARM_InfHybridCap & payOff):ARM_InfHybridPayOffValue( payOff ){}
	virtual ~ARM_InfHybridCapValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double &  );

};


class	ARM_InfPayOffVisitor:	public ARM_AcyclicVisitor,
								public ARM_Visitor<ARM_InfHybridPayOff>,
								public ARM_Visitor<ARM_InfHybridCap>{

public:
	ARM_InfPayOffVisitor(){}
	virtual ~ARM_InfPayOffVisitor(){}

	ARM_InfPayOffValuePtr GetInfPayOffValue() const{ return itsInfPayOffValue; }

	void Visit( const ARM_InfHybridPayOff & );
	void Visit( const ARM_InfHybridCap & );

protected:
	ARM_InfPayOffValuePtr itsInfPayOffValue;
};

inline 	void ARM_InfPayOffVisitor::Visit( const ARM_InfHybridPayOff & payOff){
		ARM_InfHybridPayOffValue* tmp = new ARM_InfHybridPayOffValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
}; 

inline 	void ARM_InfPayOffVisitor::Visit( const ARM_InfHybridCap & payOff){
		ARM_InfHybridCapValue* tmp = new ARM_InfHybridCapValue(payOff);
		itsInfPayOffValue	= ARM_InfPayOffValuePtr ( tmp );	
};

CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















