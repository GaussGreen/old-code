
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Option Payoff Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFOPTIONSPREADVISITOR_H
#define _INGPINFLATION_INFOPTIONSPREADVISITOR_H

#pragma warning(disable : 4786) 


#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation\infvisitor.h"

CC_USING_NS( std, map )
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )

//========> ARM_InfHybridCap
class ARM_InfHybridCap:  public ARM_InfHybridPayOff{

public:
	ARM_InfHybridCap	( ):ARM_InfHybridPayOff( ){};
	ARM_InfHybridCap	(	const ARM_MAP_Curve & );
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

//========> ARM_InfHybridDigital
class ARM_InfHybridDigit:  public ARM_InfHybridPayOff{

public:
	ARM_InfHybridDigit	( ):ARM_InfHybridPayOff( ){};
	ARM_InfHybridDigit	(	const ARM_MAP_Curve &, const  ARM_GP_CurvePtr &);
	ARM_InfHybridDigit	(	const ARM_InfHybridDigit & rhs):ARM_InfHybridPayOff(rhs){}
	virtual ~ARM_InfHybridDigit(){};
	
	ASSIGN_OPERATOR		( ARM_InfHybridDigit		);

	virtual ARM_Object* Clone	() const { return new ARM_InfHybridDigit(*this); }
 	virtual void		Accept ( ARM_AcyclicVisitor &  );

};


inline void ARM_InfHybridDigit::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfHybridDigit>* tmp =dynamic_cast< ARM_Visitor<ARM_InfHybridDigit>* > (&visitor);
	tmp->Visit( *this ) ;
}


//========> ARM_InfHybridCapValue
class ARM_InfHybridCapValue:public ARM_InfHybridPayOffValue{

public:
	ARM_InfHybridCapValue(): ARM_InfHybridPayOffValue(){}
	ARM_InfHybridCapValue( const ARM_InfHybridCap & payOff):ARM_InfHybridPayOffValue( payOff ){}
	virtual ~ARM_InfHybridCapValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double &  );

};

//========> ARM_InfHybridDigitValue
class ARM_InfHybridDigitValue:public ARM_InfHybridPayOffValue{

public:
	ARM_InfHybridDigitValue(): ARM_InfHybridPayOffValue(){}
	ARM_InfHybridDigitValue( const ARM_InfHybridDigit & payOff):ARM_InfHybridPayOffValue( payOff ){}
	virtual ~ARM_InfHybridDigitValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double &  );
	
};

CC_END_NAMESPACE()

#endif

/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















