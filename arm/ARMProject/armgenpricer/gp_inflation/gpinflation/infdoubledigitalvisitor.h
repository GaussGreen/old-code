
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Double Digital Visitor												 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFDOUBLEDIGITALVISITOR_H
#define _INGPINFLATION_INFDOUBLEDIGITALVISITOR_H


#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation\infvisitor.h"


CC_BEGIN_NAMESPACE( ARM )


//========> ARM_InfDoubleDigit
class ARM_InfDoubleDigit:  public ARM_InfHybridPayOff{

public:
	ARM_InfDoubleDigit	( ):ARM_InfHybridPayOff( ){};
	ARM_InfDoubleDigit	(	const ARM_MAP_Curve &, const ARM_MAP_Curve &, const ARM_MAP_Curve &  );
	ARM_InfDoubleDigit	(	const ARM_InfDoubleDigit & rhs);
	
	virtual string		toString(	const string& indent="", const string& nextIndent="") const;
	virtual ~ARM_InfDoubleDigit(){};
	
	ASSIGN_OPERATOR		( ARM_InfDoubleDigit		);

	virtual ARM_Object* Clone	( ) const { return new ARM_InfDoubleDigit(*this); }
 	virtual void		Accept  ( ARM_AcyclicVisitor &  );
	virtual void		ValidateSecurity ( const vector<string> & );

	ARM_MAP_Curve		GetM_OptCurve () const{ return itsM_OptCurve; }
	ARM_MAP_Curve		GetS_OptCurve () const{ return itsS_OptCurve; }

protected:

	ARM_MAP_Curve		itsM_OptCurve;
	ARM_MAP_Curve		itsS_OptCurve;

};


inline void ARM_InfDoubleDigit::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfDoubleDigit>* tmp =dynamic_cast< ARM_Visitor<ARM_InfDoubleDigit>* > (&visitor);
	tmp->Visit( *this ) ;
}

//========> ARM_InfDoubleDigitValue
class ARM_InfDoubleDigitValue:public ARM_InfHybridPayOffValue{

public:
	ARM_InfDoubleDigitValue(): ARM_InfHybridPayOffValue(){}
	ARM_InfDoubleDigitValue( const ARM_InfDoubleDigit & payOff):ARM_InfHybridPayOffValue( ){itsPayOff = ARM_InfHybridPayOffPtr(new ARM_InfDoubleDigit(payOff) );}
	virtual ~ARM_InfDoubleDigitValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double &  );

};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















