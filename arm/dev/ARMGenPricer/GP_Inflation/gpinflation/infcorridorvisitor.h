
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Corridor Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: July, 5th 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFCORRIDORVISITOR_H
#define _INGPINFLATION_INFCORRIDORVISITOR_H

#pragma warning(disable : 4786) 


#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation\infvisitor.h"

CC_USING_NS( std, map )
CC_USING_NS( std, string )

CC_BEGIN_NAMESPACE( ARM )


//========> ARM_InfCorridor
class ARM_InfCorridor:  public ARM_InfHybridPayOff{

public:
	ARM_InfCorridor	( ):ARM_InfHybridPayOff( ){};
	ARM_InfCorridor	(	const ARM_MAP_Curve &, const ARM_MAP_Curve &  );
	ARM_InfCorridor	(	const ARM_InfCorridor & );
	
	virtual string		toString(	const string& indent="", const string& nextIndent="") const;
	virtual ~ARM_InfCorridor(){};
	
	ASSIGN_OPERATOR		( ARM_InfCorridor		);

	virtual ARM_Object* Clone	( ) const { return new ARM_InfCorridor(*this); }
 	virtual void		Accept  ( ARM_AcyclicVisitor &  );

//	void ValidateSecurity ( const vector<string> & );

};


inline void ARM_InfCorridor::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfCorridor>* tmp =dynamic_cast< ARM_Visitor<ARM_InfCorridor>* > (&visitor);
	tmp->Visit( *this ) ;
}

//========> ARM_InfDoubleDigitValue
class ARM_InfCorridorValue:public ARM_InfHybridPayOffValue{

public:
	ARM_InfCorridorValue(): ARM_InfHybridPayOffValue(){}
	ARM_InfCorridorValue( const ARM_InfCorridor & payOff):ARM_InfHybridPayOffValue( ){itsPayOff = ARM_InfHybridPayOffPtr(new ARM_InfCorridor(payOff) );}
	virtual ~ARM_InfCorridorValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double &  );

};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















