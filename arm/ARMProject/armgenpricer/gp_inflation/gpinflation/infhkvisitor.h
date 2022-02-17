
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Hk Model Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFHKVISITOR_H
#define _INGPINFLATION_INFHKVISITOR_H


#include "gpinflation\infmodelvisitor.h"
#include "gpinflation\hybridinfindex.h"

#pragma warning(disable : 4786) 


CC_BEGIN_NAMESPACE( ARM )



//=======> ARM_InfHK

class ARM_InfHK:  public ARM_InfIntModel{

public:
	ARM_InfHK	( ){}
	ARM_InfHK		(	const string & name, const double & dis, const double & dom, const double & eps, const double & cen);
	ARM_InfHK	(	const ARM_InfHK & rhs):ARM_InfIntModel(rhs){ }
	virtual ARM_Object* Clone() const { return new ARM_InfHK(*this); }
	virtual ~ARM_InfHK(){};
	ASSIGN_OPERATOR	( ARM_InfHK );

 	void Accept ( ARM_AcyclicVisitor &  );
};

inline void ARM_InfHK::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfHK>* tmp =dynamic_cast< ARM_Visitor<ARM_InfHK>* > (&visitor);
	tmp->Visit( *this ) ;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfModelValue
///	Action : Class for Inflation Valorisation Models

/////////////////////////////////////////////////////


//=======> ARM_InfHKValue

class ARM_InfHK;
class ARM_InfHKValue:public ARM_InfModelValue{

public:
	ARM_InfHKValue():ARM_InfModelValue(){}
	ARM_InfHKValue( const ARM_InfHK & model){	itsModel = model;	}
	virtual ~ARM_InfHKValue(){}

	virtual double operator() ( const ARM_InfPayOffValuePtr		& payOff, 	
								const double					& res,
								const ARM_MAP_Double			& fwd,
								const ARM_MAP_PairDb			& cor,								
								const ARM_MAP_VolPar			& vol );

protected:
	virtual void ValidatePayOff(const ARM_InfPayOffValuePtr		& payOff){}

private:
	ARM_InfHK itsModel;
};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















