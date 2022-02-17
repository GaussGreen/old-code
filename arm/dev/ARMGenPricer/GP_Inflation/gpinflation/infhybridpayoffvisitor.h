
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Hybrid Payoff Visitor													 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFHYBRIDPAYOFFVISITOR_H
#define _INGPINFLATION_INFHYBRIDPAYOFFVISITOR_H

#pragma warning(disable : 4786) 

#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation\infpayoffvisitor.h"
#include "gpinflation\infvisitor.h"


#include <map>
#include <string>
#include <vector>



CC_BEGIN_NAMESPACE( ARM )

CC_USING_NS( std, map )
CC_USING_NS( std, string )
CC_USING_NS( std, vector )


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

	virtual void		ValidateSecurity ( const vector<string> & );

	ARM_MAP_Curve		GetOptCurve () const{ return itsOptCurve; }
	ARM_MAP_Curve		GetCpnCurve () const{ return itsCpnCurve; }

	ARM_MAP_VolType		GetVolType  () const{ return itsVolType;  }
	void				SetVolType  (  const ARM_MAP_VolType & map) { itsVolType = map;  }

protected:

	ARM_MAP_Curve	itsOptCurve;
	ARM_MAP_Curve	itsCpnCurve;
	ARM_MAP_VolType itsVolType;
};

inline void ARM_InfHybridPayOff::Accept ( ARM_AcyclicVisitor & visitor ){   
		ARM_Visitor<ARM_InfHybridPayOff>* tmp =dynamic_cast< ARM_Visitor<ARM_InfHybridPayOff>* > (&visitor);
		tmp->Visit( *this ) ;
}

typedef ARM_CountedPtr<ARM_InfHybridPayOff> ARM_InfHybridPayOffPtr;

class ARM_InfHybridPayOffValue:public ARM_InfPayOffValue{

public:
	ARM_InfHybridPayOffValue():ARM_InfPayOffValue(){}
	ARM_InfHybridPayOffValue( const ARM_InfHybridPayOff & payOff){	itsPayOff = ARM_InfHybridPayOffPtr(new ARM_InfHybridPayOff(payOff) );	}
	virtual ~ARM_InfHybridPayOffValue(){}

	virtual double operator()	(	const double &, const ARM_MAP_Double &, const ARM_MAP_Double & 	);
	virtual double CptOpt		(	const double &, const ARM_MAP_Double &  );

	virtual VolType operator[]	(	const string & idx){ return itsPayOff->GetVolType()[idx]; }

	ARM_InfHybridPayOffPtr	GetPayOff() const { return itsPayOff; }

protected:
	ARM_InfHybridPayOffPtr itsPayOff;
};

ARM_GP_Vector MergeCurve ( ARM_GP_Vector abs, ARM_GP_CurvePtr cur);


CC_END_NAMESPACE()

#endif


/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















