
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Bi Log Model Visitor													 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFBILOGVISITOR_H
#define _INGPINFLATION_INFBILOGVISITOR_H


#include "gpinflation\infmodelvisitor.h"

#pragma warning(disable : 4786) 


CC_BEGIN_NAMESPACE( ARM )


//=======> ARM_InfNumBiLog

class ARM_InfNumBiLog:  public ARM_InfIntModel{

public:
	ARM_InfNumBiLog	( ):ARM_InfIntModel(){}
	ARM_InfNumBiLog	(	const string & name, const double & dis, const double & dom, const double & eps, const double & cen):ARM_InfIntModel(name, dis, dom, eps, cen){}
	ARM_InfNumBiLog	(	const ARM_InfNumBiLog & rhs):ARM_InfIntModel(rhs){}
	virtual ARM_Object* Clone() const { return new ARM_InfNumBiLog(*this); }
	virtual ~ARM_InfNumBiLog(){};
	ASSIGN_OPERATOR	( ARM_InfNumBiLog );

 	void Accept ( ARM_AcyclicVisitor &  );
};

inline void ARM_InfNumBiLog::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfNumBiLog>* tmp =dynamic_cast< ARM_Visitor<ARM_InfNumBiLog>* > (&visitor);
	tmp->Visit( *this ) ;
}

//=======> ARM_InfBiLog

class ARM_InfBiLog:  public ARM_InfNumBiLog{

public:
	ARM_InfBiLog	( ):ARM_InfNumBiLog(){}
	ARM_InfBiLog	(	const string & name, const double & dis, const double & dom, const double & eps, const double & cen):ARM_InfNumBiLog(name, dis, dom, eps, cen){}
	ARM_InfBiLog	(	const ARM_InfBiLog & rhs):ARM_InfNumBiLog(rhs){}
	virtual ARM_Object* Clone() const { return new ARM_InfBiLog(*this); }
	virtual ~ARM_InfBiLog(){};
	ASSIGN_OPERATOR	( ARM_InfBiLog );

 	void Accept ( ARM_AcyclicVisitor &  );
};

inline void ARM_InfBiLog::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfBiLog>* tmp =dynamic_cast< ARM_Visitor<ARM_InfBiLog>* > (&visitor);
	tmp->Visit( *this ) ;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfModelValue
///	Action : Class for Inflation Valorisation Models

/////////////////////////////////////////////////////


//=======> ARM_InfNumBiLogValue

class ARM_InfNumBiLogValue:public ARM_InfModelValue{

public:
	ARM_InfNumBiLogValue():ARM_InfModelValue(){}
	ARM_InfNumBiLogValue( const ARM_InfNumBiLog & model){	itsModel = model;	}
	virtual ~ARM_InfNumBiLogValue(){}

	virtual double operator() (	const  ARM_InfPayOffValuePtr	& payOff, 	
								const  double					& res,
								const  ARM_MAP_Double			& fwd,
								const  ARM_MAP_PairDb			& cor,								
								const  ARM_MAP_VolPar			& vol );	

protected:
	virtual void ValidatePayOff( const ARM_InfPayOffValuePtr	& payOff){}

	virtual double CptVol	 (	const int						& cri,		// critère
								const ARM_InfIrIndex			& idx,
								const double					& res, 
								const double					& ten, 
								const double					& fwd, 
								const double					& str,
								const double					& eps);
	
	virtual double	CptStrike(	const string & str, 
								const double & lag, 
								const ARM_InfPayOffValuePtr & payOff, 
								const ARM_MAP_Double & fwd);

private:
	ARM_InfNumBiLog itsModel;	
};

//=======> ARM_InfBiLogValue

class ARM_InfBiLogValue:public ARM_InfNumBiLogValue{

public:
	ARM_InfBiLogValue():ARM_InfNumBiLogValue(){}
	ARM_InfBiLogValue( const ARM_InfBiLog & model){	itsModel = model;	}
	virtual ~ARM_InfBiLogValue(){}

	virtual double operator() (	const ARM_InfPayOffValuePtr		& payOff, 	
								const double					& res,
								const ARM_MAP_Double			& fwd,
								const ARM_MAP_PairDb			& cor,								
								const ARM_MAP_VolPar			& vol );	
protected:
	virtual void ValidatePayOff( const ARM_InfPayOffValuePtr	& payOff);

private:
	ARM_InfBiLog itsModel;	
};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















