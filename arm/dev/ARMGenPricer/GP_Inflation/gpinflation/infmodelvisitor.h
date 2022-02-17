
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Model Visitor															 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFMODELVISITOR_H
#define _INGPINFLATION_INFMODELVISITOR_H


#include "gpinflation\infvisitor.h"
#include "gpinflation\hybridinfindex.h"

#pragma warning(disable : 4786) 


CC_BEGIN_NAMESPACE( ARM )

struct InfIrIndex;
class ARM_InfPayOffValue;

typedef map<string, InfIrIndex >				ARM_MAP_Index;
typedef map<string, ARM_VolParam>				ARM_MAP_VolPar;
typedef map<pair<string, string>, double>		ARM_MAP_PairDb;

/////////////////////////////////////////////////////

///	Class  : ARM_InfModelValue
///	Action : Generic Class for Inflation Valorisation Models

/////////////////////////////////////////////////////

//=======> ARM_InfModelValue

class ARM_InfModelValue{
public:
	ARM_InfModelValue() {}
	virtual ~ARM_InfModelValue(){}
	virtual double operator() (	const ARM_InfPayOffValuePtr	& payOff, 	
								const double					& res,
								const ARM_MAP_Double			& fwd,	
								const ARM_MAP_PairDb			& cor,								
								const ARM_MAP_VolPar			& vol ) =0;

	virtual void ValidatePayOff(const ARM_InfPayOffValuePtr		& )=0;
};


typedef ARM_CountedPtr<ARM_InfModelValue>		ARM_InfModelValuePtr;



/////////////////////////////////////////////////////

///	Class  : ARM_InfModelValue
///	Action : Generic Class for Inflation Models (contain information)

/////////////////////////////////////////////////////

//=======> ARM_InfIntModel

class ARM_InfIntModel:  public ARM_InfModel{

public:
	ARM_InfIntModel ( ){ };
	ARM_InfIntModel	(	const string & name, const double & dis, const double & dom, const double & eps, const double & cen);
	ARM_InfIntModel	(	const ARM_InfIntModel & rhs);
	virtual ~ARM_InfIntModel(){};
 	virtual string toString(	const string& indent="", const string& nextIndent="") const;

	ARM_GP_VectorPtr	GetPosit() const { return itsPosit; }
	ARM_GP_VectorPtr	GetWeigt() const { return itsWeigt; }

	double GetDiscretization()	const { return itsDis;	}		// Number of integral discretization
	double GetDomaine()			const { return itsDom;	}
	double GetEpsilon()			const { return itsEps;	}
	double GetCenter()			const { return itsCen;	}

protected:
	double itsDis;		// Number of integral discretization
	double itsDom;		// Domain of integration
	double itsEps;		// epsilon for call spread
	double itsCen;		// Center of the domain

	ARM_GP_VectorPtr	itsPosit;
	ARM_GP_VectorPtr	itsWeigt;

};

//=======> ARM_InfNoModel

class ARM_InfNoModel:  public ARM_InfModel{

public:
	ARM_InfNoModel	(	):ARM_InfModel("NoModel"){}
	ARM_InfNoModel	(	const ARM_InfNoModel & rhs):ARM_InfModel(rhs){}
	virtual ~ARM_InfNoModel(){};
 	virtual string toString(	const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_InfNoModel(*this); }
	ASSIGN_OPERATOR	( ARM_InfNoModel );

 	void Accept ( ARM_AcyclicVisitor &  );
};

inline void ARM_InfNoModel::Accept ( ARM_AcyclicVisitor & visitor ){   
	ARM_Visitor<ARM_InfNoModel>* tmp =dynamic_cast< ARM_Visitor<ARM_InfNoModel>* > (&visitor);
	tmp->Visit( *this ) ;
}



/////////////////////////////////////////////////////

///	Class  : ARM_InfModelValue
///	Action : Class for Inflation Valorisation Models

/////////////////////////////////////////////////////

//=======> ARM_InfNoValue

class ARM_InfNoModelValue:public ARM_InfModelValue{

public:
	ARM_InfNoModelValue():ARM_InfModelValue(){}
	ARM_InfNoModelValue( const ARM_InfNoModel & model){	itsModel = model;	}
	virtual ~ARM_InfNoModelValue(){}

	virtual double operator() (	const ARM_InfPayOffValuePtr		& payOff, 	
								const  double					& res,
								const  ARM_MAP_Double			& fwd,
								const  ARM_MAP_PairDb			& cor,								
								const  ARM_MAP_VolPar			& vol );

	virtual void ValidatePayOff(const ARM_InfPayOffValuePtr		& ){ }

private:
	ARM_InfNoModel itsModel;
};




/////////////////////////////////////////////////////

///	Class  : ARM_InfModelVisitor
///	Action : Class for Inflation Models Visitor

/////////////////////////////////////////////////////

class	ARM_InfNumBiLog;
class	ARM_InfBiLog;
class	ARM_InfHK;

class	ARM_InfModelVisitor:	public ARM_AcyclicVisitor,
								public ARM_Visitor<ARM_InfNumBiLog>,
								public ARM_Visitor<ARM_InfBiLog>,
								public ARM_Visitor<ARM_InfHK>,
								public ARM_Visitor<ARM_InfNoModel>{

public:
	ARM_InfModelVisitor(){}
	virtual ~ARM_InfModelVisitor(){}

	ARM_InfModelValuePtr GetInfModelValue() const{ return itsInfModelValue; }

	void Visit( const ARM_InfNumBiLog	& );
	void Visit( const ARM_InfHK			& );
	void Visit( const ARM_InfNoModel	& );
	void Visit( const ARM_InfBiLog		& );

protected:
	ARM_InfModelValuePtr itsInfModelValue;
};


CC_END_NAMESPACE()

#endif





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















