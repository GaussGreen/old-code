
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Visitor																 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFVISITOR_H
#define _INGPINFLATION_INFVISITOR_H

#include <gpbase\rootobject.h>
#include <gpbase\typedef.h>
#include <gpbase\gpvector.h>

#include <gpinfra/gramfunctorargdict.h>
#include <gpinfra/gramfunctorarg.h>

#include <util/patterns/visitor.hpp>

CC_BEGIN_NAMESPACE( ARM )

class ARM_InfPricer: public ARM_RootObject{
public:
	ARM_InfPricer(){}
	virtual ~ARM_InfPricer(){};
	virtual void	Compute()=0;

	ARM_GramFctorArgDict GetFunctor() const { return itsFunctor;}

protected:

	ARM_GramFctorArgDict	itsFunctor;
};


class ARM_InfType: public ARM_RootObject{
public:
	ARM_InfType():ARM_RootObject() {}
	virtual ~ARM_InfType(){}
};


class ARM_InfPayOff: public ARM_InfType{

public:
	ARM_InfPayOff ( ):ARM_InfType(){}
	ARM_InfPayOff ( const ARM_InfPayOff & ){}
	virtual ~ARM_InfPayOff(){}
	
	virtual void	Accept ( ARM_AcyclicVisitor &  )=0;
	virtual string	ExportShortName() const { return "LINFP";}
 	virtual string	toString(	const string& indent="", const string& nextIndent="") const =0;
	virtual void	ValidateSecurity ( const vector<string> & )=0;
};

class ARM_InfModel: public ARM_InfType{

public:
	ARM_InfModel ( ){}
	ARM_InfModel ( const string & name):ARM_InfType(),itsName(name) {}
	ARM_InfModel ( const ARM_InfModel & rhs ):itsName(rhs.itsName){}
	virtual ~ARM_InfModel(){}
	
	virtual void	Accept ( ARM_AcyclicVisitor &  )=0;
	virtual string	GetName() const { return itsName; }
	virtual string	ExportShortName() const { return "LINFM";}
 	virtual string	toString(	const string& indent="", const string& nextIndent="") const =0;

protected:
	string itsName;
};


CC_END_NAMESPACE()

#endif














/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















