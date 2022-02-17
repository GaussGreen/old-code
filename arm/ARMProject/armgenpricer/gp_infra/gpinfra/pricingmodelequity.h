/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodelequity.h
 *
 *  \brief equity model version of the pricing model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */

#ifndef _INGPINFRA_PRICINGMODELEQUITY_H
#define _INGPINFRA_PRICINGMODELEQUITY_H
 

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "pricingmodel.h"
#include "pricingfunctionequity.h"


CC_BEGIN_NAMESPACE( ARM )
 /// macro for namespace ... define namespace only if supported


///////////////////////////////////////////////////////
/// \class ARM_PricingModelEquity
/// \brief
/// This abstract class is the standard
/// interface for equity pricing models
/// Derived from the ARM_PricingModel
///////////////////////////////////////////////////////
class ARM_PricingModelEquity :  public ARM_PricingModel,
                                public ARM_PricingFunctionEquity
{
public:
	ARM_PricingModelEquity(const ARM_ZeroCurvePtr& zc=ARM_ZeroCurvePtr(NULL), const ARM_ModelParams* params=NULL  ); 
	ARM_PricingModelEquity(const ARM_PricingModelEquity& rhs);
	virtual ~ARM_PricingModelEquity();
    ARM_PricingModelEquity& operator = (const ARM_PricingModelEquity& rhs);

	/// equity type model
	virtual int GetType() const;

    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_PricingModelEquity : abstract class !"); }
};


CC_END_NAMESPACE()

#endif

