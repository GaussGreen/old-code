/*!
 *
 * Copyright (c) CDC IXIS CM June 2004 Paris
 *
 *	\file pricingmodelequity.cpp
 *  \brief equity version of a pricing model!
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */


#include "gpinfra/pricingmodelequity.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/pricingmodeltype.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelEquity::ARM_PricingModelEquity( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params )
:	ARM_PricingModel(zc, params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelEquity::ARM_PricingModelEquity(const ARM_PricingModelEquity& rhs)
:	ARM_PricingModel(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PricingModelEquity::~ARM_PricingModelEquity()
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PricingModelEquity& ARM_PricingModelEquity::operator=(const ARM_PricingModelEquity& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModel::operator=(rhs);
        // Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelEquity
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_PricingModelEquity::GetType() const
{
	return MT_EQUITY_MODEL;
}


CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
