/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file SFRMPricingStates.cpp
 *
 *  \brief SFRMPricingStates
 *	\author  R.Guillemot, A.Schauly
 *	\version 1.0
 *	\date January 2005
 */

#include "gpbase/env.h"
#include "gpmodels/SFRMPricingStates.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SFRMPricingStatesContext::ARM_SFRMPricingStatesContext()
{
	itsFwdMinIndex = 0;
	itsNumeraireTimeIndex = 0;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: Constructor
///	Returns: 
///	Action : COPYConstructor
////////////////////////////////////////////////////
ARM_SFRMPricingStatesContext::ARM_SFRMPricingStatesContext( const ARM_SFRMPricingStatesContext& rhs ) : itsFwdMinIndex(rhs.itsFwdMinIndex),
			itsNumeraireTimeIndex(rhs.itsNumeraireTimeIndex)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: Operator=
///	Returns: 
///	Action : Operator=
////////////////////////////////////////////////////

ARM_SFRMPricingStatesContext& ARM_SFRMPricingStatesContext::operator=( const ARM_SFRMPricingStatesContext& rhs )
{
	if( this!= &rhs )
	{
		ARM_PricingStatesContext::operator = (rhs);
		itsFwdMinIndex = rhs.itsFwdMinIndex;
		itsNumeraireTimeIndex = rhs.itsNumeraireTimeIndex;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: Clone()
///	Returns: ARM_SFRMPricingStatesContext*
///	Action : Clone
////////////////////////////////////////////////////

ARM_Object* ARM_SFRMPricingStatesContext::Clone() const
{
	return new ARM_SFRMPricingStatesContext(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_SFRMPricingStatesContext::~ARM_SFRMPricingStatesContext()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRMPricingStatesContext
///	Routine: ToSFRMPricingStatesContext
///	Returns: ARM_SFRMPricingStatesContext *
///	Action : casts into ARM_SFRMPricingStatesContext *
////////////////////////////////////////////////////

ARM_SFRMPricingStatesContext * ARM_SFRMPricingStatesContext::ToSFRMPricingStatesContext() 
{ 
	return this;
}

CC_END_NAMESPACE()