/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file Q1FPricingContext.h
 *
 *  \brief Q1FPricingContext
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPMODELS_Q1FPRINCINCONTEXT_H
#define _INGPMODELS_Q1FPRINCINCONTEXT_H

#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_Q1FPricingStatesContext : public ARM_PricingStatesContext
{
private: 
	ARM_GP_VectorPtr itsForwardValue;
public: 
	ARM_Q1FPricingStatesContext()
		:	ARM_PricingStatesContext(), itsForwardValue(NULL) {}
	ARM_Q1FPricingStatesContext( const ARM_Q1FPricingStatesContext& rhs )
		:	ARM_PricingStatesContext(rhs), itsForwardValue(rhs.itsForwardValue) {}
	ARM_Q1FPricingStatesContext& operator=( const ARM_Q1FPricingStatesContext& rhs )
	{
		if( this != &rhs )
		{
			ARM_PricingStatesContext::operator =( rhs );
			itsForwardValue = rhs.itsForwardValue;
		}
		return *this;
	}

	virtual ARM_Q1FPricingStatesContext * ToQ1FPricingStatesContext() { return this; }
	virtual ARM_Object* Clone() const { return new ARM_Q1FPricingStatesContext(*this); }
	virtual ~ARM_Q1FPricingStatesContext(){};

	inline ARM_GP_VectorPtr GetFwdValue() const { return itsForwardValue; }
	inline void SetFwdValue( const ARM_GP_VectorPtr& forwardValue ) { itsForwardValue = forwardValue; }
	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_Q1FPricingStatesContext"; }
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/