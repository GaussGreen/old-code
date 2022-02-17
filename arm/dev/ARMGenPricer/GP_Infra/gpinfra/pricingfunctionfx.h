/*!
 *
 * Copyright (c) IXIS CIB 2003 Paris
 *
 *	\file pricingfunctionfx.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2005
 */

#ifndef _INGPINFRA_PRICINGFUNCTIONFX_H
#define _INGPINFRA_PRICINGFUNCTIONFX_H
 

#include "gpbase/env.h"
#include "typedef.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )
/// macro for namespace ... define namespace only if supported

///////////////////////////////////////////////////////
/// \class ARM_PricingFunctionFx
/// \brief
/// This abstract class is the interface for fx model
///////////////////////////////////////////////////////
class ARM_PricingFunctionFx : public ARM_PricingFunctionEquity
{
private:
	ARM_PricingModelPtr itsIRForModel;

public:
	ARM_PricingFunctionFx():ARM_PricingFunctionEquity(), itsIRModel(NULL) {};
	
	ARM_PricingFunctionFx(const ARM_PricingFunctionFx& rhs)
	: ARM_PricingFunctionEquity(rhs), itsIRForModel(rhs.itsIRForModel){};
	
	ARM_PricingFunctionFx& operator=( const ARM_PricingFunctionFx& rhs)
	{
		if( this!= & rhs )
		{
			ARM_PricingFunctionEquity::operator =(rhs);
			itsIRForModel =rhs.itsIRForModel;
		}
		return *this;
	}
	virtual ~ARM_PricingFunctionFx() {};

	/// for multi asset type fx model!
	virtual void SetIRForeignModel( const ARM_PricingModelPtr& irModel )
	{ ARM_THROW( ERR_INVALID_ARGUMENT, ": unimplemented method 'SetIRForeignModel'" ); }
	inline ARM_PricingModelPtr GetIRForeignModel() const { return itsIRForModel; }
};


CC_END_NAMESPACE()

#endif

