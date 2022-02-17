/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillaequityoption.h
 *
 *  \brief vanilla equity options
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPCALIB_VANILLAEQUITYOPTION_H
#define _INGPCALIB_VANILLAEQUITYOPTION_H

#include "gpbase/env.h"
#include "vanillaeqfxoption.h"

CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaEqOption
/// \brief common vanilla arg for fx option
///////////////////////////////////////////////////////////////
struct ARM_VanillaEqOption: public ARM_VanillaEqFxOption
{
	ARM_VanillaEqOption(
		const string& curveName,
		double maturityTime,
		double fwdTime,
        double payTime,
		double strike,
        int callPut )
	: ARM_VanillaEqFxOption(curveName,maturityTime,fwdTime,payTime,strike,callPut) {}
    ARM_VanillaEqOption(const ARM_VanillaEqOption& rhs)
		:	ARM_VanillaEqFxOption( rhs ) {}
    ARM_VanillaEqOption& operator=(const ARM_VanillaEqOption& rhs)
	{
		if( this != &rhs )
		{
			ARM_VanillaEqFxOption::operator =( rhs);
		}
		return *this;
	}
	virtual ~ARM_VanillaEqOption(){};

	virtual double Price(ARM_PricingModel* model ) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const ;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_EQUITYOPTION; }


	/// standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaEqOption(*this);}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
