/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillafxoption.h
 *
 *  \brief vanilla fx options
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPCALIB_VANILLAFXOPTION_H
#define _INGPCALIB_VANILLAFXOPTION_H

#include "gpbase/env.h"
#include "vanillaeqfxoption.h"

CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaFxOption
/// \brief common vanilla arg for fx option
///////////////////////////////////////////////////////////////
struct ARM_VanillaFxOption: public ARM_VanillaEqFxOption
{
	ARM_VanillaFxOption(const string& curveName,
		                double maturityTime,
		                double fwdTime,
		                double payTime,
		                double strike,
                        int callPut,
                        double fixedFx = -1.0)
	: ARM_VanillaEqFxOption(curveName,maturityTime,fwdTime,payTime,strike,callPut),
	itsPayFgnCcy(false) 
    {
        itsFixedFx = fixedFx;
    }

    ARM_VanillaFxOption(const ARM_VanillaFxOption& rhs)
	:	ARM_VanillaEqFxOption( rhs ),
	itsPayFgnCcy(rhs.itsPayFgnCcy) {}
    
    ARM_VanillaFxOption& operator=(const ARM_VanillaFxOption& rhs)
	{
		if( this != &rhs )
		{
			ARM_VanillaEqFxOption::operator =( rhs);
		}

        itsPayFgnCcy = rhs.itsPayFgnCcy;

        itsFixedFx   = rhs.itsFixedFx;

        itsAccountingPricingFlag = rhs.itsAccountingPricingFlag;

		return *this;
	}

	virtual ~ARM_VanillaFxOption(){};

	bool GetPayFgnCcy() const { return itsPayFgnCcy; }
	void SetPayFgnCcy(bool payFgnCcy ) { itsPayFgnCcy = payFgnCcy; }

    void SetFixedFx(double fx)
    {
        itsFixedFx = fx;
    }

    double GetFixedFx(void) const
    {
        return(itsFixedFx);
    }

    void SetDiscPricingMode(int discPricingMode)
    {
        itsAccountingPricingFlag = discPricingMode;
    }

    int GetDiscPricingMode(void)
    {
        return(itsAccountingPricingFlag);
    }

	virtual double Price(ARM_PricingModel* model ) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_FXOPTION; }

	/// standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaFxOption(*this);}

private:

	bool   itsPayFgnCcy;
    
    double itsFixedFx;

    int   itsAccountingPricingFlag;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
