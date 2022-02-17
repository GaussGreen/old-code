/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillafxoption.h
 *
 *  \brief vanilla fx options
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPCALIB_VANILLAEQFXOPTION_H
#define _INGPCALIB_VANILLAEQFXOPTION_H

#include "gpbase/env.h"
#include "vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaEqFxOption
/// \brief common vanilla arg for fx option
///////////////////////////////////////////////////////////////
struct ARM_VanillaEqFxOption: public ARM_VanillaArg
{
	ARM_VanillaEqFxOption(
		const string& curveName,
		double maturityTime,
		double fwdTime,
		double payTime,
		double strike,
        int callPut )
	:	
		ARM_VanillaArg( curveName, 0.0, callPut, maturityTime ),
		itsStrike( strike ),
		itsFwdTime( fwdTime ),
		itsPayTime( payTime )
	{}
    ARM_VanillaEqFxOption(const ARM_VanillaEqFxOption& rhs)
	:	ARM_VanillaArg( rhs ), itsStrike( rhs.itsStrike ), itsFwdTime( rhs.itsFwdTime),itsPayTime( rhs.itsPayTime) {}
    ARM_VanillaEqFxOption& operator=(const ARM_VanillaEqFxOption& rhs)
	{
		if( this != &rhs )
		{
			ARM_VanillaArg::operator =( rhs);
			itsStrike	= rhs.itsStrike;
			itsFwdTime	= rhs.itsFwdTime;
			itsPayTime	= rhs.itsPayTime;
		}
		return *this;
	}
	virtual ~ARM_VanillaEqFxOption(){};

    /// Accessors
    double GetStrike() const { return itsStrike; }
	double GetFwdTime() const { return itsFwdTime; }
	double GetPayTime() const { return itsPayTime; }

	virtual double Price(ARM_PricingModel* model ) const = 0;
    virtual double ImpliedVol(ARM_PricingModel* model) const = 0;
	virtual VanillaArgType GetType() const = 0;

	/// standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const = 0;
	virtual ARM_Object* Clone() const = 0;

private:
	double itsStrike;
	double itsFwdTime;
	double itsPayTime;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
