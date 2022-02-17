/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillafxswaption.h
 *
 *  \brief vanilla ir fx swaption
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date March 2006
 */


#ifndef _INGPCALIB_VANILLAIRFXSWAPTION_H
#define _INGPCALIB_VANILLAIRFXSWAPTION_H

#include "gpbase/env.h"
#include "gpbase/countedptr.h"
#include "gpbase/gpvector.h"

#include "vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_VanillaSwaptionArg;
struct ARM_VanillaFxOption;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaIrFxSwaption
/// \brief Vanilla instrument for hybride IR/FX basket option
///		   We describe a basket made of :
///				1) a forward forex strip with variable nominals
///				2) a variable nominal domestic swaption
///		   and the intrinsic value of the option is the positive
///		   part of the sum (nominals will manage call or put feature)
///		   Strike of the option is implicit in the fixed leg of the VNS
///////////////////////////////////////////////////////////////
struct ARM_VanillaIrFxSwaption:public ARM_VanillaArg
{
	ARM_VanillaIrFxSwaption(const string& curveName,
							double expiryTime,
							int callPut,
							double strike,
							const ARM_GP_Vector& fxResetTimes,
							const ARM_GP_Vector& fxSettlementTimes,
							const ARM_GP_Vector& fxPayTimes,
							const ARM_GP_Vector& fxNominals,
							const ARM_CountedPtr<ARM_VanillaSwaptionArg>& irSwaption)
		: ARM_VanillaArg( curveName, 0.0, callPut, expiryTime ),itsStrike(strike),
		itsFxResetTimes(fxResetTimes), itsFxSettlementTimes(fxSettlementTimes),
		itsFxPayTimes(fxPayTimes), itsFxNominals(fxNominals), itsIrSwaption(irSwaption)
		{}

    ARM_VanillaIrFxSwaption(const ARM_VanillaIrFxSwaption& rhs)
		: ARM_VanillaArg(rhs), itsStrike(rhs.itsStrike),itsFxResetTimes(rhs.itsFxResetTimes), itsFxSettlementTimes(rhs.itsFxSettlementTimes),
		itsFxPayTimes(rhs.itsFxPayTimes), itsFxNominals(rhs.itsFxNominals), itsIrSwaption(rhs.itsIrSwaption)
		{}
		
    ARM_VanillaIrFxSwaption& operator=(const ARM_VanillaIrFxSwaption& rhs)
	{
		if (&rhs != this)
		{ 
			this->~ARM_VanillaIrFxSwaption();
			new (this) ARM_VanillaIrFxSwaption (rhs);
		}
		return *this;
	}
	virtual ~ARM_VanillaIrFxSwaption(); /// empty code is in .cpp to avoid including here .h of ARM_CountedPtr targets

	double GetStrike() { return itsStrike; }
	const ARM_GP_Vector& GetFxResetTimes() const		{ return itsFxResetTimes; }
	const ARM_GP_Vector& GetFxSettlementTimes() const	{ return  itsFxSettlementTimes; }
	const ARM_GP_Vector& GetFxPayTimes() const			{ return  itsFxPayTimes; }
	const ARM_GP_Vector& GetFxNominals() const			{ return  itsFxNominals; }

	const ARM_VanillaSwaptionArg& GetIrSwaption() const { return *itsIrSwaption; }

	virtual double Price(ARM_PricingModel* model ) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_IRFXSWAPTION; }

	/// standard ARM Support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaIrFxSwaption(*this);}


private:

	double itsStrike;

	ARM_GP_Vector itsFxResetTimes;
	ARM_GP_Vector itsFxSettlementTimes;
	ARM_GP_Vector itsFxPayTimes;
	ARM_GP_Vector itsFxNominals;

	ARM_CountedPtr<ARM_VanillaSwaptionArg>	itsIrSwaption;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
