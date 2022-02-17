/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanilladigital.h
 *
 *  \brief vanilla digital
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_VANILLADIGITAL_H
#define _INGPCALIB_VANILLADIGITAL_H

#include "vanillacap.h"


CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaDigitalArg
/// \brief cap/caplet argument simple struct
///////////////////////////////////////////////////////////////

struct ARM_VanillaDigitalArg: public ARM_VanillaCapDigitalArg
{
    ARM_VanillaDigitalArg( 
		const string& curveName,
		double evalTime,
		int CallPut,
		ARM_GP_Vector* nominals,
		ARM_GP_Vector* resetTimes,
		ARM_GP_Vector* startTimes,
		ARM_GP_Vector* endTimes,
		ARM_GP_Vector* strikes,
		ARM_GP_Vector* payTimes,
		ARM_GP_Vector* payPeriods )
	:
		ARM_VanillaCapDigitalArg( 
			curveName,
			evalTime,
			CallPut,
			nominals,
			resetTimes,
			startTimes,
			endTimes,
			strikes,
			payTimes, 
			payPeriods )
	{}

	ARM_VanillaDigitalArg(const ARM_VanillaDigitalArg& rhs) : ARM_VanillaCapDigitalArg(rhs) {};
    ARM_VanillaDigitalArg& operator=(const ARM_VanillaDigitalArg& rhs)
	{
		if( this != &rhs )
			ARM_VanillaCapDigitalArg::operator =(rhs);
		return *this;
	}
    virtual ~ARM_VanillaDigitalArg() {}

    virtual double PriceOplet(
		ARM_PricingFunctionIR* model,
		const string& curveName, 
		double evalTime,
		double payTime,
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        double strike,
        int capFloor,
        const ARM_PricingStatesPtr& states ) const;

    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaDigitalArg(*this);}
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_DIGITAL; }
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
