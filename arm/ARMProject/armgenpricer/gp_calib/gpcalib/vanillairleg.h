/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairleg.h
 *
 *  \brief vanilla interest rates leg
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPCALIB_VANILLAIRLEG_H
#define _INGPCALIB_VANILLAIRLEG_H

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "vanillaarg.h"


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaIRSwaplegArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaIRSwaplegArg:public ARM_VanillaArg
{
    ARM_VanillaIRSwaplegArg(
		const string& curveName,
		double evalTime,
		int CallPut,
		ARM_GP_Vector* nominals,
		ARM_GP_Vector* resetTimes,
		ARM_GP_Vector* startTimes,
		ARM_GP_Vector* endTimes,
		ARM_GP_Vector* strikes,
		ARM_GP_Vector* payTimes,
		ARM_GP_Vector* payPeriods,
		double expiry	= 0.0,
		double tenor	= 0.0,
		double Spread = 0.,
		int DecompFreq = 1)
	:
		ARM_VanillaArg(curveName, evalTime, CallPut, expiry, 0, Spread, DecompFreq ),
		itsNominals(nominals),
		itsResetTimes(resetTimes),
		itsStartTimes(startTimes),
		itsEndTimes(endTimes),
		itsPayTimes(payTimes),
		itsPayPeriods(payPeriods)

	{}	
		
    ARM_VanillaIRSwaplegArg(const ARM_VanillaIRSwaplegArg& arg);
    ARM_VanillaIRSwaplegArg& operator=(const ARM_VanillaIRSwaplegArg& rhs);
    virtual ~ARM_VanillaIRSwaplegArg();
    virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const ;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_IRLEG; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors (const and non const version)
	inline const ARM_GP_Vector* GetResetTimes() const { return itsResetTimes; }
	inline const ARM_GP_Vector* GetPayTimes() const { return itsPayTimes; }
	inline const ARM_GP_Vector* GetStartTimes() const { return itsStartTimes; }
	inline const ARM_GP_Vector* GetEndTimes() const { return itsEndTimes; }
	inline const ARM_GP_Vector* GetPayPeriods() const { return itsPayPeriods; }

	inline ARM_GP_Vector* GetResetTimes()  { return itsResetTimes; }
	inline ARM_GP_Vector* GetPayTimes()  { return itsPayTimes; }
	inline ARM_GP_Vector* GetStartTimes()  { return itsStartTimes; }
	inline ARM_GP_Vector* GetEndTimes()  { return itsEndTimes; }
	inline ARM_GP_Vector* GetPayPeriods()  { return itsPayPeriods; }

	
    virtual double PriceOplet(
		ARM_PricingModelIR* model,
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

private :
	void CopyNoCleanUp(const ARM_VanillaIRSwaplegArg& rhs);
	void CleanUp();

    ARM_GP_Vector* itsNominals;
	ARM_GP_Vector* itsResetTimes; 
	ARM_GP_Vector* itsStartTimes;
    ARM_GP_Vector* itsEndTimes;
    ARM_GP_Vector* itsPayTimes; 
	ARM_GP_Vector* itsPayPeriods;

};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
