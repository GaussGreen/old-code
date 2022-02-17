/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillacap.h
 *
 *  \brief vanilla cap
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_VANILLACAP_H
#define _INGPCALIB_VANILLACAP_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/typedef.h"

/// gpinfra
#include "gpbase/port.h"
#include "vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingFunctionIR;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaCapDigitalArg
/// \brief common vanilla arg for cap and digital cap
///////////////////////////////////////////////////////////////
struct ARM_VanillaCapDigitalArg: public ARM_VanillaArg
{

    ARM_VanillaCapDigitalArg( 
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
		double Spread = 0.,
		int DecompFreq = 1,
		bool atTheMoneyFlag = false)
	:
		ARM_VanillaArg(curveName, evalTime, CallPut, Spread, DecompFreq ),
		itsNominals(nominals),
		itsResetTimes(resetTimes),
		itsStartTimes(startTimes),
		itsEndTimes(endTimes),
		itsStrikes(strikes),
		itsPayTimes(payTimes),
		itsPayPeriods(payPeriods),
		itsAtTheMoneyFlag(atTheMoneyFlag)
	{}

    ARM_VanillaCapDigitalArg(const ARM_VanillaCapDigitalArg& arg);
    ARM_VanillaCapDigitalArg& operator=(const ARM_VanillaCapDigitalArg& rhs);
    virtual ~ARM_VanillaCapDigitalArg();
    
	virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
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
        const ARM_PricingStatesPtr& states ) const= 0;

	/// Accessors (const and non const version)
	inline const ARM_GP_Vector* GetResetTimes() const   { return itsResetTimes;     }
	inline const ARM_GP_Vector* GetPayTimes() const     { return itsPayTimes;       }
	inline const ARM_GP_Vector* GetStartTimes() const   { return itsStartTimes;     }
	inline const ARM_GP_Vector* GetEndTimes() const     { return itsEndTimes;       }
	inline const ARM_GP_Vector* GetPayPeriods() const   { return itsPayPeriods;     }
    inline const ARM_GP_Vector* GetStrikes() const      { return itsStrikes;        }
    inline const ARM_GP_Vector* GetNotionals() const    { return itsNominals;       }
	inline const bool GetAtTheMoneyFlag() const			{ return itsAtTheMoneyFlag; }

    

	inline ARM_GP_Vector* GetResetTimes()   { return itsResetTimes;     }
	inline ARM_GP_Vector* GetPayTimes()     { return itsPayTimes;       }
	inline ARM_GP_Vector* GetStartTimes()   { return itsStartTimes;     }
	inline ARM_GP_Vector* GetEndTimes()     { return itsEndTimes;       }
	inline ARM_GP_Vector* GetPayPeriods()   { return itsPayPeriods;     }
	inline ARM_GP_Vector* GetStrikes()      { return itsStrikes;        }
    inline ARM_GP_Vector* GetNotionals()    { return itsNominals;       }
	inline bool GetAtTheMoneyFlag() 		{ return itsAtTheMoneyFlag; }



private :
	void CopyNoCleanUp(const ARM_VanillaCapDigitalArg& rhs);
	void CleanUp();

    ARM_GP_Vector* itsNominals;
	ARM_GP_Vector* itsResetTimes; 
	ARM_GP_Vector* itsStartTimes;
    ARM_GP_Vector* itsEndTimes;
    ARM_GP_Vector* itsStrikes;
    ARM_GP_Vector* itsPayTimes; 
	ARM_GP_Vector* itsPayPeriods;
	bool           itsAtTheMoneyFlag;
};



///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaCapArg
/// \brief cap/caplet argument simple struct
///////////////////////////////////////////////////////////////

struct ARM_VanillaCapArg: public ARM_VanillaCapDigitalArg
{
    ARM_VanillaCapArg( 
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
		double Spread,
		int DecompFreq,
		bool atTheMoneyFlag = false)
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
			payPeriods,
			Spread,
			DecompFreq,
			atTheMoneyFlag)
	{}

	ARM_VanillaCapArg(const ARM_VanillaCapArg& rhs) : ARM_VanillaCapDigitalArg(rhs) {};
    ARM_VanillaCapArg& operator=(const ARM_VanillaCapArg& rhs)
	{
		if( this != &rhs )
			ARM_VanillaCapDigitalArg::operator =(rhs);
		return *this;
	}
    virtual ~ARM_VanillaCapArg() {}

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
	virtual ARM_Object* Clone() const { return new ARM_VanillaCapArg(*this);}
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_CAP; }
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
