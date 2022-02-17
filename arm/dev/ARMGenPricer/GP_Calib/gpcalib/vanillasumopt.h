/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairswap.h
 *
 *  \brief vanilla sum option
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2005
 */

#ifndef _INGPCALIB_VANILLASUMOPT_H
#define _INGPCALIB_VANILLASUMOPT_H

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "vanillaarg.h"


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaIRSwapArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaSumOptArg:public ARM_VanillaArg
{
    ARM_VanillaSumOptArg(
		const string& curveName,
		double evalTime,
		int capFloor,
		ARM_GP_VectorPtr coeffs,
		ARM_GP_VectorPtr fwdResetTimes,
		ARM_GP_VectorPtr fwdStartTimes,
		ARM_GP_VectorPtr fwdEndTimes,
		double payTime,
		ARM_GP_VectorPtr fwdPeriods,
		double strike,
		double volatilityRatio)
	:
		ARM_VanillaArg(curveName, evalTime, capFloor),
		itsCoeffs(coeffs),
		itsFwdResetTimes(fwdResetTimes),
		itsFwdStartTimes(fwdStartTimes),
		itsFwdEndTimes(fwdEndTimes),
		itsPayTime(payTime),
		itsFwdPeriods(fwdPeriods),
		itsStrike(strike),
		itsVolatilityRatio(volatilityRatio)
	{}		
    ARM_VanillaSumOptArg(const ARM_VanillaSumOptArg& arg);
    ARM_VanillaSumOptArg& operator=(const ARM_VanillaSumOptArg& rhs);
    virtual ~ARM_VanillaSumOptArg();
    virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_SUMOPT; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors (const and non const version)
    inline const ARM_GP_VectorPtr GetCoeffs() const			{ return itsCoeffs;        }
	inline const ARM_GP_VectorPtr GetFwdResetTimes() const	{ return itsFwdResetTimes;    }
	inline const ARM_GP_VectorPtr GetFwdStartTimes() const	{ return itsFwdStartTimes;    }
	inline const ARM_GP_VectorPtr GetFwdEndTimes() const		{ return itsFwdEndTimes;    }
	inline const ARM_GP_VectorPtr GetFwdPeriods() const		{ return itsFwdPeriods;      }
    

    inline const ARM_GP_VectorPtr GetCoeffs()					{ return itsCoeffs;			}
	inline const ARM_GP_VectorPtr GetFwdResetTimes()			{ return itsFwdResetTimes;	}
	inline const ARM_GP_VectorPtr GetFwdStartTimes()			{ return itsFwdStartTimes;	}
	inline const ARM_GP_VectorPtr GetFwdEndTimes()			{ return itsFwdEndTimes;	}
	inline const ARM_GP_VectorPtr GetFwdPeriods()				{ return itsFwdPeriods;		}

	inline int GetCapFloor() const							{ return itsCapFloor;	}
	inline double GetPayTime() const						{ return itsPayTime;    }
	inline double GetStrike() const							{ return itsStrike;		}
	inline double GetVolatilityRatio() const				{ return itsVolatilityRatio; }
	inline double GetSumFwd() const							{ return itsSumFwd; }
	inline double GetSumVol() const							{ return itsSumVol; }

private:
    ARM_GP_VectorPtr		itsCoeffs;
    ARM_GP_VectorPtr		itsFwdResetTimes;
	ARM_GP_VectorPtr		itsFwdStartTimes;
	ARM_GP_VectorPtr		itsFwdEndTimes;
	ARM_GP_VectorPtr		itsFwdPeriods;

	int itsCapFloor;
    double itsPayTime;
	double itsStrike;
	double itsVolatilityRatio;

	double itsSumFwd;
	double itsSumVol;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
