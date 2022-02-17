/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairswap.h
 *
 *  \brief vanilla interest rates swap
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPCALIB_VANILLAIRSWAP_H
#define _INGPCALIB_VANILLAIRSWAP_H

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "vanillaarg.h"


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaIRSwapArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaIRSwapArg:public ARM_VanillaArg
{
    ARM_VanillaIRSwapArg(
		const string& curveName,
		double evalTime,
		int CallPut,
		double nominal,
		double expiryTime,
		double startTime,
		double endTime,
		ARM_GP_Vector* strikes,	
		ARM_GP_Vector* fixPayTimes,
		ARM_GP_Vector* fixPayPeriods,
		ARM_GP_Vector* floatResetTimes,
		ARM_GP_Vector* floatStartTimes,
		ARM_GP_Vector* floatEndTimes,
		ARM_GP_Vector* floatIntTerms,
        int FixFrequency			)
	:
		ARM_VanillaArg(curveName, evalTime, CallPut, expiryTime),
		itsNominal(nominal),
		itsStartTime(startTime),
		itsEndTime(endTime),		
		itsStrikes(strikes),
		itsFixPayTimes(fixPayTimes),
		itsFixPayPeriods(fixPayPeriods ),
		itsFloatResetTimes(floatResetTimes),
		itsFloatStartTimes(floatStartTimes),
		itsFloatEndTimes(floatEndTimes),
		itsFloatIntTerms(floatIntTerms),
        itsFixFrequency(FixFrequency)
	{}		
    ARM_VanillaIRSwapArg(const ARM_VanillaIRSwapArg& arg);
    ARM_VanillaIRSwapArg& operator=(const ARM_VanillaIRSwapArg& rhs);
    virtual ~ARM_VanillaIRSwapArg();
    virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_IRSWAP; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors (const and non const version)
    inline const ARM_GP_Vector* GetFixPayTimes() const      { return itsFixPayTimes;        }
	inline const ARM_GP_Vector* GetFixPayPeriods() const    { return itsFixPayPeriods;      }
	inline const ARM_GP_Vector* GetFloatResetTimes() const  { return itsFloatResetTimes;    }
	inline const ARM_GP_Vector* GetFloatStartTimes() const  { return itsFloatStartTimes;    }
	inline const ARM_GP_Vector* GetFloatEndTimes() const    { return itsFloatEndTimes;      }
	inline const ARM_GP_Vector* GetFloatIntTerms() const    { return itsFloatIntTerms;      }
    inline ARM_GP_Vector* GetStrikes() const     { return itsStrikes;            }

    inline ARM_GP_Vector* GetFixPayTimes()      { return itsFixPayTimes;        }
	inline ARM_GP_Vector* GetFixPayPeriods()    { return itsFixPayPeriods;      }
	inline ARM_GP_Vector* GetFloatResetTimes()  { return itsFloatResetTimes;    }
	inline ARM_GP_Vector* GetFloatStartTimes()  { return itsFloatStartTimes;    }
	inline ARM_GP_Vector* GetFloatEndTimes()    { return itsFloatEndTimes;      }
	inline ARM_GP_Vector* GetFloatIntTerms()    { return itsFloatIntTerms;      }
    inline ARM_GP_Vector* GetStrikes(){ return itsStrikes;            }
    inline int GetFixFrequency( ) const         { return itsFixFrequency;       }

private:
    double itsNominal;
	double itsStartTime;
    double itsEndTime;
    ARM_GP_Vector*	itsStrikes;
    ARM_GP_Vector*		itsFixPayTimes;
	ARM_GP_Vector*		itsFixPayPeriods;
	
	/// non strictly necessary for model
	/// but useful for market models!
	ARM_GP_Vector*	itsFloatResetTimes;
	ARM_GP_Vector* itsFloatStartTimes;
	ARM_GP_Vector* itsFloatEndTimes;
	ARM_GP_Vector* itsFloatIntTerms;

    ///Get fixed frequency
    int itsFixFrequency;
    void CopyNoCleanUp(const ARM_VanillaIRSwapArg& rhs);
	void CleanUp();
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
