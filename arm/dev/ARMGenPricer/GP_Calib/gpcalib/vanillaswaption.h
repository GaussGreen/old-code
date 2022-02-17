/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillaswaption.h
 *
 *  \brief vanilla swpation
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

#ifndef _INGPCALIB_VANILLASWAPTION_H
#define _INGPCALIB_VANILLASWAPTION_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"

#include "vanillaarg.h"



CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSwaptionArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaSwaptionArg:public ARM_VanillaArg
{
	static const string SwaptionColNamesTable [];
	enum SwaptionColAlias
    {
        ResetDate=0,
        StartDate,
        PayDate,
        Swap,
		Option
    };


    ARM_VanillaSwaptionArg(
		const string& curveName,
		double evalTime,
		int CallPut,
		ARM_GP_Vector* FundingSpread,
		ARM_GP_Vector* FixNominal,
		ARM_GP_Vector* FloatNominal,
		double resetTime,
		double startTime,
		double endTime,
		ARM_GP_Vector* strikes,	
		ARM_GP_Vector* fixPayTimes,
		ARM_GP_Vector* fixPayPeriods,
		ARM_GP_Vector* floatResetTimes,
		ARM_GP_Vector* floatStartTimes,
		ARM_GP_Vector* floatEndTimes,
		ARM_GP_Vector* floatIntTerms,
        int FixFrequency,
		int FloatFrequency,
		int FixDayCount,
		int FloatDayCount,
		double asofDate = 0,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true,
		string* genSecString = NULL,
		bool atTheMoneyFlag = false)
	:
		ARM_VanillaArg(curveName, evalTime, CallPut, resetTime),
		itsFixNominal(FixNominal),
		itsFloatNominal(FloatNominal),
		itsResetTime(resetTime),
		itsStartTime(startTime),
		itsEndTime(endTime),		
		itsStrikes(strikes),
		itsFixPayTimes(fixPayTimes),
		itsFixPayPeriods(fixPayPeriods ),
		itsFloatResetTimes(floatResetTimes),
		itsFloatStartTimes(floatStartTimes),
		itsFloatEndTimes(floatEndTimes),
		itsFloatIntTerms(floatIntTerms),
        itsFixFrequency(FixFrequency),
		itsFloatFrequency(FloatFrequency),
		itsFixDayCount(FixDayCount),
		itsFloatDayCount(FloatDayCount),
		itsIsConstantNotional(isConstantNotional),
		itsIsConstantSpread(isConstantSpread),
		itsIsConstantStrike(isConstantStrike),
		itsAsofDate(asofDate),
		itsGenSecString (genSecString),
		itsFundingSpread(FundingSpread),
		itsAtTheMoneyFlag(atTheMoneyFlag)
	{}	
		ARM_VanillaSwaptionArg(
		const string& curveName,
		double evalTime,
		int CallPut,
		double nominal,
		double resetTime,
		double startTime,
		double endTime,
		ARM_GP_Vector* strikes,	
		ARM_GP_Vector* fixPayTimes,
		ARM_GP_Vector* fixPayPeriods,
		ARM_GP_Vector* floatResetTimes,
		ARM_GP_Vector* floatStartTimes,
		ARM_GP_Vector* floatEndTimes,
		ARM_GP_Vector* floatIntTerms,
        int FixFrequency,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true,
		bool atTheMoneyFlag = false)
	:
		ARM_VanillaArg(curveName, evalTime, CallPut, resetTime),
		itsFixNominal(new ARM_GP_Vector(1,nominal)),
		itsFloatNominal(new ARM_GP_Vector(1,nominal)),
		itsResetTime(resetTime),
		itsStartTime(startTime),
		itsEndTime(endTime),		
		itsStrikes(strikes),
		itsFixPayTimes(fixPayTimes),
		itsFixPayPeriods(fixPayPeriods ),
		itsFloatResetTimes(floatResetTimes),
		itsFloatStartTimes(floatStartTimes),
		itsFloatEndTimes(floatEndTimes),
		itsFloatIntTerms(floatIntTerms),
        itsFixFrequency(FixFrequency),
		itsFloatFrequency(K_SEMIANNUAL),
		itsFixDayCount(K30_360),
		itsFloatDayCount(KACTUAL_360),
		itsIsConstantNotional(isConstantNotional),
		itsIsConstantSpread(isConstantSpread),
		itsIsConstantStrike(isConstantStrike),
		itsAsofDate(0),
		itsGenSecString(NULL),
		itsFundingSpread(new ARM_GP_Vector(1,0.0)),
		itsAtTheMoneyFlag(atTheMoneyFlag)
	{}	
    ARM_VanillaSwaptionArg(const ARM_VanillaSwaptionArg& arg);
    ARM_VanillaSwaptionArg& operator=(const ARM_VanillaSwaptionArg& rhs);
    virtual ~ARM_VanillaSwaptionArg();
	virtual ARM_GenSecurityPtr VanillaSwaptionToGenSec() const;
    virtual double Price(ARM_PricingModel* model) const;
	double Vega(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_SWAPTION; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

    /// Accessors
    inline double GetStartTime() const {return itsStartTime;}
    inline double GetEndTime() const {return itsEndTime;}
    inline double GetResetTime() const {return itsResetTime;}
	inline double GetAsofDate() const {return itsAsofDate;}

	/// Accessors (const and non const version)
    inline const ARM_GP_Vector* GetFixPayTimes() const      { return itsFixPayTimes;        }
	inline const ARM_GP_Vector* GetFixPayPeriods() const    { return itsFixPayPeriods;      }
	inline const ARM_GP_Vector* GetFloatResetTimes() const  { return itsFloatResetTimes;    }
	inline const ARM_GP_Vector* GetFloatStartTimes() const  { return itsFloatStartTimes;    }
	inline const ARM_GP_Vector* GetFloatEndTimes() const    { return itsFloatEndTimes;      }
	inline const ARM_GP_Vector* GetFloatIntTerms() const    { return itsFloatIntTerms;      }
    inline const ARM_GP_Vector* GetStrikes() const          { return itsStrikes;            }
    inline const ARM_GP_Vector*  GetFixNotional() const     { return itsFixNominal;         }
	inline const ARM_GP_Vector*  GetFloatNotional() const   { return itsFloatNominal;       }
    inline  double  GetNotional() const						{ return (*itsFixNominal)[0];	}
	inline  bool  GetIsConstantNotional() const				{ return itsIsConstantNotional;	}
	inline  bool  GetIsConstantSpread() const				{ return itsIsConstantSpread;	}
	inline  bool  GetIsConstantStrike() const				{ return itsIsConstantStrike;	}


    inline ARM_GP_Vector* GetFixPayTimes()      { return itsFixPayTimes;        }
	inline ARM_GP_Vector* GetFixPayPeriods()    { return itsFixPayPeriods;      }
	inline ARM_GP_Vector* GetFloatResetTimes()  { return itsFloatResetTimes;    }
	inline ARM_GP_Vector* GetFloatStartTimes()  { return itsFloatStartTimes;    }
	inline ARM_GP_Vector* GetFloatEndTimes()    { return itsFloatEndTimes;      }
	inline ARM_GP_Vector* GetFloatIntTerms()    { return itsFloatIntTerms;      }
    inline ARM_GP_Vector* GetStrikes()          { return itsStrikes;            }
	inline ARM_GP_Vector* GetFundingSpread() const    { return itsFundingSpread;      }
    inline int GetFixFrequency( ) const         { return itsFixFrequency;       }
	inline int GetFloatFrequency( ) const         { return itsFloatFrequency;       }	
	inline int GetFixDayCount( ) const         { return itsFixDayCount;       }
	inline int GetFloatDayCount( ) const         { return itsFloatDayCount;       }
	inline bool GetAtTheMoneyFlag() const		 { return itsAtTheMoneyFlag;}
	inline void SetAtTheMoneyFlag(bool flag)	{ itsAtTheMoneyFlag = flag;}


protected:
	ARM_GP_Vector* itsFundingSpread;
	ARM_GP_Vector* itsFixNominal;
	ARM_GP_Vector* itsFloatNominal;
	bool itsIsConstantNotional;
	bool itsIsConstantSpread;
	bool itsIsConstantStrike;
	double itsAsofDate;
	double itsResetTime;
	double itsStartTime;
    double itsEndTime;
    ARM_GP_Vector*	itsStrikes;
    ARM_GP_Vector*	itsFixPayTimes;
	ARM_GP_Vector*	itsFixPayPeriods;
	string* itsGenSecString; /// Pointer to the GenSecString of the Swaption
							/// The Pointer will be deleted only when we destruct the swaption
	
	/// non strictly necessary for model
	/// but useful for market models!
	ARM_GP_Vector*	itsFloatResetTimes;
	ARM_GP_Vector* itsFloatStartTimes;
	ARM_GP_Vector* itsFloatEndTimes;
	ARM_GP_Vector* itsFloatIntTerms;

    ///Get fixed frequency
    int itsFixFrequency;
	int itsFloatFrequency;
	int itsFixDayCount;
	int itsFloatDayCount;

	mutable bool itsAtTheMoneyFlag;

    void CopyNoCleanUp(const ARM_VanillaSwaptionArg& rhs);
	void CleanUp();
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
