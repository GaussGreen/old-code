/*!
 *
 * Copyright (c) CDC IXIS CM February 2005 Paris
 *
 *	\file vanillaspreadoption.h
 *
 *  \brief vanilla spreadoption
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date February 2005
 */

#ifndef _INGPCALIB_VANILLASPREADOPTION_H
#define _INGPCALIB_VANILLASPREADOPTION_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "vanillaarg.h"
#include "gpinfra/nummethod.h"


CC_BEGIN_NAMESPACE( ARM )
///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSpreadOptionArg
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaSpreadOptionArg:public ARM_VanillaArg
{
	static const string SpreadOptionColNamesTable [];
	static const double PayDateToleranceInDays;

public:
    enum SpreadOptionColAlias
    {
		EventDate=0,	
		PayDate,
		RateLongStartDate,
		RateLongEndDate,
		RateShortStartDate,
		RateShortEndDate,
		RateLong,
		RateShort,	
		IT,	   
		Notional,
		Discount,	 
		CoeffLong,
		CoeffShort,
		Strike,	
		SpreadOption
    };
private:

public:
		
	ARM_VanillaSpreadOptionArg(
		const string& curveName,
		double evalTime,
		int CallPut,
		double expiryTime,
		double startTime,
		double endTime,
		ARM_GP_Vector* resetTimes,
		ARM_GP_Vector* payTimes,
		ARM_GP_Vector* payPeriods,
		ARM_GP_Vector* notional,
		ARM_GP_Vector* coeffLong,
		ARM_GP_Vector* coeffShort,
		ARM_GP_Vector* strikes,
		ARM_GP_Vector* swapLongFloatStartTime,
		ARM_GP_Vector* swapLongFloatEndTime,
		const ARM_VectorVector& swapLongFixPayTimes,
		const ARM_VectorVector& swapLongFixPayPeriods,
		ARM_GP_Vector* swapShortFloatStartTime,
		ARM_GP_Vector* swapShortFloatEndTime,
		const ARM_VectorVector& swapShortFixPayTimes,
		const ARM_VectorVector& swapShortFixPayPeriods,
		bool longIndexIsLibor  = false,
		bool shortIndexIsLibor = false,
		bool isDigital = false,
		ARM_GP_Vector* fixValues = NULL,
		double spread1 = -0.1,
		double spread2 = 0.1,
		int payIndexType					 = K_FIXED,
		ARM_GP_Vector* payIndexLeverages	 = NULL,
		ARM_GP_Vector* swapPayFloatStartTime = NULL,
		ARM_GP_Vector* swapPayFloatEndTime   = NULL,
		const ARM_VectorVector& swapPayFixPayTimes = ARM_VectorVector(),
		const ARM_VectorVector& swapPayFixPayPeriods = ARM_VectorVector(),
		ARM_GP_Vector* payIndexResetTimes = NULL,
		ARM_IntVector* periodIndex		  = NULL)

	:
		ARM_VanillaArg(curveName, evalTime, CallPut,expiryTime ),
		itsStartTime(startTime),
		itsEndTime(endTime),
		itsResetTimes(resetTimes),
		itsPayTimes(payTimes),
		itsPayPeriods(payPeriods),
		itsNotional(notional),
		itsCoeffLong(coeffLong),
		itsCoeffShort(coeffShort),
		itsStrikes(strikes),
		itsSwapLongFloatStartTime(swapLongFloatStartTime),
		itsSwapLongFloatEndTime(swapLongFloatEndTime),
		itsSwapLongFixPayTimes(swapLongFixPayTimes),
		itsSwapLongFixPayPeriods(swapLongFixPayPeriods),
		itsSwapShortFloatStartTime(swapShortFloatStartTime),
		itsSwapShortFloatEndTime(swapShortFloatEndTime),
		itsSwapShortFixPayTimes(swapShortFixPayTimes),
		itsSwapShortFixPayPeriods(swapShortFixPayPeriods),
		itsLongIndexIsLibor  (longIndexIsLibor),
		itsShortIndexIsLibor (shortIndexIsLibor),
		itsIsDigital(isDigital),
		itsFixValues(fixValues),
		itsPayIndexLeverages(payIndexLeverages),
		itsSpread1(spread1),
		itsSpread2(spread2),
		itsPayIndexType(payIndexType),
		itsSwapPayFloatStartTime(swapPayFloatStartTime),
		itsSwapPayFloatEndTime(swapPayFloatEndTime),
		itsSwapPayFixPayTimes(swapPayFixPayTimes),
		itsSwapPayFixPayPeriods(swapPayFixPayPeriods),
		itsPayIndexResetTimes (payIndexResetTimes),
		itsPeriodIndex(periodIndex)
	{
	}	
		
    ARM_VanillaSpreadOptionArg(const ARM_VanillaSpreadOptionArg& arg);
    ARM_VanillaSpreadOptionArg& operator=(const ARM_VanillaSpreadOptionArg& rhs);
    virtual ~ARM_VanillaSpreadOptionArg();

	virtual ARM_GenSecurityPtr VanillaSpreadOptionToGenSec() const;
    virtual double Price(ARM_PricingModel* model) const;
	virtual double ImpliedVol(ARM_PricingModel* model) const;
	
	/// adjusts start & end dates and computes index fix pay times & periods
	virtual void ComputeIndexSchedulesAndAdjustDates (ARM_Currency* ccy, double asOfDate) ;

    virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_SPREADOPTION; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// Accessors (const and non const version)
    inline const ARM_GP_Vector* GetResetTimes()					const	{ return itsResetTimes;					}
	inline const ARM_GP_Vector* GetPayTimes()					const	{ return itsPayTimes;					}
	inline const ARM_GP_Vector* GetSwapLongFloatStartTime()		const	{ return itsSwapLongFloatStartTime;		}
	inline const ARM_GP_Vector* GetSwapLongFloatEndTime()		const	{ return itsSwapLongFloatEndTime;		}
	inline const ARM_GP_Vector* GetSwapShortFloatStartTime()	const	{ return itsSwapShortFloatStartTime;	}
	inline const ARM_GP_Vector* GetSwapShortFloatEndTime()		const	{ return itsSwapShortFloatEndTime;		}	
	inline const ARM_GP_Vector* GetSwapPayFloatStartTime()		const	{ return itsSwapPayFloatStartTime;		}
	inline const ARM_GP_Vector* GetSwapPayFloatEndTime()		const	{ return itsSwapPayFloatEndTime;		}
	inline const ARM_GP_Vector* GetPayPeriods()					const	{ return itsPayPeriods;					}
	inline const ARM_GP_Vector* GetNotional()					const	{ return itsNotional;					}
	inline const ARM_GP_Vector* GetCoeffLong()					const	{ return itsCoeffLong;					}
	inline const ARM_GP_Vector* GetCoeffShort()					const	{ return itsCoeffShort;					}
    inline const ARM_GP_Vector* GetStrikes()					const	{ return itsStrikes;					}
	inline const ARM_VectorVector& GetSwapLongFixPayTimes()		const   { return itsSwapLongFixPayTimes;			}
	inline const ARM_VectorVector& GetSwapLongFixPayPeriods()	const   { return itsSwapLongFixPayPeriods;		}
	inline const ARM_VectorVector& GetSwapShortFixPayTimes()	const   { return itsSwapShortFixPayTimes;		}
	inline const ARM_VectorVector& GetSwapShortFixPayPeriods()	const   { return itsSwapShortFixPayPeriods;		}
	inline const ARM_VectorVector& GetSwapPayFixPayTimes()		const   { return itsSwapPayFixPayTimes;			}
	inline const ARM_VectorVector& GetSwapPayFixPayPeriods()	const   { return itsSwapPayFixPayPeriods;		}
	inline const ARM_IntVector*	 GetPeriodIndex()				const	{ return itsPeriodIndex;				}
	inline const ARM_GP_Vector*	 GetPayIndexResetTimes()		const   { return itsPayIndexResetTimes;			}
	inline const size_t	 GetNumPeriods()						const	{ return itsPayIndexResetTimes->size(); }


    inline ARM_GP_Vector* GetResetTimes()					{ return itsResetTimes;					}
	inline ARM_GP_Vector* GetPayTimes()						{ return itsPayTimes;					}
	inline ARM_GP_Vector* GetSwapLongFloatStartTime()		{ return itsSwapLongFloatStartTime;		}
	inline ARM_GP_Vector* GetSwapLongFloatEndTime()			{ return itsSwapLongFloatEndTime;		}
	inline ARM_GP_Vector* GetSwapShortFloatStartTime()		{ return itsSwapShortFloatStartTime;	}
	inline ARM_GP_Vector* GetSwapShortFloatEndTime()		{ return itsSwapShortFloatEndTime;		}
	inline ARM_GP_Vector* GetSwapPayFloatStartTime()		{ return itsSwapPayFloatStartTime;		}
	inline ARM_GP_Vector* GetSwapPayFloatEndTime()			{ return itsSwapPayFloatEndTime;		}
	inline ARM_GP_Vector* GetPayPeriods()					{ return itsPayPeriods;					}
	inline ARM_GP_Vector* GetNotional()						{ return itsNotional;					}
	inline ARM_GP_Vector* GetCoeffLong()					{ return itsCoeffLong;					}
	inline ARM_GP_Vector* GetCoeffShort()					{ return itsCoeffShort;					}
	inline ARM_GP_Vector* GetStrikes()						{ return itsStrikes;					}
	inline ARM_VectorVector& GetSwapLongFixPayTimes()		{ return itsSwapLongFixPayTimes;		}
	inline ARM_VectorVector& GetSwapLongFixPayPeriods()		{ return itsSwapLongFixPayPeriods;		}
	inline ARM_VectorVector& GetSwapShortFixPayTimes()		{ return itsSwapShortFixPayTimes;		}
	inline ARM_VectorVector& GetSwapShortFixPayPeriods()	{ return itsSwapShortFixPayPeriods;		}
	inline ARM_VectorVector& GetSwapPayFixPayTimes()		{ return itsSwapPayFixPayTimes;			}
	inline ARM_VectorVector& GetSwapPayFixPayPeriods()		{ return itsSwapPayFixPayPeriods;		}
	inline ARM_IntVector*	 GetPeriodIndex()				{ return itsPeriodIndex;				}
	inline ARM_GP_Vector*	 GetPayIndexResetTimes()		{ return itsPayIndexResetTimes;			}
	inline size_t			 GetNumPeriods()				{ return itsPayIndexResetTimes->size(); }



	inline bool IsDigital() const							{ return itsIsDigital; }
	inline double GetSpread1() const						{ return itsSpread1; }
	inline double GetSpread2() const						{ return itsSpread2; }
	inline ARM_GP_Vector* GetFixValues() const				{ return itsFixValues; }
	inline ARM_GP_Vector* GetPayIndexLeverages() const		{ return itsPayIndexLeverages; }
	inline int GetPayIndexType() const						{ return itsPayIndexType;}
	
	inline  void SetSwapLongFloatStartTime(ARM_GP_Vector* vect){ itsSwapLongFloatStartTime =vect;}
	inline  void SetSwapLongFloatEndTime(ARM_GP_Vector* vect){ itsSwapLongFloatEndTime =vect;}
	inline  void SetSwapShortFloatStartTime(ARM_GP_Vector* vect){itsSwapShortFloatStartTime =vect;}
	inline  void SetSwapShortFloatEndTime(ARM_GP_Vector* vect){ itsSwapShortFloatEndTime =vect;}
	inline  void SetSwapPayFloatStartTime(ARM_GP_Vector* vect){itsSwapPayFloatStartTime =vect;}
	inline  void SetSwapPayFloatEndTime(ARM_GP_Vector* vect){ itsSwapPayFloatEndTime =vect;}
	
	inline void SetSwapLongFixPayTimes(int i, ARM_GP_Vector* vect)		{itsSwapLongFixPayTimes[i]=vect;		}
	inline void SetSwapLongFixPayPeriods(int i, ARM_GP_Vector* vect)	{itsSwapLongFixPayPeriods[i]=vect;	}
	inline void SetSwapShortFixPayTimes(int i, ARM_GP_Vector* vect)		{itsSwapShortFixPayTimes[i]=vect;	}
	inline void SetSwapShortFixPayPeriods(int i, ARM_GP_Vector* vect)	{itsSwapShortFixPayPeriods[i]=vect;	}
	inline void SetSwapPayFixPayTimes(int i, ARM_GP_Vector* vect)		{itsSwapPayFixPayTimes[i]=vect;	}
	inline void SetSwapPayFixPayPeriods(int i, ARM_GP_Vector* vect)	{itsSwapPayFixPayPeriods[i]=vect;	}


private:
	double itsStartTime;
	double itsEndTime;

	// to know if index are libor
	bool itsLongIndexIsLibor;
	bool itsShortIndexIsLibor;
			
	ARM_GP_Vector* itsResetTimes;
	ARM_GP_Vector* itsPayTimes;
	ARM_GP_Vector* itsSwapLongFloatStartTime;
	ARM_GP_Vector* itsSwapLongFloatEndTime;
	ARM_GP_Vector* itsSwapShortFloatStartTime;
	ARM_GP_Vector* itsSwapShortFloatEndTime;
	ARM_GP_Vector* itsPayPeriods;
	ARM_GP_Vector* itsNotional;
	ARM_GP_Vector* itsCoeffLong;
	ARM_GP_Vector* itsCoeffShort;
	ARM_GP_Vector* itsStrikes;
	ARM_VectorVector itsSwapLongFixPayTimes;
	ARM_VectorVector itsSwapLongFixPayPeriods;
	ARM_VectorVector itsSwapShortFixPayTimes;
	ARM_VectorVector itsSwapShortFixPayPeriods;
	
	/// index of payment period
	ARM_IntVector* itsPeriodIndex;
	
	bool itsIsDigital;
	// For Corridor Spread
	ARM_GP_Vector* itsFixValues;
	double itsSpread1;
	double itsSpread2;
	
	/// all following vectors are of size = nb of flows 
	int itsPayIndexType;
	ARM_GP_Vector* itsPayIndexLeverages;
	ARM_GP_Vector* itsPayIndexResetTimes;
	ARM_GP_Vector* itsSwapPayFloatStartTime;
	ARM_GP_Vector* itsSwapPayFloatEndTime;
	ARM_VectorVector itsSwapPayFixPayTimes;
	ARM_VectorVector itsSwapPayFixPayPeriods;

	void CopyNoCleanUp(const ARM_VanillaSpreadOptionArg& rhs);
	void CleanUp();
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
