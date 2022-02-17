/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file irrate.cpp
 *
 *  \brief irindex defines structures which represent
 *  rates
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2005
 */

#include "gpinfra/irrate.h"

/// gpbase
#include "gpbase/datestrip.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SwapRate
///	Routine: CreateSwapRate
///	Returns: ARM_SwapRatePtr
///	Action : Create a swap rate structure
///			 /// caution : input dates are absolute dates !!!
////////////////////////////////////////////////////
ARM_SwapRatePtr ARM_SwapRate::CreateSwapRate(
double asOfDate,
double startDate,
double endDate,
int fixDayCount,
int fixFreq,
char* fixCalendar)
{
	ARM_DateStrip dateStrip(
		startDate, 
		endDate, 
		fixFreq, 
		fixDayCount, 
		fixCalendar,
		K_MOD_FOLLOWING, 
		K_ADJUSTED, 
		K_SHORTSTART, 
		GETDEFAULTVALUE, 
		fixFreq, 
		GETDEFAULTVALUE,
		fixCalendar );

	ARM_SwapRate* swapRate = new ARM_SwapRate;

	int size = dateStrip.GetFlowStartDates()->size();

	swapRate->fixPayTimes.resize(size);
	swapRate->fixPayPeriods.resize(size);

	for (size_t i = 0; i < size; ++i)
	{
		swapRate->fixPayTimes[i] = (*dateStrip.GetPaymentDates())[i] - asOfDate;
		swapRate->fixPayPeriods[i] = (*dateStrip.GetInterestTerms())[i];
	}

	swapRate->resetTime = (*dateStrip.GetResetDates())[0]-asOfDate;
	swapRate->floatStartTime = (*dateStrip.GetFlowStartDates())[0]-asOfDate;
	swapRate->floatEndTime = (*dateStrip.GetFlowEndDates())[size-1]-asOfDate;

	return ARM_SwapRatePtr(swapRate);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/