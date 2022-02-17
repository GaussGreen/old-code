/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file irrate.h
 *
 *  \brief irindex defines structures which represent
 *  rates
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2005
 */


#ifndef _INGPINFRA_IRINDEX_H
#define _INGPINFRA_IRINDEX_H

// gp base
#include "gpbase/gpvector.h"

// gp infra
#include "gpinfra/typedef.h"

CC_BEGIN_NAMESPACE( ARM )



///////////////////////////////////////////////////////
/// \class ARM_SwapRate
/// \brief
/// Structure to represent a swap rate
/// or a CMS index
///////////////////////////////////////////////////////

struct ARM_SwapRate
{
public:
	// First Reset Time of the swap
	double resetTime;
	// Start Time of the floating leg
	double floatStartTime;
	// End Time of the floating leg
	double floatEndTime;
	// Pay Times of the fix leg
	ARM_GP_Vector fixPayTimes;
	// Pay Periods of the fix leg
	ARM_GP_Vector fixPayPeriods;
	
	//////////////////////////////////////////////////////////////////
	/// constructor #1: inputs are times (nb of days to asof)
	//////////////////////////////////////////////////////////////////
	ARM_SwapRate (){};
	ARM_SwapRate (	double					_resetTime, 
					double					_floatStartTime, 
					double					_floatEndTime, 
					const ARM_GP_Vector&	_fixPayTimes, 
					const ARM_GP_Vector&	_fixPayPeriods)
		:	resetTime		(_resetTime),
			floatStartTime	(_floatStartTime),
			floatEndTime	(_floatEndTime),
			fixPayTimes		(_fixPayTimes),
			fixPayPeriods	(_fixPayPeriods) {};

	//////////////////////////////////////////////////////////////////
	/// constructor #2: inputs are absolute dates
	//////////////////////////////////////////////////////////////////
	ARM_SwapRate (	double					_asOfDate,
					double					_resetDate, 
					double					_floatStartDate, 
					double					_floatEndDate, 
					const ARM_GP_Vector&	_fixPayDates, 
					const ARM_GP_Vector&	_fixPayPeriods)
		:	resetTime		(_resetDate - _asOfDate),
			floatStartTime	(_floatStartDate - _asOfDate),
			floatEndTime	(_floatEndDate - _asOfDate),
			fixPayTimes		(_fixPayDates - _asOfDate),
			fixPayPeriods	(_fixPayPeriods) {};
	
			

	/// static builder
	/// caution : input dates are absolute dates !!!
	static ARM_SwapRatePtr CreateSwapRate(
		double asOfDate,
		double startDate,
		double endDate,
		int fixDayCount,
		int fixFreq,
		char* fixCalendar);
};

CC_END_NAMESPACE()

#endif