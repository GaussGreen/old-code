#ifndef _INGPINFRA_IRINDEX_H
#define _INGPINFRA_IRINDEX_H

#include "gpbase/gpvector.h"
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
	std::vector<double> fixPayTimes;
	// Pay Periods of the fix leg
	std::vector<double> fixPayPeriods;
	
	//////////////////////////////////////////////////////////////////
	/// constructor #1: inputs are times (nb of days to asof)
	//////////////////////////////////////////////////////////////////
	ARM_SwapRate (){};
	ARM_SwapRate (	double					_resetTime, 
					double					_floatStartTime, 
					double					_floatEndTime, 
					const std::vector<double>&	_fixPayTimes, 
					const std::vector<double>&	_fixPayPeriods)
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
					const std::vector<double>&	_fixPayDates, 
					const std::vector<double>&	_fixPayPeriods)
		:	resetTime		(_resetDate - _asOfDate),
			floatStartTime	(_floatStartDate - _asOfDate),
			floatEndTime	(_floatEndDate - _asOfDate),
			fixPayTimes		(_fixPayDates),
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