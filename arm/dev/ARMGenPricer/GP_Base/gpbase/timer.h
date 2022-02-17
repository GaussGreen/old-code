/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file timer.h
 *
 *  \brief header for timing
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPBASE_TIMER_H
#define _INGPBASE_TIMER_H

#include "port.h"
#include <ctime>		/// for timer



CC_BEGIN_NAMESPACE( ARM )

class ARM_Timer
{
public:
	static void Engage();
	static void Disengage();
	static bool itsIsEngaged;
	static const double NOT_CALC_DURATION;
private:
	/// ----------------- timer part
	clock_t itsStartTime;
	clock_t itsEndTime;
	double itsDuration;

public:
	ARM_Timer( clock_t startTime=0, clock_t endTime=0, double duration = NOT_CALC_DURATION)
	:  itsStartTime(startTime), itsEndTime(endTime), itsDuration(duration) {}

	/// timer part
	void ClockStartTime();
	void ClockEndTime();
	inline double GetDuration() const { return itsDuration;}
	inline void SetDuration( double duration ) { itsDuration = duration; }
	inline void SetStartTime( clock_t startTime ) { itsStartTime = startTime; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

