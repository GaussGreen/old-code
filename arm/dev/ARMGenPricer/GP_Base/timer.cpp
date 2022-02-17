/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file stringmanip.cpp
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2005
 */

#include "gpbase/timer.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Timer
/// Static variables
////////////////////////////////////////////////////

bool ARM_Timer::itsIsEngaged = true;
const double ARM_Timer::NOT_CALC_DURATION  = -1.00;

////////////////////////////////////////////////////
///	Class  : ARM_Timer
///	Routine: Engage
///	Returns: void
///	Action : function to engage the timer
////////////////////////////////////////////////////

void ARM_Timer::Engage()
{
	ARM_Timer::itsIsEngaged = true;
}

////////////////////////////////////////////////////
///	Class  : ARM_Timer
///	Routine: Disengage
///	Returns: void
///	Action : function to disengage the timer
/////////////////////////////////////////////////////////////////
void ARM_Timer::Disengage()
{
	ARM_Timer::itsIsEngaged = false;
}

////////////////////////////////////////////////////
///	Class  : ARM_Timer
///	Routine: ClockStartTime
///	Returns: void
///	Action : function to start the timing
/////////////////////////////////////////////////////////////////
void ARM_Timer::ClockStartTime()
{ 
	if (ARM_Timer::itsIsEngaged )
		itsStartTime = clock(); 

}

////////////////////////////////////////////////////
///	Class  : ARM_Timer
///	Routine: ClockEndTime
///	Returns: void
///	Action : function to stop the timing
/////////////////////////////////////////////////////////////////

void ARM_Timer::ClockEndTime()
{ 
	if (ARM_Timer::itsIsEngaged )
	{
		itsEndTime = clock(); 
		itsDuration = ((double)(itsEndTime-itsStartTime))/ ((double)CLOCKS_PER_SEC); 
	}
	else
		itsDuration = ARM_Timer::NOT_CALC_DURATION;
}

CC_END_NAMESPACE()