/*!
 *
 * Copyright (c) IXIS CI January 2005 Paris
 *
 *	\file holidays.cpp
 *
 *  \brief file for the holidays object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPBASE_HOLIDAYS_H
#define _INGPBASE_HOLIDAYS_H
 
#include "gpbase/port.h"
#include <map>
#include <string>
CC_USING_NS(std,string)
CC_USING_NS(std,map)
CC_USING_NS(std,less)

#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// a class for all the calendar

class ARM_Holidays
{
	/// nested structure
	struct ARM_VacationData
	{
		ARM_GP_Vector	itsVacationsDates;
		bool			itsIncludeWeekEnd;
		ARM_VacationData( const ARM_GP_Vector& vacationsDates, bool includeWeekEnd )
		:	itsVacationsDates(vacationsDates), itsIncludeWeekEnd(includeWeekEnd) {}
	};
	typedef map< string, ARM_VacationData, less<string> > DataMap;
	DataMap itsData;

public:
	ARM_Holidays( const string& vacationsFile);
	bool IsBusinessDay( const string& isoName, double julianDate) const;
};


extern ARM_Holidays VacationCalendar;

CC_END_NAMESPACE()


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
