/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file fdates.cpp
 *
 *  \brief file for the fast and flexible dates
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#if defined(__USE_BASE_LIBRARY)

#include "gpbase/fdates.h"
#include "gpbase/env.h"
#include <cmath>

/// ARM Kernel
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_FDate
///	Routine: General statics
////////////////////////////////////////////////////

#define MAXMONTHNB		12
#define MAXWEEKDAYNB	07

	const int ARM_FDate::MINMONTH	= 1;
	const int ARM_FDate::MAXMONTH	= 12;
	const int ARM_FDate::MINDAY		= 1;
	const int ARM_FDate::MAXDAY		= 31;
	const int ARM_FDate::MINWEEKDAY	= 0;
	const int ARM_FDate::MAXWEEKDAY	= 6;
	const int ARM_FDate::MINDATE     = 4749L;      /// 1 January 4700 bc
	const int ARM_FDate::MAXDATE     = 10852487L;  /// 31 December 25000
	const int ARM_FDate::MINYEAR		= 1960;
	const int ARM_FDate::MAXYEAR		= 5000;



	const int ARM_FDate::MonthDays[MAXMONTHNB] = { 
		31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	const string ARM_FDate::ShortMonths[MAXMONTHNB] = 
		{"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

	const string ARM_FDate::ShortMonthsCapitalLetter[MAXMONTHNB] = 
		{"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

	const string LongMonths[MAXMONTHNB] = {
		"January", "February", "March", "April", "May", "June", 
		"July", "August", "September", "October", "November", "December"};

	const string ShortWeekDays[MAXWEEKDAYNB] = 
		{ "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};

	const string ARM_FDate::LongWeekDays[MAXWEEKDAYNB]= 
		{"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};


#undef MAXMONTHNB
#undef MAXWEEKDAYNB


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
///	Routine: JulianToDMY
///	Returns: the corresponding day,month and year
////////////////////////////////////////////////////
ARM_SimpleDate ARM_SimpleDate::JulianToDMY(double Julian)
{
    long Ja, JAlpha, Jb, Jc, Jd, Je;
	int day,month,year;
	const int GREGOR  = 2299161L;
 
    if(  Julian>=ARM_FDate::MINDATE && Julian<=ARM_FDate::MAXDATE )
    {
       if ( Julian >= GREGOR )
       {
          JAlpha = (long) (((double)(Julian - 1867216L) - 0.25) / 36524.25);
          Ja = (long) (Julian + 1 + JAlpha - (long)(0.25*JAlpha));
       }
       else
          Ja = (long) Julian;
 
       Jb	= Ja + 1524;
       Jc	= (long) (6680.0 + ((double)(Jb - 2439870L) - 122.1) / 365.25);
       Jd	= (long) (365*Jc + (0.25*Jc));
       Je	= (long)((Jb - Jd)/30.6001);
       day	= (int)(Jb - Jd - (int)(30.6001*Je));
       month= (int)Je - 1;
 
       if (month > 12)
          month -= 12;
 
       year	= (int)(Jc - 4715);
       if (month > 2)
          --year;
       if ( year <= 0)
          --year;
		
	   ARM_SimpleDate date = { day,month,year };
	   return date;

    }
    else
       throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
             "Invalid Date : conversion from Julianian Date to dd/mm/yyyy failed");
}



////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
///	Routine: IsLeap
///	Returns: bool
///	Action : returns true for a leap year
////////////////////////////////////////////////////
bool ARM_SimpleDate::IsLeap( int y)
{
	return y%4==0 && y%4000!=0 && y%100!=0 ||y%400==0;
}


////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
///	Routine: DaysInMonth
///	Returns: int
///	Action : returns nb of days in a month
////////////////////////////////////////////////////
int ARM_SimpleDate::DaysInMonth(int m, int y)
{
	if( m==2 && IsLeap(y) )
		return 29;
	else
		return ARM_FDate::MonthDays[m-1];
}

////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
///	Routine: IsAValidDate
///	Returns: bool
///	Action : returns true if it is a valid date, otherwise false
////////////////////////////////////////////////////

bool ARM_SimpleDate::IsAValidDate(int d, int m, int y)
{
    if (   y< ARM_FDate::MINYEAR   || y > ARM_FDate::MAXYEAR 
        || m < ARM_FDate::MINMONTH || m > ARM_FDate::MAXMONTH 
		|| d > ARM_SimpleDate::DaysInMonth(m,y) )
		return false;
	else
		return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_SimpleDate
///	Routine: JulianToDMY
///	Returns: the corresponding day,month and year
////////////////////////////////////////////////////
double ARM_SimpleDate::DMYToJulian(const ARM_SimpleDate& date )
{
	int d=date.itsDay,
		m=date.itsMonth,
		y=date.itsYear;

    if( IsAValidDate(d,m,y) )
	{
		long Ja, Jm, Jy;
		double Jul;
		if ( y < 0 )
		   y++;

		if ( m > 2 )
		{
		   Jy = y;
		   Jm = m + 1;
		}
		else
		{
		   Jy = y - 1;
		   Jm = m + 13;
		}
 
		static const double IGREG = (15+31L*(10+12L*1582));
		Jul = (floor(365.25*(double)Jy)+floor(30.6001*(double)Jm)+d+1720995L);

		if ( d+31L*(m+12L*y) >= IGREG )
		{
		   Ja = (long) (0.01*Jy);
		   Jul += (2-Ja+(long)(0.25*Ja));
		}
		return(Jul);
	}
	else
		ARM_THROW(ERR_INVALID_DATA, "Invalid Date" );
}


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

