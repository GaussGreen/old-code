/*!
 *
 * Copyright (c) IXIS CI January 2005 Paris
 *
 * \file fdates.h
 *
 *  \brief file for the fast and flexible date object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPBASE_FDATES_H
#define _INGPBASE_FDATES_H
 
#include "gpbase/port.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

#define MAXMONTHNB		12
#define MAXWEEKDAYNB	07

extern const char* ARM_DEFAULT_COUNTRY;

/// simple date for conversion between julian and std date!
struct ARM_SimpleDate
{
	int itsDay;
	int itsMonth;
	int itsYear;
	static bool IsAValidDate(int d, int m, int y);
	static ARM_SimpleDate JulianToDMY( double Julian );
	static double DMYToJulian( const ARM_SimpleDate& rhs );
	static int DaysInMonth(int m, int y);
	static bool IsLeap( int y);

};

/// std class ARM_FDate
class ARM_FDate
{
public:
	enum Month { Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec };
	enum Day { Sun, Mon, Tue, Wed, Thu, Fri, Sat};

	static const int MINMONTH;
	static const int MAXMONTH;
	static const int MINDAY;
	static const int MAXDAY;
	static const int MINWEEKDAY;
	static const int MAXWEEKDAY;
	static const int MINDATE;
	static const int MAXDATE;
	static const int MINYEAR;
	static const int MAXYEAR;

	static const string ShortMonths[MAXMONTHNB];
	static const string ShortMonthsCapitalLetter[MAXMONTHNB];
	static const string LongMonths[MAXMONTHNB];
	static const string ShortWeekDays[MAXWEEKDAYNB];
	static const string LongWeekDays[MAXWEEKDAYNB];
	static const int MonthDays[MAXMONTHNB];


private:
	/// nested class to avoid conflict
	double itsJulian;
	ARM_SimpleDate itsDate;
	Day itsDayOfWeek;

public:
	/// constructor
	ARM_FDate( int day, int month, int year );
	ARM_FDate( double julian );
	ARM_FDate( const string& stringDate, const string& format = "DD.MM.YYYY" );

	/// copy constructor
    ARM_FDate(const ARM_FDate& d);
	ARM_FDate& operator=(const ARM_FDate& rhs );
    virtual ~ARM_FDate() {}

	/// acccessors
    inline double	GetJulian() const { return itsJulian;}
	inline int GetDay() const { return itsDate.itsDay;}
    inline Month GetMonth() const { return (Month) itsDate.itsMonth; }
    inline int GetYear() const { return itsDate.itsYear; }
    inline Day GetDayOfWeek() const { return itsDayOfWeek; }

	/// ARM support
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    
	/// Operators
    ARM_FDate& AddDays(int n);
    ARM_FDate& AddMonths(int months, int GOTO_END_OF_MONTH = 0);
    ARM_FDate& AddYears(int);
    ARM_FDate& AddPeriod(int, const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& AddPeriod(const string& s, const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& AddPeriodMult(int, int mult = 1, const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& NextNonWeekEndDay(int d);
    ARM_FDate& AdjustToBusDate(const char* ccy, long rule);
    ARM_FDate& NextBusinessDay(const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& NextBusinessDays(int nb, const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& PreviousBusinessDay( const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& PreviousBusinessDay(int nb, const char* ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& GoodBusinessDay(int rule=0, const char *ccy = ARM_DEFAULT_COUNTRY);
    ARM_FDate& GapBusinessDay(int gap, const char* ccy = ARM_DEFAULT_COUNTRY);

    bool IsOrdinaryYear() const;
    bool IsWeekEndDay() const;
    bool IsLeapYear() const;
	bool IsBusinessDay(const char* Country) const;
};


#undef MAXMONTHNB
#undef  MAXWEEKDAYNB

CC_END_NAMESPACE()


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
