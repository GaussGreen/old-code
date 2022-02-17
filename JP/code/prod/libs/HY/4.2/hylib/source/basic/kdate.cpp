#include "kdate.h"
#include "kexception.h"
#include "ksymboltable.h"
#include "kstring.h"
#include <time.h>

using namespace std;

KString KDate::_KDATE_HOLIDAY_FILE = "NONE";
int KDate::_daysInMonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };


void KDate::setHolidayFile (const char* file)
{
	_KDATE_HOLIDAY_FILE = CheckPathDelimiters(theSymbolTable.get(file));
}

bool KDate::isBusinessDay(KString& holidayFile)
{
	// Gto is so fucking stupid.  Any other library takes one line to do this.
	TBoolean answer;

	char* t = const_cast<char*> (holidayFile.c_str());
	if (GtoIsBusinessDay(_date, t, &answer) IS FAILURE)
		throw KException ("Failure in isBusinesDay");

	return (answer == TRUE);
}



KDate KDate::today()
{
	struct tm *theTime;
	time_t long_time;
	
	time( &long_time );  
	theTime = localtime( &long_time );
	
	long month = 1+ theTime->tm_mon;
	long day   = theTime->tm_mday;
	long year  = 1900 + theTime->tm_year;
	
	if (year < 1980) year += 100;
	
	return KDate(month,day,year);
}


KDate::KDate (int month, int day, int year)
{
	// 1980 is our pivot date.  So, it's hard coded.
	// By 2080, I'll be dead, so I don't care!
	if (year < 100 && year >= 80) year = year + 1900;
	if (year < 100 && year < 80) year = year + 2000;
	
	TMonthDayYear	mdy;
	
	mdy.month = month;
	mdy.day = day;
	mdy.year = year;

	if (GtoMDYToDate(&mdy, &_date) IS FAILURE)
		throw KException ("Bad Date Input.  Bad.  Real Bad.");
	_fracDay = 0.5;

}

KDate::KDate(const KString& dateString)
{
	TDate date;
	char* t = const_cast<char*> (dateString.c_str());
	if (GtoStringToDate(t, &date) IS FAILURE)
		throw KException ("Failed to convert KString to date");
	_date = date;
	_fracDay = 0.5;
}

KDate& KDate::set_day (int day)
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed in set_Year");

	int maxDays = GtoDaysInMonth(mdy.year, mdy.month);
	if (day > maxDays) day = maxDays;

	mdy.day = day;

	if (GtoMDYToDate(&mdy, &_date) IS FAILURE)
		throw KException ("Illegal set day: ") << day;

	return *this;
}

KDate& KDate::set_month (int month)
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed in set_Year");

	int maxDays = GtoDaysInMonth(mdy.year, month);
	if (mdy.day > maxDays) mdy.day = maxDays;

	mdy.month = month;

	if (GtoMDYToDate(&mdy, &_date) IS FAILURE)
		throw KException ("Illegal set month: ") << month;

	return *this;
}

KDate& KDate::set_year (int year)
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed in set_Year");

	int maxDays = GtoDaysInMonth(year, mdy.month);
	if (mdy.day > maxDays) mdy.day = maxDays;

	mdy.year = year;

	if (GtoMDYToDate(&mdy, &_date) IS FAILURE)
		throw KException ("Illegal set year: ") << year;

	return *this;
}

int KDate::get_day()
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed to get MDY from TDate");
	return mdy.day;
}

int KDate::get_month()
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed to get MDY from TDate");
	return mdy.month;
}

int KDate::get_year()
{
	TMonthDayYear	mdy;
	if (GtoDateToMDY(_date, &mdy) IS FAILURE)
		throw KException ("Failed to get MDY from TDate");
	return mdy.year;
}



KDate& KDate::operator=(TDate date)
{
	_date = date;
	_fracDay = 0.5;
	return *this;
}


KDate KDate::operator+(const KDateIn& dateInt) const
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	// God, this is so ugly.
										// Patch to make work with GtoDtFwdAny

	if (GtoDtFwdAdj (_date, &gtoDateInt, &temp) IS FAILURE)
		throw KException("Failed in operator+ of KDate");

	KDate ans (temp);
	ans._fracDay = _fracDay;
	return ans;
}


KDate& KDate::operator+= (const KDateIn& dateInt)
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	

	if (GtoDtFwdAdj (_date, &gtoDateInt, &temp) IS FAILURE)
		throw KException("Failed in operator+ of KDate");

	_date = temp;
	return *this;
}



KDate KDate::operator- (const KDateIn& dateInt) const
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;
	gtoDateInt.interval.prd = -gtoDateInt.interval.prd;

	if (GtoDtFwdAdj (_date, &gtoDateInt, &temp) IS FAILURE)
		throw KException("Failed in operator- of KDate");

	KDate ans (temp);
	ans._fracDay = _fracDay;
	return ans;
}

KDate& KDate::operator-= (const KDateIn& dateInt)
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	
	gtoDateInt.interval.prd = -gtoDateInt.interval.prd;

	if (GtoDtFwdAdj (_date, &gtoDateInt, &temp) IS FAILURE)
		throw KException("Failed in operator- of KDate");

	_date = temp;

	return *this;
}

double GetYearDiff (const KDate date1, const KDate& date2, int dayCount)  
{
	double	yearFrac;
	if (GtoDayCountFraction (date1, date2, dayCount, &yearFrac) != SUCCESS)
		throw KException ("Day Count Fraction failed");

	return yearFrac + (date2._fracDay - date1._fracDay)/365.;
}
int GetDaysDiff (const KDate date1, const KDate& date2, int dayCount)  
{
	long dayDiff;
	if (GtoDaysDiff(date1, date2, dayCount, &dayDiff) != SUCCESS)
		throw KException ("Day Count Fraction failed");

	return dayDiff;
}

KDate	KDate::dFwd(double t)const
{
	double	x = _date + _fracDay + t;
	KDate	ans;
	ans._date = long(x);
	ans._fracDay = x - ans._date;
	return ans;
}

ostream& operator<<(ostream& s, const KDate& a)
{
	KString temp = GtoFormatDate (a._date);
	s << temp;
	return s;
}

ostream& operator<<(ostream& s, const KDateIn& a)
{
	s << "(" << a.interval.prd << "," ;
	s << a.interval.prd_typ << "," ;
	s << a.isBusDays << "," ;
	s << a.holidayFile << ")";
	return s;
}

double	&KTimePoint::operator[](const KString &s)
{
	KMap(KString, double)::iterator	iter = _m.find(s); 
	if(iter == _m.end())
	{
		pair< KMap(KString, double)::iterator, bool> tmp;
		tmp = _m.insert(KMap(KString, double)::value_type(s, -1E10));
		if(tmp.second == false)
			KException("KTimePoint::value failed!");
		iter = tmp.first;
	}
	return iter->second;
}



