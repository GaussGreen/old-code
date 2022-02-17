#include "drdate.h"
#include "drexception.h"
#include "drutils.h"
#include "drsymboltable.h"
#include <time.h>

DRString DRDate::DRDATE_HOLIDAY_FILE = "NONE";
int DRDate::DaysInMonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };


void DRDate::SetHolidayFile (const char* file)
{
	DRDATE_HOLIDAY_FILE = CheckPathDelimiters(theSymbolTable.get(file));
}

bool DRDate::IsBusinessDay(DRString& holidayFile)
{
	// Gto is so fucking stupid.  Any other library takes one line to do this.
	TBoolean answer;

	char* t = const_cast<char*> (holidayFile.c_str());
	if (GtoIsBusinessDay(m_date, t, &answer) IS FAILURE)
		throw DRException ("Failure in IsBusinesDay");

	return (answer == TRUE);
}

int DRDate::GetDaysInMonth(int month)
{
	return DaysInMonth[month];
}


DRDate DRDate::today()
{
	struct tm *theTime;
	time_t long_time;
	
	time( &long_time );  
	theTime = localtime( &long_time );
	
	long month = 1+ theTime->tm_mon;
	long day   = theTime->tm_mday;
	long year  = 1900 + theTime->tm_year;
	
	if (year < 1980) year += 100;
	
	return DRDate(month,day,year);
}

DRDate DRDate::MakeDateFromMonthYear (int monthYear)
{
	int year = (monthYear - 1) / 12;
	int month = monthYear - 12 * year;
	return DRDate(month, 1, year);
}

DRDate::DRDate (int month, int day, int year)
{
	// 1980 is our pivot date.  So, it's hard coded.
	// By 2080, I'll be dead, so I don't care!
	if (year < 100 && year >= 80) year = year + 1900;
	if (year < 100 && year < 80) year = year + 2000;
	
	m_mdy.month = month;
	m_mdy.day = day;
	m_mdy.year = year;

	if (GtoMDYToDate(&m_mdy, &m_date) IS FAILURE)
		throw DRException ("Bad Date Input.  Bad.  Real Bad.");

	m_mdyGood = true;
}

DRDate& DRDate::SetDay (int day)
{
	if (!m_mdyGood) ComputeMDY ();

	int maxDays = GtoDaysInMonth(m_mdy.year, m_mdy.month);
	if (day > maxDays) day = maxDays;

	m_mdy.day = day;

	if (GtoMDYToDate(&m_mdy, &m_date) IS FAILURE)
		throw DRException ("Illegal set day: ") << day;

	return *this;
}

DRDate& DRDate::SetMonth (int month)
{
	if (!m_mdyGood) ComputeMDY ();

	int maxDays = GtoDaysInMonth(m_mdy.year, month);
	if (m_mdy.day > maxDays) m_mdy.day = maxDays;

	m_mdy.month = month;

	if (GtoMDYToDate(&m_mdy, &m_date) IS FAILURE)
		throw DRException ("Illegal set month: ") << month;

	return *this;
}

DRDate& DRDate::SetYear (int year)
{
	if (!m_mdyGood) ComputeMDY ();

	int maxDays = GtoDaysInMonth(year, m_mdy.month);
	if (m_mdy.day > maxDays) m_mdy.day = maxDays;

	m_mdy.year = year;

	if (GtoMDYToDate(&m_mdy, &m_date) IS FAILURE)
		throw DRException ("Illegal set year: ") << year;

	return *this;
}

void DRDate::ComputeMDY ()
{
	if (GtoDateToMDY(m_date, &m_mdy) IS FAILURE)
		throw DRException ("Failed to get MDY from TDate");

	m_mdyGood = true;
}

int DRDate::GetDay()
{
	if (!m_mdyGood) ComputeMDY();
	return m_mdy.day;
}

int DRDate::GetMonth()
{
	if (!m_mdyGood) ComputeMDY();
	return m_mdy.month;
}

int DRDate::GetYear()
{
	if (!m_mdyGood) ComputeMDY();
	return m_mdy.year;
}

int DRDate::GetMonthYear() 
{
	if (!m_mdyGood) ComputeMDY();
	return m_mdy.year * 12 + m_mdy.month;
}

DRDate& DRDate::operator=(TDate date)
{
	m_date = date;
	m_mdyGood = false;
	
	return *this;
}


DRDate DRDate::operator+(const DRDateIn& dateInt) const
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	// God, this is so ugly.
										// Patch to make work with GtoDtFwdAny

	if (GtoDtFwdAdj (m_date, &gtoDateInt, &temp) IS FAILURE)
		throw DRException("Failed in operator+ of DRDate");

	DRDate ans (temp);
	return ans;
}


DRDate& DRDate::operator+= (const DRDateIn& dateInt)
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	

	if (GtoDtFwdAdj (m_date, &gtoDateInt, &temp) IS FAILURE)
		throw DRException("Failed in operator+ of DRDate");

	m_date = temp;
	m_mdyGood = false;

	return *this;
}


DRDate operator+(const DRDateIn& dateInt, const DRDate& date)
{
	return date + dateInt;
}

DRDate operator-(const DRDateIn& dateInt, const DRDate& date)
{
	return date - dateInt;
}

DRDate DRDate::operator- (const DRDateIn& dateInt) const
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;
	gtoDateInt.interval.prd = -gtoDateInt.interval.prd;

	if (GtoDtFwdAdj (m_date, &gtoDateInt, &temp) IS FAILURE)
		throw DRException("Failed in operator- of DRDate");

	DRDate ans (temp);
	return ans;
}

DRDate& DRDate::operator-= (const DRDateIn& dateInt)
{
	TDate temp;
	TDateAdjIntvl gtoDateInt = dateInt;	
	gtoDateInt.interval.prd = -gtoDateInt.interval.prd;

	if (GtoDtFwdAdj (m_date, &gtoDateInt, &temp) IS FAILURE)
		throw DRException("Failed in operator- of DRDate");

	m_date = temp;
	m_mdyGood = false;

	return *this;
}

double GetYearFrac (const DRDate date1, const DRDate& date2, int dayCount)  
{
	double yearFrac;

	if (GtoDayCountFraction (date1, date2, dayCount, &yearFrac) ISNT SUCCESS)
		throw DRException ("Day Count Fraction failed");

	return yearFrac;
}

long GetNumDays (const DRDate date1, const DRDate& date2, int dayCount)  
{
	long numDays;

	if (GtoDaysDiff (date1, date2, dayCount, &numDays) ISNT SUCCESS)
		throw DRException ("Days Diff failed");

	return numDays;
}

DRDate toDate(long date)
{
  int year = date / 10000;
  int month = (date / 100) % 100;
  int day = date % 100;

  return DRDate(month, day, year);
}

DRDate toDate(DRString& dateString)
{
	if (dateString.at(0) == 'T') {

		if (dateString.size() > 2) {
			if (dateString.at(dateString.size() - 1) == 'B') {
				int numBDays = atoi(left_substring(right_substring(dateString, -1), -1).c_str());
				return DRDate::today() + numBDays * ONE_BUS_DAY;
			}
			else {
				int numDays = atoi (right_substring(dateString, -1).c_str());
				return DRDate::today() + numDays * ONE_DAY;
			}
		}
		else 
			return DRDate::today();
	}

	else if (dateString == "KAPITALDATE" || dateString == "KAPITALVALUEDATE") 
	{
		ifstream inFile ("disczero.dat");
		if (!inFile)
			throw DRException ("Failed to find disczero.dat in order to get KapitalDate");

		char junk[20];
		inFile >> junk;
		inFile >> junk;
		inFile >> junk;
		
		int temp;
		inFile >> temp;

		int year, day, month;
		month = (temp / 100) % 100;
		day = temp % 100;
		year = temp / 10000;

		return DRDate(month, day, year);
	}

	else if (dateString == "KAPITALTODAYDATE") 
	{
		ifstream inFile ("today.dat");
		if (!inFile)
			throw DRException ("Failed to find today.dat in order to get KapitalDate");

		int temp;
		inFile >> temp;

		int year, day, month;
		month = (temp / 100) % 100;
		day = temp % 100;
		year = temp / 10000;

		return DRDate(month, day, year);
	}

	else if (atof(dateString.c_str()) == -999) return DRDate::today();

	else if (atof (dateString.c_str()) == 0) return DRDate(0);

	else {
		TDate date;
		char* t = const_cast<char*> (dateString.c_str());
		if (GtoStringToDate(t, &date) IS FAILURE)
			throw DRException ("Failed to convert string to date");
		
		return DRDate(date);
	}
}

ostream& operator<<(ostream& s, const DRDate& a)
{
	DRString temp = GtoFormatDate (a.m_date);
	s << temp;
	return s;
}

DRString toString(const DRDate& date) 
{return DRString(GtoFormatDate((TDate) date));}




