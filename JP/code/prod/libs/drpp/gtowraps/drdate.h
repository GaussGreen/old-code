#ifndef __drdate__h
#define __drdate__h

extern "C" {
#include "ldate.h"
#include "dateconv.h"
#include "convert.h"
}

#include "drdatein.h"
#include "drstring.h"

//	This class wraps a TDate
//	They should cast freely into each other.

//  I like DRDate cause it gives me C++ printing and the ability to overload
//	operator+ and operator- so I can add and subtract DRDateIn's.

//	For example, I can write neat stuff like:
//
//	DRDate XmasDay(12,25,1997);
//	cout << "Ho! Ho! Ho! The day two months after Xmas is " << XMasDay + 2 * ONE_MONTH << endl;


class DRDate {
public:
	static DRDate today ();
	static DRDate MakeDateFromMonthYear (int);
	static void SetHolidayFile (const char*);
	static int GetDaysInMonth(int month);

	DRDate(TDate date = 0);
	DRDate(int month, int day, int year);

	int GetDay();
	int GetMonth();
	int GetYear();

	DRDate& SetDay (int);
	DRDate& SetMonth(int);
	DRDate& SetYear(int);

	int GetMonthYear();

	bool IsBusinessDay(DRString& holidayFile = DRDate::DRDATE_HOLIDAY_FILE);

	DRDate& operator=(TDate);

	operator TDate();
	operator TDate() const;

	DRDate operator+(const DRDateIn&) const;
	friend DRDate operator+(const DRDateIn&, const DRDate&);
	DRDate& operator+=(const DRDateIn&);

	DRDate operator-(const DRDateIn&) const;
	friend DRDate operator-(const DRDateIn&, const DRDate&);
	DRDate& operator-= (const DRDateIn&);

	friend double GetYearFrac(const DRDate, const DRDate&, int dayCount = GTO_ACT_365F);
	friend long GetNumDays(const DRDate, const DRDate&, int dayCount = GTO_ACT_365F);

	friend ostream& operator<<(ostream&, const DRDate&);

private:
	static DRString DRDATE_HOLIDAY_FILE;
	static int DaysInMonth[];

	TDate m_date;

	TMonthDayYear m_mdy;
	bool m_mdyGood;

	void ComputeMDY();
};

DRDate toDate(long); // input in YYYYMMDD format

DRString toString(const DRDate& date);

DRDate toDate(DRString&);

inline DRDate::operator TDate() {return m_date;}

inline DRDate::operator TDate() const {return m_date;}

inline DRDate::DRDate(TDate date) {m_date = date; m_mdyGood = false;}

#endif


