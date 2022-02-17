#ifndef __KDATE__H
#define __KDATE__H

#include "kplatform.h"

extern "C" {
#include "ldate.h"
#include "dateconv.h"
#include "convert.h"
#include "bastypes.h"
#include "busday.h"
}

#include "kstring.h"
#include "kplatdep.h"
#include "kexception.h"
#include "kcurve.h"
#include <map>
#include IOSTREAM_H



class KDateIn : public TDateAdjIntvl{
public:
	KDateIn(int numPeriods, char prdType, TBoolean isBusDaysP, char* holidayFileP = (char*) "NONE");
	friend KDateIn operator*(int, const KDateIn&);
	friend std::ostream& operator<<(std::ostream&, const KDateIn&);
};

// Some useful constants

const KDateIn ONE_MONTH = KDateIn(1, 'M', FALSE);
const KDateIn ONE_YEAR = KDateIn(1, 'A', FALSE);
const KDateIn ONE_DAY = KDateIn(1, 'D', FALSE);
const KDateIn ONE_SEMI = KDateIn (1, 'S', FALSE);
const KDateIn ONE_QUARTER = KDateIn(1, 'Q', FALSE);
const KDateIn ONE_BUS_DAY = KDateIn(1, 'D', TRUE);

//day count methods
const	int	B30_360   = GTO_B30_360;
const	int	ACT_ACT   = GTO_ACT_365;
const	int	ACT_360   = GTO_ACT_360;
const	int	ACT_365   = GTO_ACT_365F;
const	int	ACT_365J  = GTO_ACT_365FJ;
const	int	B30E_360  = GTO_B30E_360;
const	int	B30E_360I = GTO_B30E_360I;
const	int	B30EP_360 = GTO_B30EP_360;

inline int	StringToDayCountConv(const KString &strConv)
{
	int	conv;
	upper_string str = strConv;
	if(str == "B30_360")
		conv = B30_360;
	else if(str == "ACT_ACT")
		conv = ACT_ACT;
	else if(str == "ACT_360")
		conv = ACT_360;
	else if(str == "ACT_365")
		conv = ACT_365;
	else if(str == "ACT_365J")
		conv = ACT_365J;
	else if(str == "B30E_360")
		conv = B30E_360;
	else if(str == "B30E_360I")
		conv = B30E_360I;
	else if(str == "B30EP_360")
		conv = B30EP_360;
	else
		KException("Wrong day count convention!");
	return conv;
}
inline KString	DayCountConvToString(int conv)
{
	KString	str;
	if(conv == B30_360)
		str = "B30_360";
	else if(conv == ACT_ACT)
		str = "ACT_ACT";
	else if(conv == ACT_360)
		str = "ACT_360";
	else if(conv == ACT_365)
		str = "ACT_365";
	else if(conv == ACT_365J)
		str = "ACT_365J";
	else if(conv == B30E_360)
		str = "B30E_360";
	else if(conv == B30E_360I)
		str = "B30E_360I";
	else if(conv == B30EP_360)
		str = "B30EP_360";
	else
		KException("Wrong day count convention!");
	return str;
}
inline KDateIn::KDateIn (int numPeriods, char prdType, TBoolean isBusDaysP, char* holidayFileP)
{
	SET_TDATE_INTERVAL (interval, numPeriods, prdType);
	isBusDays = isBusDaysP;
	holidayFile = holidayFileP;
    badDayConv = GTO_BAD_DAY_NONE;
}

inline KDateIn operator*(int num, const KDateIn& a)
{
	int newPrds = num * a.interval.prd;

	return KDateIn(newPrds, a.interval.prd_typ, a.isBusDays, a.holidayFile);
}

//	This class wraps a TDate
//	They should cast freely into each other.
class KDate {
public:
	static KDate today ();
	static void setHolidayFile (const char*);
	
	KDate(TDate date = 0){_date = date; _fracDay = .5;}
	KDate(int month, int day, int year);
	KDate(const KString&); //input in "MM/DD/YYYY" or "YYYYMMDD"format
	KDate& operator=(TDate);
	operator TDate();
	operator TDate() const;

	int get_day();
	int get_month();
	int get_year();

	KDate& set_day (int);
	KDate& set_month(int);
	KDate& set_year(int);

	bool isBusinessDay(KString& holidayFile = KDate::_KDATE_HOLIDAY_FILE);

	KDate operator+(const KDateIn&) const;
	KDate& operator+=(const KDateIn&);

	KDate operator-(const KDateIn&) const;
	KDate& operator-= (const KDateIn&);
	
	KDate	dFwd(double)const;
	KString	str()const{return KString(GtoFormatDate(_date));}
	friend double	GetYearDiff(const KDate, const KDate&, int dayCount = ACT_365);
	friend int		GetDaysDiff(const KDate, const KDate&, int dayCount = ACT_365);
	friend std::ostream& operator<<(std::ostream&, const KDate&);

private:
	static KString _KDATE_HOLIDAY_FILE;
	static int _daysInMonth[];
	TDate _date;
	double	_fracDay; //[0, 1)
};

inline KDate::operator TDate() {return _date;}
inline KDate::operator TDate() const {return _date;}


//class type for special values and types for a date
typedef	KCurve<KString, double>	KTimePoint;
//class type for datelist
typedef	KCurve<KDate, KTimePoint>	KTimeLine;

#endif


