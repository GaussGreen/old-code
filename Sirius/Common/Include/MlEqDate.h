//	MlEqDate.h :			 Date handle class. (Essentially a GDA wrapper.)
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQDATE_H_
#define _MLEQDATE_H_

#include "smart.h"

typedef std::map<int, double>			MlEqDateYearFractionCache;

class MlEqDate : public RCObject
{
public:		
	static unsigned long				ExcelDateToJulianDate(long nExcelDate);	
	static long							GetCurrentDate(void);
	static long							GetNumBusinessDays(long nStartDate, long nEndDate, const std::string& szCalendarOpt);
	static long							JulianDateToExcelDate(unsigned long nJulianDate);
	static std::string					SiriusToGDADayCountConvention(const std::string& sz);
	
	// The lack of a default constructor is deliberate.
	explicit MlEqDate(long nDate);
	explicit MlEqDate(long nDate, DayCountConventionEnum dcc);
	explicit MlEqDate(const MlEqDate& date);
	virtual								~MlEqDate();

	MlEqDateHandle						Add(long nDays, long nMonths, long nYears, BusinessDayConventionEnum bdc, const std::string& szCalendar, bool bPreserveEoM) const;
	static long							AddBusinessDays(long nExcelDate, long nDays, const std::string& szCalendar);
	MlEqDateHandle						AddBusinessDays(long nDays, const std::string& szCalendar);
	void								AddMonth(long nMonths);
	MlEqDateHandle						AddTenor(const std::string& szTenor, BusinessDayConventionEnum bdc = NoBusinessDayConvention, const std::string& szCalendar = "", bool bPreserveEoM = false) const;	
	std::string							GetTenor(long nEndDate, BusinessDayConventionEnum bdc = NoBusinessDayConvention, const std::string& szCalendar = "", bool bPreserveEoM = false) const;
	long								GetDate(void) const;
	DayCountConventionEnum				GetDayCountConvention(void) const;
	estring								GetDayCountConventionStr(void) const;
	static long							GetDayOfMonth(long nExcelDate);
	long								GetDayOfMonth(void) const;	
	WeekdayEnum							GetDayOfWeek(void) const;
	long								GetDaysInYear(void) const;
	MonthEnum							GetMonth(void) const;
	std::string							GetString(void) const;
	std::string							GetStringYMD(void) const;
	long								GetYear(void) const;	
	double								GetYearFraction(long nDate) const;	
	static bool							IsBusinessDay(long nExcelDate, const std::string& szCalendarOpt);
	void								PutDate(long nDate);
	void								PutDayCountConvention(DayCountConventionEnum dcc);	
	void								PutDayOfMonth(long nDay);
	MlEqDateHandle						SetToBusinessDay(BusinessDayConventionEnum bdc, std::string& szCalendar) const;

protected:
	static std::string					GetString(SYSTEMTIME lt);
	bool								IsEqualDayOfMonth(SYSTEMTIME stStart, SYSTEMTIME stEnd) const;
		
	explicit MlEqDate(unsigned long nJulianDate, DayCountConventionEnum dcc);
	mutable MlEqDateYearFractionCache	m_mapYearFractions;	
	static const std::string			s_aszMonths[];	
	long								m_nExcelDate;
	unsigned long						m_nJulianDate;
	DayCountConventionEnum				m_dcc;	
};

#endif