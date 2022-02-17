//	MlEqDate.cpp :             Implementation of the MlEqDate class
//
//	Author :				   David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqdate.h"

/*static*/ const std::string			MlEqDate::s_aszMonths[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };


/*explicit*/ MlEqDate::MlEqDate(long nDate)
{
	m_dcc = ActualActual;
	PutDate(nDate);
}

/*explicit*/ MlEqDate::MlEqDate(long nDate, DayCountConventionEnum dcc)
{
	m_dcc = dcc;
	PutDate(nDate);
}

/*explicit*/ MlEqDate::MlEqDate(unsigned long nJulianDate, DayCountConventionEnum dcc)
{	
	m_dcc = dcc;
	m_nJulianDate = nJulianDate;
	m_nExcelDate = JulianDateToExcelDate(nJulianDate);
}

/*explicit*/ MlEqDate::MlEqDate(const MlEqDate& date)
{	
	m_mapYearFractions = date.m_mapYearFractions;
	m_nExcelDate = date.m_nExcelDate;
	m_nJulianDate = date.m_nJulianDate;
	m_dcc = date.m_dcc;
}

MlEqDate::~MlEqDate()
{
}

MlEqDateHandle MlEqDate::Add(long nDays, long nMonths, long nYears, BusinessDayConventionEnum bdc, const std::string& szCalendar, bool bPreserveEoM) const
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Date", GDA::HDElement(m_nJulianDate));
	if (nDays) hdeParams.insert("Days", GDA::HDElement(nDays));
	if (nMonths) hdeParams.insert("Months", GDA::HDElement(nMonths));
	if (nYears) hdeParams.insert("Years", GDA::HDElement(nYears));
	if (bdc){
		std::string szBdc;
		CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdc, &szBdc);
		hdeParams.insert("BDC", GDA::HDElement(szBdc.data()));
	}
	if (szCalendar.size()) hdeParams.insert("Calendar", GDA::HDElement(szCalendar.data()));
	hdeParams.insert("PreserveEoM", GDA::HDElement(bPreserveEoM));		
	return new MlEqDate(GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.AddYMD", hdeParams, GDA::HDElement("ComputedDate"), GDA::GetDefaultContext()).lookup("ComputedDate").asJulian(), m_dcc);
}


/*static*/ long MlEqDate::AddBusinessDays(long nExcelDate, long nDays, const std::string& szCalendar)
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Date", GDA::HDElement(ExcelDateToJulianDate(nExcelDate)));
	hdeParams.insert("DaysToAdd", GDA::HDElement(nDays));
	if (szCalendar.size()) hdeParams.insert("Calendar", GDA::HDElement(szCalendar.data()));	
	long nJulianDate = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.AddBusDays", hdeParams, GDA::HDElement("ComputedDate"), GDA::GetDefaultContext()).lookup("ComputedDate").asJulian();
	return JulianDateToExcelDate(nJulianDate);	
}

MlEqDateHandle MlEqDate::AddBusinessDays(long nDays, const std::string& szCalendar)
{	
	long n = AddBusinessDays(m_nExcelDate, nDays, szCalendar);	
	return new MlEqDate(n, m_dcc);
}

void MlEqDate::AddMonth(long nMonths)
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Date", GDA::HDElement(m_nJulianDate));
	if (nMonths) hdeParams.insert("Months", GDA::HDElement(nMonths));
	unsigned long nJulianDate = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.AddYMD", hdeParams, GDA::HDElement("ComputedDate"), GDA::GetDefaultContext()).lookup("ComputedDate").asJulian();
	long nExcelDate = JulianDateToExcelDate(nJulianDate);
	PutDate(nExcelDate);
}

MlEqDateHandle MlEqDate::AddTenor(const std::string& szTenor, BusinessDayConventionEnum bdc /*=NoBusinessDayConvention*/, const std::string& szCalendar /*=""*/, bool bPreserveEoM /*=false*/) const
{	
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Tenor", GDA::HDElement(szTenor.data()));
	hdeParams.insert("Date", GDA::HDElement(m_nJulianDate));		
	if (bdc != NoBusinessDayConvention){
		std::string szBdc;
		CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdc, &szBdc);
		hdeParams.insert("BDC", GDA::HDElement(szBdc.data()));
	}
	if (szCalendar.size()) hdeParams.insert("Calendar", GDA::HDElement(szCalendar.data()));
	if (bPreserveEoM) hdeParams.insert("PreserveEoM", GDA::HDElement(bPreserveEoM));		
	
	GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.AddTenor", hdeParams, GDA::HDElement("ComputedDate"), GDA::GetDefaultContext()).lookup("ComputedDate");
	if (hdeResult.isException()){
		hdeResult.raise();
	}		
	return new MlEqDate(hdeResult.asJulian(), m_dcc);
}

//	Returns the number of business days BETWEEN this date and an input date.

/*static*/ long MlEqDate::GetNumBusinessDays(long nStartDate, long nEndDate, const std::string& szCalendarOpt)
{
	if (nStartDate == nEndDate) return 0;

	long nSign = 1L;

	if (nStartDate > nEndDate){
		std::swap(nStartDate, nEndDate);
		nSign = -1L;
	}
			
	long nUpper = nEndDate - nStartDate;	// Largest possible number of business days.
	long nLower = 0L;						// Smallest possible number of business days.
	long nDays = 0L;	

	for (long nPreviousDays = ~nDays; nDays != nPreviousDays;){
		nPreviousDays = nDays;
		nDays = (nUpper + nLower) / 2L;
		long nTrialDate = AddBusinessDays(nStartDate, nDays, szCalendarOpt);
		if (nTrialDate <= nEndDate){
			// More business days could be required
			nLower = nDays;
		} else {
			// Less business days are required
			nUpper = nDays;
		}
	}
	return nSign * nDays;
}

//	returns the current system date
/*static*/ long MlEqDate::GetCurrentDate(void)
{	
	try {
		GDA::Context_ref	ctx = GDA::GetLibraryInstance()->getDefaultNamespace()->getContext();
		GDA::HDElement		hdeCurrentDate = ctx->resolve("SystemDate");
		return JulianDateToExcelDate(hdeCurrentDate.asJulian());
	} catch (...){
		// Gda system date not available
		SYSTEMTIME		lt;	
		DATE			date;
		::GetLocalTime(&lt);	
		SystemTimeToVariantTime(&lt, &date);		
		return date;			
	}
}

long MlEqDate::GetDate(void) const
{
	return m_nExcelDate;
}

DayCountConventionEnum MlEqDate::GetDayCountConvention(void) const
{
	return m_dcc;
}

estring MlEqDate::GetDayCountConventionStr(void) const
{
	return CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, GetDayCountConvention());
}


/*static*/ long MlEqDate::GetDayOfMonth(long nExcelDate)
{
	SYSTEMTIME							lt;
	::VariantTimeToSystemTime(nExcelDate, &lt);
	return lt.wDay;
}

long MlEqDate::GetDayOfMonth(void) const
{
	return GetDayOfMonth(m_nExcelDate);	
}

WeekdayEnum MlEqDate::GetDayOfWeek(void) const
{
	std::string							szWeekday;
	WeekdayEnum							weekday;

	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Date", GDA::HDElement(m_nJulianDate));
	szWeekday = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.DayOfWeek", hdeParams, GDA::HDElement("DayOfWeek"), GDA::GetDefaultContext()).lookup("DayOfWeek").asCStr();
	if (CEnumMap::GetEnum("WeekdayEnum", LIBID_Sirius, szWeekday, &weekday)) throw "MlEqDate::GetDayOfWeek failed";
	return weekday;
}

long MlEqDate::GetDaysInYear(void) const
{		
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Date", GDA::HDElement(m_nJulianDate));			
	return GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.DaysInYear", hdeParams, GDA::HDElement("Days"), GDA::GetDefaultContext()).lookup("Days").asInteger();
}

MonthEnum MlEqDate::GetMonth(void) const
{
	SYSTEMTIME							lt;
	::VariantTimeToSystemTime(m_nExcelDate, &lt);
	return (MonthEnum)lt.wMonth;
}

/*static*/ std::string MlEqDate::GetString(SYSTEMTIME lt)
{
	// This function must support invalid lt.
	std::stringstream					ssOut;
	if (lt.wMonth < 1 || lt.wMonth > 12){
		ssOut << lt.wDay << '-' << lt.wMonth << '-' << lt.wYear;
	} else {
		ssOut << lt.wDay << '-' << s_aszMonths[lt.wMonth - 1] << '-' << lt.wYear;
	}
	return ssOut.str();
}

std::string MlEqDate::GetString(void) const
{	
	SYSTEMTIME							lt;	
	::VariantTimeToSystemTime(m_nExcelDate, &lt);
	return GetString(lt);
}

std::string MlEqDate::GetStringYMD(void) const
{
	SYSTEMTIME							lt;
	std::stringstream					ssOut;
		
	::VariantTimeToSystemTime(m_nExcelDate, &lt);
	ssOut << lt.wYear << '-' << s_aszMonths[lt.wMonth - 1] << '-' << lt.wDay;
    return ssOut.str();
}

//	Performs the reverse operation of AddTenor
std::string MlEqDate::GetTenor(long nEndDate, BusinessDayConventionEnum bdc /*= NoBusinessDayConvention*/, const std::string& szCalendar /*= ""*/, bool bPreserveEoM /*= false*/) const
{
	long								nStartDate = m_nExcelDate;	
	estring								szTenor;
	SYSTEMTIME							stStart, stEnd;
		
	if (bdc != NoBusinessDayConvention) throw "MlEqDate::GetTenor() does not support a specified business day convention";
	if (szCalendar.size()) throw "MlEqDate::GetTenor() does not support a specified calendar";
	if (bPreserveEoM) throw "MlEqDate::GetTenor() does not support a true value for bPreserveEoM";
		
	if (nStartDate == nEndDate) return "0d";
			
	::VariantTimeToSystemTime(nStartDate, &stStart);
	::VariantTimeToSystemTime(nEndDate, &stEnd);
	if (stStart.wDay == stEnd.wDay && stStart.wMonth == stEnd.wMonth){
		// integral years
		return estring(stEnd.wYear - stStart.wYear) + "y";
	} else if (IsEqualDayOfMonth(stStart, stEnd)){
		// integral months
		return estring((stEnd.wYear - stStart.wYear) * 12 + (stEnd.wMonth - stStart.wMonth)) + "m";
	} else if (stStart.wDayOfWeek == stEnd.wDayOfWeek){
		// integral weeks		
		return estring((nEndDate - nStartDate) / 7) + "w";
	}
	
	// default to integral days		
	return estring(nEndDate - nStartDate) + "d";	
}

long MlEqDate::GetYear(void) const
{
	SYSTEMTIME							lt;	
	::VariantTimeToSystemTime(m_nExcelDate, &lt);
	return lt.wYear;
}

double MlEqDate::GetYearFraction(long nDate) const
{
	MlEqDateYearFractionCache::const_iterator it = m_mapYearFractions.find(nDate);
	if (it != m_mapYearFractions.end()){
		double f = it->second; 
		return f;
	} else {
		std::string szDcc;
		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		hdeParams.insert("StartDate", GDA::HDElement(m_nJulianDate));
		hdeParams.insert("EndDate", GDA::HDElement(ExcelDateToJulianDate(nDate)));
		CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, m_dcc, &szDcc);
		szDcc.assign(SiriusToGDADayCountConvention(szDcc));		
		hdeParams.insert("DCM", GDA::HDElement(szDcc.data()));		
		GDA::HDElement		hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.YearFrac", hdeParams, GDA::HDElement("YearFrac"), GDA::GetDefaultContext()).lookup("YearFrac");
		double				f = hdeResult.asDouble();
		m_mapYearFractions[nDate] = f;
		return f;
	}
}

/*static*/ unsigned long MlEqDate::ExcelDateToJulianDate(long nExcelDate)
{
	if (nExcelDate <= 60){
		// Microsoft thinks that 1900 was a leap year! Day 60 is 29-Feb-1900.
		return nExcelDate + 2415019 + 1;		
	} else {	
		return nExcelDate + 2415019;
	}
}

/*static*/ bool MlEqDate::IsBusinessDay(long nExcelDate, const std::string& szCalendarOpt)
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);	
	hdeParams.insert("Date", GDA::HDElement(ExcelDateToJulianDate(nExcelDate)));
	if (szCalendarOpt.size()){
		hdeParams.insert("Calendar", GDA::HDElement(szCalendarOpt.data()));
	}
	GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("Date.TestDate", hdeParams, GDA::HDElement("IsGood"), GDA::GetDefaultContext()).lookup("IsGood");
	return hdeResult.asBoolean();
}

bool MlEqDate::IsEqualDayOfMonth(SYSTEMTIME stStart, SYSTEMTIME stEnd) const
{
	if (stStart.wDay == stEnd.wDay){
		return true;
	} else if (stStart.wDay < stEnd.wDay){
		// e.g. 29-Feb-2004 to 31-Mar-2004 should return true since this is one month.
		// Return true if every increment of stStart.wDay up to and including the valud of stEnd.wDay invalidates stStart.
		FILETIME ftDummy;
		while (++stStart.wDay <= stEnd.wDay){
			if (::SystemTimeToFileTime(&stStart, &ftDummy)) return false;
		}		
		return true;
	} else /*stStart.wDay > stEnd.wDay*/ {
		// e.g. 31-Jan-2001 to 28-Feb-2001 should return true since this is one month.
		// Return true if every increment of stEnd.wDay up to and including the value of stStart.wDay invalidates stEnd.
		FILETIME ftDummy;
		while (++stEnd.wDay <= stStart.wDay){
			if (::SystemTimeToFileTime(&stEnd, &ftDummy)) return false;
		}
		return true;
	}
	return false;
}

/*static*/ long MlEqDate::JulianDateToExcelDate(unsigned long nJulianDate)
{
	if (nJulianDate <= 2415019 + 60){
		// Microsoft thinks that 1900 was a leap year! Day 60 is 29-Feb-1900.
		return nJulianDate - 2415019 - 1;
	} else {		
		return nJulianDate - 2415019;
	}
}

void MlEqDate::PutDate(long nDate)
{
	m_mapYearFractions.clear();
	m_nExcelDate = nDate;
	m_nJulianDate = ExcelDateToJulianDate(m_nExcelDate);
}

void MlEqDate::PutDayCountConvention(DayCountConventionEnum dcc)
{
	m_mapYearFractions.clear();
	m_dcc = dcc;
}

void MlEqDate::PutDayOfMonth(long nDay)
{
	SYSTEMTIME							lt;
	double								fNewDate;	
	::VariantTimeToSystemTime(m_nExcelDate, &lt);
	lt.wDay = nDay;
	if (!::SystemTimeToVariantTime(&lt, &fNewDate)){
		// Failure
		throw "Invalid date '" + GetString(lt) + "'";
	}
	PutDate(fNewDate);
}

// Ensures the date is a business day
MlEqDateHandle MlEqDate::SetToBusinessDay(BusinessDayConventionEnum bdc, std::string& szCalendar) const
{
	return Add(0, 0, 0, bdc, szCalendar, true);
}

//
//	Transforms Sirius-Style day count conventions that GDA does not
//	recognise to GDA equivalents.
//
/*static*/ std::string MlEqDate::SiriusToGDADayCountConvention(const std::string& sz)
{
	if (sz == "Thirty360"){
		return "30/360";
	} else if (sz == "Act365F"){
		return "Act/365F";
	}
	return sz;
}