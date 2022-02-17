//	MlEqDividendSchedule.cpp : Implementation of the Dividend Schedule
//
//	Authors :				   David Cuin / Romain Barc
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqdividendschedule.h"
#include "mleqzerocurve.h"

/*static*/ const long					MlEqDividendSchedule::s_nMaxDate(2958465L);			// corresponds to 31-Dec-9999

MlEqDividendSchedule::MlEqDividendSchedule(void) : Discrete(m_discrete), Continuous(m_continuous), m_dcc(ActualActual), m_ds(NoDataSource)
{
}

/*explicit*/ MlEqDividendSchedule::MlEqDividendSchedule(const std::string& szName, DataSourceEnum ds, long nDate, DayCountConventionEnum dcc, MlEqZeroCurveHandle hForeignCurve, const std::map<long, double> mapDiscrete, const std::map<long, double> mapContinuous) :
	Discrete(m_discrete),
	Continuous(m_continuous), 	
	m_szName(szName),
	m_ds(ds),
	m_dcc(dcc),
	m_hForeignCurve(hForeignCurve),
	m_discrete(mapDiscrete, nDate),
	m_continuous(mapContinuous, nDate),
	m_discrete_original(mapDiscrete),
	m_continuous_original(mapContinuous)
{
}

//	Extrapolates existing dividends up to m_nDate + szContinuousYieldAfter using an input growth
//  rate. Then, any discrete terms after szContinuousYieldAfter are removed.
void MlEqDividendSchedule::ApplyGrowth(double fGrowth, const std::string& szContinuousYieldAfter)
{	
	if (m_continuous.size()) throw "MlEqDividendSchedule::ApplyGrowth cannot be applied to dividend schedules with continuous yield terms";
	if (m_discrete.size() <= 1) return;		// Nothing to apply. We need at least two points to compute a growth rate.
				
	long								nSource_Last = m_discrete.GetLastDate();	
	long								nSource_Penultimate = m_discrete.GetPenultimateDate();
	long								nPeriod = nSource_Last - nSource_Penultimate;
	long								nFrom = GetDate();
	long								nTo = MlEqDate(nFrom, m_dcc).AddTenor(szContinuousYieldAfter)->GetDate() - 1;
	std::map<long, double>::iterator	it;
		
	// Extrapolate dividends
	for (std::map<long, double>::const_iterator it = m_discrete.begin(); it != m_discrete.end(); ++it){
		// Only add dividend terms if their dates are in the rance nSource_Last + nPeriod to nTo
		long nAt = MlEqDate(it->first, m_dcc).AddTenor("1y")->GetDate();
		if (nAt > nTo) break;
		if (nAt <= nSource_Last + nPeriod) continue;
		double fValue = it->second * (1.0 + fGrowth);
		m_discrete[nAt] = fValue;		
	}
				
	// Delete any discrete dividends after nTo
	while (m_discrete.size()){
		it = m_discrete.end();
		--it;
		if (it->first <= nTo) break;
		m_discrete.erase(it);
	}	
}

// Remove all the continuous yields, and insert a single final yield growth term.
void MlEqDividendSchedule::CalibrateYield(double fSpot, MlEqZeroCurveHandle hzc, const std::string& szTenor, const std::string& szAddAt)
//	szTenor - this is the range we consider when deriving (enter as 1y etc).
//  szAddAt - this is the tenor (w.r.t. the dividend schedule date) where the yield term is added (typically "50y")
{	
	if (!m_continuous.size() && !m_discrete.size()) return;		// Nothing to calibrate
	
	if (!m_discrete.size()) throw "The dividend schedule '" + m_szName + "' does not contain any discrete dividends";
		
	long								nYieldDate = MlEqDate(GetDate(), m_dcc).AddTenor(szAddAt)->GetDate();				
	MlEqDateHandle						hDateTo = new MlEqDate(m_discrete.GetLastDate(), m_dcc);
	MlEqDateHandle						hDateFrom = hDateTo->AddTenor("-" + szTenor);

	if (nYieldDate <= m_discrete.GetLastDate()) throw "The AddAt tenor value is too small - discrete dividends occur this date";
	
	double								t = hDateFrom->GetYearFraction(hDateTo->GetDate());
	double								fPV1 = GetPresentValue(hzc, GetDate(), hDateFrom->GetDate());
	double								fPV2 = GetPresentValue(hzc, GetDate(), hDateTo->GetDate());	
	double								fYield = -log((1.0 - fPV2 / fSpot) / (1.0 - fPV1 / fSpot)) / t;

	// Add this value to the point defined by GetDate() + szAddAt
	m_continuous.clear();
	PutContinuousValueAt(m_discrete.GetLastDate(), 0.0);	
	PutContinuousValueAt(nYieldDate, fYield);
}

//
//  Here we compute the largest integer less than or equal to the year fraction between
//  the dividend schedule date and each discrete dividend point. Then we multiply it 
//  by (1 + nYearFraction * Amount).
//
//  The trick here is to compute the year fractions a minimum number of times.
//
void MlEqDividendSchedule::DiscreteShift(double fGrowth, const std::string& szTenor)
{	
	MlEqDateHandle						hDate = new MlEqDate(GetDate(), m_dcc);
	long								nYearFraction = 0;
	MlEqDateHandle						hDate_NextPeriod = hDate->AddTenor(szTenor);
		
	for (std::map<long, double>::iterator it = m_discrete.begin(); it != m_discrete.end(); ++it){
		if (it->first >= hDate_NextPeriod->GetDate()){
			hDate = hDate_NextPeriod;
			hDate_NextPeriod = hDate->AddTenor(szTenor);
			++nYearFraction;
		}
		it->second *= (1.0 + (double)nYearFraction * fGrowth);
	}
}

DataSourceEnum MlEqDividendSchedule::GetDataSource(void) const
{
	return m_ds;
}

long MlEqDividendSchedule::GetDate(void) const
{	
	if (m_discrete.GetDate() != m_continuous.GetDate()) throw "Inconsistent dates encountered in the dividend schedule '" + m_szName + "'";
	return m_discrete.GetDate();		
}

DayCountConventionEnum MlEqDividendSchedule::GetDayCountConvention(void) const
{
	return m_dcc;
}

estring MlEqDividendSchedule::GetDayCountConventionStr(void) const
{
	return CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, GetDayCountConvention());
}

std::string MlEqDividendSchedule::GetName(void) const
{
	return m_szName;
}

//	GetForwad - This is the lowest level forward function in Sirius
/*static*/ double MlEqDividendSchedule::GetForward(long nStart, long nMaturity, double forwardStartSpot, MlEqZeroCurveHandle hzc, MlEqDividendScheduleHandle hDivs, MlEqZeroCurveHandle hzcFXYield)
// This is the lowest level implementer of GetNaturalForward and GetQuantoForward functions; all overloads must call into this
{	
	double								fPresentValue(0.0);
	double								fYield(0.0);
	double								fDiscountFactor(1.0);
	double								fFXYieldDiscountFactor(1.0);
	double								fYearFraction(0.0);

	if (nMaturity < nStart){
		throw "Error in returning the forward. The input start date (" + MlEqDate(nStart).GetString() + ") is more recent than the maturity (" + MlEqDate(nMaturity).GetString() + ")";
	}
	
	// Check for consistency between hzc, hDivs and hzcFXYield
	if (!!hzc && !!hDivs && hzc->GetReferenceDate() != hDivs->GetDate()){
		throw "Cannot obtain the forward value since the zero curve date ('" + hzc->GetName() + "', " + MlEqDate(hzc->GetReferenceDate()).GetString() + ") differs from the dividend schedule date ('" + hDivs->GetName() + "', " + MlEqDate(hDivs->GetDate()).GetString() + ")";
	}
	if (!!hzc && !!hzcFXYield && hzc->GetReferenceDate() != hzcFXYield->GetReferenceDate()){
		throw "Cannot obtain the forward value since the zero curve date ('" + hzc->GetName() + "', " + MlEqDate(hzc->GetReferenceDate()).GetString() + ") differs from the FX zero curve date ('" + hzcFXYield->GetName() + "', " + MlEqDate(hzcFXYield->GetReferenceDate()).GetString() + ")";
	}

	if (!!hzc) fDiscountFactor = hzc->GetDiscountFactor(nStart, nMaturity);
	if (!!hzcFXYield) fFXYieldDiscountFactor = hzcFXYield->GetDiscountFactor(nStart, nMaturity);
	if (!!hDivs){
		fPresentValue = hDivs->GetPresentValue(hzc, nStart, nMaturity);
		fYield = hDivs->GetYield(nStart, nMaturity);		
		MlEqDate date(nStart,  hDivs->GetDayCountConvention());
		fYearFraction = date.GetYearFraction(nMaturity);
	}
	
	double f = (forwardStartSpot - fPresentValue) / fDiscountFactor * fFXYieldDiscountFactor * exp(-fYearFraction * fYield);
	return f;
}

//  returns the present value of the discrete dividends
double MlEqDividendSchedule::GetPresentValue(MlEqZeroCurveHandle hzc, long nFrom, long nTo, long nDiscountFromDate)
{		
	double								fPresentValue(0.0);
	std::vector<long>					v;

	if (!nDiscountFromDate) nDiscountFromDate = nFrom;
	m_discrete.GetDates(v);

	// Discrete part:	Advance to the first dividend in the range nFrom < date <= nTo. We do NOT include a dividend payable on nDiscountFromDate in the result.
	std::map<long, double>::const_iterator it = m_discrete.begin();
	for (; it != m_discrete.end() && it->first <= nFrom; ++it);
	for (; it != m_discrete.end() && it->first <= nTo; ++it){
		fPresentValue += (it->second * (!!hzc ? hzc->GetDiscountFactor(nDiscountFromDate, it->first) : 1.0));		// ToDo - this is slow due to ? : 
	}
	return fPresentValue;
}

// this function returns a local dividend rate (just the number in the dividend schedule, not the average yield)
double MlEqDividendSchedule::GetRate(long nDate)
{
	std::map<long, double>::const_iterator it = m_continuous.begin();
	long itBegin = it->first;
	for (; it != m_continuous.end() && it->first < nDate; ++it);
		
	if (it == m_continuous.end()){
		return 0.0;
	}
	double f = it->second;
	return f;
}

//	returns the average yield of the continuous divideds
double MlEqDividendSchedule::GetYield(long nFrom, long nTo)
{
	if (nFrom == nTo) return 0.0;
	
	// Continuous part:	Advance to the first yield in the range From < date <= To. We do not need any yield term occuring on nFrom.
	std::map<long, double>::const_iterator it = m_continuous.begin();
	long itBegin = it->first;
	for (; it != m_continuous.end() && it->first <= nFrom; ++it);
	if (it == m_continuous.end()) return 0.0;	

	double fYield;
				
	// Starting contribution.
	if (it->first > nTo){		
		// First dividend term is beyond nTo.
		fYield = it->second * (nTo - nFrom);
	} else {
		// First contribution is from nFrom to it->first with yield it->second.		
		fYield = it->second * (it->first - nFrom);
	}
	
	// Middle contribution(s)	
	long nPrevious = it->first;
	for (; it != m_continuous.end() && it->first <= nTo; ++it){
		int itfirst=it->first;
		double itsecond=it->second;
		fYield += (it->second * (it->first - nPrevious));
		nPrevious = it->first;
	}
	
	
	// Final contribution
	if (it->first > nTo && nTo > nPrevious){		
		fYield += (it->second * (nTo - nPrevious));
	}

	// compute and return the average yield
	double f = fYield / (nTo - nFrom);
	return f;
}

MlEqZeroCurveHandle MlEqDividendSchedule::GetForeignCurve(void) const
{
	return m_hForeignCurve;
}

void MlEqDividendSchedule::PutContinuous(const std::map<long, double>& mapContinuous)
{	
	m_continuous.PutSchedule(mapContinuous);
	m_continuous_original = mapContinuous;
}

void MlEqDividendSchedule::PutContinuousValueAt(long nDate, double fValue)
{
	m_continuous_original[nDate] = m_continuous[nDate] = fValue;
}

void MlEqDividendSchedule::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

void MlEqDividendSchedule::PutDate(long nDate)
{
	m_discrete.PutDate(nDate);
	m_continuous.PutDate(nDate);	
}

void MlEqDividendSchedule::PutDayCountConvention(DayCountConventionEnum dcc)
{
	m_dcc = dcc;
}

void MlEqDividendSchedule::PutDiscrete(const std::map<long, double>& mapDiscrete)
{
	m_discrete.PutSchedule(mapDiscrete);
	m_discrete_original = mapDiscrete;
}

void MlEqDividendSchedule::PutDiscreteValueAt(long nDate, double fValue)
{
	m_discrete_original[nDate] = m_discrete[nDate] = fValue;
}

void MlEqDividendSchedule::PutName(const std::string& szName)
{
	m_szName.assign(szName);
}

void MlEqDividendSchedule::PutForeignCurve(MlEqZeroCurveHandle hzc)
{
	m_hForeignCurve = hzc;
}

void MlEqDividendSchedule::Reset(void)
{
	m_discrete.PutSchedule(m_discrete_original);
	m_continuous.PutSchedule(m_continuous_original);
}

// Add the fYield value to all the continuous yield terms.
void MlEqDividendSchedule::Shift(double fYield)
{	
	for (std::map<long, double>::iterator it = m_continuous.begin(); it != m_continuous.end(); ++it){
		it->second += fYield;
	}
	// Now add the yield value term to m_nMaxDate.
	if (!m_continuous.size()){
		m_continuous[s_nMaxDate] = fYield;
	} else {
		std::map<long, double>::reverse_iterator	it = m_continuous.rbegin();
		long										nDate = it->first;
		if (nDate >= s_nMaxDate){
			// We have already added the contribution
		} else {
			m_continuous[s_nMaxDate] = fYield;
		}
	}
}

//	Subsequent resets on this object will revert to this state.
void MlEqDividendSchedule::Stick(void)
{	
	m_discrete_original = m_discrete;				
	m_continuous_original = m_continuous;	
}