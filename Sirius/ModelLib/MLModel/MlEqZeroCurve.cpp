//	MlEqZeroCurve.cpp : Implementation of MlEqZeroCurve
//
//	Author :            David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqzerocurve.h"
#include "utility.h"
#include "mleqdate.h"
#include "MlEqInterpolator.h"

MlEqZeroCurve::MlEqZeroCurve(void) : m_ds(NoDataSource), m_nCachedReferenceDate(0L)
{
}

/*explicit*/ MlEqZeroCurve::MlEqZeroCurve(GDA::Functor_const_ref hYieldCurve) : m_ds(NoDataSource), m_nCachedReferenceDate(0L), m_hYieldCurve(hYieldCurve)
{
}

/*virtual*/ MlEqZeroCurve::~MlEqZeroCurve(void)
{
}

// ToDo - remove all this once we are calling the ocean libraries directly
void MlEqZeroCurve::ClearCache(void) const
{
	m_mapDiscountFactors.clear();	
	m_nCachedReferenceDate = 0L;	
}

/*static*/ std::string MlEqZeroCurve::Default(const std::string& szCurrency, const std::string& szRequest)
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("Currency", GDA::HDElement(szCurrency.data()));
	GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveDefault", hdeParams, GDA::HDElement(szRequest.data()), GDA::GetDefaultContext()).lookup(szRequest.c_str());
	if (hdeResult.isException()) hdeResult.raise();
	return hdeResult.asCStr();
}

double MlEqZeroCurve::GetCashRate(const std::string& szTermOrDate) const
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	hdeParams.insert("Term", GDA::HDElement(szTermOrDate.data()));
	return GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.CashRate", hdeParams, GDA::HDElement("CashRate"), GDA::GetDefaultContext()).lookup("CashRate").asDouble();
}

DataSourceEnum MlEqZeroCurve::GetDataSource(void) const
{
	return m_ds;
}

double MlEqZeroCurve::GetDiscountFactor(long nMaturity) const
{					
	MlEqZeroCurveDiscountFactorCache::const_iterator it = m_mapDiscountFactors.find(nMaturity);
	if (it != m_mapDiscountFactors.end()){
		double res = it->second; 
		return res;
	} else {
		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
		hdeParams.insert("Date", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nMaturity)));
		GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.DF", hdeParams, GDA::HDElement("DF"), GDA::GetDefaultContext()).lookup("DF");
		if (hdeResult.isException()) hdeResult.raise();
		double f = hdeResult.asDouble();
		m_mapDiscountFactors[nMaturity] = f;
		return f;
	}
}

double MlEqZeroCurve::GetDiscountFactor(long nFrom, long nTo) const
{
	double f = GetDiscountFactor(nTo) / GetDiscountFactor(nFrom);
	return f;
}

double MlEqZeroCurve::GetFxSpot(void) const
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	return GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveProperty", hdeParams, GDA::HDElement("Fx"), GDA::GetDefaultContext()).lookup("Fx").asDouble();
}

std::string MlEqZeroCurve::GetInterpolator(void) const
{
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));	
	return GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveInfo", hdeParams, GDA::HDElement("Interpolator"), GDA::GetDefaultContext()).lookup("Interpolator").asCStr();
}

std::string MlEqZeroCurve::GetName(void) const
{	
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	return GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveInfo", hdeParams, GDA::HDElement("Currency"), GDA::GetDefaultContext()).lookup("Currency").asCStr();
}

long MlEqZeroCurve::GetReferenceDate(void) const
{
	if (!m_nCachedReferenceDate){
		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));	
		m_nCachedReferenceDate = MlEqDate::JulianDateToExcelDate(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveInfo", hdeParams, GDA::HDElement("CurveDate"), GDA::GetDefaultContext()).lookup("CurveDate").asJulian());
	}
	return m_nCachedReferenceDate;
}

GDA::HDElement MlEqZeroCurve::GetYieldCurve(void) const
{
	return GDA::HDElement(m_hYieldCurve);
}

void MlEqZeroCurve::PutFxSpot(double fFxSpot)
{
	ClearCache();
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	hdeParams.insert("Fx", GDA::HDElement(fFxSpot));	
	m_hYieldCurve = GDA::Functor_const_ref(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.ModifyCurve", hdeParams, GDA::HDElement("YieldCurve"), GDA::GetDefaultContext()).lookup("YieldCurve").asConstObject());	
}

void MlEqZeroCurve::PutInterpolator(const std::string& szInterpolator)
{
	ClearCache();
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	hdeParams.insert("Interpolator", GDA::HDElement(szInterpolator.data()));
	m_hYieldCurve = GDA::Functor_const_ref(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.ModifyCurve", hdeParams, GDA::HDElement("YieldCurve"), GDA::GetDefaultContext()).lookup("YieldCurve").asConstObject());	
}

void MlEqZeroCurve::PutDataSource(DataSourceEnum ds)
{	
	m_ds = ds;
}

void MlEqZeroCurve::PutReferenceDate(long nDate)
{
	ClearCache();
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	hdeParams.insert("CurveDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nDate)));	
	m_hYieldCurve = GDA::Functor_const_ref(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.ModifyCurve", hdeParams, GDA::HDElement("YieldCurve"), GDA::GetDefaultContext()).lookup("YieldCurve").asConstObject());	
}

void MlEqZeroCurve::PutShift(double fShift)
{
	ClearCache();
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
	hdeParams.insert("Scenario", GDA::HDElement("ShiftAll"));
	hdeParams.insert("Shock", GDA::HDElement(fShift));
	m_hYieldCurve = GDA::Functor_const_ref(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.ShockCurve", hdeParams, GDA::HDElement("YieldCurve"), GDA::GetDefaultContext()).lookup("YieldCurve").asConstObject());	
}

void MlEqZeroCurve::PutShift(const MlEqInterpolatorHandle& hInterpolator)
{		
	ClearCache();
	long					nInstruments = 0L;	
	long					nEndDate = 0L; 
	double					fShockVal(0.0); 					
	GDA::HDElement			hdeCurveInfoParams(GDA::HDElement::Dictionary);
	GDA::HDElement			hdeCurveInstDatesParams(GDA::HDElement::Dictionary);

	if (!hInterpolator) throw "No interpolator supplied for void MlEqZeroCurve::PutShift";
	hdeCurveInfoParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
		
	GDA::HDElement  hdeCurveInstDates(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.YieldCurveInfo", hdeCurveInfoParams, GDA::HDElement("CurveInstDates"), GDA::GetDefaultContext()));

	if (hdeCurveInstDates.isException()) hdeCurveInstDates.raise();

	hdeCurveInstDates = GDA::HDElement::dictToArray(hdeCurveInstDates);
		
	nInstruments = hdeCurveInstDates.rows();
	
	// Loop through all the instruments and shock them accordingly with respect to the Interpolator
	for (long idx = 0; idx < nInstruments; ++idx){

		nEndDate = MlEqDate::JulianDateToExcelDate(hdeCurveInstDates(idx, 1).asJulian());
		fShockVal = hInterpolator->getValue(nEndDate);

		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		hdeParams.insert("YieldCurve", GDA::HDElement(m_hYieldCurve));
		hdeParams.insert("Scenario", GDA::HDElement("Add"));
		hdeParams.insert("Shock", GDA::HDElement(fShockVal));
		hdeParams.insert("Desc", GDA::HDElement(idx));

		m_hYieldCurve = GDA::Functor_const_ref(GDA::GetLibraryInstance()->getFactory()->createAndApply("Curve.ShockCurve", hdeParams, GDA::HDElement("YieldCurve"), GDA::GetDefaultContext()).lookup("YieldCurve").asConstObject());
	}	
}

void MlEqZeroCurve::PutYieldCurve(GDA::Functor_const_ref hYieldCurve)
{
	ClearCache();
	m_hYieldCurve_original = m_hYieldCurve = hYieldCurve;	
	#ifdef _DEBUG
		m_szName.assign(GetName());
	#endif
}

void MlEqZeroCurve::Reset(void)
{
	ClearCache();
	m_hYieldCurve = m_hYieldCurve_original;	
}

void MlEqZeroCurve::Stick(void)
{
	ClearCache();
	m_hYieldCurve_original = m_hYieldCurve;
}