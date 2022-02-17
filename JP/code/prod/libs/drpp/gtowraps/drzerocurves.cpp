// drzerocurves.cpp: implementation of the drzerocurves class.
//
//////////////////////////////////////////////////////////////////////

#include "drzerocurves.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const DRString KAPITAL_DISC_ZERO	= "disczero.dat";	// libor
const DRString KAPITAL_INDEX_ZERO		= "zero.dat";		// cmt

DRZeroCurves::DRZeroCurves(FileType type, DRDate& envDate, DRString cmtFile, 
						   DRString yldFile)
						   :m_parRatesPtr(NULL), m_usingParRates(false), 
						   m_indexCurves(1)
{
	if (type == SWAP) {
		m_usingParRates = true;
		m_parRatesPtr = DRParRatesPtr(new DRParRates(envDate, cmtFile, yldFile));
		SetZerosFromParRates(envDate);
	}
	else {
		m_discountCurve	= DRCurve(KAPITAL_DISC_ZERO);
		m_indexCurves.operator[](0)= DRCurve(KAPITAL_INDEX_ZERO);
	}
}

DRZeroCurves::DRZeroCurves(TCurve* liborCurve, TCurve* cmtCurve)
						   :m_parRatesPtr(NULL), m_usingParRates(false),
						   m_indexCurves(1)
{
	m_discountCurve	= DRCurve(liborCurve);
	m_indexCurves.operator[](0)= DRCurve(cmtCurve);
}


DRZeroCurves::DRZeroCurves(DRDate& envDate, MbsIO& mbsIO)
:m_usingParRates(false), m_parRatesPtr(NULL), m_indexCurves(1)
{
	int i;
	
	MbsIOMap& mbsIOMap = mbsIO.get_map();
	
	MbsIOVector& liborDateVec = mbsIOMap["LIBORDATES"].get_vector();
	MbsIOVector& liborRateVec = mbsIOMap["LIBORRATES"].get_vector();
	
	if (liborDateVec.size() != liborRateVec.size()) 
		throw DRException ("Size mismatch");
	
	DRDateList liborDates;
	DArray liborRates (liborDateVec.size());
	
	for (i = 0; i < liborDateVec.size(); i++) {
		liborDates.AddDate (toDate(liborDateVec[i].get_string()));
		liborRates[i] = toNumber(liborRateVec[i].get_string()) / 100;
	}
	
	double liborBasis = toNumber (mbsIOMap.get("LiborBasis"));
	
	long liborDayCountConv = toDayCount (mbsIOMap.get("LiborDayCount"));
	
	m_discountCurve	= 	DRCurve(liborDates, liborRates, envDate, liborBasis, 
		liborDayCountConv);

	MbsIOVector& cmtDateVec = mbsIOMap["CMTDATES"].get_vector();
	MbsIOVector& cmtRateVec = mbsIOMap["CMTRATES"].get_vector();
	
	if (cmtDateVec.size() != cmtRateVec.size()) 
		throw DRException ("Size mismatch");
	
	DRDateList cmtDates;
	DArray cmtRates (cmtDateVec.size());
	
	for (i = 0; i < cmtDateVec.size(); i++) {
		cmtDates.AddDate (toDate(cmtDateVec[i].get_string()));
		cmtRates[i] = toNumber(cmtRateVec[i].get_string()) / 100;
	}
	
	double cmtBasis = toNumber (mbsIOMap.get("CMTBasis"));
	
	long cmtDayCountConv = toDayCount (mbsIOMap.get("CMTDayCount"));
	
	m_indexCurves.operator[](0)	= 	DRCurve(cmtDates, cmtRates, envDate, cmtBasis, 
		cmtDayCountConv);	
}


bool DRZeroCurves::operator==(const DRZeroCurves& a) const
{
	bool ans = (m_indexCurves.size() == a.m_indexCurves.size()) && 
		(m_discountCurve == a.m_discountCurve);
	
	for (int i = 0; ans && i < m_indexCurves.size(); i++) {
		ans &= m_indexCurves[i] == a.m_indexCurves[i];
	}

	return ans;
}


void DRZeroCurves::SetZerosFromParRates(DRDate& envDate)
{
	int size = m_parRatesPtr->parRates().size();
	DArray cmtRates(size), liborRates(size);
	
	cmtRates	= m_parRatesPtr->parRates() + m_parRatesPtr->cmt10Yhumps();
	liborRates	= m_parRatesPtr->parRates() + m_parRatesPtr->swapSpreads();
	
	DRDateList dates = m_parRatesPtr->dates();
	m_discountCurve	= DRCurve(envDate, liborRates, dates);
	m_indexCurves.operator[](0)= DRCurve(envDate, cmtRates, dates);
}

void DRZeroCurves::TweakParRate(Tweak tweak, double bpTweak)
{
	if (!m_usingParRates) 
		throw DRException("Can't tweak par rates, cause not using par rates");

	m_parRatesPtr->TweakParRate(tweak, bpTweak);
	DRDate envDate = m_discountCurve.baseDate();
	
	SetZerosFromParRates(envDate);
}

double DRZeroCurves::GetForwardRate (DRRate& rate, DRDate& date)
{
	double value;
	
	if (GtoFloatRateArrayValue (rate, date, m_discountCurve, m_indexCurves, 
		1, &value) ISNT SUCCESS)
		throw DRException("Some kinda of failure in figuring out rate");
	
	return value;
}

ostream& operator<< (ostream& s, DRZeroCurves& a)
{
	s << "Discount Curve" << endl << a.m_discountCurve;

	s << "Index Curves" << endl << a.m_indexCurves;

	s << "Using Par Rates:\t" << a.m_usingParRates << endl;

	return s;
}


