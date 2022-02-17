//	MlEqDateSchedule.cpp :		Implementation of the date schedule class
//
//	Author :					David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqdateschedule.h"

#undef max

void MlEqDateSchedule::Create(long nStartDate, long nEndDate, const std::string& szFrequency, const std::string& szCalendar, BusinessDayConventionEnum bdc)
{
	std::string							szBdc;
	GDA::Functor_const_ref				hDateSchedule;
	
	CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdc, &szBdc);	
	GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
	hdeParams.insert("StartDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nStartDate)));
	hdeParams.insert("EndDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nEndDate)));
	hdeParams.insert("Frequency", GDA::HDElement(szFrequency.data()));
	hdeParams.insert("Calendar", GDA::HDElement(szCalendar.data()));
	hdeParams.insert("BDC", GDA::HDElement(szBdc.data()));
	hDateSchedule = GDA::GetLibraryInstance()->getFactory()->create("Date.Schedule", hdeParams, GDA::GetDefaultContext());
	GDA::HDElement hdeResult = hDateSchedule->apply(GDA::HDElement("Strip")).lookup("Strip");
	
	// build the schedule from the hdeResult function
	clear();
	for (unsigned int nRow = 0; nRow < hdeResult.rows(); nRow++){
		long n = MlEqDate::JulianDateToExcelDate(hdeResult(nRow, 0).asJulian());
		insert(std::pair<long, double>(n, 0.0));
	}
}

// returns an array corresponding to a column in the date schedule
void MlEqDateSchedule::GetColumnBeforeToday(long nToday, long nColumn, std::vector<double>& af) const
{	
	af.resize(size());	// This is the maximum size of the output
	long nElements = 0;
	for (std::map<long, std::vector<double> >::const_iterator it = begin(); it != end(); ++it){
		if (it->first >= nToday){
			af.resize(nElements);
			return;
		}		
		if (nColumn >= it->second.size()) throw "No value is defined at " + MlEqDate(it->first).GetString() + " in column " + estring(nColumn + 1);
		af[nElements++] = (it->second)[nColumn];
	}
}

void MlEqDateSchedule::GetColumnBeforeToday(long nToday, long nColumn, CMatrix& m) const
//	m - we size this to one column
{
	std::vector<double>					af;
	GetColumnBeforeToday(nToday, nColumn, af);
	m.resize(1, af.size(), 0.0);
	for (long n = 0; n < af.size(); n++){
		m[0][(int)n] = af[n];
	}
}

void MlEqDateSchedule::GetColumnsBeforeToday(long nToday, CMatrix& m) const
{	
/*	long								nElementsConsidered;
	long								nWidth = GetWidth(nToday, &nElementsConsidered);
	int									nRow = 0;

	ATLASSERT(false);	// ToDo - test
	m.resize(nElementsConsidered, nWidth);
	for (std::map<long, std::vector<double> >::const_iterator it = begin(); it != end(); ++it){
		for (long nCol = 0; nCol < it->second.size(); nCol++){
			m[nRow][(int)nCol] = it->second[nCol];
		}
		nRow++;
	}
*/
	long								nElementsConsidered;
	long								nWidth = GetWidth(nToday, &nElementsConsidered);
	int									nRow = 0;

//	ATLASSERT(false);	// ToDo - test
	m.resize(nWidth, nElementsConsidered);

	for (std::map<long, std::vector<double> >::const_iterator it = begin(); it != end(); ++it){
		for (long nCol = 0; nCol < it->second.size(); nCol++){
			m[(int)nCol][nRow] = it->second[nCol];
		}
		nRow++;
	}

}

// returns the length of the longest std::vector<double> in the schedule
long MlEqDateSchedule::GetWidth(long nIncludeBefore /*= std::numeric_limits<long>::max()*/, long* pnElementsConsidered /*= NULL*/) const
//	nIncludeBefore - we only consider elements with dates before this value
//	pnElementsConsidered (returned, nullable) - returns the nummber of elements in the schedule before nIncludeBefore
{
	long								nWidth = 0;
	if (pnElementsConsidered) *pnElementsConsidered = 0;
	for (std::map<long, std::vector<double> >::const_iterator it = begin(); it != end(); ++it){		
		if (it->first >= nIncludeBefore) return nWidth;
		if (pnElementsConsidered) (*pnElementsConsidered)++;
		long nSize = it->second.size();
		nWidth = std::max(nWidth, nSize);
	}
	return nWidth;
}