// drcurve.cpp: implementation of the drcurve class.
//
//////////////////////////////////////////////////////////////////////


#include "drcurve.h"
#include "drsymboltable.h"

extern "C" {
#include "gtonpi.h"
}
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DRCurve::DRCurve()
{
	m_tcurve = NULL;
}

DRCurve::~DRCurve()
{
	if (m_tcurve) 
		GtoFreeTCurve (m_tcurve);
}

DRCurve::DRCurve (TCurve* tcurve)
{
	if (tcurve)
		m_tcurve = GtoCopyCurve (tcurve);
	else 
		m_tcurve = NULL;
}


DRCurve::DRCurve(const DRCurve& rhs)
{
	if (rhs.m_tcurve == NULL) m_tcurve = NULL;	
	else m_tcurve = GtoCopyCurve (rhs.m_tcurve);
}

DRCurve& DRCurve::operator=(const DRCurve& rhs)
{
	if (m_tcurve) GtoFreeTCurve (m_tcurve);
	
	if (rhs.m_tcurve)
		m_tcurve = GtoCopyCurve (rhs.m_tcurve);
	else 
		m_tcurve = NULL;

	return *this;
}

DRCurve::DRCurve(DRDateList& dates, DArray& values, DRDate& baseDate, double basis, 
				 long dayCountConv)
{
	m_tcurve = GtoMakeTCurve(baseDate, dates, values, dates.size(), basis, 
		dayCountConv);
}

DRDateList DRCurve::GetDates() const
{
	int size = m_tcurve->fNumItems;

	LArray temp (size);

	for (int i = 0; i < size; i++) {
		temp[i] = m_tcurve->fArray[i].fDate;
	}

	DRDateList dates (temp, size);

	return dates;
}

DArray DRCurve::GetRates() const
{
	int size = m_tcurve->fNumItems;

	DArray temp (size);

	for (int i = 0; i < size; i++) {
		temp[i] = m_tcurve->fArray[i].fRate;
	}

	return temp;
}

bool DRCurve::operator==(const DRCurve& a) const
{
	bool ans;

	ans = (m_tcurve->fNumItems == a.m_tcurve->fNumItems) &&
		(m_tcurve->fBaseDate == a.m_tcurve->fBaseDate) &&
		(m_tcurve->fBasis == a.m_tcurve->fBasis) &&
		(m_tcurve->fDayCountConv == a.m_tcurve->fDayCountConv);

	for (int i = 0; ans && (i<m_tcurve->fNumItems); i++) {	
		ans &= (m_tcurve->fArray[i].fDate == a.m_tcurve->fArray[i].fDate);
		ans &= (m_tcurve->fArray[i].fRate == a.m_tcurve->fArray[i].fRate);
	}

	return ans;
}

/*
void DRCurve::GetDatesAndRates(DRDateList& dateArray, DArray& values)
{
	int size = m_tcurve->fNumItems;

	LArray temp (size);
	values.resize(size);

	for (int i = 0; i < size; i++) {
		temp[i] = m_tcurve->fArray[i].fDate;
		values[i] = m_tcurve->fArray[i].fRate;
	}
	dateArray = DRDateList (temp, size);
}
*/

void DRCurve::Print(ostream& stream) const
{
	stream << "NumPoints<" << m_tcurve->fNumItems;
	stream << "> Base Date<" << GtoFormatDate(m_tcurve->fBaseDate);
	stream << "> Basis<" << m_tcurve->fBasis;
	stream << "> YearFrac Day Count Conv: " << m_tcurve->fDayCountConv << endl;

	int i;
	for (i = 0; i<m_tcurve->fNumItems; i++) 
 		stream << GtoFormatDate(m_tcurve->fArray[i].fDate) << " " ;
	stream << endl;

	for (i = 0; i<m_tcurve->fNumItems; i++) 
		stream << m_tcurve->fArray[i].fRate << " " ;
	stream << endl;
}

DRCurve::DRCurve (DRDate& envDate, DArray& parRates, DRDateList& dates)
{
	if (parRates.size() != 25 || dates.size() != 25)
		throw DRException("Wrong size arrays -- must use standard long par rates");

	static char names[] = "MMMMMMMMMMSSSSSSSSSSSSSSS";
	DArray prices (25);
	prices = 1;

	m_tcurve = GtoNPiZC(envDate,			// (I) Value date                      
						360,				// (I) Denominator for cash instruments
						2,					// (I) Swap coupon frequency           
						GTO_B30_360,		// (I) See GtoDayCountConvention       
						parRates,			// (I) Array of cash/swap rates        
						dates,				// (I) Array of cash/swap dates      
						prices,				// (I) Array of cash/swap prices      
						names,				// (I) Array of instrument names       
						GtoSETINSTR,	// (I) Default flag: dates,names,prices
						25,					// (I) Number of cash/swap instruments 
						NULL,				// (I) Array of futures prices         
						NULL,				// (I) Array of futures dates          
						0,					// (I) Number of futures to use        
						0,					// (I) Number of FRAs 
						NULL,				// (I) End dates of FRAs 
						0,					// (I) Vol for futures to cash adj     
						0,					// (I) Futures stub method to use      
						NULL,				// (I) Futures stub data               
						0,					// (I) Number of elements in stub data 
						0,					// (I) Interp method for synthetic swaps
						0,					// (I) Interp method for overall zero cv
						'0',				// (I) For fwd smoothing, length of fwds
						0);					// (I) Iff swap dates bad-bus-day adjust
	if (!m_tcurve)
		throw DRException("GtoNPiZC Failure");

}

DRCurve::DRCurve (DRString filename)
{
	ifstream inFile(CheckPathDelimiters(theSymbolTable.get(filename)).c_str());
	if (!inFile) 
		throw DRException ("Failure to Open Kapital File ") << filename;

	char temp [20], trash[255];

	inFile.getline(trash, 255);
	TDate startDate;
	long startDateTemp;
	inFile >> startDateTemp;
	startDate = toDate(startDateTemp);

	int i;
	for (i = 0; i < 7; i++) {
		inFile >> temp;
		inFile.getline(trash, 255);
	}

	int numRates;
	inFile >> numRates;

	inFile >> temp;
	inFile.getline(trash, 255);

	MbsTDateArray dates (numRates);
	DArray rates (numRates);

	long tempDate;
	for (i = 0 ; i < numRates ; i++) {
		inFile >> tempDate;
		dates[i] = toDate(tempDate);
		inFile >> rates[i];
	}

	rates /= 100.;
	m_tcurve = GtoMakeTCurve (startDate, dates, rates, numRates, 1, GTO_ACT_365F);
}


ostream& operator<<(ostream& s, const DRCurve& a)
{
	a.Print(s);
	return s;
}
