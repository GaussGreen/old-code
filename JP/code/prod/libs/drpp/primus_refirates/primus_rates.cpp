#include "drplatdep.h"
#include FSTREAM_H
#include "primus_rates.h"
#include "drsymboltable.h"
#include "drplatdep.h"

PrimusHistRates::PrimusHistRates (DRDate date, const char* ccFile, const char* rateFile)
	: m_date(date)
{
	DRString tempCCFile = CheckPathDelimiters(theSymbolTable.get (ccFile));

	DRString tempRateFile = CheckPathDelimiters(theSymbolTable.get (rateFile));

	ifstream in (tempCCFile.c_str());
	if (!in) 
		throw DRException ("Missing Primus rate file ") << tempCCFile;

	ifstream in2 (tempRateFile.c_str());
	if (!in2) 
		throw DRException ("Missing Primus rate file ") << tempRateFile;
	
	DRDate tempDate = date;

	m_histRates = PpmRegressHistRatesNew();

	char* t1 = const_cast<char*> (tempCCFile.c_str());
	char* t2 = const_cast<char*> (tempRateFile.c_str());

    if (ppmRegressReadHistFile(t1, t2, m_histRates, 
		tempDate.GetYear(), tempDate.GetMonth(), tempDate.GetDay())
		ISNT SUCCESS)
		throw DRException ("Error in reading ppm_regress_historical files");
}

PrimusHistRates::~PrimusHistRates ()
{
	PpmRegressHistRatesFree(m_histRates);
}

DArray PrimusHistRates::Get30YrMtgRates (DRDate date, int numRates)
{
	int startIndex, numHistRates;

	CalculateIndicies (date, numRates, startIndex, numHistRates);
	
	DArray ans (numRates);
	int i;
	for (i = 0; i < numHistRates; i++) {
		ans[i] = m_histRates->mtgRates30yr[i];
	}

	double lastPoint = m_histRates->mtgRates30yrPts[numHistRates - 1];
	for (i = numHistRates; i < numRates; i++) {
		ans[i+numHistRates] = lastPoint;
	}

	return ans;
}

DRDate PrimusHistRates::GetLastHistRateDate()
{
	if (m_histRates->isLastHistRateIncomplete) {		
		int month = m_histRates->partialMonthDateMonth;
		int year = m_histRates->partialMonthDateYear;
		int day = m_histRates->partialMonthDateDay;
		return DRDate (month, day, year);
	}
	else {
		DRDate date (m_histRates->firstDateMonth, 1, m_histRates->firstDateYear);
		date += (m_histRates->nRates - 1) * ONE_MONTH;
		date.SetDay(DRDate::GetDaysInMonth(date.GetMonth()));
		return date;
	}
}

int PrimusHistRates::CalculateOffset (DRDate& date)
{
	int desiredMY = date.GetMonthYear();
	int histStartMY = m_histRates->firstDateYear * 12 + m_histRates->firstDateMonth;
	int offset = desiredMY - histStartMY;

	if (offset >= m_histRates->nRates)
		throw DRException ("Not enough Primus historical rates");

	return offset;
}

double PrimusHistRates::GetMtgRate(DRDate date, int term)
{
	int offset = CalculateOffset(date);
	double ans = (term == 30) ? m_histRates->mtgRates30yr[offset] : m_histRates->mtgRates15yr[offset];
	ans /= 100;

	return ans;
}

double PrimusHistRates::GetPoint(DRDate date, int term)
{
	int offset = CalculateOffset(date);
	double ans = (term == 30) ? m_histRates->mtgRates30yrPts[offset] : m_histRates->mtgRates15yrPts[offset];
	ans /= 100;

	return ans;
}

double PrimusHistRates::GetLatestMtgPts(int term)
{
	double ans = (term == 30) ? m_histRates->latestMtgRate30yrPts : m_histRates->latestMtgRate15yrPts;
	ans /= 100;
	return ans;
}

DArray PrimusHistRates::Get30YrPoints (DRDate date, int numRates)
{
	int startIndex, numHistRates;

	CalculateIndicies (date, numRates, startIndex, numHistRates);
	
	DArray ans (numRates);
	int i;
	for (i = 0; i < numHistRates; i++) {
		ans[i] = m_histRates->mtgRates30yrPts[i];
	}

	double lastPoint = m_histRates->mtgRates30yrPts[numHistRates - 1];
	for (i = numHistRates; i < numRates; i++) {
		ans[i+numHistRates] = lastPoint;
	}

	return ans;
}	

void PrimusHistRates::CalculateIndicies (DRDate date, int num, int& startIndex, int& numHistRates)
{
	int startMY = date.GetMonthYear();
	
	int histStartMY;
	histStartMY = m_histRates->firstDateYear * 12 + m_histRates->firstDateMonth;

	startIndex = startMY - histStartMY;
	numHistRates = m_histRates->nRates - startIndex;
}

double PrimusHistRates::GetCMSBlend(double percent1Y)
{
	double libor10yrForSpread = m_histRates->tsy10yrForSpread + m_histRates->ssprd10yrForSpread/100;
	
	double blend = m_histRates->lib1yrForSpread * percent1Y 
		+ libor10yrForSpread * (1 - percent1Y);

	return blend;
}

double PrimusHistRates::GetFN30YrServicing()
{
	return m_histRates->fn30yrServicing;
}

double PrimusHistRates::GetFN30YrNetCpn()
{
	return m_histRates->fn30yrNetCurrCpn;
}

double PrimusHistRates::GetCMTSpread(int term, double percent1Y)
{
	double blend = m_histRates->tsy1yrForSpread * percent1Y 
		+ m_histRates->tsy10yrForSpread * (1 - percent1Y);

	double parRate = GetParRate(term);

	double spread = parRate-blend;
	spread += (term == 30) ? m_histRates->offRunMtgSpread30YrAdj : m_histRates->offRunMtgSpread15YrAdj;
	spread *= 100;	// return value is in bp

	return spread;
}

double PrimusHistRates::GetParRate (int term)
{
	double parRate;

	if (term == 30) {
		parRate = m_histRates->fn30yrServicing + m_histRates->fn30yrNetCurrCpn;
	}
	else if (term == 15) {
		parRate = m_histRates->fn15yrServicing + m_histRates->fn15yrNetCurrCpn;
	}
	else 
		throw DRException ("Bad term ") << term;

	return parRate;
}


double PrimusHistRates::GetCMSSpread(int term, double percent1Y)
{
	double libor10yrForSpread = m_histRates->tsy10yrForSpread + m_histRates->ssprd10yrForSpread/100;
	
	double blend = m_histRates->lib1yrForSpread * percent1Y 
		+ libor10yrForSpread * (1 - percent1Y);

	double parRate = GetParRate(term);

	double spread = (parRate - blend);
	spread *= 100;	// return value is in bp

	return spread;
}

void PrimusHistRates::Get30YrRates (DRDate& date, DArray& rates, 
									DArray& points, double& latest30YrPts)
{
	DRDate tempDate (m_histRates->firstDateMonth, m_histRates->firstDateDay, 
		m_histRates->firstDateYear);

	date = tempDate;

	rates.resize(m_histRates->nRates);
	points.resize(m_histRates->nRates);

	for (int i = 0; i < m_histRates->nRates; i++) {
		rates[i] = m_histRates->mtgRates30yr[i] / 100.;
		points[i] = m_histRates->mtgRates30yrPts[i] / 100.;
	}
	latest30YrPts = m_histRates->latestMtgRate30yrPts / 100;
}
