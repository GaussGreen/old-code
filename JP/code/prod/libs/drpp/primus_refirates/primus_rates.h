#ifndef __primus_rates_h
#define __primus_rates_h

#include "ppmregresshistorical.h"
#include "drdate.h"
#include "drvalarray.h"
#include "drptr.h"

const char DefaultPrimusHistCCFile [] = "$HIST_REFIRATES_DIR/ppm_regress_historical_cc_rates.dat";
const char DefaultPrimusHistRatesFile [] = "$HIST_REFIRATES_DIR/ppm_regress_historical_rates.dat";

class PrimusHistRates {
public:
	PrimusHistRates (DRDate, const char* ccFile = DefaultPrimusHistCCFile, 
		const char* rateFile = DefaultPrimusHistRatesFile);

	~PrimusHistRates ();
	DArray Get30YrMtgRates (DRDate date, int numRates);
	DArray Get30YrPoints (DRDate date, int numRates);

	double GetCMTSpread(int term = 30, double percent1Y = .215);
	double GetCMSSpread(int term = 30, double percent1Y = .215);
	double GetParRate (int term = 30);

	double GetCMSBlend (double percent1Y = .215);

	double GetOffRun30YrAdjustment () {return m_histRates->offRunMtgSpread30YrAdj;}
	// Warning: the OffRun30YrAdjustment assumes the primus
	// ratios of 21.5/78.5 1/10yr blend

	void Get30YrRates (DRDate& date, DArray& rates, DArray& points, double& latest30YrPts);

	DRDate GetLastHistRateDate();
	double GetMtgRate(DRDate, int term = 30);
	double GetPoint(DRDate, int term = 30);
	double GetLatestMtgPts(int term = 30); 
	double GetFN30YrServicing();
	double GetFN30YrNetCpn();

private:
	void CalculateIndicies (DRDate date, int num, int& startIndex, int& numHistRates);
	int CalculateOffset (DRDate&);

	PpmRegressHistoricalRates* m_histRates;
	DRDate m_date;

	PrimusHistRates(const PrimusHistRates&);	// make these private, 
	PrimusHistRates& operator=(const PrimusHistRates&);
};

typedef DRPtr<PrimusHistRates> PrimusHistRatesPtr;
#endif
