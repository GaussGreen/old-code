//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GMinInvBenefit.cpp
//
//   Author      : Jan Cieslikiewicz / Francois Lu
//
//   Description   Guaranteed Minimum Investment Benefit payoff
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SampleList.hpp"
DRLIB_BEGIN_NAMESPACE

//////// instrument class //////////
class GMinInvBenefit: public Generic1Factor, 
                        virtual public LastSensDate,
                        virtual public IMCIntoProduct{
protected:
    /// fields ////////
    SampleListSP            spotSamples;    // all historic spots
	bool					priceFees;				// true=price fees, false=price structure
	int                     typeGMIB; 
	int                     typeGMDB;
	double					feeGMIB;
	double					initialWaitPeriod;
	double					initialAmount;
	int                     limitAgeBenefitIncrease; // age limit for MAV death benefit calculation
	double					adjustmentTreasurySettlement;			    
	DoubleArray             mortalityRate;      // mortality rate as num of years
    DoubleArray             surrenderBaseRate;      // mortality rate as num of years
    IntArray                ageAtStart;         // age group when contract start
	DoubleArray             groupSizeAtStart;   // percentage of people in the corresponding age group
	DoubleMatrix			CurrRates;			// Current Settlement Rates
	DoubleMatrix			GuarRates;			// Guaranteed Settlement Rates
	double					swap10;				// 10-year swap rate - not really sued in current implementation
	double					rollUpGMIB;			// percentage by which Roll Up GMIB increases each period
	double					rollUpGMDB;			// percentage by which Roll Up GMDB increases each period
	double					capGMIB;			// cap on Roll Up GMIB in relation to premium
	double					minRefCR;			// min rate used to reference Current Rate
	double					maxRefCR;			// max rate used to reference Current Rate
	int						treasuryRateUsed;	// treasury rate used to read tables
	int						surrenderPattern;	// surrender rate pattern used
	int						terminPattern;		// termination "pattern" uses - at this year people change behavior
	bool					terminationGMIB;	// true - allow termination of GMIB at a calculated rate											
												// constants for formulae
	double surrender1, surrender2, surrender3, surrender4, surrender5;
	double termin1, termin2, termin3, termin4, termin5, termin6, termin7;
	double annuit1, annuit2, annuit3;
	double with1, with2, with3, with4, with5;


public:
    static CClassConstSP const TYPE;
    friend class GMinInvBenefitProd;

    virtual void Validate(){
        static const string routine("GMinInvBenefit::Validate");
        // NB don't call parent's validate - issues with fwdStart
        if (fwdStarting){
            throw ModelException(routine, "Fwd starting flag on "
                                 "Generic1Factor is not used");
        }
        if (oneContract){
            throw ModelException(routine, "oneContract flag on "
                                 "Generic1Factor is not used");
        }

        if (ageAtStart[0] <= 0) {
            throw ModelException(routine, "ageAtStart and the array size must be >0.");
        }

        if (ageAtStart.size() != groupSizeAtStart.size()) {
            throw ModelException(routine, "ageAtStart and groupSizeAtStart arrays must be the same size.");
        }

		double totGroup = 0;
		int i;
		for (i = 0; i < groupSizeAtStart.size(); i++) totGroup += groupSizeAtStart[i];
		if (totGroup != 1) {
			throw ModelException(routine, "sizes of all age groups must add up to 100%");
		}

        if (spotSamples->getDates().size() == 0) {
            throw ModelException(routine, "there must be at least one sample date.");
        }

        if (spotSamples->getFirstDate() > valueDate) {
            throw ModelException(routine, "there must be one sample date in the past (or = valueDate).");
        }

		if (surrenderBaseRate.size() < spotSamples->getDates().size()) {
			throw ModelException(routine, "there must be a Surrender Base Rate for each sample date.");
		}

		if (mortalityRate[mortalityRate.size() - 1] != 1) {
			throw ModelException(routine, "table with Mortality Rates must end with an age for which the rate is 100%");
		}

		if (typeGMIB < 0 || typeGMIB > 2 || typeGMDB < 0 || typeGMDB > 2)  {
			throw ModelException(routine, "GMIB or GMDB type incorrect.");
		}

		if (minRefCR > maxRefCR) {
			throw ModelException(routine, "minRefCR cannot be greater than maxRefCR.");
		}
 
		if (initialAmount <= 0 || limitAgeBenefitIncrease <= 0 || capGMIB <= 0 || minRefCR <= 0 || 
			maxRefCR <= 0 || treasuryRateUsed <= 0) {
			throw ModelException(routine,
				"one of InitialAmount, BenefitIncreaseAgeLimit, GMIB cap, minRefCR, or maxRefCR less or equal to 0");
		}

		if (initialWaitPeriod < 0 || surrenderPattern < 0) {
			throw ModelException(routine, "Waiting Period or Surrender Pattern less than 0");
		}

		if (terminPattern < 2) {
			throw ModelException(routine, "Termination Rate Pattern must be at least 2");
		}

		for (i = 2; i < CurrRates.numRows(); i++) {
			if (CurrRates[0][i] != CurrRates[0][i - 1] + 1) {
				throw ModelException(routine, "missing age in Current Settlement Rates table");
			}
		}

		for (i = 2; i < GuarRates.numRows(); i++) {
			if (GuarRates[0][i] != GuarRates[0][i - 1] + 1) {
				throw ModelException(routine, "missing age in Guaranteed Settlement Rates table");
			}
		}
		
		for (i = 2; i < GuarRates.numCols(); i++) {
			if (GuarRates[i][0] != GuarRates[i - 1][0] + 1) {
				throw ModelException(routine, "missing calendar year in Guaranteed Settlement Rates table");
			}
		}

        AssetUtil::assetCrossValidate(asset.get(),
                                      false, //fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity*) const{
        DateTime end = valueDate.rollDate(365*(mortalityRate.size() - ageAtStart[0]));
        return end;
    }

private:
    GMinInvBenefit(): Generic1Factor(TYPE) {}

    // for reflection
    GMinInvBenefit(const GMinInvBenefit& rhs); // not implemented
    GMinInvBenefit& operator=(const GMinInvBenefit& rhs); // not implemented

    static IObject* defaultGMinInvBenefit(){
        return new GMinInvBenefit();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta){
        // use valueDate before it changes
        spotSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
                         asset.get()); // then roll our past values
        Generic1Factor::sensShift(theta); // and then call parent's method
        return true; // continue to tweak components which implement Theta
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GMinInvBenefit, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultGMinInvBenefit);
        FIELD(spotSamples, "historical values");
		FIELD(priceFees,    "priceFees");
		FIELD(typeGMIB,    "typeGMIB");
		FIELD(typeGMDB,    "typeGMDB");
		FIELD(feeGMIB,    "feeGMIB");
		FIELD(initialWaitPeriod,    "initialWaitPeriod");
		FIELD(initialAmount,    "initialAmount");
		FIELD(limitAgeBenefitIncrease,    "limitAgeBenefitIncrease");
		FIELD(adjustmentTreasurySettlement,    "adjustmentTreasurySettlement");
		FIELD(mortalityRate,    "mortalityRate");
		FIELD(surrenderBaseRate,    "surrenderBaseRate");
		FIELD(ageAtStart,    "ageAtStart");
		FIELD(groupSizeAtStart,    "groupSizeAtStart");
		FIELD(CurrRates,    "CurrRates");
		FIELD(GuarRates,    "GuarRates");
		FIELD(swap10,    "swap10");
		FIELD(rollUpGMIB,    "rollUpGMIB");
		FIELD(rollUpGMDB,    "rollUpGMDB");
		FIELD(capGMIB,    "capGMIB");
		FIELD(minRefCR,    "minRefCR");
		FIELD(maxRefCR,    "maxRefCR");
		FIELD(treasuryRateUsed,    "treasuryRateUsed");
		FIELD(surrenderPattern,    "surrenderPattern");
		FIELD(terminPattern,    "terminPattern");
		FIELD(surrender1,    "surrender1");
		FIELD(surrender2,    "surrender2");
		FIELD(surrender3,    "surrender3");
		FIELD(surrender4,    "surrender4");
		FIELD(surrender5,    "surrender5");
		FIELD(termin1,    "termin1");
		FIELD(termin2,    "termin2");
		FIELD(termin3,    "termin3");
		FIELD(termin4,    "termin4");
		FIELD(termin5,    "termin5");
		FIELD(termin6,    "termin6");
		FIELD(termin7,    "termin7");
		FIELD(annuit1,    "annuit1");
		FIELD(annuit2,    "annuit2");
		FIELD(annuit3,    "annuit3");
		FIELD(with1,    "with1");
		FIELD(with2,    "with2");
		FIELD(with3,    "with3");
		FIELD(with4,    "with4");
		FIELD(with5,    "with5");
		FIELD(terminationGMIB,    "terminationGMIB");

        clazz->setPublic(); // make visible to EAS/spreadsheet

    }
};

/* MC product class for super rainbow */
//////// product class //////////
class GMinInvBenefitProd : public IMCProduct, virtual public IMCProductLN{
private:
    const GMinInvBenefit*   inst; // reference to original instrument
    DateTimeArray             simulationDates;

    int             yrIndex;    // just needed for year until initialWaitPeriod years
    DoubleArray     histSpots;
    DateTimeArray   histDates;
    DoubleArray     growthFactors; // growth factor fom payout date to last date


public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    GMinInvBenefitProd(const GMinInvBenefit*         inst,
                IRefLevelSP refLevel,
                const SimSeriesSP&             simSeries) :
                IMCProduct(inst->asset.get(),
                          inst->valueDate,
                          inst->discount.get(),
                          refLevel,
                          simSeries,
                          inst->spotSamples,
                          inst->instSettle.get(),
                          simSeries->getLastDate()),
                          inst(inst){
                    
                    simulationDates = simSeries->getAllDates();

                    histDates = inst->spotSamples->getDates();
                    for (yrIndex=1; yrIndex < histDates.size(); yrIndex++)
                    {
                        if (histDates[yrIndex] > inst->valueDate)
                            break;
                    }
		           	histSpots = inst->spotSamples->getValues();

                    // compute df for payment dates
                    growthFactors.resize(histDates.size());
					DateTime lastDate = histDates[histDates.size() -1]; // this is the last date
                    for (int i=0; i<histDates.size(); i++)
                    {
						if(inst->valueDate <= histDates[i])
							growthFactors[i] = 1/inst->discount->pv(histDates[i], lastDate);
						else
							growthFactors[i] = 99999999; // this number should never be used, or a bug
                    }
                }

    /** print extra output **/
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const;
    /** vol interp for LN */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;


    void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

    double calcDeathBenefit(double AV, double MAVDB, int index);

	double PriceOption(int startAge, const double *path);
	double UpdateAV(int year, double AV, const double *path);
	double UpdateGMIB(int year, int typeGMIB, double AV, int startAge, double premium, vector <double> &GMIBs);
	double UpdateGMDB(int year, int typeGMDB, double AV, double GMDB, int startAge, double premium);
	double CalcGMIBITM(double GMIB, double AV, int year, int startAge);
	double CalcAnnuitRate(int year, double itmGMIB);
	double CalcSurrRate(int year, double itmGMIB);
	double CalcTerminRate(int year, double itmGMIB);
	double CalcMortalRate(int year, int startAge);
	double CalcPayoff(int year, int startAge, double GMIB, double AV);
	double CalcWithdrawRate(double itmGMIB);
	double WithdrawGMDB(double withdrawRate, double GMDB);
	double WithdrawGMIB(double withdrawRate, double AV, int typeGMIB, vector <double> &GMIBs);
	double WithdrawPremium(double withdrawRate, double premium);
	double WithdrawAV(double withdrawRate, double AV);
	double UpdatePayoff(int year, double payoffYear, double payoffTotal);
	double UpdateFee(int year, double feeYear, double feeTotal);
	double GetGR(int year, int startAge);
	double GetCR(int year, int startAge);
	double GrowthFactor(DateTime grownDate, DateTime grownToDate);

};
// end of class GMinInvBenefitProd




double GMinInvBenefitProd::calcDeathBenefit(double AV, double MAVDB,int index)
{
	return 0;
}




double GMinInvBenefitProd::PriceOption(int startAge, const double *path)
{
	// prices the option for a particular age group

	int year;				// periods after first value date
	double AV, GMIB, GMDB;	// account value, guaranteed minimum investment benefit, guaranteed minimum death benefit
	double premium;			// initial investment, but adjusted pro rate for withdrawals
	vector <double> GMIBs;	// ponter to two types of GMIB benefit - needed becasue witdrawal updates them differently
	double itmGMIB;			// percentage in the money of GMIB with respect to AV
	double annuitRate, surrRate, terminRate, mortalRate, withdrawRate;	// different rates
	double payoffYear, feeYear, payoffTotal, feeTotal;	// payoff and fee from the option for particular year and whole time
	double proportionLeft;	// proportion of people left at current time, equal to 1 at time 0	
		
		
	// initialize variables
	AV = GMIB = GMDB = premium = inst->initialAmount;
	GMIBs.resize(2);
	GMIBs[0] = inst->initialAmount;
	GMIBs[1] = inst->initialAmount;
	proportionLeft = 1;
	payoffTotal = feeTotal = 0;


	// get payoff for each year
	for (year = 1; year < histDates.size(); year++) {

		// update things with new value of the spot
		AV = UpdateAV(year, AV, path);
		GMIB = UpdateGMIB(year, inst->typeGMIB, AV, startAge, premium, GMIBs);
		GMDB = UpdateGMDB(year, inst->typeGMDB, AV, GMDB, startAge, premium);
		itmGMIB = CalcGMIBITM(GMIB, AV, year, startAge);	//change GR, CR functions

		// calculate how many people annuitized
		annuitRate = CalcAnnuitRate(year, itmGMIB);  // change waiting period

		// calculate payoff and fees for the year
		payoffYear = annuitRate * proportionLeft * CalcPayoff(year, startAge, GMIB, AV);
		feeYear = inst->feeGMIB * GMIB * proportionLeft;
		
		// add payoffs and fees from this year to total values for this simulation
		payoffTotal = UpdatePayoff(year, payoffYear, payoffTotal);
		feeTotal = UpdateFee(year, feeYear, feeTotal);

		// calculate how many people drop the policy or die
		surrRate = CalcSurrRate(year, itmGMIB);
		terminRate = CalcTerminRate(year, itmGMIB);
		mortalRate = CalcMortalRate(year, startAge);

		// calculate how much has been withdrawed and reduce things according to withdrawal rules
		withdrawRate = CalcWithdrawRate(itmGMIB);
		GMDB = WithdrawGMDB(withdrawRate, GMDB);
		GMIB = WithdrawGMIB(withdrawRate, AV, inst->typeGMIB, GMIBs);
		premium = WithdrawPremium(withdrawRate, premium);
		AV = WithdrawAV(withdrawRate, AV);
		
		// recaluclate how many people left for next year
		proportionLeft = (1 - annuitRate) * (1 - surrRate) * (1 - terminRate) * (1 - mortalRate) * proportionLeft;		
		if (proportionLeft == 0) break;  // no one left, end, no more payoffs to come
	}

	if (inst->priceFees) return feeTotal; //this needs to be changed when calculating both fee and payoff at once
	return payoffTotal; 	
}


double GMinInvBenefitProd::UpdateAV(int year, double AV, const double *path) {
	// calculates new value of account value - depending on the return of the asset
	double newAV = AV * path[year] / path[year - 1];
	return newAV;
}

double GMinInvBenefitProd::UpdateGMIB(int year, int typeGMIB, double AV, int startAge, double premium, vector <double> &GMIBs) {
	
	// calculates new GMIB, using new account value and type of GMIB
	// types: 0 - MAV, 1 - 5% RollUp, 2 - max(MAX, 5% RollUp)
	// cap on GMIB 5% type equal 2 times the premium
	// GMIB doesn't increase when the person reaches 80 (limitAgeBenefitIncrease) years
	
	
	int yearDiff = (int) histDates[0].yearFrac(histDates[year]);
	
	if (startAge + yearDiff <= inst->limitAgeBenefitIncrease) {
		GMIBs[0] = Maths::max(AV, GMIBs[0]);
		GMIBs[1] = GMIBs[1] * (1 + inst->rollUpGMIB);
	}
	if (GMIBs[1] > inst->capGMIB * premium) {
		GMIBs[1] = inst->capGMIB * premium;
	}

	double GMIB;
	if (typeGMIB == 0) {
		GMIB = GMIBs[0];
	} else if (typeGMIB == 1) {
		GMIB = GMIBs[1];
	} else {
		GMIB = Maths::max(GMIBs[0], GMIBs[1]);
	}
	return GMIB;
	
}

double GMinInvBenefitProd::UpdateGMDB(int year, int typeGMDB, double AV, double GMDB, int startAge, double premium) {
	
	// calculates new GMDB, using new account value and type of GMDB
	// types: 0 - MAV, 1 - 5% RollUp, 2 - max(MAX, 5% RollUp)
	// GMDB doesn't increase when the person reaches 80 years

	double newGMDB;

	int yearDiff = (int) histDates[0].yearFrac(histDates[year]);
	if (startAge + yearDiff > inst->limitAgeBenefitIncrease) {
		return GMDB;
	}

	if (typeGMDB == 0) {
		newGMDB = Maths::max(AV, GMDB);
	} else if (typeGMDB == 1) {
		newGMDB = premium * ::pow((1 + inst->rollUpGMDB), year);
	} else {
		newGMDB = Maths::max(Maths::max(AV, GMDB), premium * ::pow((1 + inst->rollUpGMDB), year));
	}

	return newGMDB;

}

double GMinInvBenefitProd::CalcGMIBITM(double GMIB, double AV, int year, int startAge) {
	
	// calculates percentage in the money of GMIB
	double guaranteedRate = GetGR(year, startAge);
	double currentRate = GetCR(year, startAge);
	double itm = (GMIB * guaranteedRate) / (AV * currentRate);
	return itm;
}


double GMinInvBenefitProd::CalcAnnuitRate(int year, double itmGMIB) {
	
	// calculates fraction of people that convert to annuity at any point in time
	int yearDiff = (int) histDates[0].yearFrac(histDates[year]);
	double rate;
	if (yearDiff <= inst->initialWaitPeriod) {
		rate = inst->annuit1;
	} else {
		rate = Maths::min(inst->annuit2, Maths::max(inst->annuit3, itmGMIB - 1));
	}
	return rate;
}


double GMinInvBenefitProd::CalcSurrRate(int year, double itmGMIB) {
	
	// calculates fraction of people that give up their benefits at any time. Surrender pattern has to be given in periods matching simulation dates
	double rate;
	if (year > inst->surrenderPattern) {
		rate = inst->surrenderBaseRate[year] * Maths::max(inst->surrender1, inst->surrender2 - inst->surrender3 * 
												Maths::max(itmGMIB - inst->surrender4, inst->surrender5));
	} else {
		rate = inst->surrenderBaseRate[year];

	}
	return rate;
}

double GMinInvBenefitProd::CalcTerminRate(int year, double itmGMIB) {

	// calculates fraction of people that terminate their policy
	double rate = 0;

	if (!inst->terminationGMIB) return rate;
	
	if (itmGMIB < 1) {
		
		if (year == 1) rate = inst->termin1;
		if (year == inst->terminPattern) rate = Maths::min(inst->termin2, 1 - itmGMIB);
		if (year > inst->terminPattern) rate = Maths::min(inst->termin3, 1 - itmGMIB);

	} else {

		if (year == 1) rate = inst->termin4;
		if (year == inst->terminPattern) rate = Maths::max(inst->termin5, inst->termin6 - inst->termin7 * (itmGMIB - 1));
	}
	return rate;
}


double GMinInvBenefitProd::CalcMortalRate(int year, int startAge) {
	// obtains the right mortality fate from the table. (mortality from begining of previous period used)
	double yearDiff = histDates[0].yearFrac(histDates[year]);
	double yearDiffPrev = histDates[0].yearFrac(histDates[year - 1]);
	int prevAge = startAge + (int) yearDiffPrev;
	double rateBeginPeriod = inst->mortalityRate[prevAge];
	double ratePeriod = (yearDiff - yearDiffPrev < 1 ? 
						rateBeginPeriod * (yearDiff - yearDiffPrev): 1 - pow(1 - rateBeginPeriod, yearDiff - yearDiffPrev));
	return ratePeriod;

}

double GMinInvBenefitProd::CalcPayoff(int year, int startAge, double GMIB, double AV) {
	// calculate payoff from the benefit
	double guaranteedRate = GetGR(year, startAge);
	double currentRate = GetCR(year, startAge);
	return Maths::max(GMIB * guaranteedRate/1000 - AV * currentRate/1000, 0.0) * 1000 / currentRate;
	
}

double GMinInvBenefitProd::CalcWithdrawRate(double itmGMIB) {
	// calculate partial withdrawal rate
	double rate = inst->with1 + (inst->with2 - inst->with3) * Maths::min(inst->with4, Maths::max(inst->with5, itmGMIB - 1));
	return rate;
}

double GMinInvBenefitProd::WithdrawGMDB(double withdrawRate, double GMDB) {
	// calculate new level of GMDB after withdrawal
	return GMDB * (1 - withdrawRate);
}

double GMinInvBenefitProd::WithdrawGMIB(double withdrawRate, double AV, int typeGMIB, vector <double> &GMIBs) {
	// calculate new level of GMIB after withdrawal
	double GMIB;
	
	GMIBs[0] = (1 - withdrawRate) * GMIBs[0];
	GMIBs[1] = Maths::max(0.0, GMIBs[1] - withdrawRate * AV); // GMIB for 5% RollUp reduced dollar for dollar, but can't be negative

	if (typeGMIB == 0) {
		GMIB = GMIBs[0]; 
	} else if (typeGMIB == 1) {
		GMIB = GMIBs[1];
	} else {
		GMIB = Maths::max(GMIBs[0], GMIBs[1]);
	}
	return GMIB;
}

double GMinInvBenefitProd::WithdrawPremium(double withdrawRate, double premium) {
	// calculate new level of premium after withdrawal
	return premium * (1 - withdrawRate);
}

double GMinInvBenefitProd::WithdrawAV(double withdrawRate, double AV) {
	// calculate new level of AV after withdrawal
	return AV * (1 - withdrawRate);
}

double GMinInvBenefitProd::UpdatePayoff(int year, double payoffYear, double payoffTotal) {
	// adds payoff for this year to the total payoffs, forwards it to the last simulation date
	return payoffTotal + payoffYear * growthFactors[year];
}

double GMinInvBenefitProd::UpdateFee(int year, double feeYear, double feeTotal) {
	// adds fee for this year to the total fees, forwards it to the last simulation date
	return feeTotal + feeYear * growthFactors[year];
}

double GMinInvBenefitProd::GetGR(int year, int startAge) {
	
	// calculates guaranteed settlement rate for the year
	int startYear = histDates[0].getYear();
	int yearDiff = (int) histDates[0].yearFrac(histDates[year]);
	int currentCalendarYear = startYear + yearDiff;
	int currentAge = startAge + yearDiff;
	
	int i;
	int row = 0;
	int col = 0;
	//find right row and column in data matrix
	for (i = 1; i < inst->GuarRates.numCols(); i++) {
		if (currentCalendarYear == inst->GuarRates[i][0]) col = i;
	}
	for (i = 1; i < inst->GuarRates.numRows(); i++) {
		if (currentAge == inst->GuarRates[0][i]) row = i;
	}
	if (col == 0) col = inst->GuarRates.numCols() - 1;
	if (row == 0) row = inst->GuarRates.numRows() - 1;

	double rate = inst->GuarRates[col][row];
	return rate;
}

double GMinInvBenefitProd::GetCR(int year, int startAge) {
	
	// calculates current settlement rate for the year
	int yearDiff = (int) histDates[0].yearFrac(histDates[year]);
	int currentAge = startAge + yearDiff;
	
	double currentContinous10Rate = log(GrowthFactor(histDates[year], histDates[year].rollDate(inst->treasuryRateUsed * 365)))
																							/ inst->treasuryRateUsed;
	double currentSwap10 = ::pow(exp(currentContinous10Rate), 1) - 1;
	double referenceRate = currentSwap10 - inst->adjustmentTreasurySettlement;
	//double currentContinous10Rate = log(growthFactors[year] / growthFactors[year + 10]) / 10;
	//double startContinous10Rate = log(growthFactors[0] / growthFactors[10]) / 10;
	//double startDiffRates = startContinous10Rate - inst->swap10;
	//double currentSwap10 = currentContinous10Rate - startDiffRates;
	//currentSwap10 = currentSwap10 - startDiffRates;
	

	// "referenceRate" has to be between 2% and 8% for this particular policy
	referenceRate = Maths::max(inst->minRefCR, referenceRate);
	referenceRate = Maths::min(inst->maxRefCR, referenceRate);

	// find right row in data matrix
	int row = 0;
	int i;
	for (i = 1; i < inst->CurrRates.numRows(); i++) {
		if (currentAge == inst->CurrRates[0][i]) row = i;
	}
	if (row == 0) row = inst->CurrRates.numRows() - 1;
	
	//find right column in data matrix and return right rate
	double rate;
	int col = 0;
	// find whether in between any interest rates in the table
	for (i = 1; i < inst->GuarRates.numCols() - 1; i++) {
		if (referenceRate >= inst->CurrRates[i][0] && referenceRate < inst->CurrRates[i + 1][0]) {
			col = i;
		}
	}
	if (referenceRate <= 0) { // rate came out negative, just take value on left side of table
		col = 1;
		rate = inst->CurrRates[col][row];
	} else if (col == 0) { // rate highre than highest rate, take value on right side of table
		col = inst->CurrRates.numCols() - 1;
		rate = inst->CurrRates[col][row];
	} else { // rate in between 2 values, look for the closest and take that one
		rate = inst->CurrRates[col][row] + (inst->CurrRates[col+ 1][row] - inst->CurrRates[col][row]) /
			(inst->CurrRates[col + 1][0] - inst->CurrRates[col][0]) * (referenceRate - inst->CurrRates[col][0]);
	}
	return rate;
}

double GMinInvBenefitProd::GrowthFactor(DateTime grownDate, DateTime grownToDate) {
	// calculates growth factor for particular year - 0 rate between year and last simulation date
	double growthFactor = 1/inst->discount->pv(grownDate, grownToDate);
	return growthFactor;
}


void GMinInvBenefitProd::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {

    const double* path = pathGen->Path(0,0); // access path
	int startAge;
	double value;  // payoff for a particular age
	double pathValue; // payoff for all ages combined, weighted by size of the age group


	if (pathGen->doingPast())
	{
		// compute realised ITM, Deaths, past lapses
	}
	else
	{
	 // compute payoff for each age group
		pathValue = 0;
		
		//double x = GrowthFactor(histDates[0]);
		
		for (int ageGroup =0; ageGroup < inst->ageAtStart.size(); ageGroup++)
		{
			startAge = inst->ageAtStart[ageGroup];
			value = PriceOption(startAge, path);
			pathValue += value * inst->groupSizeAtStart[ageGroup];
		}
	}

	if (!pathGen->doingPast())
	{
		prices.add(pathValue);
	
	}
}







// for the LogNormal path generator
CVolRequestLNArray GMinInvBenefitProd::getVolInterp(const IMCPathGenerator* pathGen,
                                int                     iAsset) const {
    DateTime imntStartDate = inst->fwdStarting? 
                             inst->startDate: inst->valueDate;
    CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

    reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
        inst->asset->getSpot(),  //ATM
        imntStartDate, 
        inst->endDate(0),
        inst->fwdStarting));
    
    return reqarr;
}


// control variate is done here
void GMinInvBenefitProd::recordExtraOutput(Control*, Results*, const IMCPrices& ) const
{
	/*
	OutputNameConstSP valueout(new OutputName("payment"));
    results->storeGreek(IObjectSP(CDouble::create(paymentStorage)), Results::DEBUG_PACKET, valueout);
	paymentStorage=0;
	Seed=-10;
	*/
}


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* GMinInvBenefit::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(spotSamples->getAllDates());

    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(spotSamples->getDates()[0]));
    return new GMinInvBenefitProd(this, refLevel, simSeries);
}

CClassConstSP const GMinInvBenefit::TYPE = CClass::registerClassLoadMethod(
    "GMinInvBenefit", typeid(GMinInvBenefit), GMinInvBenefit::load);

// force linker to include this file (avoid having header file) */
extern bool GMinInvBenefitLoad()
{
    return true && GMinInvBenefit::TYPE;
}

DRLIB_END_NAMESPACE
