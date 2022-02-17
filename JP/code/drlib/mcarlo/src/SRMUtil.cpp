
//----------------------------------------------------------------------------
//
//   Group       : QR cross asset
//
//   Filename    : SRMUtil.hpp
//
//   Description : 
//
//   Author      : Henrik Rasmussen
//
//   Date        : 16 Dec 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/SRMConstants.hpp"

DRLIB_BEGIN_NAMESPACE

int SRMRound(double x)
{
	int y = (int) ( 0.5 * (1. + floor(2. * x)) );

	return (y);
}

double SRMExp(double exponent) 
{
    return (Maths::isPositive(exponent - SRMConstants::EXP_CUTOFF)) ? 
        exp(SRMConstants::EXP_CUTOFF) : exp(exponent); 
}

// XXX this is a fudge
// this problem needs to be addressed at the DateTime level
// basically if dateFrom.time == EOD and dateTo.time < SOD,
// e.g. dateTo.time == BEFORE_EX_DIV_TIME, then year frac between
// the two dates is negative !!!
double SRMYearFrac(const DateTime &dateFrom,
                   const DateTime &dateTo)
{
    if (dateFrom == dateTo)
	{
        return 0.0;
    }
    
	bool increasing = dateFrom < dateTo;
    
	double dt = dateFrom.yearFrac(dateTo);
    
	if (increasing && Maths::isNegative(dt))
	{
        throw ModelException("SRMYearFrac",
                             "Negative year fraction ("
                             + Format::toString(dt)
                             + ") between the two following increasing dates: "
                             + dateFrom.toString()
                             + ", "
                             + dateTo.toString());
    }

    if (!increasing && Maths::isPositive(dt))
	{
        throw ModelException("SRMYearFrac",
                             "Positive year fraction ("
                             + Format::toString(dt)
                             + ") between the two following decreasing dates: "
                             + dateFrom.toString()
                             + ", "
                             + dateTo.toString());
    }
    
	if (Maths::isZero(dt))
	{
        dt = 0.0;
    }
    
	return dt;
}

//////////////////////////////////////////////////////////////////////////////
//  From eqvol.c:SimpleFwdCurve
//
//   Given NbTP points, it outputs (NbTP - 1) "Forward rates * TimePeriod",
//     where FwdRate[i] spans (TPDate[i],TPDate[i+1]).
//
//  Note : this is a general utility, so we should remove versions from other
//         fiels as soon as convenient
//////////////////////////////////////////////////////////////////////////////

vector<double> simpleFwdCurve(
    IYieldCurveConstSP   yc,        // (I)
    const DateTimeArray& TPDate)    // (I) date of time point 
{
    static const string method("SimpleFwdCurve");
    
	if (TPDate.empty()) 
	{
        throw ModelException(method, "Zero dates supplied");
    }
    
	vector<double> FwdRate(TPDate.size()-1); // output
   
	// inverse of discount factor at start of period
	// assumes flat extrapolation 
    double Discount1 = yc->pv(TPDate[0]); 
    
	for (unsigned int i = 0; i < FwdRate.size(); i++) 
	{
        // inverse of discount factor at end of period   
		// assumes flat extrapolation
        double Discount2 = yc->pv(TPDate[i+1]);
        FwdRate[i] = Discount1 / Discount2 - 1.0;
        Discount1 = Discount2;
    }
    return FwdRate;
}// SimpleFwdCurve 


/*****  volExtend **************************************************/
/*
*      Extends a curve using the flat interp method and modifies
*      the date indexing of the extended curve. It copes with cases where
*      a period in the extended curve covers multiple periods in the source
*      curve.
*
*      On input:
*      ---------
*      SrcCurve[t] is applicable between InputDates[t-1] to
*      InputDates[t]. (when t=0, it will be between Today and InputDates[0])
*
*      Note: InputDates must be in strictly ascending order
*
*      On output:
*      ----------
*      ExtCurve[t] is the spot vol applicable between ExtendDates[t] and 
*      ExtendDates[t+1].
*
*      Note: 1) We assume ExtendDates[0] to be Today
*            2) Size of ExtSpotVol must be the same as that of ExtendDatess
*            3) We assume ExtSpotVol[NbExtendDatess-1] is never used because
*               it is meaningless.
*            4) ExtendDatess must be in strictly ascending order
*
*
*  From util_s::ExtendSpotVol 
*/
vector<double> SRMUtil::extendVol(
 	const DateTimeArray&  InputDates, // vol expiry dates
	const vector<double>& SrcCurve,
	const DateTimeArray&  ExtendDates) 
{
    static const string method("SRMRatesUtil::extendSpotVol");
    vector<double> ExtSpotVol(ExtendDates.size()-1);
    /* Loop through each period ExtendDates[i] to ExtendDates[i+1] */
    for (int i = 0; i < ExtendDates.size()-1; i++) {
        /* find smallest expiry date >= ExtendDates[i] */
        int LeftExpOfs = ExtendDates[i].findUpper(InputDates);
        if (LeftExpOfs == InputDates.size()) {
            /* if all expiries are < ExtendDates[i], then use last spot vol   */
            ExtSpotVol[i] = SrcCurve.back();
        } else if (ExtendDates[i+1] <= InputDates[LeftExpOfs]) {
            /* if this exp covers the whole ExtendDates prd, 
               then use its vol */
            ExtSpotVol[i] = SrcCurve[LeftExpOfs];
        } else {
            /* find largest expiry date <= ExtendDates[i+1] */
            int RightExpOfs = ExtendDates[i+1].findLower(InputDates);
            if (RightExpOfs < 0){
                throw ModelException(method, "Internal error");
            }
            /* integrate from ExtendDates[i] to InputDates[LeftExpOfs] */
            double var = Maths::square(SrcCurve[LeftExpOfs])*
                SRMYearFrac(ExtendDates[i], InputDates[LeftExpOfs]);
            
            /* integrate from InputDates[LeftExpOfs] to
               InputDates[RightExpOfs] */
            for (int j = LeftExpOfs + 1; j <= RightExpOfs; j++) {
                var += SrcCurve[j] * SrcCurve[j] *
                    SRMYearFrac(InputDates[j-1L], InputDates[j]);
            }
            
            /* integrate from InputDates[RightExpOfs] to ExtendDates[i+1] */
            if (RightExpOfs == (int)SrcCurve.size()-1){
                var += Maths::square(SrcCurve[RightExpOfs]) *
                    SRMYearFrac(InputDates[RightExpOfs], ExtendDates[i+1]);
            } else {
                var += Maths::square(SrcCurve[RightExpOfs+1L]) *
                    SRMYearFrac(InputDates[RightExpOfs], ExtendDates[i+1]);
            }
            
            double yearFraction = SRMYearFrac(ExtendDates[i], ExtendDates[i+1]);
            /* we only have to check zeroYearFraction, since strange cases are dealt
               with in SRMYearFrac */
            if (Maths::isZero(yearFraction)) {
                if (i>0) {
                    ExtSpotVol[i] = ExtSpotVol[i-1];
                } else {
                    ExtSpotVol[0] = 0.0; // not tested, more a theoretical case? 
                }
            } else {
                ExtSpotVol[i] = var / yearFraction;
            }

			// this can cause platform differences
            if (ExtSpotVol[i] < 0.) {
                throw ModelException(method, "Negative variance between "
                             + ExtendDates[i].toString()
                             + " and "
                             + ExtendDates[i+1].toString());
            }
            ExtSpotVol[i] = sqrt(ExtSpotVol[i]);
        }
    } /* for each i */
    return ExtSpotVol;
} /* ExtendSpotVol */

bool SRMUtil::SRMExplosion(double logIr) 
{
    return (Maths::isPositive(fabs(logIr) - SRMConstants::EXP_CUTOFF)) ? 
        true : false;
}

// TO DO - move this into a common utility file
//// the integral of exp(-b*t) between t and t+dt
double SRMUtil::GFAC(double t, double dt, double b)
{
    return (b<=SRMConstants::SRM_TINY ? 
            dt: ((1.0 - exp(-b*dt)) * exp(-b*t) /b));
}

double SRMUtil::GFACSimple(double t, double b)
{
    return b<=SRMConstants::SRM_TINY ? 
            t: (1.0 - exp(-b*t))/b;
}

// Return merged today plus all future dates inside allDates;
DateTimeArray SRMUtil::calcDiffusionDates(const DateTime& today, const DateTimeArray& allDates)
{
    return DateTime::merge(DateTimeArray(1, today), today.getFutureDates(allDates));
}
// Return the size of DateTimeArray returned by calcDiffusionDates()
int SRMUtil::getNumSimDates(const DateTime& today, const DateTimeArray& allDates)
{
    // return calcDiffusionDates(today, allDates).size(); // slow
    return 1 + today.numFutureDates(allDates); // O(log n)
}

/** populates sqrtYearFracs with sqrt of year fracs between sim dates using
    trading time */
vector<double> SRMUtil::computeSqrtYearFrac(const DateTimeArray& dates)
{
    vector<double> sqrtYearFracs(dates.size()-1);
    for (int i = 0; i+1 < dates.size(); ++i)
	{
        sqrtYearFracs[i] = sqrt(SRMYearFrac(dates[i], dates[i+1]));
    }
    return sqrtYearFracs;
}

	// year fractions, offset from today
vector<double> SRMUtil::computeYearFrac(const DateTime& today, const DateTimeArray& dates)
{
    vector<double> yearFracs(dates.size());
    for (int i = 0; i < dates.size(); ++i)
	{
        yearFracs[i] = SRMYearFrac(today, dates[i]);
    }
    return yearFracs;
}

/** From the list of dates requested add additional points as specified
    by timePointsPerYear */
DateTimeArraySP SRMUtil::fillInTimeLine(
        const DateTime&            start,      
        const DateTimeArray&       allDates, // dates in order
        double                     timePointsPerYear) 
{
    DateTimeArray inputDates = start.getNotPastDates(allDates);
    if (!inputDates.empty() && inputDates[0]!=start)
        inputDates.insert(inputDates.begin(), start); 
    // so we have the inputDates timeline starting and including "start" date

    // Populate timeline. Algorithm from timeline.c:GenerateTimeLine
    /* "extend the compact input date list according to the increasing */
    /* time-step extension algorithm"                                  */
    int Intval, NbStep, Step, Step0, Rest, Rest0, i, j;
    int MaxStep = (int)floor (372.0/timePointsPerYear);
    int LastStep = 0;

    // start with current date
    DateTimeArraySP simDates(new DateTimeArray(1, start));
    simDates->reserve(inputDates.size());
    for (i = 0; i < inputDates.size()-1; i++){
        /* Calculate the number of time points within each interval */
        Intval = inputDates[i+1].daysDiff(inputDates[i]);
        NbStep = 0;
        Step   = LastStep;
        do {
            if (Step < MaxStep) {
                Step++;
            }
            NbStep++;
            Intval -= Step;
        }while (Intval > 0);

        /* place the timepoints */
        Rest  = -Intval; // Days left over in current interval     
        Rest0 = Rest / NbStep;
        Rest -= Rest0 * NbStep;
        
        Step = LastStep;
        for (j = 0; j < NbStep; j++) {
            if (Step < MaxStep) {
                Step++;
            }
            if (j < NbStep-1) {
                /* Length of current time step after adjust   */
                Step0 = Step - Rest0 - (j >= NbStep - Rest);
                simDates->push_back(DateTime(simDates->back().rollDate(Step0)));
            } 
            else {
                // manually add the explicit date - as it might
                // have a different time of day to the date
                // created in above loop if we went to NbStep
                simDates->push_back(inputDates[i+1]);
            }
        }
        LastStep = Step;

    }/* for i */
    return simDates;
}



DRLIB_END_NAMESPACE
