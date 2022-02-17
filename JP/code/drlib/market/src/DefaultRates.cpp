//----------------------------------------------------------------------------
//
//   Group       : EDR
//
//   Filename    : DefaultRates.cpp
//
//   Description : Essentially a clean spread curve with explicit dates
//                 Moved here from CDSHelper
//
//   Author      : André Segger
//
//   Date        : 24 August 2004
//
//

#include "edginc/config.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"


DRLIB_BEGIN_NAMESPACE

DefaultRates::~DefaultRates(){}

/** calculates the probability of a default due to a jump-to-zero event.
    Default implementation here returns 1 */
double DefaultRates::calcJumpPV(
    const DateTime&      fromDate,
    const DateTime&      toDate) const
{
    return 1.0;
}

/** A bit undocumented - seems to combine calcDefaultPV and calcJumpPV */
double DefaultRates::calcTotalDefaultPV(
    const DateTime&     fromDate,
    const DateTime&     toDate) const
{
    double defaultPV = calcDefaultPV(fromDate, toDate);
    double jumpPV    = calcJumpPV(fromDate, toDate);
    double koProb    = (1. - defaultPV) + defaultPV * (1. - jumpPV);
    return 1. - koProb;
    //return (drCalcDefaultPV(fromDate, toDate));//1. - koProb;
}
                                               
/** Calculates an implied par spread for a spot or forward CDS */
double DefaultRates::cdsParSpread(
    const DateTime&                effDate,
    int                            frequency,
    double                         notional,
    double                         recovery,
    bool                           accrueFee,
    const DateTime&                accruedEffectiveDate,
    const DateTime&                maturity,
    YieldCurveConstSP              discount,
    bool                           isE2C,
    const DayCountConvention*      dcc) const  // should this be in DefaultRates
{
    return CDSHelper::CDSParSpread(getValueDate(),
                                   effDate,
                                   frequency,
                                   notional,
                                   recovery,
                                   accrueFee,
                                   accruedEffectiveDate,
                                   maturity,
                                   discount,
                                   isE2C,
                                   this,
                                   dcc);
}

/** Calculate spot rates
 * (implicitly from forward rates inside DefaultRates).
 * It would probably be preferable to explicitly have a
 * getFwdRates() and a getSpotRates() method (or a more 
 * generic getRates(RateType type) methods). */
CashFlowArraySP DefaultRates::convertToSpot(
    const DateTime& valueDate,
    const DayCountConvention& dcc) const
{
    static const string method = "DefaultRates::convertForwardToSpot";
    try {
        DateTime prevDate;
        double yearFrac;
        double yearToBM;
        double nonFwdRate;
        double annNonFwdRate;
        double discountFactor;
        const DateTimeArray defaultDates = getDates();
        const CDoubleArray defaultRates = getRates();
        int nbDates = defaultRates.size();

        CashFlowArraySP annDefaultRates(new CashFlowArray(nbDates));
        
        for (int i = 0; i < nbDates; i++) {
            (*annDefaultRates)[i].date = defaultDates[i];
            if (i == 0) { // the first default rate is already a 'non-fwd' rate
                yearToBM = dcc.years(valueDate, defaultDates[i]);
                nonFwdRate = defaultRates[i];
            }
            else {
                nonFwdRate = 0.0;
                prevDate = valueDate;
                
                // Rq: this could be optimised !
                for (int j=0; j<i+1; j++) {
                    yearFrac = dcc.years(prevDate, defaultDates[j]);
                    nonFwdRate += defaultRates[j] * yearFrac;
                    prevDate = defaultDates[j];
                }
                yearToBM = dcc.years(valueDate, defaultDates[i]);
                nonFwdRate = nonFwdRate / yearToBM;
            }
            
            // convert to an annualised rate
            discountFactor = RateConversion::rateToDiscountYearFrac(
                nonFwdRate,
                yearToBM,
                CompoundBasis::CONTINUOUS);
            
            annNonFwdRate = RateConversion::discountToRateYearFrac(
                discountFactor,
                yearToBM,
                CompoundBasis::ANNUAL);
                
            (*annDefaultRates)[i].amount = annNonFwdRate;
        }
        return annDefaultRates;
    }
    catch (exception &e) {
        throw ModelException(&e, method);
    }    
}
class DefaultRates::DefaultLogOfPVKey: public virtual IKey{
public:
    DefaultLogOfPVKey(const DefaultRates* defRates):defRates(defRates){}
    /** Returns the log of the default PV between the two dates */
    virtual double calc(const DateTime&  loDate,
                        const DateTime&  hiDate){
        return (log(defRates->calcDefaultPV(loDate, hiDate)));
    }
private:
    const DefaultRates* defRates;
};

/** Returns a key used to optimise repeated calculations of
    default probabilities. The calc method for this key returns
    the natural logarithm of 1 - the default probability between
    the two dates (or equivalently the product of the forward
    default rate (continuous, Act/365F) and the negative year
    fraction (Act/365F) betweeen the two dates. Or another way of
    saying this is to say that it returns minus the integral of
    the forward default rate. The default implementation has no
    performance improvements. */
DefaultRates::IKey* DefaultRates::logOfDefaultPVKey() const{
    return new DefaultLogOfPVKey(this);
}

DRLIB_END_NAMESPACE
