//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : Computes the conditional survival probability at a given 
//               timepoint. Essentially multiplies the "inner names" 
//               conditional survival probabilities - And then integrates
//               them across market factors, to compute the fee leg price.
//               This class is very FtD specific.
//
// Date        : Oct 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/ConditionalFtDLossIntegrand.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"


DRLIB_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
//                               Static method
////////////////////////////////////////////////////////////////////////////////
DefaultRatesSP getDefaultRatesFromSurvivalProb(const DoubleArray& survivalProb,
                                               const DateTimeArray& timeLine,
                                               const DateTime& valueDate)
{
    // Default rates are continuous forwards; survival probs are discount factors.
    // Need Act365F to convert one to the other (JLHP - double check this)
    DayCountConventionConstSP act365f(new Actual365F());
    // Compute the default rates
    const int nbDates = timeLine.size();
    DoubleArray spotDefaultRate(nbDates, 1.0);
    for (int t=0; t < nbDates; ++t) {
        if (!Maths::isZero(survivalProb[t])) {
            spotDefaultRate[t] = 
                RateConversion::discountToRate(survivalProb[t],
                                               valueDate, 
                                               timeLine[t], 
                                               act365f.get(), 
                                               CompoundBasis::CONTINUOUS);
        }
    }
    
    // convert from spots to forwards
    DoubleArraySP fwdDefaultRates = 
        RateConversion::spotsToForwards(valueDate, 
                                        timeLine,
                                        spotDefaultRate, 
                                        act365f.get());

    // build the  default rates object
    DefaultRatesSP defaultRates(new CDSHelper::CParSpreadDefaultRates(
        valueDate,
        timeLine, 
        *fwdDefaultRates));

    return defaultRates;
}


////////////////////////////////////////////////////////////////////////////////


// destructor
ConditionalFtDLossIntegrand::~ConditionalFtDLossIntegrand() 
{}

// constructor
ConditionalFtDLossIntegrand::ConditionalFtDLossIntegrand(
        const vector<ICondLossDistributionsGenKeyArraySP> keysByDate,
        DateTimeArraySP timeLine,
        DateTimeArrayConstSP integrationDates,
        CDoubleArrayConstSP namesLoss,
        int mfDim,
        const DateTime& valueDate,
        const double recoveryDelay,
        IDiscountCurveSP discount) :
    FunctionNDDouble(mfDim, *InfiniteRange::createInfiniteRangeArray(mfDim)),
    integrationDates(integrationDates),
    keysByDate(keysByDate),
    timeLine(timeLine),
    namesLoss(namesLoss),
    valueDate(valueDate),
    recoveryDelay(recoveryDelay),
    discount(discount)
{}



// function
double ConditionalFtDLossIntegrand::operator() (
    const CDoubleArray& marketFactor) const
{
    static const string method("ConditionalFtDLossIntegrand::operator()");

    int nbIntegrationDates = integrationDates->size();
    //check if the deal has matured
    if (nbIntegrationDates <= 1) { // typically 1, the valuationDate
        return 0.0;
    }

    if (!timeLine) {
        throw ModelException(method, "Timeline is null");
    }
    if (!timeLine) {
        throw ModelException(method, "timeLine is null");
    }
    const int nbDates = keysByDate.size();
    if (nbDates == 0) {
        throw ModelException(method, "No keysByDate provided");
    }
    const int nbTimeLineDates = timeLine->size();
    if (nbTimeLineDates == 0) {
        throw ModelException(method, "No time line dates provided");
    }

    const int timeLineOffset = ((*timeLine)[0] == valueDate ? 1 : 0);
    if (nbDates + timeLineOffset != nbTimeLineDates) {
        throw ModelException(method, "The number of elements in keysByDate (" +
                             Format::toString(nbDates) + ") is not consistent "
                             "with the number of elements in the timeline (" +
                             Format::toString(nbTimeLineDates) + ").");
    }    


    DateTimeArray reducedTimeLine(timeLine->begin() + timeLineOffset, 
                                  timeLine->end());

    IMarketFactorValueSP mfValue(new DoubleArrayMarketFactor(marketFactor));

    // Get all names' survival probabilities
    const int nbNames = keysByDate[0]->size();
    DoubleArray jointSurvivalProb(nbDates, 1.0);
    DoubleArrayArray nameSurvProbByName(nbNames, DoubleArray(nbDates));
    for (int i=0; i < nbNames; ++i) {
        for (int t=0; t < nbDates; ++t) {
            // Obtain the names and joint survival probabilities for this time and name
            nameSurvProbByName[i][t] = 
                (*(keysByDate[t]))[i]->conditionalSurvProb(mfValue);
            jointSurvivalProb[t] *= nameSurvProbByName[i][t];
        }
    }

    // Compute the joint default rates
    DefaultRatesSP jointDefaultRates = 
        getDefaultRatesFromSurvivalProb(jointSurvivalProb, 
                                        reducedTimeLine, 
                                        valueDate);

    IDecretionCurveConstSP prepaySP(new NullDecretionCurve("NullPrepayCurve"));

    // Now, integrate for each name
    double pv = 0.0;
    for (int i=0; i < nbNames; ++i) {
        // Compute the name default rates
        DefaultRatesSP nameDefaultRates = 
            getDefaultRatesFromSurvivalProb(nameSurvProbByName[i], 
                                            reducedTimeLine, 
                                            valueDate);
        

        // Build the integrator
        CDSPricer::FtDIntegrator integrator(valueDate,
                                            integrationDates,
                                            discount,
                                            nameDefaultRates,
                                            jointDefaultRates,
                                            prepaySP,
                                            CDSPricer::DELAYED_DISCOUNTING,
                                            recoveryDelay);

        // And integrate
        if (!(integrator.empty())) {
            double namesIntegratedPV(0.0);
            while (integrator.nextStep()) {
                double integralValue = integrator.integralOfFtDPVTimesProbDensity();
                double namesPV = (*namesLoss)[i] * integralValue;
                namesIntegratedPV += namesPV;
            }
            pv += namesIntegratedPV;
        }
    }
    return pv;
}

DRLIB_END_NAMESPACE
