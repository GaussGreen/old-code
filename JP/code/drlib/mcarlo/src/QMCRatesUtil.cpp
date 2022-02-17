//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QMCRatesUtil.cpp
//
//   Description : Helper for SRM - used for holding intermediate data
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/QMCRatesUtil.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include <cassert>

DRLIB_BEGIN_NAMESPACE


///// constructs and populates SRMRatesUtil object
QMCRatesUtil::QMCRatesUtil(
            const DateTime&      baseDate,
            IYieldCurveConstSP   discYC,
            IYieldCurveConstSP   diffYC): // eg 6M
        baseDate(baseDate), 
        discYC(discYC), 
        diffYC(diffYC), 
//        processedVol(_processedVol),
//        momentMatching(false),
        initialized(false)
{
}


/** Calculates a deterministic discount factor from date to today. Uses
    diffused curve if useDiffusedCurve is true else uses discount curve */
double QMCRatesUtil::pv(const DateTime& date, bool useDiffusedCurve) const
{
    return (useDiffusedCurve? diffYC: discYC)->pv(date);
}

/** Calculates a deterministic discount factor between the 2 dates. Uses
    diffused curve if useDiffusedCurve is true else uses discount curve */
double QMCRatesUtil::pv(
    const DateTime& loDate, 
    const DateTime& hiDate,
    bool useDiffusedCurve) const
{
    return (useDiffusedCurve? diffYC: discYC)->pv(loDate, hiDate);
}
    
//// computes log of DF [determinstic] from discount zero curve between dates
//// on extended time line and
//// populates logDiscFactor (and sets it to right length)
void QMCRatesUtil::computeLogDiscFactor(vector<double>& logDiscFactor) const
{
    logDiscFactor.resize((*dates).size()-1);
    auto_ptr<IYieldCurve::IKey> key(discYC->logOfDiscFactorKey());
    for (size_t i = 0; i < logDiscFactor.size(); i++){
        logDiscFactor[i] = key->calc((*dates)[i], (*dates)[i+1]);
    }
}

void QMCRatesUtil::computeLogDiscFactor(const DateTimeArray& myDates,
                                     vector<double> &  logDiscFactor,
                                     IYieldCurveConstSP yc) const
{
    logDiscFactor.resize(myDates.size());
    auto_ptr<IYieldCurve::IKey> key(yc->logOfDiscFactorKey()); // The state variables use
    for (size_t i = 0; i < logDiscFactor.size(); i++){
            logDiscFactor[i] = key->calc(baseDate, myDates[i]);
    }
}

//// populates extendedTimeLine field. Requires dates. 
//// From SRM3:ExtendFullTimeLine. I think I'm assuming that we only need to
//// deal with one zero curve here
// FIXME: EDFs may use different YC, shouldn't we collect all future zeroDates() ?
void QMCRatesUtil::calcExtendedTimeLine()
{
    // get the zero curve dates
    const DateTimeArray& zeroDates = diffYC->zeroDates();
    // identify where we switch from dates to zeroDates
    int numZeroDatesToUse = (*dates).back().numFutureDates(zeroDates);
    extendedTimeLine.reserve((*dates).size() + numZeroDatesToUse);
    extendedTimeLine = (*dates);
    extendedTimeLine.insert(extendedTimeLine.end(),
                            zeroDates.end()-numZeroDatesToUse,
                            zeroDates.end());
}

void QMCRatesUtil::assertInitialized(void) const
{
	assert(initialized);
}
DRLIB_END_NAMESPACE
