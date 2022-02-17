//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditUtil.cpp
//
//   Description : Base class for SRM Credit util classes- calibration etc
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/CRCalib.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
SRMCreditUtil::SRMCreditUtil(
    const DateTime& baseDate,
    const string& smileParamsKey,
    ICDSParSpreadsConstSP stochCDSCurve) :
    baseDate(baseDate), 
    stochCDSCurve(stochCDSCurve),
    momentMatching(false),
    initialized(false)
{
}

/** destructor */
SRMCreditUtil::~SRMCreditUtil() {
}

/** returns log of the deterministic survival probability between two dates
    i.e. returns log [S(0,T) / S(0,t)] = - int_t^T lambda(0,u) du */
double SRMCreditUtil::logSurvProbRatio(const DateTime& loDate, 
                                       const DateTime& hiDate) const {
    // get access to DefaultRates
    DefaultRatesSP defRates(stochCDSCurve->defaultRates());

    // HACK -- assumption that baseDate can be used here ... is this really the case ... ?!?
    double temp1 = log(defRates->calcDefaultPV(baseDate, loDate));
    double temp2 = log(defRates->calcDefaultPV(baseDate, hiDate));
    return (temp2 - temp1);
}

/** returns log of the deterministic survival probability from today to a given date, 
    i.e. returns log S(0,t) = -int_0^t lambda(0,u) du */
double SRMCreditUtil::logSurvProbRatio(const DateTime& date) const {
    return logSurvProbRatio(baseDate, date);
}

/** computes and stores log of deterministic ratio of survival probabilities from first date in a 
    vector of given dates to each date (populating logFwdProbSimple and setting it to right length)
    result_i = log [S(0,T_i) / S(0,T_0)] = - int_{T_0}^{T_i} lambda(0,u) du */
void SRMCreditUtil::computeLogFwdProbSimple(const DateTimeArray& myDates,
                                            vector<double>& fwdLnProbRatios) const {
    fwdLnProbRatios.resize(myDates.size());
    if (myDates.empty())
        return;
    
    for (size_t i = 0; i < fwdLnProbRatios.size(); ++i) {
        // fwdLnProbRatios[i] = logSurvProbRatio(myDates[0], myDates[i]);
         fwdLnProbRatios[i] = logSurvProbRatio(baseDate, myDates[i]);
    }
}

/** computes and stores log of the deterministic ratio of consecutive survival probabilities
    on extended time line (populating logFwdProbSimple and setting it to right length)
    result_i = log [S(0,T_{i+1}) / S(0,T_{i})] = - int_{T_{i}}^{T_{i+1}} lambda(0,u) du */
void SRMCreditUtil::computeLogFwdProbSimple(vector<double>& logFwdProbSimple) const {
    logFwdProbSimple.resize(dates->size()-1);
    for (unsigned int i = 0; i < logFwdProbSimple.size(); i++) {
        logFwdProbSimple[i] = logSurvProbRatio((*dates)[i], (*dates)[i+1]);
    }
}

/** returns the CDS curve that is diffused */
ICDSParSpreadsConstSP SRMCreditUtil::getCdsCurve() const { 
    return stochCDSCurve; 
}

/** returns all the simulation dates excluding today */
const DateTimeArray& SRMCreditUtil::getSimDates() const { 
    assert(initialized); 
    return *dates; 
}

/** returns CDS dates merged with sim date */
const DateTimeArray& SRMCreditUtil::getExtendedTimeLine() const { 
    assert(initialized); 
    return extendedTimeLine;
}

/** returns today */
const DateTime& SRMCreditUtil::getBaseDate() const { 
    return baseDate; 
}

/** populates extendedTimeLine field 
    dates (date of simStart plus dates after simStart) + extra dates from input 
    along the lines of SRM3:ExtendFullTimeLine */
void SRMCreditUtil::calcExtendedTimeLine() {
    // get the CDS expiry dates 
    DefaultRatesSP defRates(stochCDSCurve->defaultRates());
    DateTimeArray cdsExpiries = defRates->getDates();
    // identify where we switch from dates to zeroDates
    int numCdsExpiriesToUse = dates->back().numFutureDates(cdsExpiries);
    extendedTimeLine.reserve(dates->size() + numCdsExpiriesToUse);
    extendedTimeLine = *dates;
    extendedTimeLine.insert(extendedTimeLine.end(),
                            cdsExpiries.end()-numCdsExpiriesToUse,
                            cdsExpiries.end());
}

DRLIB_END_NAMESPACE
