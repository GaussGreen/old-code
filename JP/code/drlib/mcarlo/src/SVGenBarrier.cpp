//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenBarrier.cpp
//
//   Description : Barrier state variable
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenBarrier.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/MCPathConfigSRMGenSV.hpp"


DRLIB_BEGIN_NAMESPACE

const double EXTREME = 1.0e+99;

BarrierPerAssetData::BarrierPerAssetData(const string& monitorType,
                                         ScheduleSP levels,
                                         ScheduleSP economicLevels,
                                         DateTimeArraySP monitoringDates,
                                         const DateTime& smoothDate,
                                         const DateTime& valueDate,
                                         bool isOut,
                                         bool isUp,
                                         bool isHit):
levels(levels), economicLevels(economicLevels), monitoringDates(monitoringDates),
isOut(isOut), isUp(isUp), isHit(isHit), barrierLevels(monitoringDates->size(), isUp ? EXTREME : -EXTREME) {
    static const string routine = "BarrierPerAssetData::BarrierPerAssetData";

    try {
        if(valueDate < monitoringDates->front() && isHit) {
            throw ModelException(
                "Input isHit flag is set to TRUE but first monitoring date " +
                monitoringDates->front().toString() +
                " is in the future");
        }

        DateTimeArray barDates = levels->getDates();
        // monitoringDates is not always barrierDates, especially for INTERP_NONE case.
        int numSteps = monitoringDates->size();
        if (levels->getInterp() == Schedule::INTERP_NONE){
            bool dummy;
            barMap = DateTime::createMapping(*monitoringDates, barDates, dummy);
        }
        else
            barMap.resize(numSteps+1, 0); // all days are sampling date

        // Create the interpolated barriers for the asset
        bool hasEconomic = economicLevels.get() && !(economicLevels->getDates().empty());
        int iStep = 0;
        for (iStep = barMap[iStep]; iStep < numSteps; iStep++, iStep+=barMap[iStep]) {
            // exclude the period which doesn't have barrier.
            if (barDates[0] <= (*monitoringDates)[iStep] && (*monitoringDates)[iStep] <= barDates.back()){
                if (hasEconomic && (*monitoringDates)[iStep]<=valueDate)    // use economic barrier for past.
                    barrierLevels[iStep] = economicLevels->interpolate((*monitoringDates)[iStep]);
                else
                    barrierLevels[iStep] = levels->interpolate((*monitoringDates)[iStep]);
            }
        }

        barrierSmoothDate = levels->interpolate(smoothDate);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


DateTimeArray BarrierPerAssetData::getFirstFutureDailyMonDate(const CAsset& asset,
                                                              const DateTime& valueDate,
                                                              const string& monitorType) const {
    static const string routine = "BarrierPerAssetData::getFirstFutureMonDate";

    try {
        DateTimeArray result;

        if(CString::equalsIgnoreCase(monitorType, BarrierData::DAILY_MONITORING)) {
            const DateTime& firstAssetMonDate = monitoringDates->front();
            const DateTime& lastAssetMonDate = monitoringDates->back();
            if(firstAssetMonDate < valueDate && valueDate <= lastAssetMonDate) {
                DateTime firstFutureMonDate;
                int barrierTime = firstAssetMonDate.getTime();
                if(valueDate.getTime() != barrierTime) {
                    // Find the date x such that
                    // 1) x > valueDate
                    // 2) x.getTime() = valueDate.getTime()
                    DateTime test(valueDate.getDate(), barrierTime);
#if 0
                    // ADD BUSINESS DATE
                    HolidayConstSP hols = AssetUtil::getHoliday(&asset);
                    if(test > valueDate) {
                        if(hols->isHoliday(test)) {
                            // That is the case if valueDate is a holiday
                            firstFutureMonDate = hols->addBusinessDays(test, 1);
                        } else {
                            // Correct date
                            firstFutureMonDate = test;
                        }
                    } else {
                        // Add a business date to test
                        firstFutureMonDate = hols->addBusinessDays(test, 1);
                    }
#else
                    // ADD CALENDAR DATE
                    if(test > valueDate) {
                        firstFutureMonDate = test;
                    } else {
                        firstFutureMonDate = test.rollDate(1);
                    }
#endif
                }

                // Check that this is a valid barrier date
                if(firstFutureMonDate > lastAssetMonDate) {
                    // We can't do that so ignore it
                    firstFutureMonDate = DateTime();
                }

                if(firstFutureMonDate > valueDate) {
                    result.push_back(firstFutureMonDate);
                }
            }
        }

        return result;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


//////////////////////////////////////////////////////////////////////


SmoothingData::SmoothingData(bool smoothing, double lowSpread, double highSpread,
                             const DateTime& smoothDate):
smoothing(smoothing), lowSpread(lowSpread), highSpread(highSpread),
smoothDate(smoothDate) {}


//////////////////////////////////////////////////////////////////////


BarrierData::BarrierData(const string& monitorType,
                         const string& closedFormDates,
                         const string& closedFormMethod,
                         const string& maxNumDivsPerYearString,
                         const vector<BarrierPerAssetDataSP>& barrierPerAsset,
                         int numHits,
                         SmoothingDataSP smoothingData,
                         const DateTime& hitDate,
                         const DateTime& valueDate):
monitorType(monitorType), closedFormDates(closedFormDates),
closedFormMethod(closedFormMethod), barrierPerAsset(barrierPerAsset),
numHits(numHits), smoothingData(smoothingData), hitDate(hitDate),
valueDate(valueDate) {
    static const string routine = "BarrierData::BarrierData";

    try {

        numAssets = barrierPerAsset.size();

        // Find if the barrier is in the past
        DateTime lastMonDate;
        int iAsset;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTime& lastAssetMonDate = barrierPerAsset[iAsset]->
                monitoringDates->back();
            lastMonDate = lastAssetMonDate < lastMonDate ? lastMonDate : lastAssetMonDate;
        }
        isBarrierInPast = lastMonDate <= valueDate ? true : false;

        amendAllAssetsHit();

        if(!european()) {
            // Make sure that all assets have more than 1 monitoring dates
            for(iAsset = 0; iAsset < numAssets; iAsset++) {
                if(barrierPerAsset[iAsset]->monitoringDates->size() < 2) {
                    throw ModelException(routine,
                        "Need at least 2 monitoring dates for monitor type " + monitorType +
                        ". Check asset number " +
                        Format::toString(iAsset + 1));
                }
            }
        }

        // Figure out number of divs per year
        if (CString::equalsIgnoreCase(maxNumDivsPerYearString, USE_ALL_DIVS)) {
            // All dividends
            maxNumDivsPerYear = -1;
        } else if(CString::equalsIgnoreCase(maxNumDivsPerYearString, NO_DIVS)) {
            // No dividends
            maxNumDivsPerYear = 0;
        } else {
            // Figure out the number of dividends
            char* endPtr;
            maxNumDivsPerYear = (int)strtol(maxNumDivsPerYearString.c_str(), &endPtr, 10 /* base 10 */);
            if (*endPtr != '\0') {
                throw ModelException(routine,
                    "Failed to convert string "+ string(maxNumDivsPerYearString)+ " to a int");
            }
            if(maxNumDivsPerYear < 0) {
                throw ModelException(routine,
                    "Maximum number of divs must be non-negative: " + string(maxNumDivsPerYearString));
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


BarrierPerAssetDataSP BarrierData::getAssetBarrier(int iAsset) const {
    return barrierPerAsset[iAsset];
}


SmoothingDataSP BarrierData::getSmoothingData() const {
    return smoothingData;
}


bool BarrierData::european() const {
    return CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING);
}


void BarrierData::amendIsHit(bool isAssetHit, int iAsset) {
    static const string routine = "BarrierData::amendIsHit";

    BarrierPerAssetDataSP assetBarrier = barrierPerAsset[iAsset];
    bool isAssetHitInput = assetBarrier->isHit;
    if(!isAssetHitInput && isAssetHit) {
        assetBarrier->isHit = isAssetHit;
        amendAllAssetsHit();
    } else {
        throw ModelException(routine,
            "Internal error: asset " + Format::toString(iAsset+1) +
            " is already hit. Forbidden to amend isHit flag.");
    }
}


void BarrierData::amendAllAssetsHit() {
    // Find if all assets hit in past for CF
    bool allAssetsHitTmp;
    if(european()) {
        allAssetsHitTmp = false;   // no meaning in case of european
    } else {
        allAssetsHitTmp = true;
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            allAssetsHitTmp = allAssetsHitTmp && barrierPerAsset[iAsset]->isHit;
        }
    }

    allAssetsHit = allAssetsHitTmp;
}



// Monitoring strings
const string BarrierData::EUROPEAN_MONITORING   = "E";
const string BarrierData::DAILY_MONITORING      = "D";
const string BarrierData::CONTINUOUS_MONITORING = "C";

const string BarrierData::DEFAULT   = "DEFAULT";

// Dates specific
const string BarrierData::INPUT     = "INPUT";
const string BarrierData::BENCHMARK = "BENCHMARK";
const string BarrierData::ENDPOINTS = "ENDPOINTS";
const string BarrierData::AUTO      = "AUTO";

// Choice of methodology
const string BarrierData::BROWNIAN_BRIDGE    = "BB";
const string BarrierData::BARRIER_ADJUSTMENT = "ADJ";

// Dividend treatment for Brownian Bridges
const string BarrierData::DEFAULT_MAX_DIVS_PER_YEAR = "None";
const string BarrierData::USE_ALL_DIVS = "All";
const string BarrierData::NO_DIVS = "None";
const int BarrierData::MIN_DIVS_PER_PERIOD = 1;


void BarrierData::validateMethodology(const string& monitorType, bool useStateVars,
                                      string& method, string& dates) {
    static const string routine = "BarrierData::validateMethodology";

    try {
        // Validate choice of dates
        if( !CString::equalsIgnoreCase(dates, INPUT) &&
            !CString::equalsIgnoreCase(dates, BENCHMARK) &&
            !CString::equalsIgnoreCase(dates, ENDPOINTS) &&
            !CString::equalsIgnoreCase(dates, AUTO) &&
            !CString::equalsIgnoreCase(dates, DEFAULT)) {
            throw ModelException("Unrecognized methodology for monitoring dates: " +
                                 dates +
                                 ". It has to be one of: " +
                                 INPUT + ", " + BENCHMARK + ", " + ENDPOINTS + ", " + AUTO + " or " + DEFAULT);
        }

        // Validate choice of method
        if( !CString::equalsIgnoreCase(method, BROWNIAN_BRIDGE) &&
            !CString::equalsIgnoreCase(method, BARRIER_ADJUSTMENT) &&
            !CString::equalsIgnoreCase(method, DEFAULT)) {
            throw ModelException("Unrecognized monitoring methodology: " +
                                 method +
                                 ". It has to be one of: " +
                                 BROWNIAN_BRIDGE + ", " +  BARRIER_ADJUSTMENT + " or " + DEFAULT);
        }

        // Validate choice of monitor type
        if(!CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING) &&
           !CString::equalsIgnoreCase(monitorType, DAILY_MONITORING) &&
           !CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) {
            throw ModelException("Unrecognized monitoring type: " +
                                 monitorType +
                                 ". It has to be one of: " +
                                 EUROPEAN_MONITORING + ", " +  DAILY_MONITORING + " or " + CONTINUOUS_MONITORING);
        }

        // Validate combination of flags
        if(CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
            // Only default is supported
            if(!CString::equalsIgnoreCase(dates, DEFAULT)) {
                throw ModelException("ClosedFormDates can only be " + DEFAULT +
                                     " for european monitoring. Here: " + dates);
            }
            if(!CString::equalsIgnoreCase(method, DEFAULT)) {
                throw ModelException("ClosedFormMethod can only be " + DEFAULT +
                                     " for european monitoring. Here: " + method);
            }
            return;
        }

        if(!useStateVars) {
            // Only default is supported
            if(!CString::equalsIgnoreCase(dates, DEFAULT)) {
                throw ModelException("ClosedFormDates can only be " + DEFAULT +
                                     " when not using state variables. Here: " + dates);
            }
            if(!CString::equalsIgnoreCase(method, DEFAULT)) {
                throw ModelException("ClosedFormMethod can only be " + DEFAULT +
                                     " when not using state variables. Here: " + method);
            }
            return;
        }

        // UseStateVars is true. Map DEFAULT values to actual methodology

        // Default method is Brownian Bridge
        if(CString::equalsIgnoreCase(method, DEFAULT)) {
            method = BROWNIAN_BRIDGE;
        }

        // Validate Barrier Adjustment
        if(CString::equalsIgnoreCase(method, BARRIER_ADJUSTMENT)) {
            // Default dates is INPUT
            if(CString::equalsIgnoreCase(dates, DEFAULT)) {
                dates = INPUT;
            }

            // Adjustment can only be INPUT
            if(!CString::equalsIgnoreCase(dates, INPUT)) {
                throw ModelException("ClosedFormDates can only be " + INPUT + " or " + DEFAULT +
                    " when using state variables with barrier adjustment. Here: " + dates);
            }
        }

        // Validate Barrier Adjustment
        if(CString::equalsIgnoreCase(method, BROWNIAN_BRIDGE)) {
            // Default dates is AUTO
            if(CString::equalsIgnoreCase(dates, DEFAULT)) {
                dates = AUTO;
            }
        }

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

//////////////////////////////////////////////////////////////////////


/** Basic implementation of IMCBarrierHVStateVar */
class SVGenBarrierHV::MCBarrierHVStateVarTrivial: virtual public IStateVar {
public:
    /** Constructor */
    MCBarrierHVStateVarTrivial(BarrierDataSP data,
                               const IMultiFactors* mAsset,
                               IStateVariableGen::IStateGen* pathGen,
                               IRefLevel::IStateVarGenSP refLevelGen,
                               SVGenSpotSP spotPathGen) {
        int numAssets = data->numAssets;
        hitNoHitValues.resize(numAssets);
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            BarrierPerAssetDataConstSP assetBarrier = data->getAssetBarrier(iAsset);
            bool isAssetHitInput = assetBarrier->isHit;
            bool isAssetHit = isAssetHitInput;
            bool isUp = assetBarrier->isUp;
            if(pathGen->doingPast() && !isAssetHit) {
                // Get refLevel state variable
                IRefLevel::IStateVarSP refLevelSV;
                refLevelSV = refLevelGen->getRefLevelSV(refLevelSV,pathGen);
                double refLevel = refLevelSV->refLevel(iAsset);
                ScheduleSP schedule = assetBarrier->levels;

                // Compare past values to barrier
                if(spotPathGen.get()) {
                    // Get state variables
                    SVGenSpot::IStateVarSP spotPathSV = spotPathGen->getSpotSV(pathGen);
                    const DateTimeArray& monDates = spotPathGen->getSimSeries()->getDates(iAsset);
                    const SVPath& path = spotPathSV->path(iAsset);

                    // See if asset hits
                    for(int iStep = path.begin(); !isAssetHit && iStep < path.end(); iStep++) {
                        double barrier = refLevel * schedule->interpolate(monDates[iStep]);
                        double spot = path[iStep];
                        if( ( isUp && spot >= barrier) ||
                            (!isUp && spot <= barrier) ) {
                            isAssetHit = true;
                        }
                    }
                }

                // Process current spot against barrier
                if(!isAssetHit) {
                    const DateTime& today = data->valueDate;
                    const DateTimeArray& monitoringDates = *assetBarrier->monitoringDates;
                    if( monitoringDates.front() <= today && today <= monitoringDates.back() ) {
                        // Check if current today is monitoring day
                        bool isDailyMonitoring = CString::equalsIgnoreCase(
                            data->monitorType, BarrierData::DAILY_MONITORING);
                        bool isContMonitoring = CString::equalsIgnoreCase(
                            data->monitorType, BarrierData::CONTINUOUS_MONITORING);

                        if( (isDailyMonitoring && today.getTime() ==  monitoringDates.front().getTime() ) ||
                             isContMonitoring) {
                            // Compare against current spot
                            double barrier = refLevel * schedule->interpolate(today);
                            double spot = mAsset->assetGetSpot(iAsset);
                            if( ( isUp && spot >= barrier) ||
                                (!isUp && spot <= barrier) ) {
                                isAssetHit = true;
                            }
                        }
                    }
                }

                // Amend the input flag
                if(isAssetHit) {
                    data->amendIsHit(isAssetHit, iAsset);
                }
            }

            // Use the amended flags for the state variable
            double binPayOnHit = assetBarrier->isOut ? 0.0 : 1.0;
            double binPayNoHit = 1.0 - binPayOnHit;
            double payoff;
            if (isAssetHit) {
                payoff = binPayOnHit;
            } else {
                payoff = binPayNoHit;
            }

            hitNoHitValues[iAsset] = payoff;
        }
    }

    /** Hit value per asset */
    virtual double hitValue(int iAsset) const {
        return hitNoHitValues[iAsset];
    }

    /** Indicates past Vs future */
    virtual bool doingPast() const {
        return true;
    }

    /** Records BARRIER_LEVEL output requests in results */
    virtual void recordBarrierLevels(CControl* control,
                                     Results*  results,
                                     const IMultiFactors* mAsset) const {
        // No reporting since this is for past only
    }

private:
    DoubleArray     hitNoHitValues;     //!< Hit not hit indicators
};


//////////////////////////////////////////////////////////////////////


double SVGenBarrierHV::hitValuePayoff(const SmoothingDataConstSP& smoothingData,
                                      const BarrierPerAssetDataConstSP& assetBarrier,
                                      const SVPath& spotSmoothPath,
                                      double refLevel,
                                      bool isAssetHit,
                                      bool isAssetHitInPast,
                                      bool doingPast) {
    double binPayOnHit = assetBarrier->isOut ? 0.0 : 1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    double payoff;
    /* Do not smooth if
           1) no smoothing is required OR
           2) asset has hit in past OR
           3) if we are processing the past. The smoothing date
              might be in the future and the memory has not been allocated. Note
              that the value will be ignored anyway if there is future. */
    if (!smoothingData->smoothing ||
        isAssetHitInPast ||
        (spotSmoothPath.begin() == spotSmoothPath.end())) {
        if (isAssetHit) {
            payoff = binPayOnHit;
        } else {
            payoff = binPayNoHit;
        }
    } else {
        double highSpread = smoothingData->highSpread;
        double lowSpread = smoothingData->lowSpread;
        double dSpread = highSpread - lowSpread;
        double OoI = assetBarrier->isOut? -1.0 : 1.0;
        double UoD = assetBarrier->isUp? 1.0 : 0.0;

        double spotAtSmooth = spotSmoothPath[spotSmoothPath.begin()];

        // Fwd
        double fwd = (spotAtSmooth / refLevel) - assetBarrier->barrierSmoothDate;

        // Put butterfly
        double pb = ( Maths::max(0., highSpread - fwd) +
                      Maths::max(0., -highSpread - fwd) -
                      2.0 * Maths::max(0., -fwd) ) / dSpread;

        // Call butterfly
        double cb = (Maths::max(0., fwd - lowSpread) +
                     Maths::max(0., fwd + lowSpread) -
                     2.0 * Maths::max(0., fwd) ) / dSpread;

        // Finally compute payoff
        double smoothFactor = UoD * (cb - pb);
        if (isAssetHit) {
            payoff = binPayOnHit + OoI * (smoothFactor - cb);
        } else {
            payoff = binPayNoHit + OoI * (smoothFactor + pb);
        }
    }
    return payoff;
}


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVEur::StateVar::StateVar(
    BarrierDataConstSP data,
    IRefLevel::IStateVarGenSP refLevelGen,
    SVGenSpotSP spotMonGen,
    SVGenSpotSP spotSmoothGen):
data(data), refLevelGen(refLevelGen), spotMonGen(spotMonGen), spotSmoothGen(spotSmoothGen),
hitNoHit(data->numAssets), hitNoHitInPast(data->numAssets),
hitNoHitSoFar(data->numAssets) {
    // Initialize hit flags to input data
    for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
        hitNoHit[iAsset]       = data->getAssetBarrier(iAsset)->isHit;
        hitNoHitInPast[iAsset] = hitNoHit[iAsset];
        hitNoHitSoFar[iAsset]  = hitNoHit[iAsset];
    }
}


void SVGenBarrierHVEur::StateVar::update(IStateVariableGen::IStateGen* pathGen) {
    static const string routine = "SVGenBarrierHVEur::StateVar::update";

    try{
        // Pull the state variables out of the PathGen:
        isPast = pathGen->doingPast();
        monitoringPath = spotMonGen->getSpotSV(pathGen);
        smoothPath = spotSmoothGen->getSpotSV(pathGen);
        refLevel = refLevelGen->getRefLevelSV(refLevel, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double SVGenBarrierHVEur::StateVar::hitValue(int iAsset) const {
    // Copy over from preservation to work
    bool isAssetHit = hitNoHitSoFar[iAsset];

    // Obtain the path and the barrier information for the Asset
    const SVPath& path = monitoringPath->path(iAsset);
    const BarrierPerAssetDataConstSP& assetBarrier = data->getAssetBarrier(iAsset);
    bool isUp = assetBarrier->isUp;
    const DoubleArray& barrierLevels = assetBarrier->barrierLevels;

    // Compare against barrier
    double ref = refLevel->refLevel(iAsset);
    for (int iStep = path.begin(); !isAssetHit && iStep < path.end(); iStep++) {
        if ( (isUp && path[iStep] >= barrierLevels[iStep] * ref) ||
            (!isUp && path[iStep] <= barrierLevels[iStep] * ref)) {
            isAssetHit = true;
        }
    }

    if(isPast) {
        // Preserve value for next call
        hitNoHitSoFar[iAsset]  = isAssetHit;
        hitNoHitInPast[iAsset] = isAssetHit;
    }

    // Compute payoff from hitting value i.e. apply smoothing
    double payoff = hitValuePayoff(data->getSmoothingData(),
                                   assetBarrier,
                                   smoothPath->path(iAsset),
                                   ref,
                                   isAssetHit,
                                   hitNoHitInPast[iAsset],
                                   isPast);
    return payoff;
}


bool SVGenBarrierHVEur::StateVar::doingPast() const {
    return isPast;
}


void SVGenBarrierHVEur::StateVar::recordBarrierLevels(CControl* control,
                                                      Results*  results,
                                                      const IMultiFactors* mAsset) const {

    static const string method = "SVGenBarrierHVEur::StateVar::recordBarrierLevels";

    try {
        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request) {
            for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
                BarrierPerAssetDataConstSP assetData = data->getAssetBarrier(iAsset);
                DateTimeArraySP monitoringDates = assetData->monitoringDates;
                const DateTime& valueDate = data->valueDate;

                // Only consider assets for which reference level has been determined
                // and monitoring has not been completed (this includes assets for which
                // the reference level is determined but monitoring has not started
                if(!hitNoHitInPast[iAsset] &&
                    refLevelGen->getDates(iAsset).back() <= valueDate &&
                    valueDate <= monitoringDates->back()) {

                    // Get the economic barrier if it exists
                    const Schedule* schedule = assetData->economicLevels.get() &&
                                                !(assetData->economicLevels->getDates().empty()) ?
                        assetData->economicLevels.get() : assetData->levels.get();

                    // Since we are in the middle of monitoring, the refLevel has been determined
                    double ref = refLevel->refLevel(iAsset);

                    // Construct reports for barrier from now to some upper date
                    DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
                    DoubleArray values(monitoringDates->size());
                    int iDate;
                    for(iDate = 0; iDate< values.size(); iDate++) {
                        values[iDate] = schedule->interpolate((*monitoringDates)[iDate]);
                    }
                    Schedule complete(*monitoringDates.get(),
                                      values,
                                      Schedule::INTERP_NONE);
                    CashFlowArraySP subset(complete.subset(valueDate, upperDate));
                    BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));
                    for (iDate = 0; iDate < subset->size(); iDate++) {
                        // needs to be absolute
                        BarrierLevel bl(assetData->isUp,
                                        (*subset)[iDate].date,(*subset)[iDate].amount * ref,
                                        false);
                        reportLevels->push_back(bl);
                    }

                    // Store them in results
                    if (!reportLevels->empty()) {
                        OutputRequestUtil::recordBarrierLevels(
                            control, results,
                            mAsset->assetGetTrueName(iAsset),
                            reportLevels.get());
                    }
                }
            }
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}



//////////////////////////////////////////////////////////////////////


SVGenBarrierHVEur::SVGenBarrierHVEur(BarrierDataSP data,
                                     IRefLevel::IStateVarGenSP refLevelGen):
data(data), refLevelGen(refLevelGen) {
    static const string routine = "SVGenBarrierHVEur::SVGenBarrierHVEur";

    try {
#if 0
        // Removed this in order to allow MC LV with inserted daily dates to operate
        // The only call seems to be conditioned safely anyway, so should be safe.
        // Sanity check
        if(!data->european()) {
            throw ModelException("Inconsistent european barrier flags.");
        }
#endif

        int numAssets = data->numAssets;

        // Validate that reference level is before all monitoring dates
        int iAsset;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTimeArraySP& assetDates = data->getAssetBarrier(iAsset)->
                monitoringDates;
            const DateTime& firstAssetMonDate = assetDates->front();
            const DateTime& lastAssetRefDate  = refLevelGen->getDates(iAsset).back();
            if(firstAssetMonDate< lastAssetRefDate) {
                throw ModelException("First monitoring date " +
                                     firstAssetMonDate.toString() +
                                     " must be after last reference level date " +
                                     lastAssetRefDate.toString() +
                                     " for asset " + Format::toString(iAsset+1));
            }
        }

        // Create an SVGenSpot for the smoothing date
        const DateTime& smoothDate = data->getSmoothingData()->smoothDate;
        if(smoothDate > data->valueDate) {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(numAssets, smoothDate));
        } else {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(SimSeriesSP(new SimSeries(numAssets))));
        }

        // Create an SVGenSpot for the monitoring dates
        // Stuff dates in a SimSeries object
        SimSeriesSP monSimSeries(new SimSeries(numAssets));
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            monSimSeries->addDates(iAsset, *(data->getAssetBarrier(iAsset))->monitoringDates);
        }
        spotMonGen = SVGenSpotSP(new SVGenSpot(monSimSeries));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void SVGenBarrierHVEur::collectStateVars(IStateVariableCollectorSP svCollector) const {
    static const string routine = "SVGenBarrierHVEur::collectStateVars";

    try {
        // Append smoothing and reflevel
        svCollector->append(spotSmoothGen.get());  // Smoothing date
        svCollector->append(refLevelGen.get());    // Reference level
        svCollector->append(spotMonGen.get());     // Discrete path
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


SVGenBarrierHVEur::IStateVarSP SVGenBarrierHVEur::getHitValueSV(
    IStateVarSP oldStateVar, IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVEur::getHitValueSV";

    try {
        if(oldStateVar.get()) {
            // Update european barrier with new state variables etc.
            StateVar* eurBarrier = dynamic_cast<StateVar*>(oldStateVar.get());
            if(!eurBarrier) {
                throw ModelException(
                    "Expected class of type SVGenBarrierHVEur::StateVar");
            }
            eurBarrier->update(pathGen);
            return oldStateVar;
        } else {
            // Create a new one
            StateVarSP barrierEuropean(new StateVar(
                data, refLevelGen, spotMonGen, spotSmoothGen));
            barrierEuropean->update(pathGen);
            return barrierEuropean;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IStateVariableSP SVGenBarrierHVEur::create(IStateVariableSP             oldStateVar,
                                        IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVEur::create";

    try {
        IStateVarSP oldHVStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
        return getHitValueSV(oldHVStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


BarrierDataConstSP SVGenBarrierHVEur::getBarrierData() const {
    return data;
}


IRefLevel::IStateVarGenSP SVGenBarrierHVEur::getRefLevel() const {
    return refLevelGen;
}


SVGenSpotSP SVGenBarrierHVEur::getSmoothPath() const {
    return spotSmoothGen;
}


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVBarAdj::SVGenBarrierHVBarAdj(BarrierDataSP data,
                                           IRefLevel::IStateVarGenSP refLevelGen,
                                           const IMultiFactors* mAsset):
data(data), refLevelGen(refLevelGen), mAsset(mAsset) {
    static const string routine = "SVGenBarrierHVBarAdj::SVGenBarrierHVBarAdj";

    try {
        // Sanity check
        if(data->european()) {
            throw ModelException("Inconsistent european barrier flags.");
        }

        int numAssets = data->numAssets;

        // Validate that reference level is before all monitoring dates
        int iAsset;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTimeArraySP& assetDates = data->getAssetBarrier(iAsset)->
                monitoringDates;
            const DateTime& firstAssetMonDate = assetDates->front();
            const DateTime& lastAssetRefDate  = refLevelGen->getDates(iAsset).back();
            if(firstAssetMonDate< lastAssetRefDate) {
                throw ModelException("First monitoring date " +
                                     firstAssetMonDate.toString() +
                                     " must be after last reference level date " +
                                     lastAssetRefDate.toString() +
                                     " for asset " + Format::toString(iAsset+1));
            }
        }

        // Create an SVGenSpot for the smoothing date
        const DateTime& smoothDate = data->getSmoothingData()->smoothDate;
        if(smoothDate > data->valueDate) {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(numAssets, smoothDate));
        } else {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(SimSeriesSP(new SimSeries(numAssets))));
        }

        // Create an SVGenSpot for the monitoring dates
        // Stuff dates in a SimSeries object
        SimSeriesSP monSimSeries(new SimSeries(numAssets));
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            monSimSeries->addDates(iAsset, *(data->getAssetBarrier(iAsset))->monitoringDates);
        }
        spotMonGen = SVGenSpotSP(new SVGenSpot(monSimSeries));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void SVGenBarrierHVBarAdj::collectStateVars(IStateVariableCollectorSP svCollector) const {
    static const string routine = "SVGenBarrierHVBarAdj::collectStateVars";

    try {
        // Append smoothing and reflevel
        svCollector->append(spotSmoothGen.get());      // Smoothing date
        svCollector->append(refLevelGen.get());        // Reference level
        svCollector->append(spotMonGen.get());         // Discrete path

        // CAREFUL: the statement "do not append(this) if data->allAssetsHit is true" is
        //          flawed. The reason is that allAssetsHit might be false at price and
        //          true at delta. This will result in an exception because the pathGen
        //          will configure itself at delta in a radically different way and the
        //          caches will complain that the dimensions are wrong. So we always
        //          append this except if barrier is in past (note that on theta tweaks
        //          the caches get regenerated so we are fine).
        // HOWEVER: one can still take a shortcut in getHitValueSV() and create a trivial
        //          state variable if allAssetsHit is true
        if(data->isBarrierInPast) {
            // Do not add elementary state var if barrrier is in the past
        } else {
            svCollector->appendElementary(this);       // Elementary state variable
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


SVGenBarrierHVBarAdj::IStateVarSP SVGenBarrierHVBarAdj::getHitValueSV(
    IStateVarSP oldStateVar, IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVBarAdj::getHitValueSV";

    try {
        if(pathGen->doingPast() || data->isBarrierInPast || data->allAssetsHit) {
            // Trivial monitoring in the following cases:
            //     1) past OR
            //     2) all barrier dates in past OR
            //     3) all assets hit in past
            // See comments in collectStateVars() method
            return IStateVarSP(new
                MCBarrierHVStateVarTrivial(data, mAsset, pathGen, refLevelGen, spotMonGen));
        } else {
            // Future implementation delegated to the path generator
            IStateVarSP barrierSV(&dynamic_cast<IStateVar&>(*pathGen->create(this)));
            return barrierSV;
        }

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IStateVariableSP SVGenBarrierHVBarAdj::create(IStateVariableSP             oldStateVar,
                                           IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVBarAdj::create";

    try {
        IStateVarSP oldHVStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
        return getHitValueSV(oldHVStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


BarrierDataConstSP SVGenBarrierHVBarAdj::getBarrierData() const {
    return data;
}


IRefLevel::IStateVarGenSP SVGenBarrierHVBarAdj::getRefLevel() const {
    return refLevelGen;
}


SVGenSpotSP SVGenBarrierHVBarAdj::getSpotPath() const {
    return spotMonGen;
}


SVGenSpotSP SVGenBarrierHVBarAdj::getSmoothPath() const {
    return spotSmoothGen;
}


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVBB::StateVar::StateVar(IStateVariableGen::IStateGen* pastPathGen,
                                     IStateVariableGen::IStateGen* futPathGen,
                                     BarrierDataConstSP            data,
                                     IRefLevel::IStateVarGenSP     refLevelGen,
                                     SVGenSpotSP                      spotSmoothGen,
                                     SVGenSpotSP                      spotAtMonStartGen,
                                     SVGenSpotSP                      spotAtMonEndGen,
                                     const vector<IntArray>&       offsets,
                                     const vector<IntArray>&       hitNoHitPath):
data(data), offsets(offsets), hitNoHitPath(hitNoHitPath),
refLevelGen(refLevelGen), spotSmoothGen(spotSmoothGen),
spotAtMonStartGen(spotAtMonStartGen), spotAtMonEndGen(spotAtMonEndGen) {

    static const string method = "SVGenBarrierHVBB::StateVar::StateVar";

    try {
        // Figure out if we need to do a European monitoring at start for "D" monitoring
        // and a European monitoring at the end for "S" type barriers
        int numAssets   = data->numAssets;
        monitorAtStart  = vector<int>(numAssets, 0);
        monitorAtEnd    = vector<int>(numAssets, 0);
        barrierPctStart = DoubleArray(numAssets, 0.0);
        barrierPctEnd   = DoubleArray(numAssets, 0.0);

        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            ScheduleSP schedule = data->getAssetBarrier(iAsset)->levels;

            if(CString::equalsIgnoreCase(schedule->getInterp(), Schedule::INTERP_STAIRS)) {
                const DateTimeArray& monEndDates = spotAtMonEndGen->getSimSeries()->getDates(iAsset);
                if(monEndDates.size()) {
                    monitorAtEnd[iAsset]  = 1;
                    barrierPctEnd[iAsset] = schedule->interpolate(monEndDates.front());
                }
            }

            const DateTimeArray& monStartDates = spotAtMonStartGen->getSimSeries()->getDates(iAsset);
            if(monStartDates.size()) {
                monitorAtStart[iAsset]  = 1;
                barrierPctStart[iAsset] = schedule->interpolate(monStartDates.front());
            }
        }

        // Get old state variable from past path gen and create new one
        IRefLevel::IStateVarSP oldStateVar = refLevelGen->getRefLevelSV(
            IRefLevel::IStateVarSP(   ), pastPathGen);
        refLevel = refLevelGen->getRefLevelSV(oldStateVar, futPathGen);

        // Get smooth path and endpoint state variables
        smoothPath         = spotSmoothGen->getSpotSV(futPathGen);
        spotAtMonStartPath = spotAtMonStartGen->getSpotSV(futPathGen);
        spotAtMonEndPath   = spotAtMonEndGen->getSpotSV(futPathGen);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


double SVGenBarrierHVBB::StateVar::hitValue(int iAsset) const {
    const BarrierPerAssetDataConstSP& assetBarrier = data->getAssetBarrier(iAsset);

    bool hitNoHitInPast = assetBarrier->isHit;
    bool isAssetHit = hitNoHitInPast;

    double assetRefLevel = refLevel->refLevel(iAsset);

    // Explicit monitoring at start when required
    if(monitorAtStart[iAsset] && !isAssetHit) {
        const SVPath& monStartPath = spotAtMonStartPath->path(iAsset);
        if(monStartPath.begin() != monStartPath.end()) {
            double barrierAtStart = barrierPctStart[iAsset] * assetRefLevel;
            double spotAtStart = monStartPath[monStartPath.begin()];
            bool isUp = assetBarrier->isUp;
            if ( (isUp && spotAtStart >= barrierAtStart) ||
                (!isUp && spotAtStart <= barrierAtStart)) {
                isAssetHit = true;
            }
        }
    }

    // Loop over BBs and figure out if it's a hit or no hit
    const IntArray& assetOffsets = offsets[iAsset];
    const int* assetHitNoHit     = &hitNoHitPath[iAsset][0];
    for(int iStep = 0; !isAssetHit && iStep < assetOffsets.size(); iStep++) {
        if (assetHitNoHit[assetOffsets[iStep]]) {
            isAssetHit = true;
        }
    }

    // Explicit monitoring at the end when required
    if(monitorAtEnd[iAsset] && !isAssetHit) {
        const SVPath& monEndPath = spotAtMonEndPath->path(iAsset);
        if(monEndPath.begin() != monEndPath.end()) {
            double barrierAtEnd = barrierPctEnd[iAsset] * assetRefLevel;
            double spotAtEnd = monEndPath[monEndPath.begin()];
            bool isUp = assetBarrier->isUp;
            if ( (isUp && spotAtEnd >= barrierAtEnd) ||
                (!isUp && spotAtEnd <= barrierAtEnd)) {
                isAssetHit = true;
            }
        }
    }

    const SVPath& spotSmoothPath = smoothPath->path(iAsset);
    // Compute payoff from hitting value i.e. apply smoothing
    double payoff = SVGenBarrierHV::hitValuePayoff(data->getSmoothingData(),
                                                   assetBarrier,
                                                   spotSmoothPath,
                                                   assetRefLevel,
                                                   isAssetHit,
                                                   hitNoHitInPast,
                                                   false);
    return payoff;
}


bool SVGenBarrierHVBB::StateVar::doingPast() const {
    return false;
}


void SVGenBarrierHVBB::StateVar::recordBarrierLevels(CControl* control,
                                                     Results*  results,
                                                     const IMultiFactors* mAsset) const {

    static const string method = "SVGenBarrierHVBB::StateVar::recordBarrierLevels";

    try {
        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request) {
            for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
                BarrierPerAssetDataConstSP assetData = data->getAssetBarrier(iAsset);
                DateTimeArraySP monitoringDates = assetData->monitoringDates;
                const DateTime& valueDate = data->valueDate;

                // Only consider live assets during their monitoring period
                if(!assetData->isHit &&
                    refLevelGen->getDates(iAsset).back() <= valueDate &&
                    valueDate <= monitoringDates->back()) {
                    // Get the economic barrier if it exists
                    const Schedule* schedule = assetData->economicLevels.get()  &&
                                                !(assetData->economicLevels->getDates().empty()) ?
                        assetData->economicLevels.get() : assetData->levels.get();

                    // Since we are in the middle of monitoring, the refLevel has been determined
                    double ref = refLevel->refLevel(iAsset);

                    // Construct reports for barrier from now to some upper date
                    DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
                    CashFlowArraySP subset(schedule->subset(valueDate, upperDate));
                    BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));
                    for (int iDate = 0; iDate < subset->size(); iDate++) {
                        BarrierLevel bl(assetData->isUp,
                                        (*subset)[iDate].date,
                                        (*subset)[iDate].amount * ref,
                                        CString::equalsIgnoreCase(data->monitorType, BarrierData::CONTINUOUS_MONITORING));
                        reportLevels->push_back(bl);
                    }

                    // Store them in results
                    if (!reportLevels->empty()) {
                        OutputRequestUtil::recordBarrierLevels(
                            control, results, mAsset->assetGetTrueName(iAsset), reportLevels.get());
                    }
                }
            }
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVBB::SVGenBarrierHVBB(BarrierDataSP data,
                                   IRefLevel::IStateVarGenSP refLevelGen,
                                   const IMultiFactors* mAsset):
data(data), refLevelGen(refLevelGen), mAsset(mAsset) {
    static const string routine = "SVGenBarrierHVBB::SVGenBarrierHVBB";

    try {
        // Sanity check
        if(data->european()) {
            throw ModelException("Inconsistent european barrier flags.");
        }

        int numAssets = data->numAssets;

        SimSeriesSP startMonDates(new SimSeries(numAssets));
        SimSeriesSP endMonDates(new SimSeries(numAssets));
        const DateTime& valueDate = data->valueDate;

        int iAsset;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTimeArray& assetDates = *data->getAssetBarrier(iAsset)->monitoringDates;
            const DateTime& firstAssetMonDate = assetDates.front();
            const DateTime& lastAssetRefDate  = refLevelGen->getDates(iAsset).back();
            // Validate that reference level is before all monitoring dates
            if(firstAssetMonDate < lastAssetRefDate) {
                throw ModelException("First monitoring date " +
                                     firstAssetMonDate.toString() +
                                     " must be after last reference level date " +
                                     lastAssetRefDate.toString() +
                                     " for asset " + Format::toString(iAsset+1));
            }

            // Add first future monitoring date
            DateTimeArray firstFutureMonDate = data->getAssetBarrier(iAsset)->
                getFirstFutureDailyMonDate(mAsset->getAsset(iAsset), valueDate, data->monitorType);
            startMonDates->addDates(iAsset, firstFutureMonDate);

            const DateTime& lastBarrierDate = assetDates.back();
            if(lastBarrierDate > valueDate) {
                endMonDates->addDates(iAsset, DateTimeArray(1, lastBarrierDate));
            } else {
                endMonDates->addDates(iAsset, DateTimeArray(0));
            }
        }

        // Create an SVGenSpot for the smoothing date
        const DateTime& smoothDate = data->getSmoothingData()->smoothDate;
        if(smoothDate > valueDate) {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(numAssets, smoothDate));
        } else {
            spotSmoothGen = SVGenSpotSP(new SVGenSpot(SimSeriesSP(new SimSeries(numAssets))));
        }

        // Create an SVGenSpot for the first and last monitoring dates
        spotAtMonStartGen = SVGenSpotSP(new SVGenSpot(startMonDates));
        spotAtMonEndGen = SVGenSpotSP(new SVGenSpot(endMonDates));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void SVGenBarrierHVBB::collectStateVars(IStateVariableCollectorSP svCollector) const {
    static const string routine = "SVGenBarrierHVBB::collectStateVars";

    try {
        // Append smoothing and reflevel
        svCollector->append(refLevelGen.get());        // Reference level
        svCollector->append(spotSmoothGen.get());      // Smoothing date
        svCollector->append(spotAtMonStartGen.get());  // Start monitoring date
        svCollector->append(spotAtMonEndGen.get());    // End monitoring date

        // CAREFUL: the statement "do not append(this) if data->allAssetsHit is true" is
        //          flawed. The reason is that allAssetsHit might be false at price and
        //          true at delta. This will result in an exception because the pathGen
        //          will configure itself at delta in a radically different way and the
        //          caches will complain that the dimensions are wrong. So we always
        //          append this except if barrier is in past (note that on theta tweaks
        //          the caches get regenerated so we are fine).
        // HOWEVER: one can still take a shortcut in getHitValueSV() and create a trivial
        //          state variable if allAssetsHit is true
        if(data->isBarrierInPast) {
            // Do not add elementary state var if barrrier is in the past
        } else {
            svCollector->appendElementary(this);       // Elementary state variable
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


SVGenBarrierHVBB::IStateVarSP SVGenBarrierHVBB::getHitValueSV(
    IStateVarSP oldStateVar, IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVBB::getHitValueSV";

    try {
        if(pathGen->doingPast() || data->isBarrierInPast || data->allAssetsHit) {
            // Trivial monitoring in the following cases:
            //     1) past OR
            //     2) all barrier dates in past OR
            //     3) all assets hit in past
            // See comments in collectStateVars() method
            return IStateVarSP(new
                MCBarrierHVStateVarTrivial(data, mAsset, pathGen, refLevelGen, SVGenSpotSP(   )));
        } else {
            // Future implementation delegated to the path generator
            IStateVarSP barrierSV(&dynamic_cast<IStateVar&>(*pathGen->create(this)));
            return barrierSV;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


IStateVariableSP SVGenBarrierHVBB::create(IStateVariableSP             oldStateVar,
                                       IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVBB::create";

    try {
        IStateVarSP oldHVStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
        return getHitValueSV(oldHVStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


BarrierDataConstSP SVGenBarrierHVBB::getBarrierData() const {
    return data;
}


IRefLevel::IStateVarGenSP  SVGenBarrierHVBB::getRefLevelGen() const {
    return refLevelGen;
}


SVGenSpotSP SVGenBarrierHVBB::getSpotSmoothGen() const {
    return spotSmoothGen;
}


SVGenSpotSP SVGenBarrierHVBB::getSpotAtMonStartGen() const {
    return spotAtMonStartGen;
}


SVGenSpotSP SVGenBarrierHVBB::getSpotAtMonEndGen() const {
    return spotAtMonEndGen;
}


//////////////////////////////////////////////////////////////////////


/** Constructs a driver generator for continuous barrier monitoring */
DateTimeArraySP SVGenBarrierHVBB::BarrierDatesHelper::getBarrierTimeline(
    const MCProductClient* prod,
    const vector<const SVGenBarrierHVBB*>& barrierGens,
    DateTimeArray& simDates) {

    static const string routine = "SVGenBarrierHVBB::BarrierDatesHelper::getBarrierTimeline";

    try {
        if(!barrierGens.size()) {
            return DateTimeArraySP(new DateTimeArray());
        }

        // Get the product cliquet start dates
        // Pass null path generator is it's not used in SVsbarrierHVBBGenArray);
        DateTimeArray cliquetStartDates = prod->getVolStartDates(0);
        const DateTime& today = prod->getToday();
        if(cliquetStartDates.size() > 1) {
            throw ModelException("Cannot support cliquets yet.");
        }
        DateTime refDate = cliquetStartDates.front();
        if(refDate < today) {
            refDate = today;
        }

        int numAssets = barrierGens[0]->getBarrierData()->numAssets;
        SimSeriesSP simSeries(new SimSeries(numAssets));

        for(unsigned int iBarrier = 0; iBarrier < barrierGens.size(); iBarrier++) {
            BarrierDataConstSP data = barrierGens[iBarrier]->getBarrierData();

            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                BarrierPerAssetDataConstSP assetBarrier = data->getAssetBarrier(iAsset);
                // Figure out the dates for monitoring
                DateTimeArray barrierDates;
                const DateTime& startDate = assetBarrier->monitoringDates->front();
                const string& method = data->closedFormDates;
                if(CString::equalsIgnoreCase(method, BarrierData::INPUT)) {
                    // Monitor between input dates
                    barrierDates = refDate.getFutureDates(*assetBarrier->monitoringDates);
                    if(startDate <= refDate) {
                        barrierDates.insert(barrierDates.begin(), refDate);
                    }
                } else if(CString::equalsIgnoreCase(method, BarrierData::BENCHMARK) ||
                          CString::equalsIgnoreCase(method, BarrierData::AUTO)) {
                    // Does not make sense to use interp style "N" with
                    // rolling monitoring dates
                    ScheduleSP schedule = assetBarrier->levels;
                    string interp = schedule->getInterp();
                    if(CString::equalsIgnoreCase(interp, Schedule::INTERP_NONE)) {
                        throw ModelException(
                            "Barrier schedule interpolation type " +
                             Schedule::INTERP_NONE +
                             " is inconsistent with monitorType "+
                             data->monitorType +
                             ". Use a continuous interpolation scheme.");
                    }

                    DateTimeArray benchmarkDates;
                    DateTime thisStartDate      = refDate < startDate ? startDate : refDate;
                    const DateTime& thisEndDate = assetBarrier->monitoringDates->back();
                    benchmarkDates = createDates(refDate, tenors, thisEndDate.getTime());
                    benchmarkDates = thisEndDate.getPastDates(benchmarkDates);
                    benchmarkDates = thisStartDate.getFutureDates(benchmarkDates);
                    benchmarkDates.insert(benchmarkDates.begin(), thisStartDate);
                    benchmarkDates = DateTime::merge(benchmarkDates, DateTimeArray(1, thisEndDate));

                    if(CString::equalsIgnoreCase(method, BarrierData::BENCHMARK) ||
                       schedule->isFlat()) {
                        // We are done
                        barrierDates = benchmarkDates;
                    } else {
                        // Need to add dates depending on interpolation style
                        DateTimeArray barrierCriticalDates;

                        const DoubleArray& scheduleLevels  = schedule->getValueArray();
                        const DateTimeArray& scheduleDates = schedule->getDateArray();

                        if(CString::equalsIgnoreCase(interp, Schedule::INTERP_LINEAR)) {
                            // Insert points at dates where slope changes
                            // That means we have more than 2 dates
                            if(scheduleDates.size() > 2) {
                                // According to Schedule.hpp (so much for encapsulation...)
                                // interpolation ignores the time
                                int dt1 = scheduleDates[1].getDate() - scheduleDates[0].getDate();
                                double previousSlope = dt1 ? (scheduleLevels[1] - scheduleLevels[0]) / dt1 : 0.0;

                                for(int i = 2; i < scheduleDates.size(); i++) {
                                    // According to Schedule.hpp (so much for encapsulation...)
                                    // interpolation ignores the time
                                    int dt = scheduleDates[i].getDate() - scheduleDates[i-1].getDate();
                                    double currentSlope = dt ? (scheduleLevels[i] - scheduleLevels[i-1]) / dt : 0.0;

                                    if(!Maths::equals(previousSlope, currentSlope)) {
                                        // Slope changes - add the start date of the period
                                        barrierCriticalDates.push_back(scheduleDates[i-1]);
                                    }
                                    previousSlope = currentSlope;
                                }
                            }

                            // Now insert points regularly depending on convexity of log-barrier
                        } else if(CString::equalsIgnoreCase(interp, Schedule::INTERP_STAIRS)) {
                            // Insert points at dates where barrier level changes
                            double previousLevel = scheduleLevels[0];
                            for(int i = 1; i < scheduleDates.size(); i++) {
                                double currentLevel = scheduleLevels[i];
                                if(!Maths::equals(previousLevel, currentLevel)) {
                                    // Level changes - add a date
                                    barrierCriticalDates.push_back(scheduleDates[i]);
                                }
                                previousLevel = currentLevel;
                            }
                        } else {
                            throw ModelException(
                            "Unrecognized barrier schedule interpolation type: " + interp);
                        }

                        // Use only future critical barrier dates
                        barrierCriticalDates = refDate.getFutureDates(barrierCriticalDates);

                        // Merge the benchmark dates with the critical dates
                        barrierDates = DateTime::merge(benchmarkDates, barrierCriticalDates);
                    }

                } else if(CString::equalsIgnoreCase(method, BarrierData::ENDPOINTS)) {
                    // Just use endpoints
                    const DateTime& endDate = assetBarrier->monitoringDates->back();
                    barrierDates.push_back(refDate < startDate ? startDate : refDate);
                    barrierDates.push_back(refDate < endDate ?   endDate : refDate);
                } else {
                    throw ModelException(
                        "Unrecognized methodology for inserting monitoring dates: " + method);
                }
                simSeries->addDates(iAsset, barrierDates);
            }
        }


        DateTimeArraySP driverDates(new DateTimeArray(simSeries->getAllDates()));
        return driverDates;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


DateTimeArray SVGenBarrierHVBB::BarrierDatesHelper::createDates(
    const DateTime& startDate,
    const MaturityPeriodArray& tenors,
    int time) {

    DateTimeArray maturities(tenors.size());
    for(int i = 0; i < tenors.size(); i++) {
        maturities[i] = DateTime(tenors[i]->toDate(startDate).getDate(), time);
    }
    return maturities;
}


MaturityPeriodArray SVGenBarrierHVBB::BarrierDatesHelper::createTenors() {
    MaturityPeriodArray tenors;
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("1W")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("2W")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("1M")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("2M")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("3M")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("6M")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("9M")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("1Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("2Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("3Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("4Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("5Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("7Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("10Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("20Y")));
    tenors.push_back(MaturityPeriodSP(new MaturityPeriod("30Y")));

    return tenors;
}


const MaturityPeriodArray SVGenBarrierHVBB::BarrierDatesHelper::tenors = createTenors();


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVStruct::StateVar::StateVar(SVGenBarrierHVSP barrierGen, bool isStructOut):
barrierGen(barrierGen), isStructOut(isStructOut) {
    numAssets = barrierGen->getBarrierData()->numAssets;
    hitValues = DoubleArray(numAssets);
    numHits = barrierGen->getBarrierData()->numHits;
}


double SVGenBarrierHVStruct::StateVar::hitValue() const {
    // double hitValue = 0.0;
    for(int iAsset = 0; iAsset < numAssets; iAsset++) {
        hitValues[iAsset] = barrierSV->hitValue(iAsset);
    }
    Algorithm::shellSort(hitValues);
    int idx =  isStructOut ? (numAssets - numHits) : (numHits - 1);
    return hitValues[idx];
}


void SVGenBarrierHVStruct::StateVar::update(IStateVariableGen::IStateGen* pathGen) {
    static const string routine = "MCBarrierHVStructSV::update";

    try{
        // Pull the state variables out of the PathGen:
        isPast = pathGen->doingPast();
        barrierSV = barrierGen->getHitValueSV(barrierSV, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


bool SVGenBarrierHVStruct::StateVar::doingPast() const {
    return isPast;
}


//////////////////////////////////////////////////////////////////////


SVGenBarrierHVStruct::SVGenBarrierHVStruct(SVGenBarrierHVSP barrierGen):
barrierGen(barrierGen) {
    static const string routine = "SVGenBarrierHVStruct::SVGenBarrierHVStruct";

    BarrierDataConstSP data = barrierGen->getBarrierData();
    bool isOut = data->getAssetBarrier(0)->isOut;
    for(int iAsset = 1; iAsset < data->numAssets; iAsset++) {
        if(isOut != data->getAssetBarrier(iAsset)->isOut) {
            throw ModelException(routine,
                "All assets must have common isOut flag. Check asset number " +
                Format::toString(iAsset + 1));
        }
    }

    isStructOut = isOut;
}


SVGenBarrierHVStruct::SVGenBarrierHVStruct(SVGenBarrierHVSP barrierGen, bool isStructOut):
barrierGen(barrierGen), isStructOut(isStructOut) {}


void SVGenBarrierHVStruct::collectStateVars(IStateVariableCollectorSP svCollector) const {
    svCollector->append(barrierGen.get());  // Underlying barrier per asset
}


SVGenBarrierHVStruct::StateVarSP SVGenBarrierHVStruct::getHitValueStructSV(
    StateVarSP                oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const {

    static const string routine = "SVGenBarrierHVStruct::getHitValueStrucSV";

    if(oldStateVar.get()) {
        // Update european barrier with new state variables etc.
        StateVar* perAssetBarrier = dynamic_cast<StateVar*>(oldStateVar.get());
        if(!perAssetBarrier) {
            throw ModelException(routine,
                "Expected class of type SVGenBarrierHVStruct::StateVar");
        }
        perAssetBarrier->update(pathGen);
        return oldStateVar;
    } else {
        // Create a new one
        StateVarSP barrierStruct(new StateVar(barrierGen, isStructOut));
        barrierStruct->update(pathGen);
        return barrierStruct;
    }
}



IStateVariableSP SVGenBarrierHVStruct::create(IStateVariableSP             oldStateVar,
                                           IStateVariableGen::IStateGen* pathGen) const {
    static const string routine = "SVGenBarrierHVStruct::create";

    try {
        StateVarSP oldHVStateVar(&dynamic_cast<StateVar&>(*oldStateVar));
        return getHitValueStructSV(oldHVStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


BarrierDataConstSP SVGenBarrierHVStruct::getBarrierData() const {
    return barrierGen->getBarrierData();
}


void SVGenBarrierHVStruct::StateVar::recordBarrierLevels(CControl* control,
                                                         Results*  results,
                                                         const IMultiFactors* mAsset) const {
    static const string method = "SVGenBarrierHVStruct::StateVar::recordBarrierLevels";

    try {
        // Delegate to hit value per asset state variable
        barrierSV->recordBarrierLevels(control, results, mAsset);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

void SVGenBarrierHVBarAdj::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

void SVGenBarrierHVBB::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}


//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// SVGenBarrierHVT class can generate a HitValue state variable and Time
// only for simHit & Eur Case.
////////////////////////////////////////////////////////////////////////
/** Constructor */
SVGenBarrierHVTSimEur::StateVar::StateVar(BarrierDataConstSP data,
                                    IRefLevel::IStateVarGenSP refLevelGen,
                                    SVGenSpotSP spotMonGen):
data(data),refLevelGen(refLevelGen),spotMonGen(spotMonGen),
hitNoHit(data->numAssets),metAtStepPast(0), isPastHit(false){
    static const string routine = "SVGenBarrierHVTSimEur::StateVar::StateVar";
    try{
        // Initialize hit flags to input data
        for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
            hitNoHit[iAsset]       = data->getAssetBarrier(iAsset)->isHit;
        }

        bool isOut = data->getAssetBarrier(0)->isOut;
        for(int iAsset = 1; iAsset < data->numAssets; iAsset++) {
            if(isOut != data->getAssetBarrier(iAsset)->isOut) {
                throw ModelException(routine, 
                    "All assets must have common isOut flag. Check asset number " + 
                    Format::toString(iAsset + 1));
            }
        }
        isStructOut = isOut;
    }catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Allows the barrier class to switch from past to future */
void SVGenBarrierHVTSimEur::StateVar::update(IStateVariableGen::IStateGen* pathGen){
static const string routine = "SVGenBarrierHVTSimEur::StateVar::update";
    try{ 
        // Pull the state variables out of the PathGen:
        isPast = pathGen->doingPast();
        monitoringPath = spotMonGen->getSpotSV(pathGen);
        refLevel = refLevelGen->getRefLevelSV(refLevel, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Hit value per asset */
bool SVGenBarrierHVTSimEur::StateVar::hitValueTime(double &hitValue, int &metAtStep) const{

    bool isConditionMet = isPastHit;
    // Obtain the path and the barrier information for the Asset
    // and pack into array.

    int iAsset, iStep, iBegin, iEnd;
    
    // still need a check that the path have common simDate!!
    const SVPath& path = monitoringPath->path(0);
    iBegin = path.begin();
    iEnd = path.end();

    // initialize hiValue by assuming isConditionMet = false
    hitValue = isStructOut ? 1.0 : 0.0;    

    int numRealHit = 0;
    metAtStep = metAtStepPast;

    // Compare against barrier
    for (iStep = iBegin; !isConditionMet && iStep < iEnd; iStep++) {
        for (iAsset=0; iAsset<data->numAssets; iAsset++){
            const SVPath& path = monitoringPath->path(iAsset);
            bool isUp = data->getAssetBarrier(iAsset)->isUp;
            double level = data->getAssetBarrier(iAsset)->barrierLevels[iStep];
            level *= refLevel->refLevel(iAsset);
            if ( (isUp && path[iStep] >= level) ||
                (!isUp && path[iStep] <= level)   ){
                numRealHit ++;
            }
        }
        if (numRealHit >= data->numHits){
            isConditionMet = true;
            metAtStep = iStep;
        }
    }
    // Preserve value for next call
    if(isPast) {
        isPastHit = isConditionMet;
        metAtStepPast = metAtStep;
    }
    if (isConditionMet)
        hitValue = isStructOut ? 0.0 : 1.0;    

    return isConditionMet;
}

/** Indicates past Vs future */
bool SVGenBarrierHVTSimEur::StateVar::doingPast() const{
    return isPast;    
};

/** Records BARRIER_LEVEL output requests in results */
// copy from SVGenBarrierHVEur::StateVar::recordBarrierLevels
void SVGenBarrierHVTSimEur::StateVar::recordBarrierLevels(CControl* control,
                                        Results*  results,
                                        const IMultiFactors* mAsset) const{
    
    static const string method = "SVGenBarrierHVEur::StateVar::recordBarrierLevels";
    
    try {
        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request) {
            for(int iAsset = 0; iAsset < data->numAssets; iAsset++) {
                BarrierPerAssetDataConstSP assetData = data->getAssetBarrier(iAsset);
                DateTimeArraySP monitoringDates = assetData->monitoringDates;
                const DateTime& valueDate = data->valueDate;
            
                // Only consider assets for which reference level has been determined
                // and monitoring has not been completed (this includes assets for which
                // the reference level is determined but monitoring has not started
                if(//!hitNoHitInPast[iAsset] &&
                    refLevelGen->getDates(iAsset).back() <= valueDate &&
                    valueDate <= monitoringDates->back()) {
                
                    // Get the economic barrier if it exists
                    const Schedule* schedule = assetData->economicLevels.get() && 
                                                !(assetData->economicLevels->getDates().empty()) ? 
                        assetData->economicLevels.get() : assetData->levels.get();

                    // Since we are in the middle of monitoring, the refLevel has been determined
                    double ref = refLevel->refLevel(iAsset); 
                
                    // Construct reports for barrier from now to some upper date
                    DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
                    DoubleArray values(monitoringDates->size());
                    int iDate;
                    for(iDate = 0; iDate< values.size(); iDate++) {
                        values[iDate] = schedule->interpolate((*monitoringDates)[iDate]);
                    }
                    Schedule complete(*monitoringDates.get(),
                                      values,
                                      Schedule::INTERP_NONE);
                    CashFlowArraySP subset(complete.subset(valueDate, upperDate));
                    BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));
                    for (iDate = 0; iDate < subset->size(); iDate++) {
                        // needs to be absolute
                        BarrierLevel bl(assetData->isUp,
                                        (*subset)[iDate].date,(*subset)[iDate].amount * ref,
                                        false);
                        reportLevels->push_back(bl);
                    }
                
                    // Store them in results
                    if (!reportLevels->empty()) {
                        OutputRequestUtil::recordBarrierLevels(
                            control, results,
                            mAsset->assetGetTrueName(iAsset),
                            reportLevels.get());
                    }
                }
            }
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Constuctor */
SVGenBarrierHVTSimEur::SVGenBarrierHVTSimEur(BarrierDataSP data,
                                 IRefLevel::IStateVarGenSP refLevelGen):
data(data), refLevelGen(refLevelGen) {
    static const string routine = "SVGenBarrierHVTSimEur::SVGenBarrierHVTSimEur";
    
    try {
        if(!data->european()) {
            throw ModelException("Inconsistent european barrier flags.");
        }
        
        int numAssets = data->numAssets;

        // Validate that reference level is before all monitoring dates
        int iAsset;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            const DateTimeArraySP& assetDates = data->getAssetBarrier(iAsset)->
                monitoringDates;
            const DateTime& firstAssetMonDate = assetDates->front();
            const DateTime& lastAssetRefDate  = refLevelGen->getDates(iAsset).back();
            if(firstAssetMonDate< lastAssetRefDate) {
                throw ModelException("First monitoring date " +
                                     firstAssetMonDate.toString() + 
                                     " must be after last reference level date " +
                                     lastAssetRefDate.toString() +
                                     " for asset " + Format::toString(iAsset+1));
            }
        }

        // Create an SVGenSpot for the monitoring dates
        // Stuff dates in a SimSeries object
        SimSeriesSP monSimSeries(new SimSeries(numAssets));
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            monSimSeries->addDates(iAsset, *(data->getAssetBarrier(iAsset))->monitoringDates);
        }
        spotMonGen = SVGenSpotSP(new SVGenSpot(monSimSeries));

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


/** Implementation of IStateVariableClient */
void SVGenBarrierHVTSimEur::collectStateVars(IStateVariableCollectorSP svCollector) const{
    static const string routine = "SVGenBarrierHVEur::collectStateVars";

    try {
        // Append spot and reflevel
        svCollector->append(refLevelGen.get());    // Reference level
        svCollector->append(spotMonGen.get());     // Discrete path
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}
    
/** Implementation of IStateVariableGen */
SVGenBarrierHVTSimEur::IStateVarSP SVGenBarrierHVTSimEur::getHitValueTimeSV(IStateVarSP oldStateVar,
                                IStateVariableGen::IStateGen* pathGen) const{
    static const string routine = "SVGenBarrierHVTSimEur::getHitValueSV";

    try {
        if(oldStateVar.get()) {
            // Update european barrier with new state variables etc.
            StateVar* eurBarrier = dynamic_cast<StateVar*>(oldStateVar.get());
            if(!eurBarrier) {
                throw ModelException(
                    "Expected class of type SVGenBarrierHVEur::StateVar");
            }
            eurBarrier->update(pathGen);
            return oldStateVar;
        } else {
            // Create a new one
            StateVarSP barrierEuropean(new StateVar(data, refLevelGen, spotMonGen));
            barrierEuropean->update(pathGen);
            return barrierEuropean;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}
    
/** Part of IStateVariableGen */
IStateVariableSP SVGenBarrierHVTSimEur::create(IStateVariableSP             oldStateVar,
                            IStateVariableGen::IStateGen* pathGen) const{
    static const string routine = "SVGenBarrierHVTSimEur::create";
    
    try {
        SVGenBarrierHVTSimEur::StateVarSP oldHVStateVar(&dynamic_cast<SVGenBarrierHVTSimEur::StateVar&>(*oldStateVar));
        return getHitValueTimeSV(oldHVStateVar, pathGen);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

 /** Allows access to data */
 BarrierDataConstSP SVGenBarrierHVTSimEur::getBarrierData() const{
    return data;
}

 /** Allows access to ref level */
 IRefLevel::IStateVarGenSP  SVGenBarrierHVTSimEur::getRefLevel() const{
    return refLevelGen;
}

DRLIB_END_NAMESPACE
