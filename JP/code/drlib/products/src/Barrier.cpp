//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Barrier.cpp
//
//   Description : Barrier Methods
//
//   Date        : Feb 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_BARRIER_CPP
#include "edginc/Barrier.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/RootFinder.hpp"
//#include "edginc/IDoubleArray.hpp"
#include "edginc/BarrierBend.hpp"
//#include "edginc/FunctionWrapper.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/******************************************************************/
// Barrier
/******************************************************************/

const string Barrier::EUROPEAN_MONITORING = "E";
const string Barrier::DAILY_MONITORING = "D";
const string Barrier::CONTINUOUS_MONITORING = "C";
const double Barrier::MINIMUM_SPREAD = 0.0001;
const double Barrier::ADJUST_CONSTANT = 0.5826;
// XXX Should be using trading time instead of this "hand-crafting".
const double Barrier::BUS_DAYS_IN_YEAR = 262.0;
const double Barrier::SQRT_ONE_DAY_YEAR_FRAC = sqrt(1.0/BUS_DAYS_IN_YEAR);
// XXX Should be using trading time equal to 252.0 for barrier adjustment
const double Barrier::BUS_DAYS_IN_YEAR_BA = 252.0;
const double Barrier::SQRT_ONE_DAY_YEAR_FRAC_BA = sqrt(1.0/BUS_DAYS_IN_YEAR_BA);

Barrier::~Barrier(){}

Barrier::Barrier(): CObject(TYPE)
{}

Barrier::Barrier(CClassConstSP clazz): CObject(clazz)
{}

double Barrier::digitalRoot(DigitalDiffFunc& digital,
                            double initGuess)
{
    static const string method = "Barrier::digitalRoot";
    try // XXX is this hurting performance?
    {
        /* check that lowStrike guess is positive */
        if (!Maths::isPositive(initGuess))
        {
            throw ModelException(method, "lowStrike guess is non-positive.");
        }


        /* Need to bracket the root first before we solve for it. */
        double lStrikeLower = 0.0;                 // lower bound for variance
        double lStrikeUpper = initGuess; // upper bound for variance
        ZBracPositive_bracket(digital,
                              lStrikeLower,
                              lStrikeUpper,
                              true);

        /* Root is now bracketed. Can look for it. */
        double lstrikeBound = ZBrent_solve(digital,
                                           lStrikeLower,
                                           lStrikeUpper,
                                           MINIMUM_SPREAD);

        if (!Maths::isPositive(lstrikeBound)) {
            throw ModelException(method,
                                 "ZBrent_solve returned negative adjusted barrier level.");
        }

        // We have found the low strike bound of the adjusted barrier level Badj that we are looking for.
        // The high strike will have been set in the digital structure so we want the mid value
        double adjB = (lstrikeBound + digital.highStrike) / 2.;
        if (!digital.fwdStarting) {
            // always need adjB in percentage format
            adjB /= digital.initialSpot;
        }
        return adjB;

    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/* Before the actual closed from adjustment is applied we need to create 2 arrays.
   One of the arrays contains the num bus days between mon dates, and the other
   the fwd vol between dates (interpolated at the barrier level at the end of the period)
   To keep things simple the arrays are created the same size as the interpolated barrier schedule
   and the first entry in each array is set to zero indicating no adjustment is necessary
   for the first barrier level. */
void Barrier::calculateFwdVols(const  DoubleArray&   interpBarrier,
                               const  CAsset&        asset,
                               double                refLevel,
                               const  DateTime&      startDate,
                               const  DateTimeArray& monDates,
                               const  DateTime&      today,
                               DoubleArray&          vols,
                               DoubleArray&          days,
                               DoubleArray&          spotVols,
                               DoubleArray&          spotDays)
{
   static const string method = "Barrier::calculateFwdVols";
    try
    {
        HolidaySP hols(Holiday::weekendsOnly());
        bool      fwdStarting = startDate.isGreater(today);

        // XXX This should really use methods from the ProcessedVol directly.
        for (int iStep = 0; iStep < monDates.size()-1; iStep++)
        {
            //interpLevel = pathGen->refLevel(iAsset, pathIdx);
            // No adjustment for monitoring periods in the past
            if (today >= monDates[iStep + 1])
            {
                vols[iStep + 1] = 0.;
            }
            else
            {
                double interpLevel = interpBarrier[iStep + 1];
                if (!fwdStarting) {
                    interpLevel *= refLevel;
                }
                CVolRequestLNSP req(new LinearStrikeTSVolRequest(
                    interpLevel,
                    startDate,
                    monDates.back(),
                    fwdStarting));

                CVolProcessedSP vol(asset.getProcessedVol(req.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));
                TimeMetricConstSP metric = vol->GetTimeMetric();
                double vol2 = volBS->CalcVol(today, monDates[iStep + 1]);
                double t2 = metric->yearFrac(today, monDates[iStep + 1]);
                vols[iStep + 1]  = vol2;
                spotDays[iStep + 1] = t2;
                spotVols[iStep + 1] = vol2;

                if (today < monDates[iStep])
                {
                    double t1 = metric->yearFrac(today, monDates[iStep]);
                    double vol1 = volBS->CalcVol(today, monDates[iStep]);

                    // Compute vol if there is trading time otherwise set to 0.0
                    double dt = t2 - t1;
                    if(Maths::isZero(dt)) {
                        vols[iStep + 1]  = 0.0;
                    } else {
                        vols[iStep + 1]  = sqrt((vol2*vol2*t2 - vol1*vol1*t1) / dt);
                    }
                }
            }

            // XXX not sure we need fromDate at all here?
            // if today within this monitoring period only want vol from now to next mon date
            // However, we want to include the vol from today and since busDaysDiff works from
            // [From + 1 -> To] then subtract 1
            DateTime fromDate = ((today > monDates[iStep]))? today: monDates[iStep];
            DateTime yesterday = DateTime(today.getDate()-1, today.getTime());

            days[iStep + 1] = hols->businessDaysDiff((today > monDates[iStep])? yesterday: fromDate,
                                                     monDates[iStep+1]);

        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double Barrier::DigitalDiffFunc::operator()(double lStrike)const
{
    double spread, pv;
    // Pin the highStrike to the low Strike so that we mimic changing only one variable
    if (fwdStarting)
    {
        highStrike = lStrike + MINIMUM_SPREAD;
        spread = initialSpot * (highStrike - lStrike);
    }
    else
    {
        highStrike = lStrike + (MINIMUM_SPREAD * initialSpot);
        spread = highStrike - lStrike;
    }
    pv = instSettle->pv(valueDate,
                        matDate,
                        discount,
                        asset);

    double cdf = CVanilla::priceSpread(valueDate,
                                       startDate,
                                       matDate,
                                       isCall,
                                       fwdStarting,
                                       oneContract,
                                       notional,
                                       initialSpot,
                                       lStrike,
                                       highStrike,
                                       instSettle,
                                       asset,
                                       discount);

    // Normalise and take out the pv factor
    cdf /= (pv * spread);

    cdf -= target;

    return cdf;
}

void Barrier::BarrierAdjustment(double vol, bool isUp, double& barrierLevel)
{
    double adj = exp(ADJUST_CONSTANT*vol*SQRT_ONE_DAY_YEAR_FRAC_BA);
    // if up barrier shift up / if down barrier shift down
    barrierLevel = isUp ? barrierLevel*adj : barrierLevel/adj;
}


/**  set up barrier levels for each step */
IntArray Barrier::setStepLevel(ScheduleSP in, double scale, const DateTimeArray& stepDates, CDoubleArray& out)
{
    DateTimeArray inDates = in->getDates();
    DoubleArray inValues = in->getValues();
    DateTimeArray datesHasOut;

    int i, j;
    int i_org;
    int numSteps = stepDates.size();
    // set upper barrier levels for each step
    if(in->getInterp() == Schedule::INTERP_NONE)
    {// monitor on barrier dates only
        j = stepDates[0].findUpper(inDates);  // start from the next coming barrier date.
        for (i=0; i<numSteps && j<inDates.size(); i++)
        {
            if (inDates[j].getDate() <= stepDates[i].getDate()){
                i_org = i;
                while (i<numSteps && inDates[j].equals(stepDates[i], false))
                {// do not compare time to allow multiple steps in one day
                    out[i] = scale * inValues[j];     //interpolate is no good. (It's always replaced by first barrier node on the date).
                    datesHasOut.push_back(stepDates[i]);
                    //out[i] = scale*in->interpolate(inDates[j]);
                    if (j < in->length()-1 && inDates[j].equals(stepDates[i], true))
                    {
                        if (inDates[j+1].equals(inDates[j], false))
                            j ++; //proceed one step for barrier, if next step is also on the same date.
                    }
                    i ++;
                }
                if (i>i_org)
                    i --;   // Push back to avoid check barrir on next i step.
                j ++;
                /*
                   j_org = j;
                   // InitTree() must have placed barrier dates on a time step
                   while (j<in->length() && inDates[j_org].equals(inDates[j], false))
                   j++;       // Go to next barrier step which are not on same date.
                 */                //break;
            }
        }
    }
    else
    {
        for (i=0; i<numSteps; i++)
        {
            if((in->length() == 1) && stepDates[numSteps-1] >= in->lastDate())
                out[i] = scale*in->lastValue();
            else
            {
                if (stepDates[i]>= in->firstDate()
                    && in->lastDate()>=stepDates[i])
                {// check if barrier is between start and end monitoring period, inclusive
                    out[i] = scale*in->interpolate(stepDates[i]);
                    datesHasOut.push_back(stepDates[i]);
                }
            }
        }
    }
    bool dummy;
    IntArray mapOut = DateTime::createMapping(stepDates, datesHasOut, dummy);
    return mapOut;
}

class BarrierHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPrivate();
        REGISTER(Barrier, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(Theta::Shift);
    }
};

CClassConstSP const Barrier::TYPE = CClass::registerClassLoadMethod(
    "Barrier", typeid(Barrier), BarrierHelper::load);
bool   BarrierLoad() {
    return (Barrier::TYPE != 0);
   }



/******************************************************************/
// BarrierSchedule
/******************************************************************/

BarrierSchedule::BarrierSchedule(): Barrier(TYPE),
     numHits(0), isHit(new IntArray(0)),
     smoothing(false),
     lowSpread(0.0), highSpread(0.0),
     isSimultaneousBreach(false), closedFormDates(BarrierData::DEFAULT),
     closedFormMethod(BarrierData::DEFAULT), maxNumDivsPerYear(BarrierData::DEFAULT_MAX_DIVS_PER_YEAR),
     economicLevels(0),
     numAssets(0), smoothHitValues(0),
     isHitPerAsset(new CBoolArray(0)),
     isHitPerAssetSoFar(new CBoolArray(0)),
     useAdjBarrier(false),
     forceDailySteps(false),
     numRealHitsSoFar(0),
     isConditionMetSoFar(false),
     globalPastHitTime(0),
     amendedIsHit(new CBoolArray(0)),
     doneAmendedIsHit(false),
     iBarrier(),
     isInternalMonitoringDates(false) {}

BarrierSchedule::BarrierSchedule(CClassConstSP clazz) : Barrier(clazz),
     numHits(0), isHit(new IntArray(0)),
     smoothing(false),
     lowSpread(0.0), highSpread(0.0),
     isSimultaneousBreach(false), closedFormDates(BarrierData::DEFAULT),
     closedFormMethod(BarrierData::DEFAULT), maxNumDivsPerYear(BarrierData::DEFAULT_MAX_DIVS_PER_YEAR),
     economicLevels(0),
     numAssets(0), smoothHitValues(0),
     isHitPerAsset(new CBoolArray(0)),
     isHitPerAssetSoFar(new CBoolArray(0)),
     useAdjBarrier(false),
     forceDailySteps(false),
     numRealHitsSoFar(0),
     isConditionMetSoFar(false),
     globalPastHitTime(0),
     amendedIsHit(new CBoolArray(0)),
     doneAmendedIsHit(false),
     iBarrier(),
     isInternalMonitoringDates(false) {}

BarrierSchedule::~BarrierSchedule() {
    // empty
}

BarrierSchedule::BarrierSchedule(ScheduleSP   levels,
                                 bool         isOut,
                                 bool         isUp,
                                 ScheduleSP   economicLevels,
                                 int          numHits,
                                 IntArraySP   isHit,
                                 DateTime     hitDate,
                                 string       monitorType,
                                 bool         smoothing,
                                 double       lowSpread,
                                 double       highSpread,
                                 bool         isSimultaneousBreach,
                                 bool         isInternalMonitoringDates):
    numHits(numHits),
    isHit(isHit),
    hitDate(hitDate),
    monitorType(monitorType),
    smoothing(smoothing),
    lowSpread(lowSpread),
    highSpread(highSpread),
    levels(levels),
    isOut(isOut),
    isUp(isUp),
    isSimultaneousBreach(isSimultaneousBreach),
    economicLevels(economicLevels),
    isHitPerAsset(new CBoolArray(0)),
    isHitPerAssetSoFar(new CBoolArray(0)),
    useAdjBarrier(false),
    forceDailySteps(false),
    globalPastHitTime(0),
    amendedIsHit(new CBoolArray(0)),
    doneAmendedIsHit(false),
    iBarrier(),
    isInternalMonitoringDates(isInternalMonitoringDates){}

// sort helper class
BarrierSchedule::IntSortHelper::IntSortHelper(IntArray ints): ints(ints){}

// (ie smallest to largest)
bool BarrierSchedule::IntSortHelper::operator()(int i1, int i2){
    return (ints[i1] < ints[i2]);
}

string BarrierSchedule::monTypeString() const {
    string monType;
    if(CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
        monType = BarrierBreach::EUROPEAN;
    } else if(CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) {
        monType = BarrierBreach::DAILY;
    } else if(CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) {
        monType = BarrierBreach::CONTINUOUS;
    } else {
        monType = BarrierBreach::NOT_APPLICABLE;
    }
    return monType;
}

/** Creates a SVGenBarrierHV from a BarrierSchedule */
SVGenBarrierHVSP BarrierSchedule::convertBarrier(const DateTimeArray& monitoringDates,
                                                 const DateTime& smoothDate,
                                                 const DateTime& valueDate,
                                                 IRefLevel::IStateVarGenSP refLevelGen,
                                                 const IMultiFactors* mAsset) const {

    static const string method = "BarrierSchedule::convertBarrier";

    try {

        DateTimeArraySP monDates(new DateTimeArray(monitoringDates));
        // Create barrier per asset
        vector<BarrierPerAssetDataSP> barrierPerAsset(numAssets);
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            //bool isAssetHit = (*isHit)[iAsset] ? true : false;
            bool isAssetHit = (*amendedIsHit)[iAsset];
            barrierPerAsset[iAsset] = BarrierPerAssetDataSP(new BarrierPerAssetData(
                monitorType,
                levels,
                economicLevels,
                monDates,
                smoothDate,
                valueDate,
                isOut,
                isUp,
                isAssetHit));
        }

        // Create smoothing data
        SmoothingDataSP smoothingData(new SmoothingData(
            smoothing,
            lowSpread,
            highSpread,
            smoothDate));

        // Create the big structure
        BarrierDataSP data(new BarrierData(
            monitorType,
            closedFormDates,
            closedFormMethod,
            maxNumDivsPerYear,
            barrierPerAsset,
            numHits,
            smoothingData,
            hitDate,
            valueDate));

        SVGenBarrierHVSP barrierHV;
        if(data->european() || forceDailySteps) {
            barrierHV = SVGenBarrierHVSP(new SVGenBarrierHVEur(data, refLevelGen));
        } else {
            // Figure out what methodology we want
            const string& cfMethod = data->closedFormMethod;
            if(CString::equalsIgnoreCase(cfMethod, BarrierData::BROWNIAN_BRIDGE)) {
                barrierHV = SVGenBarrierHVSP(new SVGenBarrierHVBB(data, refLevelGen, mAsset));
            } else if(CString::equalsIgnoreCase(cfMethod, BarrierData::BARRIER_ADJUSTMENT)) {
                barrierHV = SVGenBarrierHVSP(new SVGenBarrierHVBarAdj(data, refLevelGen, mAsset));
            } else {
                throw ModelException("Unknown closed form methodology: " + cfMethod);
            }
        }

        return barrierHV;

    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Creates a SVGenBarrierHV from a BarrierSchedule */
SVGenBarrierHVTSP BarrierSchedule::convertBarrier(const DateTimeArray& monitoringDates,
                                                 const DateTime& smoothDate,
                                                 const DateTime& valueDate,
                                                 IRefLevel::IStateVarGenSP refLevelGen,
                                                 const IMultiFactors* mAsset,int iDummy) const {

    static const string method = "BarrierSchedule::convertBarrier";

    try {
        
        DateTimeArraySP monDates(new DateTimeArray(monitoringDates));
        // Create barrier per asset
        vector<BarrierPerAssetDataSP> barrierPerAsset(numAssets);
        for(int iAsset = 0; iAsset < numAssets; iAsset++) {
            //bool isAssetHit = (*isHit)[iAsset] ? true : false;
            bool isAssetHit = (*amendedIsHit)[iAsset];
            barrierPerAsset[iAsset] = BarrierPerAssetDataSP(new BarrierPerAssetData(
                monitorType,
                levels,
                economicLevels,
                monDates,
                smoothDate,
                valueDate,
                isOut,
                isUp,
                isAssetHit));
        }

        // Create smoothing data
        SmoothingDataSP smoothingData(new SmoothingData(
            smoothing, 
            lowSpread,
            highSpread,
            smoothDate));

        // Create the big structure
        BarrierDataSP data(new BarrierData(
            monitorType,
            closedFormDates,
            closedFormMethod,
            maxNumDivsPerYear,
            barrierPerAsset,
            numHits,
            smoothingData,
            hitDate,
            valueDate));

        SVGenBarrierHVTSP barrierHVT;
        barrierHVT = SVGenBarrierHVTSP(new SVGenBarrierHVTSimEur(data, refLevelGen));
        return barrierHVT;

    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


void BarrierSchedule::validatePop2Object()
{
    static const string method = "BarrierSchedule::validatePop2Object";

    try
    {
        if (!isHit.get() ||
            isHit->size()<1) {
            // XXX should throw error here, but IMS requires that we build ANYTHING
            // XXX and only complain later... Do Me!
        }
        numAssets = isHit->size();
        smoothHitValues = DoubleArraySP(new DoubleArray(numAssets));
        pastHitTime = IntArray(numAssets);

        // Initialize it to input
        amendedIsHit->resize(isHit->size());
        for(int i = 0; i < isHit->size(); i++) {
            (*amendedIsHit)[i] = (*isHit)[i] ? true: false;
        }

        // In case of European monitoring assume that the input isHit is correct
        // (it is also irrelavant...)
        if(CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
            doneAmendedIsHit = true;
        }

        if ((CString::equalsIgnoreCase(monitorType, DAILY_MONITORING) ||
             CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) &&
            levels->length() < 2)
        {
            // XXX this will cause probs in IMS - need to defer the exception
            throw ModelException(method,
                                 "The barrier schedule must contain more than one value "
                                 "to allow daily or continuous monitoring.");
        }

        if (isSimultaneousBreach &&
            !CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
            throw ModelException(method,
                                 "Simultanous breach requires E-style barrier.");
        }

        if (isSimultaneousBreach && smoothing) {
            throw ModelException(method,
                                 "Smoothing is not available when using Simultanous breach.");
        }
        // need to enforce users putting in economic levels but don't want to make
        // the filed compulsory to proetct all our old tests. Validate instead to
        // stop schedules which are non-null but which have no dates
        if (economicLevels.get() && economicLevels->getDates().empty()) {
            throw ModelException(method,
                                 "Contractual barrier schedule is empty.");
        }

        // avoid alloc anywhere near simulation
        areMet = IntArray(numAssets);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

void BarrierSchedule::validate(int nbAssets)const
{
    static const string method = "BarrierSchedule::validate";
    try
    {
        // check numHits required for breach is <= num assets
        if (numHits > nbAssets || numHits <= 0)
        {
            throw ModelException(method,
                                 "num hits " +
                                 Format::toString(numHits) +
                                 " required for breach must be > 0 and <= " +
                                 Format::toString(nbAssets) +
                                 " (number of assets)");
        }

        // check isHit array is same size as asset array
        if (isHit->size() != nbAssets)
        {
            throw ModelException(method,
                                 "isHit array length " +
                                 Format::toString(isHit->size()) +
                                 " is not equal the number of assets " +
                                 Format::toString(nbAssets));
        }

        // check for a meaningful monitor type
        if (!(CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) &&
            !(CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) &&
            !(CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)))
        {
             throw ModelException(method,
                                  "Barrier monitoring must be one of 'European' E, 'Daily' D, or Continuous C.");
        }

        if (!(CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)))
        {
            useAdjBarrier = true;
        }

        if (smoothing)
        {
            // ensure lowSpread is not positive
            if (Maths::isPositive(lowSpread))
            {
                throw ModelException(method,
                                     "lowSpread ("+
                                     Format::toString(lowSpread)+
                                     ">0.0) must be 0 or -ve");
            }

            // ensure highSpread is not negative
            if (Maths::isNegative(highSpread))
            {
                throw ModelException(method,
                                     "highSpread ("+
                                     Format::toString(highSpread)+
                                     "<0.0) must be 0 or +ve");
            }

            if (!Maths::isPositive(highSpread - lowSpread))
            {
                throw ModelException(method,
                                     "highSpread-lowSpread (" + Format::toString(highSpread - lowSpread) +
                                     ") must be > 0");
            }
        }
        isHitPerAsset->resize(isHit->size());
        isHitPerAssetSoFar->resize(isHit->size());
        for (int i = 0; i < isHitPerAssetSoFar->size(); i++)
        {
            (*isHitPerAssetSoFar)[i] = (*isHit)[i]?true:false;
        }

        // check numHits required for breach is <= num assets
        if (numHits > nbAssets || numHits <= 0)
        {
            throw ModelException(method,
                                 "num hits " +
                                 Format::toString(numHits) +
                                 " required for breach must be > 0 and <= " +
                                 Format::toString(nbAssets) +
                                 " (number of assets)");
        }

        // check isHit array is same size as asset array
        if (isHit->size() != nbAssets)
        {
            throw ModelException(method,
                                 "isHit array length " +
                                 Format::toString(isHit->size()) +
                                 " is not equal the number of assets " +
                                 Format::toString(nbAssets));
        }

    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

void BarrierSchedule::validate(const MonteCarlo* model,
                               int nbAssets) {

    static const string routine = "BarrierSchedule::validate";

    try {
        // Basic validation
        validate(nbAssets);

        // Ensure closedForm strings are used meaningfully
        BarrierData::validateMethodology(monitorType, model->stateVarUsed(),
                                         closedFormMethod, closedFormDates);

        // Validate barrier schedule for "D" monitoring.
        // All dates must have common time
        if(CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) {
            int time = levels->getDateArray()[0].getTime();
            if(!levels->timesAreAll(time)) {
                throw ModelException(
                    "All monitoring dates on barrier schedule must have identical "
                    "time of day for D monitoring.");
            }
        }
        // Validate model type for C or D monitor, which is only supported by MCLN and MCImplied
        // Refine this to allow special insert daily monitoring dates when we have LV. This is
        // not ideal since we rely on the product to correctly trigger the insertion of daily dates.
        smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();
        if ((CString::equalsIgnoreCase(monitorType, DAILY_MONITORING) ||
             CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) &&
            CClass::forName("MCPathConfigLV")->isInstance(pathConfig)) {
            // (can be) ok
        } else if (CString::equalsIgnoreCase(monitorType, DAILY_MONITORING) ||
                   CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) {
            CClassConstSP MCLNTYPE = CClass::forName("MCPathConfigLN");
            CClassConstSP MCITYPE = CClass::forName("MCPathConfigImplied");
            if(!(MCLNTYPE->isInstance(pathConfig) ||
                 MCITYPE->isInstance(pathConfig))){
                throw ModelException("Model of type "+ pathConfig->getClass()->getName() +
                                     " does not support "+ monitorType +" monitoring.");
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void BarrierSchedule::amendIsHit(const IMultiFactors* mFactors,
                                 IRefLevel::IStateVarSP refLevelSV,
                                 const DateTimeArrayArray& allDates,
                                 const DateTime& today) {
    static const string routine = "BarrierSchedule::amendIsHit";


    try {
        if(!doneAmendedIsHit) {
            // amend the input flags by comparing against spot

            // If outside monitoring interval or if daily but current time is not schedule time
            // then we do not have extra information so use the input flag directly
            bool isDailyMonitoring = CString::equalsIgnoreCase(
                monitorType, BarrierData::DAILY_MONITORING);
            bool isContMonitoring = CString::equalsIgnoreCase(
                monitorType, BarrierData::CONTINUOUS_MONITORING);

            for(int iAsset = 0; iAsset < numAssets; iAsset++) {
                // See if barrier is live
                const DateTimeArray& monDates = allDates[iAsset];
                if( monDates.front() <= today && today <= monDates.back() ) {
                    bool isTodayMonitoringPoint =
                        (isDailyMonitoring && today.getTime() ==  monDates.front().getTime() ) ||
                         isContMonitoring;
                    if(isTodayMonitoringPoint) {
                        double refLevel = refLevelSV->refLevel(iAsset);
                        double spot = mFactors->assetGetSpot(iAsset);
                        double barrier = refLevel * levels->interpolate(today);
                        bool assetHitsNow = ( isUp && spot >= barrier) ||
                                            (!isUp && spot <= barrier);

                        // amend input flag
                        (*amendedIsHit)[iAsset] = (*isHit)[iAsset] || assetHitsNow;
                    } else {
                        // No more information than input flag
                        (*amendedIsHit)[iAsset] = (*isHit)[iAsset] ? true : false;
                    }
                }
            }

            // Set the flag to true to avoid doing it again
            doneAmendedIsHit = true;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void BarrierSchedule::createInterpBarrier(const DateTime& valueDate,
                                          const DateTimeArray& samples)const
{
    static const string method = "BarrierSchedule::createInterpBarrier";
    try
    {
        if (samples.size()<1) {
            throw ModelException(method,
                                 "Require at least one sample date!");
        }
        // constness forces a copy here
        // NB for dates on/before valueDate we use the legal terms if possible
		// XXX need to implement this correctly for cts/daily monitoring with AssetHistory
        barrierDates = DateTimeArraySP(new DateTimeArray(samples));
        // set up the mapping array
        isMonitorStep.resize(samples.size(),false);
        if (isInternalMonitoringDates &&
            levels->getInterp() == Schedule::INTERP_NONE){
            bool dummy;
            barMap = DateTime::createMapping(*barrierDates, levels->getDates(), dummy);
        }
        else
            barMap.resize(samples.size()+1, 0); // all days are sampling date
        if (!(interpBarrier.get()))
        {
            DoubleArraySP barLevels(new DoubleArray(samples.size(),-1.0));
            int idx = 0;
            for (idx = barMap[idx]; idx < samples.size(); idx++, idx += barMap[idx])
            {
                if (isInternalMonitoringDates &&
                    samples[idx] < levels->firstDate() || samples[idx] > levels->lastDate()){
                    isMonitorStep[idx] = false;
                }
                else{
                    isMonitorStep[idx] = true;
                    if (economicLevels.get() &&
                               samples[idx] <= valueDate ) {//use LegalTerm for Past Sampling
                        (*barLevels)[idx] = economicLevels->interpolate(samples[idx]);
                    } else {
                        (*barLevels)[idx] = levels->interpolate(samples[idx]);
                    }
                }
            }
			origDate = valueDate;
            interpBarrier = barLevels;
        }
        int nbAssets = isHit->size();
        adjBarLevels = DoubleArrayArraySP(new DoubleArrayArray(nbAssets));
        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
            (*adjBarLevels)[iAsset] = *interpBarrier; // initialize with interp barrier
        }
        // This is not appropriate to the name of this function - but then the
        // name is trash anyway.
        numRealHitsSoFar = 0;
        for (int i = 0; i < isHitPerAssetSoFar->size(); i++)
        {
            if ((*isHitPerAssetSoFar)[i]) {
                numRealHitsSoFar++;
            }
        }
        isConditionMetSoFar = (numRealHitsSoFar >= numHits);
        if (isConditionMetSoFar) {
            // then we should know when ...
            // Not concerned about performance here (it's all done)
            // Also, we are a little forgiving and pick the nearest date
            // after or on hitDate
            if (hitDate < samples[0]) {
                throw ModelException(method,
                                     "Past hit date (" + hitDate.toString() +
                                     ") cannot be before first monitoring date (" +
                                     samples[0].toString() + ")");
            }
            bool found = false;
            for(int j=0; !found && j<samples.size(); j++) {
                if (hitDate <= samples[j]) {
                    if (j==0 || hitDate > samples[j-1]) {
                        found = true;
                        if (CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
                            // require date is a match
                            if (hitDate != samples[j]) {
                                throw ModelException(method,
                                                     "For European monitoring past hit date (" + hitDate.toString() +
                                                     ") must be a monitoring date.");
                            }
                        }
                        globalPastHitTime = j;
                    }
                }
            }
            if (!found) {
                throw ModelException(method,
                                     "Past hit date (" + hitDate.toString() +
                                     ") cannot be after final monitoring date (" +
                                     samples.back().toString() + ")");
            }
        }

    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// makes sure the barrier today is the risk one for event detection
// Note if we're doing legal events the risk barrier would have
// already been overridden with the legal one so we're okay
void BarrierSchedule::overrideEventBarrier(const DateTime& valueDate)const
{
    static const string method = "BarrierSchedule::overrideEventBarrier";
    try {
		// only do anything if we've built the event barrier
        if (interpBarrier.get()) {
            for (int i = 0; i < barrierDates->size(); i++) {
                if ((*barrierDates)[i] == valueDate ) {
                    (*interpBarrier)[i] = levels->interpolate((*barrierDates)[i]);
				}
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Allow barrier to influence mon dates in inst if Daily
DateTimeArray BarrierSchedule::getFutureMonitoringDates(const DateTime&         today,
                                                        const DateTimeArray&    instMonitorDates,
                                                        smartPtr<IMCPathConfig> pathConfig) {
    // Currently condition on LV Monte Carlo model and either daily
    // or continuous monitoring. The latter is not strictly correct
    // but is not too bad in context of Archimedes (it'll give a number
    // instead of failing).
    if (CString::equalsIgnoreCase(monitorType, DAILY_MONITORING) ||
        CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING)) {
        CClassConstSP MCLVTYPE = CClass::forName("MCPathConfigLV");
        if(MCLVTYPE->isInstance(pathConfig)){
            // Class attribute in order to maintain cache of info for this decision
            forceDailySteps = true;
        }
    }

    DateTimeArray newMonitorDates;
    if (forceDailySteps && today < instMonitorDates.back()) {
        HolidaySP weekEndHols(Holiday::weekendsOnly()); // which hols I wonder?
        // add daily dates from interval (firstDate, mat] where
        // firstDate is the first future monitoring date
        int i;
        for(i=0; i<instMonitorDates.size() && instMonitorDates[i]<=today; i++) {
            // capture past dates
            newMonitorDates.push_back(instMonitorDates[i]);
        }
        DateTime nextDate = instMonitorDates.front();
        if (i>0) {
            // Have some past so start from today with daily rather than from
            // first monitoring date.
            // If we do daily from "today" the time will be AM. Rather take time from
            // the first instMonitorDates (and assume it will be the same for all).
            nextDate = DateTime(today.getDate(), nextDate.getTime());
            // See if nextDate is greater than today else add a date
            // See if today is a holiday else add a data
            if(nextDate <= today || weekEndHols->isHoliday(nextDate)) {
                nextDate = weekEndHols->addBusinessDays(nextDate, 1);
            }
        }
        for(; i<instMonitorDates.size(); i++) {
            while (nextDate<instMonitorDates[i]) {
                newMonitorDates.push_back(nextDate);
                nextDate = weekEndHols->addBusinessDays(nextDate, 1);
            }
            nextDate = instMonitorDates[i];
        }
        newMonitorDates.push_back(instMonitorDates.back());

        // Here is why this method is not const:
        // Need to make the rest of the code look at the appropriate barrier level.
        useAdjBarrier = false;

    } else {
        newMonitorDates = instMonitorDates;
    }
    return newMonitorDates;
}


ScheduleSP BarrierSchedule::getBarrierSchedule() const {
    return levels;
}

// Get the interp level from the barrier : use the last level
double BarrierSchedule::getInterpLevel(const IMCPathGenerator* pathGen,
                                       const IMCProduct* product,
                                       int iAsset)const
{
    static const string method = "BarrierSchedule::getInterpLevel";

    try
    {
        const  DateTime& today = product->getToday();
        const  DateTime&    startDate = product->getRefLevel()->getAllDates().front();
        bool   fwdStarting = startDate.isGreater(today);

        double interpLevel = interpBarrier->back();

        if (!fwdStarting)
        {
            interpLevel *= pathGen->refLevel(iAsset, 0);
        }

        return interpLevel;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// At the moment the legacy remains with a barrier class containing N assets
// and no easy way to split this up (if only it was interfaces ....). So
// for now allow 2 ways to access the "ko probs" - per asset and overall...
double BarrierSchedule::hitValue(const  IMCPathGenerator*  pathGen,
                                 int    iAsset) {
    static const string method = "BarrierSchedule::hitValue";

    if (isSimultaneousBreach) {
        throw ModelException(method, "SimultaneousBreach is not supported here");
    }

    // initialise from latest past values
    (*isHitPerAsset)[iAsset] = (*isHitPerAssetSoFar)[iAsset];

    bool isAssetHit = (*isHitPerAsset)[iAsset];
    bool wasAssetHit = isAssetHit;  // value from past
    bool doingPast = pathGen->doingPast();

    const double* path = pathGen->Path(iAsset, 0);
    int beginIdx = pathGen->begin(iAsset);
    int endIdx   = pathGen->end(iAsset);
    double refLevel = pathGen->refLevel(iAsset, 0);

    for (int iStep=beginIdx+barMap[beginIdx]; !isAssetHit && iStep<endIdx; iStep++, iStep += barMap[iStep] ) {
        double barrierLevel;
        if (isMonitorStep[iStep]){
            if (doingPast || !useAdjBarrier)
            {
                barrierLevel = refLevel * (*interpBarrier)[iStep];
            }
            else
            {
                barrierLevel = refLevel * (*adjBarLevels)[iAsset][iStep];
            }

            if ((isUp && path[iStep] >= barrierLevel) ||
                (!isUp && path[iStep] <= barrierLevel))
            {
                isAssetHit = true;
            }
        }
    }

    if (doingPast) {
        // update history
        (*isHitPerAssetSoFar)[iAsset] = isAssetHit;
    }

    // XXX The rest is not meaningful unless the path is complete, but I'm not sure how to catch that

    // XXX These should be in the class and calc'd just once - need RT class for init
    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    double payoff;
    if (!smoothing || wasAssetHit || !pathGen->Path(iAsset, 0)) { // do not smooth past hits
        if (isAssetHit) {
            payoff = binPayOnHit;
        } else {
            payoff = binPayNoHit;
        }
    } else {
        double dSpread = highSpread - lowSpread;
        double OoI = isOut? -1.0 : 1.0;
        double UoD = isUp? 1.0 : 0.0;

        double B = (useAdjBarrier)? ((*adjBarLevels)[iAsset]).back(): levels->lastValue();
        double St = pathGen->Path(iAsset, 0)[pathGen->end(iAsset)-1]; // at mat; pathIdx=0 forced here
        double S0 = pathGen->refLevel(iAsset,0); // pathIdx fixed 0 here ...
        double fwd = St/S0 - B;
        double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                     - 2.0 * Maths::max(0., -fwd))/dSpread;
        double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                     - 2.0 * Maths::max(0., fwd))/dSpread;
        double smoothFactor = UoD * (cb - pb);
        if (isAssetHit) {
            payoff = binPayOnHit + OoI * (smoothFactor - cb);
        } else {
            payoff = binPayNoHit + OoI * (smoothFactor + pb);
        }
    }
    return payoff;
}

/** This is the multi-asset version of above. Does a similar thing for
    each asset, but then combines those into a single "hit/no hit" value. */
// NBNB This only works if same dates for all assets!
double BarrierSchedule::hitValue(const  IMCPathGenerator*  pathGen) {
    static const string method = "BarrierSchedule::hitValue(multi-asset)";

    int iAsset;
    bool doingPast = pathGen->doingPast();

    if (isSimultaneousBreach) {
        throw ModelException(method, "SimultaneousBreach not supported here");
    }

    // First determine whether each asset has hit (actually)
    // Initialise from latest past values, for all assets note
    *isHitPerAsset = *isHitPerAssetSoFar;
    for(iAsset=0; iAsset<numAssets; iAsset++) {
        bool& isAssetHit = (*isHitPerAsset)[iAsset];

        for (int iStep=pathGen->begin(iAsset)+barMap[pathGen->begin(iAsset)];
             !isAssetHit && iStep<pathGen->end(iAsset); iStep++, iStep += barMap[iStep]) {
             if (isMonitorStep[iStep]){
                double barrierLevel;
                double refLevel = pathGen->refLevel(iAsset, 0);
                const double* path = pathGen->Path(iAsset, 0);

                if (doingPast || !useAdjBarrier)
                {
                    barrierLevel = refLevel * (*interpBarrier)[iStep];
                }
                else
                {
                    barrierLevel = refLevel * (*adjBarLevels)[iAsset][iStep];
                }

                if ((isUp && path[iStep] >= barrierLevel) ||
                    (!isUp && path[iStep] <= barrierLevel))
                {
                    isAssetHit = true;
                }
             }
        }

        if (doingPast) {
            // update history
            (*isHitPerAssetSoFar)[iAsset] = isAssetHit;
        }
    }

    // XXX The rest is not meaningful unless the path is complete, but I'm not sure how to catch that
    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    for(iAsset=0; iAsset<numAssets; iAsset++) {
        bool isAssetHit = (*isHitPerAsset)[iAsset];
        bool wasAssetHit = (*isHitPerAssetSoFar)[iAsset];
        double& payoff = (*smoothHitValues)[iAsset];

        if (!smoothing || wasAssetHit || !pathGen->Path(iAsset, 0)) {
            if (isAssetHit) {
                payoff = binPayOnHit;
            } else {
                payoff = binPayNoHit;
            }
        } else {
            double dSpread = highSpread - lowSpread;
            double OoI = isOut? -1.0 : 1.0;
            double UoD = isUp? 1.0 : 0.0;

            double B = (useAdjBarrier)? ((*adjBarLevels)[iAsset]).back(): levels->lastValue();
            double St = pathGen->Path(iAsset, 0)[pathGen->end(iAsset)-1]; // at mat; pathIdx=0 forced here
            double S0 = pathGen->refLevel(iAsset,0); // pathIdx fixed 0 here ...
            double fwd = St/S0 - B;
            double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                         - 2.0 * Maths::max(0., -fwd))/dSpread;
            double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                         - 2.0 * Maths::max(0., fwd))/dSpread;
            double smoothFactor = UoD * (cb - pb);
            if (isAssetHit) {
                payoff = binPayOnHit + OoI * (smoothFactor - cb);
            } else {
                payoff = binPayNoHit + OoI * (smoothFactor + pb);
            }
        }
    }

    // Sorted in ascending order
    Algorithm::shellSort(*smoothHitValues);
    // In -> use the numHits'th highest value; Out -> use the numHits'th lowest value
    int idx = isOut? (numAssets - numHits) : (numHits - 1);
    return (*smoothHitValues)[idx];
}

// Combines value methods with more info - allowing treatment "at hit"
// NBNB NOT yet tested ... should be passed through the ECO tests for instance XXX
void BarrierSchedule::hitValueAndTime(const   IMCPathGenerator*  pathGen,
                                      int     iAsset,
                                      double& value,
                                      bool&   isConditionMet,
                                      int&    metAtStep) {
    static const string method = "BarrierSchedule::hitValueAndTime";

    // initialise from latest past values
    (*isHitPerAsset)[iAsset] = (*isHitPerAssetSoFar)[iAsset];

    bool isAssetHit = (*isHitPerAsset)[iAsset];
    bool wasAssetHit = isAssetHit;  // value from past
    bool doingPast = pathGen->doingPast();

    const double* path = pathGen->Path(iAsset, 0);
    int beginIdx = pathGen->begin(iAsset);
    int endIdx   = pathGen->end(iAsset);
    double refLevel = pathGen->refLevel(iAsset, 0);

    // condition here is each asset being hit
    isConditionMet = isAssetHit;
    metAtStep = pastHitTime[iAsset];
    for (int iStep=beginIdx+barMap[beginIdx]; !isAssetHit && iStep<endIdx; iStep++, iStep += barMap[iStep]) {
        if (isMonitorStep[iStep]) {
            double barrierLevel;
            if (doingPast || !useAdjBarrier)
            {
                barrierLevel = refLevel * (*interpBarrier)[iStep];
            }
            else
            {
                barrierLevel = refLevel * (*adjBarLevels)[iAsset][iStep];
            }

            if ((isUp && path[iStep] >= barrierLevel) ||
                (!isUp && path[iStep] <= barrierLevel))
            {
                isAssetHit = true;
                isConditionMet = true;
                metAtStep = iStep;
            }
        }
    }

    if (doingPast) {
        // update history
        (*isHitPerAssetSoFar)[iAsset] = isAssetHit;
        pastHitTime[iAsset] = metAtStep;
    }

    // XXX The rest is not meaningful unless the path is complete, but I'm not sure how to catch that

    // XXX These should be in the class and calc'd just once - need RT class for init
    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    if (!smoothing || wasAssetHit || !pathGen->Path(iAsset, 0)) { // do not smooth past hits
        if (isAssetHit) {
            value = binPayOnHit;
        } else {
            value = binPayNoHit;
        }
    } else {
        double dSpread = highSpread - lowSpread;
        double OoI = isOut? -1.0 : 1.0;
        double UoD = isUp? 1.0 : 0.0;

        double B = (useAdjBarrier)? ((*adjBarLevels)[iAsset]).back(): levels->lastValue();
        double St = pathGen->Path(iAsset, 0)[pathGen->end(iAsset)-1]; // at mat; pathIdx=0 forced here
        double S0 = pathGen->refLevel(iAsset,0); // pathIdx fixed 0 here ...
        double fwd = St/S0 - B;
        double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                     - 2.0 * Maths::max(0., -fwd))/dSpread;
        double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                     - 2.0 * Maths::max(0., fwd))/dSpread;
        double smoothFactor = UoD * (cb - pb);
        if (isAssetHit) {
            value = binPayOnHit + OoI * (smoothFactor - cb);
        } else {
            value = binPayNoHit + OoI * (smoothFactor + pb);
        }
    }
}

// Have tested this one via NFBinary...
// NBNB This only works for same dates per asset. This allows the asset loop
// to be inside the timestep loop ...
void BarrierSchedule::hitValueAndTime(const  IMCPathGenerator*  pathGen,
                                      double& value,
                                      bool&   isConditionMet,
                                      int&    metAtStep) {
    static const string method = "BarrierSchedule::hitValueAndTime(multi-asset)";
    int iAsset;
    bool doingPast = pathGen->doingPast();

    // First determine whether each asset has hit (actually)
    // Initialise from latest past values, for all assets note and overall condition
    *isHitPerAsset = *isHitPerAssetSoFar;
    int numRealHits = numRealHitsSoFar; // unusually this SoFar is init'd inside AND outside the payoff
    isConditionMet = isConditionMetSoFar; // what's the point continuing if this is true?!
    metAtStep = globalPastHitTime; //XXX globalPastHitTime should be init'd from hitDate.... DO ME!

    int beginStep = pathGen->begin(0/*iAsset*/);    // see comment above ...
    int endStep = pathGen->end(0/*iAsset*/);
    for (int iStep=beginStep+barMap[beginStep]; iStep<endStep; iStep++, iStep += barMap[iStep]) {
        if (isSimultaneousBreach) {
            numRealHits = 0;
        }
        if (isMonitorStep[iStep]){
            for(iAsset=0; iAsset<numAssets; iAsset++) {
                bool& isAssetHit = (*isHitPerAsset)[iAsset];
                if (!isAssetHit || isSimultaneousBreach) {
                    // see if it hits this step
                    double refLevel = pathGen->refLevel(iAsset, 0);
                    double barrierLevel;
                    if (doingPast || !useAdjBarrier) {
                        barrierLevel = refLevel * (*interpBarrier)[iStep];
                    } else {
                        barrierLevel = refLevel * (*adjBarLevels)[iAsset][iStep];
                    }
                    double pathAtStep = pathGen->Path(iAsset, 0)[iStep];
                    if ((isUp && pathAtStep >= barrierLevel) ||
                        (!isUp && pathAtStep <= barrierLevel)) {
                        isAssetHit = true;
                        numRealHits++;

                        // Has the overall hit condition been satisfied? ...
                        // Only assigned for the earliest step at which condition met.
                        // NB! Cannot exit the asset loop early since we rely on valid info
                        // for all assets (for smoothing below). Instead we are careful when
                        // updating the "condition met & when". Note these may be set in
                        // initialisation code, if the condition is met that early.
                        if (!isConditionMet && numRealHits >= numHits) {
                            isConditionMet = true;
                            metAtStep = iStep;
                        }
                    }
                }
            }
        }
    }

    if (doingPast) {
        // update history
        numRealHitsSoFar = numRealHits; // = hits flagged by user's bool + ones detected from samples here
        *isHitPerAssetSoFar = *isHitPerAsset; // array copy
        isConditionMetSoFar = isConditionMet;
        globalPastHitTime = metAtStep;
    }

    // XXX The rest is not meaningful unless the path is complete, but I'm not sure how to catch that
    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    if (!smoothing) {
        value = isConditionMet ? binPayOnHit : binPayNoHit;
    } else {
        for(iAsset=0; iAsset<numAssets; iAsset++) {
            bool isAssetHit = (*isHitPerAsset)[iAsset];
            bool wasAssetHit = (*isHitPerAssetSoFar)[iAsset];
            double& payoff = (*smoothHitValues)[iAsset];

            if (wasAssetHit || !pathGen->Path(iAsset, 0)) { // don't smooth contribution from past hit assets
                payoff = isAssetHit ? binPayOnHit : binPayNoHit;
            } else {
                double dSpread = highSpread - lowSpread;
                double OoI = isOut? -1.0 : 1.0;
                double UoD = isUp? 1.0 : 0.0;

                double B = (useAdjBarrier)? ((*adjBarLevels)[iAsset]).back(): levels->lastValue();
                double St = pathGen->Path(iAsset, 0)[pathGen->end(iAsset)-1]; // at mat; pathIdx=0 forced here
                double S0 = pathGen->refLevel(iAsset,0); // pathIdx fixed 0 here ...
                double fwd = St/S0 - B;
                double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                             - 2.0 * Maths::max(0., -fwd))/dSpread;
                double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                             - 2.0 * Maths::max(0., fwd))/dSpread;
                double smoothFactor = UoD * (cb - pb);
                if (isAssetHit) {
                    payoff = binPayOnHit + OoI * (smoothFactor - cb);
                } else {
                    payoff = binPayNoHit + OoI * (smoothFactor + pb);
                }
            }
        }

        // Sorted in ascending order
        Algorithm::shellSort(*smoothHitValues);
        // In -> use the numHits'th highest value; Out -> use the numHits'th lowest value
        int idx = isOut? (numAssets - numHits) : (numHits - 1);
        value = (*smoothHitValues)[idx];
    }
}

// Barrier Event for RainbowRangeKO (RRK) and (untested) TriggerECO 
void BarrierSchedule::getEvents(const IMCPathGenerator* pathGen,
                                EventResults* events,
                                const string barrierName,
                                const IAggregateMakerSP& bsk) {

    static const string method = "BarrierSchedule::getEvents";
    if (smoothing) {
        throw ModelException(method, "smoothing is not supported here");
    }
    if (!CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
        throw ModelException(method, "monitorType = " + monitorType +
            " is not allowed. Only European Monitoring is supported here.");
    }
    if (!isOut) {
        throw ModelException(method, "isOut = false (Knock-In) is not allowed.");
    }
    try {
        int iAsset;
        bool doingPast = pathGen->doingPast();
        int nbAssets = isHit->size();
        
        SimpleDoubleArray assetComps(nbAssets, 0.0);
        IAggregateSP koBasket = IAggregateSP(bsk->getAggregate(&assetComps));

        int beginStep = pathGen->begin(0/*iAsset*/);
        int endStep = pathGen->end(0/*iAsset*/);

        double barrierLevel = 0.0;
        double bskLvlAtStep = 0.0;

        for (int iStep=beginStep+barMap[beginStep]; iStep<endStep; ++iStep, iStep += barMap[iStep]) {
            if (isMonitorStep[iStep]) {
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetComps[iAsset] = pathGen->Path(iAsset, 0)[iStep]/pathGen->refLevel(iAsset, 0);
                }
                if (doingPast || !useAdjBarrier) {
                    barrierLevel = (*interpBarrier)[iStep];
                } else {
                    barrierLevel = (*adjBarLevels)[iAsset][iStep];
                }
                bskLvlAtStep = koBasket->aggregate();
                if ((isUp && bskLvlAtStep >= barrierLevel) ||
                    (!isUp && bskLvlAtStep <= barrierLevel)) {

                        StringArraySP assetNames(new StringArray(0));
                        DoubleArraySP assetLevels(new DoubleArray(0));
                        DoubleArraySP barrLevels(new DoubleArray(0));

                        assetNames->push_back("Basket");
                        assetLevels->push_back(bskLvlAtStep);
                        barrLevels->push_back(barrierLevel);

                        events->addEvent(new BarrierBreach((*barrierDates)[iStep],
                            barrierName,
                            EUROPEAN_MONITORING,
                            BarrierBreach::KNOCK_OUT,
                            isUp,
                            0,
                            assetNames,
                            assetLevels,
                            barrLevels));
                        break;
                }//if ((isUp && bskLvlAtStep >= barrierLevel)
            }//if (isMonitorStep[iStep])
        }// for (int iStep=beginStep+barMap[beginStep]
    }catch (exception& e) {
        throw ModelException(e,method);
    }

}
// Have tested this one via TriggerECO...
// NBNB This only works for same dates per asset.
// In barrier, MonType = C or D are designed to work.
// No smoothing is available.
// How to catch up the step, at which bsk level are in smoothing but not touch barrier???
void BarrierSchedule::hitValueAndTimeAsBasket(const  IMCPathGenerator*  pathGen,
                                              const  IAggregateMakerSP& bsk,
                                              double& value,
                                              bool&   isConditionMet,
                                              int&    metAtStep) {
    static const string method = "BarrierSchedule::hitValueAndTimeAsBasket";
    if (smoothing) {
        throw ModelException(method, "smoothing is not supported here");
    }
    if (!CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING)) {
        throw ModelException(method, "monitorType = " + monitorType +
                            " is not allowed. Only European Monitoring is supported here.");
    }
    if (!isOut) {
        throw ModelException(method, "isOut = false (Knock-In) is not allowed.");
    }

    int iAsset;
    bool doingPast = pathGen->doingPast();
    int nbAssets = isHit->size();

    isConditionMet = isConditionMetSoFar; // what's the point continuing if this is true?!
    metAtStep = globalPastHitTime; //XXX globalPastHitTime should be init'd from hitDate.... DO ME!

    if(!isConditionMet){
        SimpleDoubleArray assetComps(nbAssets, 0.0);
        IAggregateSP koBasket = IAggregateSP(bsk->getAggregate(&assetComps));

        int beginStep = pathGen->begin(0/*iAsset*/);
        int endStep = pathGen->end(0/*iAsset*/);

        double barrierLevel = 0.0;
        double bskLvlAtStep = 0.0;
        // make reference Level
        for (int iStep=beginStep+barMap[beginStep]; iStep<endStep; iStep++, iStep += barMap[iStep]) {
            if (isMonitorStep[iStep]){
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetComps[iAsset] = pathGen->Path(iAsset, 0)[iStep]/pathGen->refLevel(iAsset, 0);
                }
                if (doingPast || !useAdjBarrier) {
                    barrierLevel = (*interpBarrier)[iStep];
                } else {
                    barrierLevel = (*adjBarLevels)[iAsset][iStep];
                }
                bskLvlAtStep = koBasket->aggregate();
                if ((isUp && bskLvlAtStep >= barrierLevel) ||
                    (!isUp && bskLvlAtStep <= barrierLevel)) {
                    isConditionMet = true;
                    metAtStep = iStep;
                    break;
                }
            }
        }
        if (doingPast) {
            // update history
            isConditionMetSoFar = isConditionMet;
            globalPastHitTime = metAtStep;
        }
    }
    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    value = isConditionMet ? binPayOnHit : binPayNoHit;
}

/*
//---Report Barrier assuming the other asset is fixed---
//Model report the barrier level for each asset, by assuming the other assets is fixed.  
//For example, basket is a plain equal weight basket of 2 assets, stock A and B.  
//Barrier is, say, up out at 110%.  Now, stock A is 120% and B is 90%.  
//For stock A, the barrier level is 130%, because the stock B is fixed at 90%.  
//For stock B, barrier level of B is 100%, because stock A is fixed at 120%.
//Another example is for worst rainbow case.  For A, basket never breach barrier by changing stock A performance.  
//thus, no barrier report for stock A.   For B, barrier level is 110%.  However, if A and B is very close, user can see only worst performer's barrier level.  
//+ ) I guess no Atlas side change is necessary.  Only works in QLIB side.
//-  ) Not full information and could not be correct.  Even when the stock cross the corresponding barrier level, the inmt wouldn't be breached the barrier when the other asset moves (a lot.)
class BarrierSchedule::BskPerfFunc{
public:
    // Full constructor 
    BskPerfFunc(const  IAggregateMakerSP& bsk,
                DoubleArray assetPerf,
                int nbAssets,
                int assetIdx):
    assetComps(nbAssets, 0.0),nbAssets(nbAssets),assetIdx(assetIdx),spots(assetPerf){
        for (int i=0;i<nbAssets;i++)
            assetComps[i] = assetPerf[i];
        koBasket = IAggregateSP(bsk->getAggregate(&assetComps));    
    };

    bool calcBar(double trgLvl, double* sol){
        static const string routine = "BarrierSchedule::BskPerfFunc::calcBar";
        try {
            // bracket the root
            double low = 0.0;
            double high = 1.e+10;
            BskFunc bskFunc(this, &BskPerfFunc::bskLvl);
            BskBarFunc bskBarFunc(bskFunc, -trgLvl);
            try{
                ZBrac_bracket(bskBarFunc,
                              low,
                              high);
                // find the root using ZBrent
                *sol = ZBrent_solve(bskBarFunc,
                                         low,
                                         high,
                                         1.e-5);
            } catch (exception& e) {
                return false; // not found
            }            
            return true;
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    double bskLvl(double level){
        // needs to refresh assetComp each time.  
        // currently, aggregate would sourt assetComp, for rainbow performer case.  
        for (int i=0;i<nbAssets;i++)
            assetComps[i] = spots[i];   
        // change the asset perf.
        assetComps[assetIdx]=level;
        double bskLvl = koBasket->aggregate();
        return bskLvl;
    };

    // Typedefs for Functors
    typedef BskPerfFunc* Ptr;
    typedef double (BskPerfFunc::* Func1DConstPtr)(double) ;
    typedef MemFuncWrapper<Ptr, Func1DConstPtr> BskFunc;
    typedef FuncShiftWrapper<BskFunc> BskBarFunc;

    SimpleDoubleArray assetComps;
    IAggregateSP koBasket;
    int nbAssets;
    int assetIdx;
    DoubleArray spots;   // no need, but seems to me better to keep the other asset's spot levels.
};

BarrierLevelArraySP BarrierSchedule::reportLevelsAsBasket(const  IAggregateMakerSP& bsk,
                                                          const DateTime& valueDate,
                                                          const DoubleArray assetPerf,
                                                          double          refLevel,
                                                          int             assetIdx) const {

    BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));
        
    if (Maths::isPositive(refLevel)) {         // wait until ref levels are available

        // first, get the usual barrier levels.

        /// ----------------  ////
        // Use economic barrier levels (instead of risk ones) if they exist
        const Schedule* myBar = economicLevels.get() ? 
                                        economicLevels.get() : levels.get();
        
        // report barrier levels over a date range
        DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

        // The barrier is defined in a rather loose way. The schedule defines a level (by date)
        // but does not indicate whether that level is relevant. For that we need barrierDates
        // and monitorType.
        // Easiest way may be to create a Schedule which contains all this info, then use subset() ...
        bool useSchedule = (!barrierDates || barrierDates->empty());
        CashFlowArraySP subset;
        if (CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING) && !useSchedule) {
            // dates are from barrierDates, levels interpolated from myBar.
            DoubleArray values(barrierDates->size());
            for(int i=0; i<values.size(); i++) {
                values[i] = myBar->interpolate((*barrierDates)[i]);
            }
            Schedule complete(*barrierDates.get(),
                              values,
                              Schedule::INTERP_NONE);
            subset = CashFlowArraySP(complete.subset(valueDate, upperDate));
        } else {
            // intervals so myBar is pretty much ok but need to make sure all days
            // are represented - but that's ok since they must be, else the pricing 
            // would fail. Well, ...., almost - but problem only if economic barrier
            // is not consistent with risk barrier, and checking that is getting a 
            // bit daft. Aha - I was wrong : N interp can get through for C/D monitorType
            // ... but I don't wish to throw exception here (it'll fail the pricing)
            // and does it really matter that much!? Ignore...
            subset = CashFlowArraySP(myBar->subset(valueDate, upperDate));
        }
        /// ----------------  ////

        // now, calculate the true levels.
        double epsilon = 1.e-10;
        int nbAssets = assetPerf.size();

        // set up function
        BskPerfFunc bskBar(bsk,assetPerf,nbAssets,assetIdx);

        double big = 1./epsilon;
        bool isContinuous = CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING);
        double solution;
        for (int i=0; i<subset->size(); i++){
            if (bskBar.calcBar((*subset)[i].amount, &solution)){
                (*subset)[i].amount = solution;
            }
            else{
                // no barrier level for this asset.  Give dummy for time being....
                (*subset)[i].amount = (isUp ? big : epsilon);
            }
            // needs to be absolute
            BarrierLevel bl(isUp,(*subset)[i].date,(*subset)[i].amount * refLevel,
                            isContinuous);
            reportLevels->push_back(bl);
        }
    }

    return reportLevels;
}
*/

/** Yet again we see weakness in the Barrier class : because it does not own the
    dates at which the barrier is active, we need to fuss with sim dates separately */
BarrierLevelArraySP BarrierSchedule::reportLevels(const DateTime& valueDate,
                                                  double          refLevel,
                                                  int             assetIdx) const {

    BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));

    if ((isSimultaneousBreach ||               // if all have to hit together OR
         !(*isHitPerAssetSoFar)[assetIdx]) &&  // ignore assets already "hit" in the past
        Maths::isPositive(refLevel)) {         // and wait until ref levels are available

        // Use economic barrier levels (instead of risk ones) if they exist
        const Schedule* myBar = economicLevels.get() ?
                                        economicLevels.get() : levels.get();

        // report barrier levels over a date range
        DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

        // The barrier is defined in a rather loose way. The schedule defines a level (by date)
        // but does not indicate whether that level is relevant. For that we need barrierDates
        // and monitorType.
        // Easiest way may be to create a Schedule which contains all this info, then use subset() ...
        bool useSchedule = (!barrierDates || barrierDates->empty());
        CashFlowArraySP subset;
        if (CString::equalsIgnoreCase(monitorType, EUROPEAN_MONITORING) && !useSchedule) {
            // dates are from barrierDates, levels interpolated from myBar.
            DoubleArray values;
            DateTimeArray barDates;
            int idx = 0;
            for (idx = barMap[idx]; idx < barrierDates->size(); idx++, idx += barMap[idx]){
                values.push_back(myBar->interpolate((*barrierDates)[idx]));
                barDates.push_back((*barrierDates)[idx]);
            }
            Schedule complete(barDates,
                              values,
                              Schedule::INTERP_NONE);
            subset = CashFlowArraySP(complete.subset(valueDate, upperDate));
        } else {
            // intervals so myBar is pretty much ok but need to make sure all days
            // are represented - but that's ok since they must be, else the pricing
            // would fail. Well, ...., almost - but problem only if economic barrier
            // is not consistent with risk barrier, and checking that is getting a
            // bit daft. Aha - I was wrong : N interp can get through for C/D monitorType
            // ... but I don't wish to throw exception here (it'll fail the pricing)
            // and does it really matter that much!? Ignore...
            subset = CashFlowArraySP(myBar->subset(valueDate, upperDate));
        }

        bool isContinuous = CString::equalsIgnoreCase(monitorType, CONTINUOUS_MONITORING);
        for (int i = 0; i < subset->size(); i++) {
            // needs to be absolute
            BarrierLevel bl(isUp,(*subset)[i].date,(*subset)[i].amount * refLevel,
                            isContinuous);
            reportLevels->push_back(bl);
        }
    }

    return reportLevels;
}

// is the barrier active at the given date
bool BarrierSchedule::isMonitoringPoint(const DateTimeArray& monDates,
                               const DateTime& today) {
    if (monitorType == Barrier::EUROPEAN_MONITORING) {
        try {
            today.find(monDates);
        } catch (exception&) {
            return false;
        }
    } else {
        if(monDates.front() <= today && today <= monDates.back()) {
            if (monitorType == Barrier::DAILY_MONITORING) {
                return (today.getTime() == monDates.front().getTime());
            }
        } else {
            return false;
        }
    }
    return true;
}

// checks if there have been barrier breaches
// and generates an event for each one
// used by NFB and RKO
void BarrierSchedule::getEvents(const  IMCPathGenerator*  pathGen,
                                EventResults* events,
                                bool useIsHitFlags,
                                const string barrierName,
                                StringArray* assNames) {
    // note we don't bother with any validation here. We assume all features
    // of BarrierUnion. This is only for events. Validation elsewhere will
    // stop the instrument pricing if things go wrong
    static const string method = "BarrierSchedule::getEvents";

    int iAsset;

    int numRealHits = 0;
    int numHitsToReport = 0;
    bool breached = false;
    int breachStep = -1;
    BoolArray isAssetHit(numAssets, false);
    BoolArray isAssetHitToReport(numAssets, false);

    // all assets share same dates
    int beginStep = pathGen->begin(0);
    int endStep = pathGen->end(0);

    // note can't report about these breaches as we don't have info
    if (useIsHitFlags) {
        for (iAsset = 0; iAsset < numAssets; iAsset++) {
            isAssetHit[iAsset] = (*isHit)[iAsset] ? true : false;
            if (isAssetHit[iAsset]) {
                numRealHits++;
            }
        }
    }
    if (numRealHits >= numHits) {
        breached = true;
    }

    if (isSimultaneousBreach) {
        numRealHits = 0;
    }

    for (int iStep=beginStep; iStep<endStep && !breached; iStep=iStep+barMap[iStep]+1) {
        numHitsToReport = 0;
        if (isMonitorStep[iStep]){
            for(iAsset=0; iAsset<numAssets; iAsset++) {
                if (!isAssetHit[iAsset] || isSimultaneousBreach) {
                    // see if it hits this step
                    double refLevel = pathGen->refLevel(iAsset, 0);
                    double barrierLevel = refLevel * (*interpBarrier)[iStep];
                    double pathAtStep = pathGen->Path(iAsset, 0)[iStep];
                    if ((isUp && pathAtStep >= barrierLevel) ||
                        (!isUp && pathAtStep <= barrierLevel)) {
                        isAssetHit[iAsset] = true;
                        numRealHits++;
                        // keep track of breaches
                        numHitsToReport++;
                        isAssetHitToReport[iAsset] = true;
                        // Can't exit asset loop early since we rely on valid info
                        // for all assets - we have to report all asset breaches
                        if (!breached && numRealHits >= numHits) {
                            breached = true;
                            breachStep = iStep;
                        }
                    }
                }
            }
        }
        // we don't report breaches if simultaneous and it didn't breach
        if (isSimultaneousBreach && !breached) {
            numHitsToReport = 0;
            numRealHits = 0;
        }

        //now unpick all the info - need to make sure we report in groups
        // so we get the number of remaining hits correct
        if (numHitsToReport > 0) {
            StringArraySP assetNames(new StringArray(numHitsToReport));
            DoubleArraySP assetLevels(new DoubleArray(numHitsToReport));
            DoubleArraySP barrLevels(new DoubleArray(numHitsToReport));
            int iBreach = 0;
            for(iAsset=0; iAsset<numAssets; iAsset++) {
                if (isAssetHitToReport[iAsset]) {
                    (*assetNames)[iBreach] = (*assNames)[iAsset];
                    double refLevel = pathGen->refLevel(iAsset, 0);
                    (*barrLevels)[iBreach] = refLevel * (*interpBarrier)[iStep];
                    (*assetLevels)[iBreach] = pathGen->Path(iAsset, 0)[iStep];
                    iBreach++;
                }
                isAssetHitToReport[iAsset] = false; //reset
            }

            events->addEvent(new BarrierBreach((*barrierDates)[iStep], barrierName,
                monTypeString(), isOut ? BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN,
                isUp, breached ? 0 : numHits - numRealHits,
                assetNames, assetLevels, barrLevels));
        }
    }
}

// generates an event if there has been a barrier breach
// called after MC has run the past - used by CalendarDrop
// a bit crap as the breach decision has been made outside but we
// need barrier levels and numHits etc
void BarrierSchedule::getEvents(EventResults* events,
                                const string barrierName,
                                int   thisPeriodIdx,
                                const IntArrayArray& metAtStep,
                                const StringArray& assNames,
                                const DoubleMatrix& assLevels,
                                const DoubleMatrix& refLevels) {
    static const string method = "BarrierSchedule::getEvents";

    // step through periods one at a time
    for (int periodIdx = 0; periodIdx <= thisPeriodIdx; periodIdx++) {
        int hitsThisPeriod = 0;
        for (int i = 0; i < numAssets; i++) {
            if ((*metAtStep[periodIdx])[i] > -1) {
                hitsThisPeriod++;
            }
        }
        if (hitsThisPeriod > 0) {
            int remainingHitsThisPeriod = numHits;
            // some breaches happened in this period - let's unpick them

            // let's get the indices of metAtStep in the order the
            // breaches occurred
            IntArray sortedMetAtStep(numAssets,0);
            for (int j = 1; j < numAssets; j++) {
                sortedMetAtStep[j] = j;
            }
            IntSortHelper sortHelper(*metAtStep[periodIdx]);
            sort(sortedMetAtStep.begin(), sortedMetAtStep.end(), sortHelper);

            int currentStep = -1;
            int startIdx = 0;
            int endIdx = 0;

            while(startIdx < numAssets && remainingHitsThisPeriod > 0) {
                // find all the breaches on one step
                while (startIdx < numAssets &&
                    ((*metAtStep[periodIdx])[sortedMetAtStep[startIdx]]
                                                          == currentStep)) {
                    startIdx++;
                }
                endIdx = startIdx;
                while (endIdx < numAssets &&
                        ((*metAtStep[periodIdx])[sortedMetAtStep[startIdx]] ==
                                        (*metAtStep[periodIdx])[sortedMetAtStep[endIdx]])) {
                    endIdx++;
                }
                currentStep = (*metAtStep[periodIdx])[sortedMetAtStep[startIdx]];
                int numHitsThisStep = endIdx - startIdx;
                remainingHitsThisPeriod -= numHitsThisStep;
                startIdx = endIdx;

                StringArraySP assetNames(new StringArray(numHitsThisStep));
                DoubleArraySP assetLevels(new DoubleArray(numHitsThisStep));
                DoubleArraySP barrLevels(new DoubleArray(numHitsThisStep));
                int iBreach = 0;
                for(int iAsset=0; iAsset<numAssets; iAsset++) {
                    if ((*metAtStep[periodIdx])[iAsset] == currentStep) {
                        (*assetNames)[iBreach] = assNames[iAsset];
                        double refLevel = refLevels[periodIdx][iAsset];
                        (*barrLevels)[iBreach] = refLevel
                                            * (*interpBarrier)[currentStep];
                        (*assetLevels)[iBreach] = assLevels[periodIdx][iAsset];
                        iBreach++;
                    }
                }

                events->addEvent(new BarrierBreach((*barrierDates)[currentStep],
                    barrierName, monTypeString(),
                    isOut ? BarrierBreach::KNOCK_OUT : BarrierBreach::KNOCK_IN,
                    isUp, Maths::max(0,remainingHitsThisPeriod),
                    assetNames, assetLevels, barrLevels));
            }
        }
    }
}

//// roll through time replacing risk barrier in past by legal one
bool BarrierSchedule::sensShift(Theta* theta) {
	// XXX need to implement this correctly for cts/daily monitoring with AssetHistory

	// if we've rolled at all we might need to do something
	// but only if we have already built the interpolated barrier
	// and also only if we've got a legal barrier
	// origDate is only stored after first call to createInterpBarrier()
    if (interpBarrier.get() && economicLevels.get()) {
		const Theta::Util thetaUtil = theta->getUtil(origDate);
		const DateTime& newDate = thetaUtil.getNewValueDate();

		if (newDate > origDate) {
			// barrierDates contains the dates
			// see if we've rolled over any of them
			// if so, amend the level there to legal one
			for (int i = 0; i < barrierDates->size(),
							(*barrierDates)[i] <= newDate; i++) {
				if ((*barrierDates)[i] > origDate) {
					(*interpBarrier)[i] =
								economicLevels->interpolate((*barrierDates)[i]);
				}
			}
		}
    }
	return true; // keep shifting
}

// Set BarrierSchedule to look like termsheet
bool BarrierSchedule::sensShift(LegalTerms* shift) {
    // 1. Replace levels with economicLevels.
    //    If there are no "economicLevels" we do not fail, simply proceed without change.
    if (economicLevels.get()) {
        levels = economicLevels;
    }
    // 2. Stop any bending.
    bendMaker = smartPtr<IBarrierBendMaker>(0);

    // 3. Stop any smoothing.
    smoothing = false;

    // 4. For now we shall leave the barrier adjusted approximations for "Daily".

    return true; // continue shifting
}

/**************************************************************************************/
// Several routines, somewhat hacky for multi-period monitoring. Pending redesign with StateVars
// SNN - I got this completely wrong. All the attributes should be in the Barrier class and
// there should be more hooks into the product : init, step, mat, if-doingPast etc
void BarrierSchedule::hitValueAndTimeAtStep(const IMCPathGenerator*  pathGen,
                                            int                      iStep,
                                            int                      iAsset,
                                            double                   refLevel,
                                            double&                  value,
                                            bool&                    isConditionMet,
                                            int&                     metAtStep) {
    static const string method = "BarrierSchedule::hitValueAndTimeAtStep";

    if (isSimultaneousBreach) {
        throw ModelException(method, "SimultaneousBreach not supported here");
    }

    if (!isConditionMet) {
        const double* path = pathGen->Path(iAsset, 0);
        bool doingPast = pathGen->doingPast();
        if (isMonitorStep[iStep]){
            double barrierLevel;
            if (doingPast || !useAdjBarrier)
            {
                barrierLevel = refLevel * (*interpBarrier)[iStep];
            }
            else
            {
                barrierLevel = refLevel * (*adjBarLevels)[iAsset][iStep];
            }

            if ((isUp && path[iStep] >= barrierLevel) ||
                (!isUp && path[iStep] <= barrierLevel))
            {
                isConditionMet = true;
                metAtStep = iStep;
            }
        }

        double binPayOnHit = isOut?0.0:1.0;
        double binPayNoHit = 1.0 - binPayOnHit;
        value = isConditionMet? binPayOnHit : binPayNoHit;

        if (doingPast) {
            // update history
            (*isHitPerAssetSoFar)[iAsset] = isConditionMet;
            pastHitTime[iAsset] = metAtStep;
        }
    }

}

// Relies on earlier setting of isConditionMet (e.g. from above) , and provides smoothed value
// We assume this is called at a "maturity" step!
void BarrierSchedule::hitValueAndTimeAtMat(const IMCPathGenerator*  pathGen,
                                           int                      iStep,
                                           int                      iAsset,
                                           double                   refLevel,
                                           bool                     isConditionMetSoFar,
                                           double&                  value,
                                           bool&                    isConditionMet,
                                           int&                     metAtStep) {
    static const string method = "BarrierSchedule::hitValueAndTimeAtMat";

    double binPayOnHit = isOut?0.0:1.0;
    double binPayNoHit = 1.0 - binPayOnHit;
    if (!smoothing || isConditionMetSoFar) { // do not smooth past hits
        if (isConditionMet) {
            value = binPayOnHit;
        } else {
            value = binPayNoHit;
        }
    } else if (!isMonitorStep[iStep]){
        value = binPayNoHit;
    } else {
        double dSpread = highSpread - lowSpread;
        double OoI = isOut? -1.0 : 1.0;
        double UoD = isUp? 1.0 : 0.0;

        double B = (useAdjBarrier)? ((*adjBarLevels)[iAsset])[iStep]: (*interpBarrier)[iStep];
        double St = pathGen->Path(iAsset, 0)[iStep]; // pathIdx=0 forced here
        double S0 = refLevel;
        double fwd = St/S0 - B;
        double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                     - 2.0 * Maths::max(0., -fwd))/dSpread;
        double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                     - 2.0 * Maths::max(0., fwd))/dSpread;
        double smoothFactor = UoD * (cb - pb);
        if (isConditionMet) {
            value = binPayOnHit + OoI * (smoothFactor - cb);
        } else {
            value = binPayNoHit + OoI * (smoothFactor + pb);
        }
    }
}

// Takes N values and provides a smoothed collective one; also gives single values from isConditionMet
// and metAtStep arrays.
void BarrierSchedule::multiHelper(DoubleArray&       hitValues,
                                  const BoolArray&   isConditionMet,
                                  const IntArray&    metAtStep,
                                  double&            overallHitValue,
                                  bool&              overallIsConditionMet,
                                  int&               overallMetAtStep) {
    // XXX This deliberately overrides the class variable since some assets may have been dropped
    // XXX Not entirely convinced by this...
    int numAssets = hitValues.size();
    Algorithm::shellSort(hitValues);
    // In -> use the numHits'th highest value; Out -> use the numHits'th lowest value
    int idx = isOut? (numAssets - numHits) : (numHits - 1);
    if (idx<0 || idx>=numAssets) {
        throw ModelException("BarrierSchedule::multiHelper",
                             "Too few assets (" + Format::toString(numAssets) +
                             ") to monitor for multiple (" +
                             Format::toString(numHits) + ") hits");
    }
    overallHitValue = hitValues[idx];

    overallIsConditionMet = false;
    overallMetAtStep = 0;
    int hitCount = 0;
    for (int i=0; i<numAssets; i++) {
        if (isConditionMet[i]) {
            hitCount++;
            if (hitCount>=numHits) {
                overallIsConditionMet = true;
                break; // no need to go on
            }
        }
    }
    if (overallIsConditionMet) {
        // how to treat overallMetAtStep? Need to find the earliest metAtStep for which
        // numHits have isConditionMet, so do something like ... filter the metAtStep
        // to have an array of metAtSteps which also have isConditionMet true, then
        // sort this and pick out the numHits'th one
        int areMetSize = 0;
        for(int j=0; j<metAtStep.size(); j++) {
            if (isConditionMet[j]) {
                areMet[areMetSize] = metAtStep[j];
                areMetSize++;
            }
        }
        if (areMetSize < numHits) {
            throw ModelException("internal error!");
        }
        sort(areMet.begin(), areMet.begin()+areMetSize);
        overallMetAtStep = areMet[numHits-1];
    }
}

// Takes N values and provides a smoothed and unsmoothed collective one
void BarrierSchedule::multiHelper2(DoubleArray&       hitValues,
                                   const BoolArray&   isConditionMet,
                                   double&            overallHitValue,
                                   double&            overallBinaryValue) {
    // XXX This deliberately overrides the class variable since some assets may have been dropped
    // XXX Not entirely convinced by this...
    int numAssets = hitValues.size();
    Algorithm::shellSort(hitValues);
    // In -> use the numHits'th highest value; Out -> use the numHits'th lowest value
    int idx = isOut? (numAssets - numHits) : (numHits - 1);
    if (idx<0 || idx>=numAssets) {
        throw ModelException("BarrierSchedule::multiHelper",
                             "Too few assets (" + Format::toString(numAssets) +
                             ") to monitor for multiple (" +
                             Format::toString(numHits) + ") hits");
    }
    overallHitValue = hitValues[idx];

    bool overallIsConditionMet = false;
    int hitCount = 0;
    for (int i=0; i<numAssets; i++) {
        if (isConditionMet[i]) {
            hitCount++;
            if (hitCount>=numHits) {
                overallIsConditionMet = true;
                break; // no need to go on
            }
        }
    }
    overallBinaryValue = (overallIsConditionMet && isOut) ||
        (!overallIsConditionMet && !isOut) ? 0.0:1.0;
}


//////////////////////////////////////////////////
//
// functionality to support IBarrierUtil interface
//
//////////////////////////////////////////////////

class BarrierUtilPerAsset : public IBarrierUtil
{
public:

    //true = hit data is from hit flag and not calculated
    virtual bool getHitDateAndStep(int& step, DateTime& hitDate) const{
        if (isAssetHit) {
            step = -1;
            hitDate = assetHitDt;
            return true;
        }
        else {
            step = hitStep;
            hitDate = (*simDates)[hitStep];
            return false;
        }
    }

    //returns true if it has an event the top level BarrierCoupler needs to process
    virtual bool getEvents(EventResults* events, bool useIsHitFlags, const string& barrierName, StringArray* assNames, const DateTime& hitDate, int hits) const{
        return true;
    }

    virtual bool getHitDataForBreachEvent(double& refLvl, double& barrierLevel, double& assetLevel, string& monitorType) const{
        if (isHit) {
            refLvl = assetRefLevel;
            barrierLevel = barrierHitLevel;
            assetLevel = assetHitLevel;
            monitorType=monType;
            return true;
        } else {
            return false;
        }
    }

    // we set the stickyThreshold to 0 if smoothing, so that for sticky barrier (ie. not simultaneous hit),
    // once hit, the hit value is remembered even if it's not 100% hit
    BarrierUtilPerAsset(int iAsset, bool isUp, const DateTimeArray barDates,
                    const DoubleArray barLevels, string interp, string monType, bool isAssetHit,
                    DateTime assetHitDt, bool smoothing, double lowSpread, double highSpread)
                    : iAsset(iAsset), isUp(isUp), dates(barDates), levels(barLevels), interp(interp),
                    monType(monType), isAssetHit(isAssetHit), assetHitDt(assetHitDt),
                    smoothing(smoothing), lowSpread(lowSpread), highSpread(highSpread),assetRefLevel(0.0),
                    barrierHitLevel(0.0),assetHitLevel(0.0),isHit(false)
    {
        // ideally should dis-allow daily monitoring. since we don't treat it differently from european
        lastLevel = levels.back();
    }
    ~BarrierUtilPerAsset(){};

    virtual void preprocess(DateTime valueDateIn, DateTimeArrayConstSP simDatesIn)
    {
        static const string method("BarrierUtilPerAsset::preprocess");
        valueDate = valueDateIn;
        simDates = simDatesIn;

        // sanity check
        if( !DateTime::isSubset(*simDates, dates) )
            throw ModelException(method, "barrier dates must be subset of sim dates");

        // assign levels to each step to avoid lookup
        relevantMonDts.clear();
        levelsPerStep.resize(simDates->size());
        if( CString::equalsIgnoreCase(monType, Barrier::CONTINUOUS_MONITORING) )
        {
            Schedule schedule(dates, levels, interp);
            for(int i=0; i<simDates->size(); i++)
            {
                if( (*simDates)[i] >= dates.front() && (*simDates)[i] <= dates.back() )
                {
                    try {
                        levelsPerStep[i] = schedule.interpolate((*simDates)[i]);
                    } catch (exception& e) {
                        throw ModelException(e, method,
                                             "Continuous barrier requires level for all dates, but failed at " +
                                             (*simDates)[i].toString() +
                                             ". May be due to schedule having interp type None.");
                    }
                    relevantMonDts.push_back((*simDates)[i]);
                }
                else
                    levelsPerStep[i] = -1.0;
            }
        } else {
            for(int i=0, j=0; i<simDates->size(); i++)
            {
                if( j<dates.size() && (*simDates)[i] == dates[j] )
                {
                    levelsPerStep[i] = levels[j++];
                    relevantMonDts.push_back((*simDates)[i]);
                }
                else
                    levelsPerStep[i] = -1.0;
            }
        }

    }

    /////////////////////////////////////////////////////////////////////////////////////////
    // this is not properly implemented now since monitorDates are still not part of barrier
    // therefore we return 'continuous' so long as it's not INTERP_NONE so that all sim dates
    // are relevant observation date if not INTERP_NONE
    /////////////////////////////////////////////////////////////////////////////////////////
    void relevantDates(DateTimeArray &relevantDates, bool &hasContinuous) const
    {
        relevantDates = dates;
        hasContinuous = CString::equalsIgnoreCase(monType, Barrier::CONTINUOUS_MONITORING);
    }

    virtual const DateTimeArray& relevantMonDates() const
    {
        return relevantMonDts;
    }

    void pathUpdated(const IMCPathGenerator* pathGenIn)
    { pathGen = pathGenIn; }

    // special instance to get the base schedule hit value
    double hitValue(int step) const
    {
        if( isAssetHit ) {
            if (assetHitDt > (*simDates)[step])
                return 0.0;
            else {
                hitStep = -1;
                return 1.0;
            }
        }

        double refLvl = pathGen->refLevel(iAsset, 0);
        double barrierLevel = refLvl * levelsPerStep[step];
        double spot = pathGen->Path(iAsset, 0)[step];
        isHit = (isUp && spot >= barrierLevel) || (!isUp && spot <= barrierLevel);

        assetRefLevel = refLvl;
        barrierHitLevel = barrierLevel;
        assetHitLevel = spot;
        hitStep = step;
        return smoothing?smoothHitValue(isHit, step):isHit;
    }

    const DoubleArray& getLevelsPerStep() const
    { return levelsPerStep; }

    // used to reset barrier level, used for barrier adjustment
    // should only be called for continuous monitoring
    virtual void setAdjFactors(const DoubleArray adjFactors)
    {
        if( adjFactors.size() != levelsPerStep.size() )
            throw ModelException("BarrierUtilPerAsset:setAdjLevels",
                                "Internal error. input level size must equal sim date size");

        // update levels on relevant steps. only those in the future
        for(int i=0; i<simDates->size(); i++)
        {
            if( (*simDates)[i] > valueDate )
                levelsPerStep[i] *= adjFactors[i];
        }
        lastLevel = levelsPerStep.back();
    }

protected:

    // only return smooth'ed hit value if isHit (unless it's last step).
    // to conform with current BarrierSchedule hitValueAndTime implementation
    inline double smoothHitValue(bool isHit, int step) const
    {
        if( !isHit && step < (levelsPerStep.size()-1) ||
            !pathGen->Path(iAsset, 0))
            return 0.0;

        double dSpread = highSpread - lowSpread;
        double St = pathGen->Path(iAsset, 0)[pathGen->end(iAsset)-1]; // at mat; pathIdx=0 forced here
        double S0 = pathGen->refLevel(iAsset,0); // pathIdx fixed 0 here ...
        double fwd = St/S0 - lastLevel;
        double pb = (Maths::max(0., highSpread - fwd) + Maths::max(0., -highSpread - fwd)
                     - 2.0 * Maths::max(0., -fwd))/dSpread;
        double cb = (Maths::max(0., fwd - lowSpread) + Maths::max(0., fwd + lowSpread)
                     - 2.0 * Maths::max(0., fwd))/dSpread;

        return isHit?(1.0 - (isUp?pb:cb)):(isUp?cb:pb);
    }

protected:

    // *** transient fields below ***
    DateTime                    valueDate;
    DateTimeArrayConstSP        simDates;
    const IMCPathGenerator*     pathGen;

    int                         iAsset;
    bool                        isUp;
    DateTimeArray               dates;
    DoubleArray                 levels;         // eco level interpolated for each sim date, maybe adjusted
    string                      interp;
    string                      monType;

    bool                        isAssetHit;     // for historical status override
    DateTime                    assetHitDt;     // date hit

    // For smoothing
    bool                        smoothing;         // Is smoothing to be used ?
    double                      lowSpread;
    double                      highSpread;

    // transient fields
    mutable double              lastLevel;      // last barrier level for smoothing
    DoubleArray                 levelsPerStep;
    DateTimeArray               relevantMonDts;

    mutable double              assetRefLevel;
    mutable double              barrierHitLevel;
    mutable double              assetHitLevel;
    mutable bool                isHit;
    mutable int                 hitStep;
};

typedef refCountPtr<BarrierUtilPerAsset> BarrierUtilPerAssetSP;

class BarrierUtilPerAssetBend : public BarrierUtilPerAsset, virtual public IBarrierUtil::ICanHaveBend
{
public:

    BarrierUtilPerAssetBend(int iAsset, bool isUp, const DateTimeArray barDates,
                    const DoubleArray barLevels, string interp, string monType, bool isAssetHit,
                    DateTime assetHitDt, bool smoothing, double lowSpread, double highSpread,
                    IBarrierBendMakerSP bendMaker)
                    : BarrierUtilPerAsset(iAsset, isUp, barDates, barLevels, interp, monType, isAssetHit,
                    assetHitDt, smoothing, lowSpread, highSpread), bendMaker(bendMaker)
    {}
    ~BarrierUtilPerAssetBend(){};

    virtual void preprocess(DateTime valueDateIn, DateTimeArrayConstSP simDatesIn)
    {
        static const string method = "BarrierUtilPerAssetBend::preprocess";

        BarrierUtilPerAsset::preprocess(valueDateIn, simDatesIn);

        bend = bendMaker->getBarrierBend(valueDate, *simDatesIn);

        // create mapping to bend period start/end dates
        DateTimeArray bendStartDts, bendEndDts;
        bendDates(bendStartDts, bendEndDts);
        if( bendEndDts.size() ==0 ) return;

        DateTime::ensureIncreasing(bendEndDts, method + ". barrrier bend end dates must in strict ascending order", true);
        if( !DateTime::isSubset(*simDates, bendEndDts) )
            throw ModelException(method, ". barrier bend period end dates must be subset of sim dates");

        bool dummy;
        bendEndMap = DateTime::createMapping(*simDates, bendEndDts, dummy);

        bendStartStep.resize(bendEndDts.size());
        bendAmt.resize(bendEndDts.size());
        int i, j;
        for(i=0, j=0; i<simDates->size(); i++)
        {
            if( bendEndMap[i] == 0 )
            {
                // look for start of bend period
                bendStartStep[j] = i;
                while( bendStartStep[j] >= 0 &&
                    (*simDates)[bendStartStep[j]] >= bendStartDts[j] )
                    bendStartStep[j]--;
                bendStartStep[j]++;

                bendEndMap[i] = j;
                bend->getBending(i, bendAmt[j]);
                j++;

                if( j== bendEndDts.size() )
                {
                    i++;
                    while(i<simDates->size())
                    {
                        bendEndMap[i++] = -1;
                    }
                }
            }
            else
                bendEndMap[i] = j;
        }

        int bendPeriodIdx = bendEndMap[levelsPerStep.size()-1];
        if( bendPeriodIdx>=0 )
            lastLevel *= bendAmt[bendPeriodIdx][levelsPerStep.size()-1-bendStartStep[bendPeriodIdx]];
    }

    // get barrier bend period starts/ends
    void bendDates(DateTimeArray &startDates, DateTimeArray &endDates) const
    {
        if( !bend )
        {
            startDates.clear();
            endDates.clear();
        }
        else
        {
            IntArray endIdx = bend->getBendingEnd();
            startDates.resize(endIdx.size());
            endDates.resize(endIdx.size());
            for(int i=0; i<endIdx.size(); i++)
            {
                endDates[i] = (*simDates)[endIdx[i]];
                startDates[i] = (*simDates)[bend->getBendingStartByEnd(endIdx[i])];
            }
        }
    }

    // calculate and return hit status for a step if not already calculated
    double hitValueBend(int step, int bendStep) const
    {
        if( isAssetHit ) {
            if (assetHitDt > (*simDates)[step])
                return 0.0;
            else {
                hitStep = -1;
                return 1.0;
            }
        }

        double refLevel = pathGen->refLevel(iAsset, 0);
        double barrierLevel = refLevel * levelsPerStep[step];
        int bendPeriodIdx = bendEndMap[bendStep];
        if( bendPeriodIdx>=0 && step>=bendStartStep[bendPeriodIdx] )
            barrierLevel *= bendAmt[bendPeriodIdx][step-bendStartStep[bendPeriodIdx]];
        double spot = pathGen->Path(iAsset, 0)[step];
        isHit = (isUp && spot >= barrierLevel) || (!isUp && spot <= barrierLevel);
        assetRefLevel = refLevel;
        barrierHitLevel = barrierLevel;
        assetHitLevel = spot;
        hitStep = step;
        return smoothing?smoothHitValue(isHit, step):isHit;
    }

    virtual void setAdjFactors(const DoubleArray adjFactors)
    {
        BarrierUtilPerAsset::setAdjFactors(adjFactors);

        int bendPeriodIdx;
        if(!bendEndMap.empty()) {
            bendPeriodIdx = bendEndMap[levelsPerStep.size()-1];
            if (bendPeriodIdx >=0) {
                lastLevel *= bendAmt[bendPeriodIdx][levelsPerStep.size()-1-bendStartStep[bendPeriodIdx]];
            }
        }
    }

private:
    IBarrierBendMakerSP         bendMaker;


    // transient fields
    mutable IntArray            bendEndMap;
    mutable IntArray            bendStartStep;
    mutable DoubleArrayArray    bendAmt;        // one array per bend period
    IBarrierBendSP              bend;
};

void BarrierSchedule::pathUpdated(const IMCPathGenerator* pathGen)
{
    if( !!iBarrier )
        iBarrier->pathUpdated(pathGen);
}

bool BarrierSchedule::createBarrierUtil(const DateTime& valueDate) const
{
    static const string method = "BarrierSchedule::createBarrierUtil";

    // smooothing is actually not supported. because current scheme focus on smoothing for
    // maturity option but may lead to inconsistency with rebate calc. hitValue based on
    // if an asset has been hit on any point on the path, however, rebate is chosen for
    // the 1st point of the path where num asset hit limit. for instance, if
    // numHit=1, asset A may hit barrier at step m, so rebate is chosen for step m
    // however another asset B hit barrier at a later step n, if B has higher hit value,
    // this is chosen (reasonable), however, rebate is from m.
    // our implementation for smoothing is the best effort. for regression reason, it's disabled
    bool hasBend = (!!bendMaker && !bendMaker->trivial());
    if( smoothing )
    {
        if( hasBend )
            throw ModelException(method, "Can not use smoothing with barrier bend");
        return false;
    }

    // we check barrierDates to see if createInterpBarrier was called
    // this is used to indicate if we want to use all the sim dates (old way)
    // or only dates from the level schedule (new way)
    bool useSchedule = (!barrierDates || barrierDates->empty());
    const DateTimeArray& barDates = useSchedule?levels->getDateArray():*barrierDates;
    const DoubleArray& barLevels = useSchedule?levels->getValueArray():*interpBarrier;
    if( useSchedule &&
        CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) {
        throw ModelException(method, "Can not use daily monitoring if get monitor dates directly from barrier schedule");
    }

    ibsPerAsset.resize(numAssets);
    vector<IBarrierUtilSP> ibs(numAssets);
    for(int iAsset=0; iAsset<numAssets; iAsset++)
    {
        if( !hasBend )
            ibsPerAsset[iAsset] = IBarrierUtilSP(new BarrierUtilPerAsset(
                iAsset, isUp, barDates, barLevels, levels->getInterp(), monitorType,
                (*isHit)[iAsset]?true:false, hitDate,
                smoothing, lowSpread, highSpread));
        else
            ibsPerAsset[iAsset] = IBarrierUtilSP(new BarrierUtilPerAssetBend(
                iAsset, isUp, barDates, barLevels, levels->getInterp(), monitorType,
                (*isHit)[iAsset]?true:false, hitDate,
                smoothing, lowSpread, highSpread, bendMaker));

        // asset barrier hit value sticky if not simultaneous
        ibs[iAsset] = isSimultaneousBreach?ibsPerAsset[iAsset]:BarrierFactory::makeSticky(ibsPerAsset[iAsset]);
    }

    iBarrier = BarrierFactory::makeCoupler(ibs, numHits, isOut, isUp);
    if( isOut ) iBarrier = BarrierFactory::makeReverser(iBarrier);

    // probably should only make it sticky if requested by payoff
    // but then should make sure we don't duplicate sticky versions of the barrier
    iBarrier = BarrierFactory::makeSticky(iBarrier);

    return true;
}

IBarrierUtilSP BarrierSchedule::getBarrierUtil() const
{
    return iBarrier;
}

/* Creates an adjusted barrier schedule for each asset using a closed form formula
   and the fwd vols and year fracs between for each monitor date interval */
void BarrierSchedule::adjustLN(const  IMCPathGenerator*  pathGen,
                               const  IMCProduct*         product,
                               const IMultiFactors*      assets,
                               int pathIdx) // XXX not sure this works for multiple paths
{
    static const string method = "BarrierSchedule::adjustLN";
    try
    {
        if (useAdjBarrier)
        {
            // only get here if daily or cts monitoring
            const DateTime& today = product->getToday();
            const SimSeries* simSeries = product->getSimSeries();

            for (int iAsset = 0; iAsset < assets->NbAssets(); iAsset++)
            {
                const DoubleArray *interpBarrier; // local interp barrier array to accommodate iBarrier
                if( !iBarrier )
                    interpBarrier = this->interpBarrier.get();
                else {
                    BarrierUtilPerAssetSP ibPerAsset = DYNAMIC_POINTER_CAST<BarrierUtilPerAsset>(ibsPerAsset[iAsset]);
                    interpBarrier = &ibPerAsset->getLevelsPerStep();
                }


                const CAsset& asset = assets->getAsset(iAsset); // XXX a shame to do this
                // XXX Note assumption that all simulation dates are monitoring dates
                const DateTimeArray& monDates = simSeries->getDates(iAsset);

                DoubleArray days(monDates.size());
                DoubleArray vols(monDates.size());
                DoubleArray spotDays(monDates.size()); // outputs not required for LN adjustment
                DoubleArray spotVols(monDates.size()); // outputs not required for LN adjustment
                DoubleArray adjFactors(monDates.size());
                days[0] = 0.;
                vols[0] = 0.; // No adjustment for the 1st barrier level

                double refLevel = pathGen->refLevel(iAsset, pathIdx);
                const  DateTime& startDate = product->getRefLevel()->getAllDates().front();
                // Get the fwd vols and year Fracs between barrier periods for this asset
                calculateFwdVols(*interpBarrier, asset, refLevel, startDate,
                                 monDates, today, vols, days, spotVols, spotDays);

                int UoD = isUp? -1 : 1;
                int i;
                // Now apply the closed form formula to each barrier level
                // XXX Wonder if we should be using vol via variance and internally handled trading time?
                // For non-European monitoring we adjust from interval (say monthly) to continuous
                for (i = 0; i < (*adjBarLevels)[iAsset].size(); i++) {
                    adjFactors[i] = exp(UoD * ADJUST_CONSTANT * vols[i] * sqrt(days[i]/BUS_DAYS_IN_YEAR));
                }
                if (CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) {
                    // for daily adjust back from adjusted-cts to daily interval
                    for (i = 0; i < (*adjBarLevels)[iAsset].size(); i++) {
                        adjFactors[i] *= exp(-UoD * ADJUST_CONSTANT * vols[i] * SQRT_ONE_DAY_YEAR_FRAC);
                    }
                }
                if( !iBarrier ) {
                    for (i = 0; i < (*adjBarLevels)[iAsset].size(); i++) {
                        (*adjBarLevels)[iAsset][i] = (*interpBarrier)[i] * adjFactors[i];
                    }
                } else {
                    BarrierUtilPerAssetSP ibPerAsset = DYNAMIC_POINTER_CAST<BarrierUtilPerAsset>(ibsPerAsset[iAsset]);
                    ibPerAsset->setAdjFactors(adjFactors);
                }
            }


        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/* Creates an adjusted barrier schedule for each asset using a closed form formula
   and the fwd vols and year fracs between for each monitor date interval */
void BarrierSchedule::adjustImplied(const  IMCPathGenerator*  pathGen,
                                    const  IMCProduct* product,
                                    const IMultiFactors*     assets,
                                    const YieldCurve* discount,
                                    InstrumentSettlement* instSettle,
                                    int pathIdx)
{
    static const string method = "BarrierSchedule::adjustImplied";
    try
    {
        if (useAdjBarrier)
        {
            const DateTime& today = product->getToday();
            const SimSeries* simSeries = product->getSimSeries();
            const DateTime& startDate = product->getRefLevel()->getAllDates().front();
            bool  fwdStarting         = startDate.isGreater(today);

            for (int iAsset = 0; iAsset < assets->NbAssets(); iAsset++)
            {
                const DoubleArray *interpBarrier; // local interp barrier array to accommodate iBarrier
                if( !iBarrier )
                    interpBarrier = this->interpBarrier.get();
                else {
                    BarrierUtilPerAssetSP ibPerAsset = DYNAMIC_POINTER_CAST<BarrierUtilPerAsset>(ibsPerAsset[iAsset]);
                    interpBarrier = &ibPerAsset->getLevelsPerStep();
                }

                const CAsset& asset = assets->getAsset(iAsset);
                const DateTimeArray& monDates = simSeries->getDates(iAsset);

                DoubleArray days(monDates.size());
                DoubleArray vols(monDates.size());
                DoubleArray spotDays(monDates.size());
                DoubleArray spotVols(monDates.size());
                DoubleArray adjFactors(monDates.size());
                days[0] = 0.;
                vols[0] = 0.; // No adjustment for the 1st barrier level

                double refLevel = pathGen->refLevel(iAsset, pathIdx);
                const  DateTime& startDate = product->getRefLevel()->getAllDates().front();

                // Get the fwd vols and year Fracs between barrier periods for this asset
                calculateFwdVols(*interpBarrier, asset, refLevel, startDate,
                                 monDates, today, vols, days, spotVols, spotDays);

                adjFactors[0] = 1.0; // No adjustment for the first barrier level

                // Now create the adjusted barrier schedule for this asset
                double lowStrike, highStrike;
                DigitalDiffFunc digitalDiff;
                // "static" across steps
                digitalDiff.valueDate = today;
                digitalDiff.startDate = startDate;
                digitalDiff.isCall = false;
                digitalDiff.fwdStarting = fwdStarting;
                digitalDiff.oneContract = true;
                digitalDiff.notional = 1.;
                digitalDiff.instSettle = instSettle;
                digitalDiff.discount = discount;
                digitalDiff.initialSpot = pathGen->refLevel(iAsset, pathIdx);
                digitalDiff.asset = &assets->getAsset(iAsset);

                for (int iStep = 1; iStep < monDates.size(); iStep++)
                {
                    // No adjustment for monitoring periods in the past
                    if (today >= monDates[iStep])
                    {
                        adjFactors[iStep] = 1.0;
                    }
                    else
                    {
                        double spread, pv;
                        // Determine the low and high strikes for this barrier level B(t)
                        if (fwdStarting)
                        {
                            lowStrike = (*interpBarrier)[iStep] - (MINIMUM_SPREAD/2.);
                            highStrike = (*interpBarrier)[iStep] + (MINIMUM_SPREAD/2.);
                            spread = refLevel * (highStrike - lowStrike);
                        }
                        else
                        {   // get absolute levels
                            lowStrike = ((*interpBarrier)[iStep] - (MINIMUM_SPREAD/2.)) * refLevel;
                            highStrike = ((*interpBarrier)[iStep] + (MINIMUM_SPREAD/2.)) * refLevel;
                            spread = highStrike - lowStrike;
                        }

                        const DateTime& matDate = monDates[iStep];

                        // find the cdf of B(t) by calculating the put spread with strikes around B(t)
                        double digital = CVanilla::priceSpread(today,
                                                               startDate,
                                                               matDate,
                                                               false,  // PUT
                                                               fwdStarting,
                                                               true,   // one contract
                                                               1.,     // notional - not used here
                                                               refLevel,
                                                               lowStrike,
                                                               highStrike,
                                                               instSettle,
                                                               &assets->getAsset(iAsset),
                                                               discount);

                        // This digital price contains a PV so we need to 'extract' this before doing the inverse
                        pv = instSettle->pv(today,
                                            matDate,
                                            discount,
                                            &assets->getAsset(iAsset));

                        digital /= (spread * pv);

                        /* Dont do the implied adjustment if probability is above 99% or below 1%.
                           Instead, at these extremes do a LN shift instead */
                        if (digital > 0.99 || digital < 0.01)
                        {
                            int UoD = isUp? -1 : 1;
                            adjFactors[iStep] =
                                exp(UoD * ADJUST_CONSTANT * vols[iStep] * sqrt(days[iStep]/BUS_DAYS_IN_YEAR));

                            if (CString::equalsIgnoreCase(monitorType, DAILY_MONITORING))
                            {
                                // for daily adjust back from adjusted-cts to daily interval
                                adjFactors[iStep] *= exp(-UoD * ADJUST_CONSTANT * vols[iStep] * SQRT_ONE_DAY_YEAR_FRAC);
                            }
                        }
                        else
                        {
                            // find the corresponding random variable from this distribution
                            double X = N1Inverse(digital);

                            double sqrtDaysYearFrac = sqrt(days[iStep]/BUS_DAYS_IN_YEAR);
                            if (CString::equalsIgnoreCase(monitorType, DAILY_MONITORING)) {
                                sqrtDaysYearFrac -= SQRT_ONE_DAY_YEAR_FRAC;
                            }
                            int UoD = isUp? -1 : +1;
                            // shift the random variable using the magic formulas
                            double Xadj = X + UoD * (ADJUST_CONSTANT * (vols[iStep] * sqrtDaysYearFrac))/
                                (spotVols[iStep] * sqrt(spotDays[iStep]));

                            // Fill in the "volatile" DigitalDiffFunc fields
                            digitalDiff.matDate = matDate;
                            digitalDiff.highStrike = highStrike;
                            digitalDiff.target = N1(Xadj);  // this is cdf of the adjusted random variable

                            // The adjusted barrier level Badj(t) is now that barrier level that satisifes
                            // N1(Xadj) = Digital(Badj(t))
                            adjFactors[iStep] = digitalRoot(digitalDiff, lowStrike)/(*interpBarrier)[iStep];
                        }

                    }
                }

                if( !iBarrier ) {
                    for (int i = 0; i < (*adjBarLevels)[iAsset].size(); i++) {
                        (*adjBarLevels)[iAsset][i] = (*interpBarrier)[i] * adjFactors[i];
                    }
                } else {
                    BarrierUtilPerAssetSP ibPerAsset = DYNAMIC_POINTER_CAST<BarrierUtilPerAsset>(ibsPerAsset[iAsset]);
                    ibPerAsset->setAdjFactors(adjFactors);
                }

            }
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

class BarrierScheduleHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(BarrierSchedule, clazz);
        SUPERCLASS(Barrier);
        EMPTY_SHELL_METHOD(defaultBarrierSchedule);
        FIELD(numHits, "number of hits");
        FIELD(isHit, "is hit ? per asset");
        FIELD(hitDate, "breach date");
        FIELD(monitorType, "monitor Type");
        FIELD(smoothing, "use smoothing ?");
        FIELD(lowSpread, "low strike spread");
        FIELD(highSpread, "high strike spread");
        FIELD_MAKE_OPTIONAL(hitDate);
        FIELD_MAKE_OPTIONAL(smoothing);
        FIELD_MAKE_OPTIONAL(lowSpread);
        FIELD_MAKE_OPTIONAL(highSpread);
        FIELD(bendMaker, "barrier bend maker");
        FIELD_MAKE_OPTIONAL(bendMaker);
        FIELD(levels, "barrier schedule");
        FIELD(economicLevels, "Economic barrier schedule");
        FIELD_MAKE_OPTIONAL(economicLevels);
        FIELD(isOut, "True=>Knock Out; False=> In");
        FIELD(isUp, "breach above ?");
        FIELD(isSimultaneousBreach, "Count only simultaneous hits");
        FIELD_MAKE_OPTIONAL(isSimultaneousBreach);
        FIELD(closedFormDates, "Closed form dates");
        FIELD_MAKE_OPTIONAL(closedFormDates);
        FIELD(closedFormMethod, "Closed form methodology");
        FIELD_MAKE_OPTIONAL(closedFormMethod);
        FIELD(maxNumDivsPerYear, "Maximum number of dividends per year for BB");
        FIELD_MAKE_OPTIONAL(maxNumDivsPerYear);
        // Transient
        FIELD(numAssets, "numAssets");
        FIELD(smoothHitValues, "smoothHitValues");
        FIELD_MAKE_TRANSIENT(smoothHitValues);
        FIELD(isHitPerAsset, "snap shot per iteration");
        FIELD(isHitPerAssetSoFar, "snap shot per iteration/hist");
        FIELD(useAdjBarrier, "use adjusted barrier");
        FIELD_MAKE_TRANSIENT(numAssets);
        FIELD_MAKE_TRANSIENT(isHitPerAsset);
        FIELD_MAKE_TRANSIENT(isHitPerAssetSoFar);
        FIELD_MAKE_TRANSIENT(useAdjBarrier);
        FIELD(pastHitTime, "hit time if in past, else undefined");
        FIELD_MAKE_TRANSIENT(pastHitTime);
        FIELD(numRealHitsSoFar, "numRealHitsSoFar");
        FIELD_MAKE_TRANSIENT(numRealHitsSoFar);
        FIELD(isConditionMetSoFar, "isConditionMetSoFar");
        FIELD_MAKE_TRANSIENT(isConditionMetSoFar);
        FIELD(globalPastHitTime, "globalPastHitTime");
        FIELD_MAKE_TRANSIENT(globalPastHitTime);
        // Transient fields for state variables
        FIELD(amendedIsHit, "amended value of input isHit");
        FIELD_MAKE_TRANSIENT(amendedIsHit);
        FIELD(doneAmendedIsHit, "Whether isHitamended has been set");
        FIELD_MAKE_TRANSIENT(doneAmendedIsHit);
        FIELD(barrierDates, "dates for monitoring");
        FIELD_MAKE_TRANSIENT(barrierDates);
        FIELD(interpBarrier, "interpolated barrier");
        FIELD_MAKE_TRANSIENT(interpBarrier);
        FIELD(origDate, "value date when interpolated barrier was built");
        FIELD_MAKE_TRANSIENT(origDate);
        FIELD(adjBarLevels, "adjusted levels per asset");
        FIELD_MAKE_TRANSIENT(adjBarLevels);
        FIELD_NO_DESC(areMet);
        FIELD_MAKE_TRANSIENT(areMet);
    }

    static IObject* defaultBarrierSchedule(){
        return new BarrierSchedule();
    }
};

CClassConstSP const BarrierSchedule::TYPE = CClass::registerClassLoadMethod(
    "BarrierSchedule", typeid(BarrierSchedule), BarrierScheduleHelper::load);


/******************************************************************/
// Viewer class to BarrierSchedule to expose BarrierBending
/******************************************************************/

class BarrierScheduleBendView : public BarrierSchedule {
public:
    static CClassConstSP const TYPE;

    BarrierScheduleBendView() : BarrierSchedule(TYPE){};
    ~BarrierScheduleBendView(){};

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(BarrierScheduleBendView, clazz);
        SUPERCLASS(BarrierSchedule);
        EMPTY_SHELL_METHOD(defaultBarrierScheduleBendView);
    }

    static IObject* defaultBarrierScheduleBendView(){
        return new BarrierScheduleBendView();
    }

};

CClassConstSP const BarrierScheduleBendView::TYPE = CClass::registerClassLoadMethod(
    "BarrierScheduleBendView", typeid(BarrierScheduleBendView), BarrierScheduleBendView::load);


/******************************************************************/
// DoubleBarrier class
/******************************************************************/

DoubleBarrier::DoubleBarrier()
: Barrier(TYPE), barrierA(0), barrierB(0), isAnd(false), isSimultaneousBreach(false)
{}

DoubleBarrier::DoubleBarrier(bool           isAnd,
                             bool           isSimultaneousBreach,
                             BarrierSP      barrierA,
                             BarrierSP      barrierB)
:  Barrier(TYPE), barrierA(barrierA), barrierB(barrierB),
isAnd(isAnd), isSimultaneousBreach(isSimultaneousBreach)
{}

DoubleBarrier::~DoubleBarrier()
{}

SVGenBarrierHVSP DoubleBarrier::convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset) const
{
    throw ModelException("DoubleBarrier::convertBarrier", "Not implemented");
}

SVGenBarrierHVTSP DoubleBarrier::convertBarrier(const DateTimeArray& monitoringDates,
                                            const DateTime& smoothDate,
                                            const DateTime& valueDate,
                                            IRefLevel::IStateVarGenSP refLevelGen,
                                            const IMultiFactors* mAsset, int iDummy) const
{ 
    throw ModelException("DoubleBarrier::convertBarrier", "Not implemented"); 
}

void DoubleBarrier::validatePop2Object()
{
    static const string method = "DoubleBarrier::validatePop2Object";

    if( !barrierA || !barrierB )
        throw ModelException(method, "barrierA and/or barrierB is missing");

    // only support BarrierSchedules
    if( !BarrierSchedule::TYPE->isInstance(*barrierA) || !BarrierSchedule::TYPE->isInstance(*barrierB) )
        throw ModelException(method, "both barrierA and barrierB need to be barrier schedule");
}

void DoubleBarrier::validate(const MonteCarlo* model,
                            int nbAssets)
{
    barrierA->validate(model, nbAssets);
    barrierB->validate(model, nbAssets);
}

void DoubleBarrier::validate(int nbAssets) const
{
    barrierA->validate(nbAssets);
    barrierB->validate(nbAssets);
}

void DoubleBarrier::amendIsHit(const IMultiFactors* mFactors,
                              IRefLevel::IStateVarSP refLevel,
                              const DateTimeArrayArray& monitoringDates,
                              const DateTime& today)
{
    throw ModelException("DoubleBarrier::amendIsHit", "Not implemented");
}

void DoubleBarrier::createInterpBarrier(const DateTime& valueDate,
                                        const DateTimeArray& samples) const
{
    barrierA->createInterpBarrier(valueDate, samples);
    barrierB->createInterpBarrier(valueDate, samples);
}

void DoubleBarrier::overrideEventBarrier(const DateTime& valueDate) const
{
    barrierA->overrideEventBarrier(valueDate);
    barrierB->overrideEventBarrier(valueDate);
}

// Allow barrier to influence mon dates in inst if Daily
DateTimeArray DoubleBarrier::getFutureMonitoringDates(const DateTime&         today,
                                                      const DateTimeArray&    instMonitorDates,
                                                      smartPtr<IMCPathConfig> pathConfig) {
    static const string method = "DoubleBarrier::getFutureMonitoringDates";

    DateTimeArray d1 = barrierA->getFutureMonitoringDates(today, instMonitorDates, pathConfig);
    DateTimeArray d2 = barrierB->getFutureMonitoringDates(today, instMonitorDates, pathConfig);

    DateTimeArray newMonitorDates(DateTime::merge(d1, d2));

    return newMonitorDates;
}

void DoubleBarrier::pathUpdated(const IMCPathGenerator* pathGen)
{
    if( !!iBarrier )
    {
        iBarrier->pathUpdated(pathGen); // will propagate to iBarrierA and iBarrierB
    }
}

bool DoubleBarrier::createBarrierUtil(const DateTime& valueDate) const
{
    static const string method = "DoubleBarrier::createBarrierUtil";

    bool statA = barrierA->createBarrierUtil(valueDate);
    bool statB = barrierB->createBarrierUtil(valueDate);
    if( !statA || !statB ) return false;

    iBarrierA = barrierA->getBarrierUtil();
    iBarrierB = barrierB->getBarrierUtil();
    IBarrierUtil::ICanHaveSticky* sa = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(iBarrierA.get());
    IBarrierUtil::ICanHaveSticky* sb = dynamic_cast<IBarrierUtil::ICanHaveSticky*>(iBarrierB.get());
    if( !sa || !sb || !sa->isSticky() || !sb->isSticky() )
        throw ModelException(method, "Both barrierA/B need to be sticky");

    // package into one IBarrier.
    vector<IBarrierUtilSP> ibs(2);
    ibs[0] = iBarrierA;
    ibs[1] = iBarrierB;

    // meaning of isAnd is as follows:
    //      if A and B same type, one KI/O or need both to KI/KO
    //      if A and B are KI and KO respectively, false=KI_CANCEL_KO, true=KI_KEEP_KO
    bool aIsIn = !dynamic_cast<BarrierSchedule*>(barrierA.get())->isKO();
    bool bIsIn = !dynamic_cast<BarrierSchedule*>(barrierB.get())->isKO();
    if( aIsIn == bIsIn )
    {
        int numHits = (isAnd==aIsIn)?2:1;
        iBarrier = BarrierFactory::makeCoupler(ibs, numHits,aIsIn,false);
        if( !isSimultaneousBreach/*sticky*/ )
            iBarrier = BarrierFactory::makeSticky(iBarrier);
    }
    else
    {
        bool kiCancelKo = !isAnd;
        iBarrier = BarrierFactory::makeCoupler(ibs, 2,false,false);
        if( kiCancelKo/*itself sticky*/ )
            iBarrier = BarrierFactory::makeSticky(iBarrier);
    }

    return true;
}

IBarrierUtilSP DoubleBarrier::getBarrierUtil() const
{
    return iBarrier;
}

ScheduleSP DoubleBarrier::getBarrierSchedule() const
{
    throw ModelException("DoubleBarrier::getBarrierSchedule", "Function not well defined. Not implemented");
}

double DoubleBarrier::getInterpLevel(const IMCPathGenerator* pathGen,
                                  const IMCProduct* product,
                                  int iAsset) const
{
    throw ModelException("DoubleBarrier::getInterpLevel", "Function not well defined. Not implemented");
}


void DoubleBarrier::adjustLN(const IMCPathGenerator* pathGen,
                          const IMCProduct*        product,
                          const IMultiFactors*    assets,
                          int                     pathIdx)
{
    barrierA->adjustLN(pathGen, product, assets, pathIdx);
    barrierB->adjustLN(pathGen, product, assets, pathIdx);
}

void DoubleBarrier::adjustImplied(const IMCPathGenerator*  pathGen,
                               const IMCProduct*         product,
                               const IMultiFactors*     assets,
                               const YieldCurve*        discount,
                               InstrumentSettlement*    instSettle,
                               int                      pathIdx)
{
    barrierA->adjustImplied(pathGen, product, assets, discount, instSettle, pathIdx);
    barrierB->adjustImplied(pathGen, product, assets, discount, instSettle, pathIdx);
}


double DoubleBarrier::hitValue(const  IMCPathGenerator*  pathGen,
                            int    assetIdx)
{
    throw ModelException("DoubleBarrier::hitValue", "Function not well defined. Not implemented");
}

double DoubleBarrier::hitValue(const  IMCPathGenerator*  pathGen)
{
    throw ModelException("DoubleBarrier::hitValue", "Not implemented");
}

void DoubleBarrier::hitValueAndTime(const   IMCPathGenerator*  pathGen,
                                 int     assetIdx,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep)
{
    throw ModelException("DoubleBarrier::hitValueAndTime", "Function not well defined. Not implemented");
}

void DoubleBarrier::hitValueAndTime(const  IMCPathGenerator*  pathGen,
                                 double& value,
                                 bool&   isConditionMet,
                                 int&    metAtStep)
{
    throw ModelException("DoubleBarrier::hitValueAndTime", "Not implemented");
}

void DoubleBarrier::hitValueAndTimeAsBasket(const  IMCPathGenerator*  pathGen,
                                         const  IAggregateMakerSP& bsk,
                                         double& value,
                                         bool&   isConditionMet,
                                         int&    metAtStep)
{
    throw ModelException("DoubleBarrier::hitValueAndTimeAsBasket", "Not implemented");
}

void DoubleBarrier::hitValueAndTimeAtStep(const IMCPathGenerator*  pathGen,
                                       int                      iStep,
                                       int                      iAsset,
                                       double                   refLevel,
                                       double&                  value,
                                       bool&                    isConditionMet,
                                       int&                     metAtStep)
{
    throw ModelException("DoubleBarrier::hitValueAndTimeAtStep", "Function not well defined. Not implemented");
}

void DoubleBarrier::hitValueAndTimeAtMat(const IMCPathGenerator*  pathGen,
                                      int                      iStep,
                                      int                      iAsset,
                                      double                   refLevel,
                                      bool                     isConditionMetSoFar,
                                      double&                  value,
                                      bool&                    isConditionMet,
                                      int&                     metAtStep)
{
    throw ModelException("DoubleBarrier::hitValueAndTimeAtMat", "Function not well defined. Not implemented");
}

void DoubleBarrier::multiHelper(DoubleArray&       hitValues,
                             const BoolArray&   isConditionMet,
                             const IntArray&    metAtStep,
                             double&            overallHitValue,
                             bool&              overallIsConditionMet,
                             int&               overallMetAtStep)
{
    throw ModelException("DoubleBarrier::multiHelper", "Not implemented");
}

void DoubleBarrier::multiHelper2(DoubleArray&       hitValues,
                              const BoolArray&   isConditionMet,
                              double&            overallHitValue,
                              double&            overallBinaryValue)
{
    throw ModelException("DoubleBarrier::multiHelpe2", "Not implemented");
}

/*
BarrierLevelArraySP DoubleBarrier::reportLevelsAsBasket(const  IAggregateMakerSP& bsk,
                                                        const DateTime& valueDate,
                                                        const DoubleArray assetPerf,
                                                        double          refLevel,
                                                        int             assetIdx) const {
    // bsk must not be same....
    throw ModelException("DoubleBarrier::::reportLevelsAsBasket", 
                         "Not implemented for DoubleBarrier");    
}*/

BarrierLevelArraySP DoubleBarrier::reportLevels(const DateTime& valueDate,
                                             double          refLevel,
                                             int             assetIdx) const
{
    BarrierLevelArraySP blaA = barrierA->reportLevels(valueDate, refLevel, assetIdx);
    BarrierLevelArraySP blaB = barrierB->reportLevels(valueDate, refLevel, assetIdx);

    // need to push blaB into blaA. no simple append function for whole array
    // this function does NOT order the dates when combining since in any case
    // we may have 2 levels on the same asset on the same date!
    for(int i=0; i<blaB->size(); i++)
        blaA->push_back((*blaB)[i]);

    return blaA;
}

void DoubleBarrier::getEvents(const  IMCPathGenerator*  pathGen,
                                EventResults* events,
                                bool useIsHitFlags,
                                const string barrierName,
                                StringArray* assetNames) {
    throw ModelException("DoubleBarrier::::getEvents",
        "Event handling not implemented for DoubleBarrier");
}

void DoubleBarrier::getEvents(EventResults* events,
                                const string barrierName,
                                int   thisPeriodIdx,
                                const IntArrayArray& metAtStep,
                                const StringArray& assNames,
                                const DoubleMatrix& assLevels,
                                const DoubleMatrix& refLevels) {
    throw ModelException("DoubleBarrier::::getEvents",
        "Event handling not implemented for DoubleBarrier");
}

bool DoubleBarrier::isMonitoringPoint(const DateTimeArray& monitoringDates,
                               const DateTime& today) {
    throw ModelException("DoubleBarrier::::isMonitoringPoint",
        "Not implemented for DoubleBarrier");
}

//// roll through time replacing risk barrier in past by legal one
bool DoubleBarrier::sensShift(Theta* theta) {
    barrierA->sensShift(theta);
    barrierB->sensShift(theta);
    return true; // continue shifting
}


// Set DoubleBarrier to look like termsheet
bool DoubleBarrier::sensShift(LegalTerms* shift) {
    barrierA->sensShift(shift);
    barrierB->sensShift(shift);
    return true; // continue shifting
}

void DoubleBarrier::hitValueAndTime(const IMCPathGenerator*  pathGen,
                                 int        beginStep,
                                 int        endStep,
                                 const double *bendAmts,
                                 double&    value,
                                 bool&      isConditionMet,
                                 int&       metAtStep,
                                 int&       numRealHits,
                                 CBoolArray *isHitPerAsset,
                                 bool&      isConditionMetSoFar,
                                 CBoolArray *isHitPerAssetSoFar,
                                 int&       numRealHitsSoFar,
                                 int&       globalPastHitTime)
{
    throw ModelException("DoubleBarrier::hitValueAndTime", "Not implemented");
}


class DoubleBarrierHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DoubleBarrier, clazz);
        SUPERCLASS(Barrier);
        EMPTY_SHELL_METHOD(defaultDoubleBarrier);
        FIELD(barrierA, "Barrier A");
        FIELD(barrierB, "Barrier B");
        FIELD(isAnd,    "True if need both barrier. In case of 1 KI 1 KO, true=KI_KEEP_KO, false=KI_CANCEL_KO");
        FIELD(isSimultaneousBreach,    "Only if isAnd=true. Not relevant in case of 1 KI and 1 KO");
        FIELD_MAKE_OPTIONAL(isSimultaneousBreach);
    }

    static IObject* defaultDoubleBarrier(){
        return new DoubleBarrier();
    }
};

CClassConstSP const DoubleBarrier::TYPE = CClass::registerClassLoadMethod(
    "DoubleBarrier", typeid(DoubleBarrier), DoubleBarrierHelper::load);



/******************************************************************/
// BarrierUnion - Wrapper stuff for IMS
/******************************************************************/

#define BARRIER_TYPE_SCHEDULE       "BarrierSchedule"
#define BARRIER_TYPE_DOUBLE         "DoubleBarrier"

// validation
void BarrierUnion::validatePop2Object(){
    static const string routine = "BarrierUnion::validatePop2Object";

    if (barrierType.empty())
    {
        throw ModelException(routine, "Blank barrier type specified!");
    }
    if (barrierType==BARRIER_TYPE_SCHEDULE)
    {
        realBarrier = barrierSchedule;
    }
    else if (barrierType==BARRIER_TYPE_DOUBLE)
    {
        realBarrier = doubleBarrier;
    }
    else
    {
        throw ModelException(routine, "Unrecognised Barrier Type " + barrierType +
                             ". Expected " +  BARRIER_TYPE_SCHEDULE + "  or " + BARRIER_TYPE_DOUBLE);
    }

    if( !realBarrier )
        throw ModelException(routine, "Expected " + barrierType + " but none supplied!");
}

// for reflection
BarrierUnion::BarrierUnion(): CObject(TYPE){}

BarrierUnion::BarrierUnion(CClassConstSP clazz) : CObject(clazz){}

//// must be in source file
BarrierUnion::~BarrierUnion(){}

class BarrierUnionHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BarrierUnion, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBarrierUnion);
        FIELD(barrierType, "Only BarrierSchedule for now");
        FIELD(barrierSchedule,  "barrierSchedule");
        FIELD_MAKE_OPTIONAL(barrierSchedule);
        FIELD(doubleBarrier,    "doubleBarrier");
        FIELD_MAKE_OPTIONAL(doubleBarrier);
        FIELD(realBarrier, "realBarrier");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(realBarrier);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultBarrierUnion(){
        return new BarrierUnion();
    }
};

CClassConstSP const BarrierUnion::TYPE = CClass::registerClassLoadMethod(
    "BarrierUnion", typeid(BarrierUnion), BarrierUnionHelper::load);

/******************************************************************/
// Wrapper to be able for alternative view of barrier union
/******************************************************************/

class BarrierUnion2 : public BarrierUnion {
public:
    static CClassConstSP const TYPE;

    BarrierUnion2() : BarrierUnion(TYPE){};
    ~BarrierUnion2(){};

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(BarrierUnion2, clazz);
        SUPERCLASS(BarrierUnion);
        EMPTY_SHELL_METHOD(defaultBarrierUnion2);
    }

    static IObject* defaultBarrierUnion2(){
        return new BarrierUnion2();
    }

};

CClassConstSP const BarrierUnion2::TYPE = CClass::registerClassLoadMethod(
    "BarrierUnion2", typeid(BarrierUnion2), BarrierUnion2::load);

DRLIB_END_NAMESPACE

