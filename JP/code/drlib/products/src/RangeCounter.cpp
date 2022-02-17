//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RangeCounter.cpp
//
//   Description : A basket of underlyings is observed. For each underlying, if
//                 its current performance (the value of the underlying in
//                 proportion to a reference value) on particular monitoring
//                 dates has moved out of a specified range, the underlying is
//                 said to have hit.
//
//                 At maturity, a payoff between the minimum redemption and the
//                 maximum redemption will be given.
//                 If penaliseStayedInRange is true, a penalty will be
//                 subtracted from the maximum redemption for every underlying
//                 that DID NOT hit (stayed within the range).
//                 Otherwise, a penalty will be subtracted from the maximum
//                 redemption for every underlying that DID hit (went outside
//                 the range).
//                 By default penaliseStayedWithinRange is false.
//
//                 The price monitoring is not continuous - the price available
//                 from the model is the price at the given monitoring point,
//                 not a high-low range required for continuous monitoring.
//
//   Date        : May 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Barrier.hpp" // for BarrierLevel objects in events
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/ObservationBuilder.hpp"


DRLIB_BEGIN_NAMESPACE

class RangeCounter: public GenericNFBase,
                virtual public IMCIntoProduct,
                virtual public LegalTerms::Shift,
                virtual public BarrierBreach::IEventHandler {
    class Range: public CObject {
    public:
        // Reflection mechanism
        static CClassConstSP const TYPE;
        void validatePop2Object() {
            static const string method = "Range::validatePop2Object";
           
            try {
                // don't allow empty ranges
                if( dates.empty() ) {
                    throw ModelException(method, "No dates provided for Range!");
                }

                // check consistency of descriptions of levels
                if( dates.size() != rangeMin.size() ) {
                    throw ModelException(method, "Size of dates ("+
                                         Format::toString(dates.size())+
                                         ") did not match size of rangeMin ("+
                                         Format::toString(rangeMin.size())+
                                         ")");
                }
                if( dates.size() != rangeMax.size() ) {
                    throw ModelException(method, "Size of dates ("+
                                         Format::toString(dates.size())+
                                         ") did not match size of rangeMax ("+
                                         Format::toString(rangeMax.size())+
                                         ")");
                }

                // check ranges
                // NB all three arrays have been established as same size
                for (int i=0; i<dates.size(); i++)
                {
                    if (!Maths::isPositive(rangeMax[i] - rangeMin[i])) {
                        throw ModelException(method, "On date '" +
                                    dates[i].toString() +
                                    "', rangeMax ("+
                                     Format::toString(rangeMax[i])+
                                     ") is less than or equal to rangeMin ("+
                                     Format::toString(rangeMin[i])+")");
                    }
                }

                // Construct the Schedule objects that represent the ranges
                minSched = ScheduleSP(new Schedule( dates, rangeMin, interp ));
                maxSched = ScheduleSP(new Schedule( dates, rangeMax, interp ));
            }
            catch (exception& e) {
                throw ModelException(&e, method);
            }
        }

        ~Range() {}

        // Range class interface
        // Expose common date information
        DateTime lastDate() const {
            return dates.back();
        }
        DateTime firstDate() const {
            return dates.front();
        }
        bool coversDateRange(const DateTime& lowDate,
                             const DateTime& highDate,
                             bool datesDefineLimit) const {
            // The schedules are constructed such that both cover the same dates
            // so we only need to check one Schedule.
            return minSched->coversDateRange(lowDate, highDate, datesDefineLimit);
        }

        // Interpolation functions for the range
        double interpolateMin(const DateTime& date) const {
            return minSched->interpolate(date);
        }
        double interpolateMax(const DateTime& date) const {
            return maxSched->interpolate(date);
        }

        // Expose internal arrays
        const DateTimeArray& getDates() const {
            return dates;
        }
        const DoubleArray& getRangeMin() const {
            return rangeMin;
        }
        const DoubleArray& getRangeMax() const {
            return rangeMax;
        }
        const string getInterp() const {
            return interp;
        }

        BarrierLevelArraySP constructBarrierLevels(const DateTime& valueDate,
                                                   double          scale,
                                                   bool isContinuous ) const
        {
            BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));

            // Report barrier levels over a date range
            DateTime endDate = BarrierLevel::barrierWindow(valueDate);

            // Construct a CashFlowArray with *relative* barrier values
            CashFlowArraySP minSet(minSched->subset(valueDate, endDate));
            CashFlowArraySP maxSet(maxSched->subset(valueDate, endDate));
            // WILL create two arrays of equal length

            // Create BarrierLevel array with scaled barrier values. This
            // is required for conversion between *relative* and *absolute*
            // measurements.
            for( int i=0; i<maxSet->size(); i++ ) {
                BarrierLevel max(true, // isUp
                                 (*maxSet)[i].date,
                                 (*maxSet)[i].amount * scale,
                                 isContinuous);
                BarrierLevel min(false, // isUp
                                 (*minSet)[i].date,
                                 (*minSet)[i].amount * scale,
                                 isContinuous);

                // Add new BarrierLevels to array
                reportLevels->push_back(max);
                reportLevels->push_back(min);
            }

            return reportLevels;
        }

    protected:
        // Fields
        DateTimeArray  dates;    // Range schedule dates
        DoubleArray    rangeMin; // Range min values
        DoubleArray    rangeMax; // Range max values
        string         interp;   // Interpolation scheme

        // Transient Fields
        ScheduleSP     minSched; // Schedule representing min
        ScheduleSP     maxSched; // Schedule representing max

    private:
        // Reflection mechanism
        Range() : CObject(TYPE) {}
        Range(const Range& rhs); // Copy constructor disabled
        Range& operator=(const Range& rhs); // Assignment disabled

        static IObject* defaultRange() {
            return new Range();
        }

        static void load(CClassSP& clazz) {
            clazz->setPublic(); // make visible to EAS/spreadsheet
            REGISTER(Range, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultRange);
            FIELD(dates, "schedule dates");
            FIELD(rangeMin, "min values");
            FIELD(rangeMax, "max values");
            FIELD(interp, "interpolation scheme");

            FIELD(minSched, "Schedule representing min" );
            FIELD_MAKE_TRANSIENT( minSched );
            FIELD(maxSched, "Schedule representing max" );
            FIELD_MAKE_TRANSIENT( maxSched );
        }

    };

    typedef smartConstPtr<Range> RangeConstSP;
    typedef smartPtr<Range> RangeSP;

protected:
    /// fields ////////
    IDateBuilderSP obsBuilder; // monitoring dates builder
    double minRedemption;        // minimum redemption on payoff
    double maxRedemption;        // maximum redemption on payoff
    double step;                 // penalty for each hit
    RangeSP range;               // scoring range
    RangeSP legalRange;          // scoring range for legal terms
    bool   penaliseStayedInRange;// penalty calculation -
                                 // true - step * underlyings within barriers
                                 // false - step * underlyings out of barriers
    /// transient fields ////////
    DateTimeArraySP monitoringDates;  // all dates after ref date to the end

public:
    static CClassConstSP const TYPE;
    friend class RangeCounterSVMC;

    // validation
    virtual void validatePop2Object(){
        static const string method = "RangeCounter::validatePop2Object";
        GenericNFBase::validatePop2Object();

        try
        {
            // check normal and legal ranges start/end dates
            if( range->firstDate() != legalRange->firstDate() ||
                range->lastDate() != legalRange->lastDate() ) {
                throw ModelException(method,
                                     "Start/end dates of range ("+
                                     range->firstDate().toString()+", "+
                                     range->lastDate().toString()+
                                     ") and legal range ("+
                                     legalRange->firstDate().toString()+", "+
                                     legalRange->lastDate().toString()+
                                     ") do not match." );
            }

            // check max redemption is above min redemption
            if (!Maths::isPositive(maxRedemption - minRedemption)) {
                throw ModelException(method, "Maximum redemption ("+
                                     Format::toString(maxRedemption)+
                                     ") should be above minimum redemption ("+
                                     Format::toString(minRedemption)+")");
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    virtual void GetMarket(const IModel*          model,
                           const CMarketDataSP    market) {
        static const string method = "RangeCounter::GetMarket";
        try {
            GenericNFBase::GetMarket(model,market);

            // Get the market for the IDateBuilder
            obsBuilder->getMarket(model, market.get());
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    virtual void Validate() {
        static const string method = "RangeCounter::Validate";
        try {
            // Date validation moved to Validate, to ensure that obsBuilder has been completely
            // built by GetMarket call. This means that we can't do validation until after
            // the market data has been populated.

            GenericNFBase::Validate();

            for( int iAsset=0; iAsset<obsMap->getSampleDates().size(); iAsset++ ) {
                const DateTimeArray& sampleDates = obsMap->getSampleDates()[iAsset];

                // validate dates are not empty
                if (sampleDates.empty()) {
                    throw ModelException(method, "No monitoring dates for "+
                                         assets->getName(iAsset));
                }

                // need to make sure the ref level averaging in has finished by
                // the time we hit monitoring dates
                DateTimeArray overlap = refLevel->getFutureDates(sampleDates[0]);
                if (!overlap.empty()) {
                    throw ModelException(method,
                                         "Monitoring dates for "+
                                         assets->getName(iAsset)+
                                         " cannot start ("+
                                         sampleDates[0].toString()+
                                         ") until the last reference level date ("+
                                         overlap.back().toString()+
                                         ") has passed");
                }

                // Check whether the ranges cover the sample dates used for each asset.
                const DateTime& first = sampleDates.front();
                const DateTime& last = sampleDates.back();

                // check barrier range dates against sampling dates

                // Compare last date with ranges
                if( last > range->lastDate() ) {
                    // XXX For the moment, fail validation, but later we recover if we can
                    throw ModelException(method, "Last sampling date ("+
                                         last.toString()+
                                         ") for "+assets->getName(iAsset)+
                                         " is after risk barrier range.");
                }

                // Compare first date with ranges
                if( first < range->firstDate() ) {
                    // XXX For the moment, fail validation, but later we recover if we can
                    throw ModelException(method, "First sampling date ("+
                                         first.toString()+
                                         ") for "+assets->getName(iAsset)+
                                         " is before risk barrier range.");
                }
                // XXX We may have to tweak a flat linear or step range, and validate against
                // a sloping range.

                if( range->getInterp() == Schedule::INTERP_NONE ) {
                    // Check that sample dates is a subset of range dates, and complain if
                    // it's not
                    const DateTimeArray& rangeDates = range->getDates();

                    int sampleIdx=0;
                    for(int rangeIdx=0; rangeIdx<rangeDates.size() && sampleIdx<sampleDates.size(); rangeIdx++) {
                        if (rangeDates[rangeIdx]==sampleDates[sampleIdx]) {
                            sampleIdx++;
                        } else if (rangeDates[rangeIdx]>sampleDates[sampleIdx]) {
                            throw ModelException(method, "Sample date "+
                                            sampleDates[sampleIdx].toString()+
                                            " for "+assets->getName(iAsset)+
                                            "is missing from risk range.");
                        }
                    }
                    // Sample dates beyond the last range date has already been checked above
                }

                // Compare last date with ranges
                if( last > legalRange->lastDate() ) {
                    // XXX For the moment, fail validation, but later we recover if we can
                    throw ModelException(method, "Last sampling date ("+
                                         last.toString()+
                                         ") for "+assets->getName(iAsset)+
                                         " is after legal barrier range.");
                }
                // Compare first date with ranges
                if( first < legalRange->firstDate() ) {
                    // XXX For the moment, fail validation, but later we recover if we can
                    throw ModelException(method, "First sampling date ("+
                                         first.toString()+
                                         ") for "+assets->getName(iAsset)+
                                         " is before legal barrier range.");
                }
                // XXX We may have to tweak a flat linear or step range, and validate against
                // a sloping range.

                if( legalRange->getInterp() == Schedule::INTERP_NONE ) {
                    // Check that sample dates is a subset of range dates, and complain if
                    // it's not
                    const DateTimeArray& rangeDates = legalRange->getDates();

                    int sampleIdx=0;
                    for(int rangeIdx=0; rangeIdx<rangeDates.size() && sampleIdx<sampleDates.size(); rangeIdx++) {
                        if (rangeDates[rangeIdx]==sampleDates[sampleIdx]) {
                            sampleIdx++;
                        } else if (rangeDates[rangeIdx]>sampleDates[sampleIdx]) {
                            throw ModelException(method, "Sample date "+
                                            sampleDates[sampleIdx].toString()+
                                            " for "+assets->getName(iAsset)+
                                            " is missing from legal range.");
                        }
                    }
                    // Sample dates beyond the last range date has already been checked above
                }
            }

        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    // XXX Be careful - these are deliberately the non-ISDAfied values.
    const DateTimeArray samplingDates() const {
        // this is required for GenericNFBase::Validate to construct the obsMap.
        return *(obsBuilder->dates().get());
    }

    // set barriers to be the economic (legal) ones
    bool sensShift(LegalTerms* shift) {
        // just replace all barriers with the corresponding economic one
        range = legalRange;

        return true; // continue shifting
    }

    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach, IModel* model,
                   const DateTime& eventDate, EventResults* events) const {
        static const string method = "RangeCounter::getEvents";

        try {
            MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
            if (mc) {
                auto_ptr<IMCProduct> prod(createProduct(mc));
                MCPathGeneratorSP pastPathGenerator(
                        mc->getPathConfig()->pastPathGenerator(prod.get()));
                // Tell the product that the generator has changed
                // Do that even if there is no past so that the product gets
                // some state variables e.g. refLevel
                prod->pathGenUpdated(pastPathGenerator.get());
                pastPathGenerator->generatePath(0); // may well do nothing
                prod->getEvents(pastPathGenerator.get(), events, eventDate);
            } else {
                throw ModelException(method,
                        "Internal error - expected Monte Carlo model for RangeCounter pricing");
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    RangeCounter(): GenericNFBase(TYPE), penaliseStayedInRange(false) {}
    RangeCounter(const RangeCounter& rhs);     // not implemented
    RangeCounter& operator=(const RangeCounter& rhs); // not implemented

    static IObject* defaultRangeCounter(){
        return new RangeCounter();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RangeCounter, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultRangeCounter);
        FIELD(obsBuilder, "monitoring dates builder");
        FIELD(minRedemption, "minimum redemption");
        FIELD(maxRedemption, "maximum redemption");
        FIELD(step, "penalty per hit");
        FIELD(range, "range" );
        FIELD(legalRange, "range for legal terms" );
        FIELD(penaliseStayedInRange,
                     "penalise for being within range for duration" );
        FIELD_MAKE_OPTIONAL(penaliseStayedInRange);
        FIELD(monitoringDates, "all dates on which levels are monitored");
        FIELD_MAKE_TRANSIENT(monitoringDates);
    }
};

/* MC product class for RangeCounterSV */
class RangeCounterSVMC: public MCProductClient,
                    virtual public IMCProductLN,
                    virtual public IMCProductImplied {
private:
    const RangeCounter*       inst;         // Instrument
    int                       nbAssets;     // convenient

    // maintain state of the instrument in the past
    BoolArray         hasHitInHistory; // Flag per asset indicating hit status
    int               hitsInHistory;   // Convenienct - count of hits in history
    DoubleArray       histRefLevels; // Reference level price, used for BARRIER_LEVEL

    BoolArray         hasHit;         // saves on memory alloc in payoff
    // can store the discrete payoff values
    DoubleArray       payoffPerHitCount;

    // precalculated barrier levels per asset - more efficient than interpolating values
    // on every MC run
    DoubleArrayArray rangeMin;
    DoubleArrayArray rangeMax;
    DoubleArrayArray legalRangeMin;
    DoubleArrayArray legalRangeMax;

    // State variables and generators
    SVGenSpotSP               spotGen;      // Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  // Generator for ref level
    SVGenDiscFactorSP         dfGen;        // Generator for discount factors
    SVGenSpot::IStateVarSP    spotSV;       // Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   // Ref level state variable
    SVDiscFactorSP            dfSV;         // Df state variable

    // To allow caching of product-level info across greeks
    class PricesRangeCounter: public MCPricesSimple{
        // cached params
        vector<BoolArraySP> cachedHasHit;  // [path idx][asset idx]
        IntArraySP          cachedNumHits; // [path idx]

        // working params
        bool                isValid;
        IntArray            changedAssets;
        bool                doingGreek;

    protected:
        IMCPrices* emptyConstructor() const{
            return new PricesRangeCounter(false, 1, 1, 1);
        }

    public:
        //// the assets whose path has changed
        const IntArray& getChangedAssets(){
            return changedAssets;
        }
        //// returns the cached data
        void cacheRead(int        pathIdx,
                       BoolArray& hasHit,
                       int&       hits) {
            hasHit = *(cachedHasHit[pathIdx]);
            hits = (*cachedNumHits)[pathIdx];
        }
        //// can we read from the cache?
        bool cacheValid(){
            return (doingGreek && isValid);
        }
        //// can we write to the cache?
        bool cacheUpdateAllowed(){
            return (!doingGreek && isValid);
        }
        //// save the supplied basketParts
        void cacheWrite(int                pathIdx,
                        const BoolArray&   hasHit,
                        int                hits){
            *(cachedHasHit[pathIdx]) = hasHit;
            (*cachedNumHits)[pathIdx] = hits;
        }
        //// number of bytes used per path in current configuration
        virtual int storagePerPath(IMCProduct* product) const {
            int n = MCPricesSimple::storagePerPath(product);
            if (isValid){
                // usually overestimates the cache size since the BoolArray
                // only uses an int on a few platforms
                n += sizeof(int) * product->getNumAssets();
                n += sizeof(int);
            }
            return n;
        }

        /** Returns a deep copy of this object */
        IMCPrices* clone() const{
            PricesRangeCounter& copy = dynamic_cast<PricesRangeCounter&>(*MCPricesSimple::clone());
            copy.cachedHasHit = cachedHasHit; // shallow copy
            copy.cachedNumHits = cachedNumHits; // shallow copy
            copy.isValid = isValid;
            copy.doingGreek = doingGreek;
            return &copy;
        }
        PricesRangeCounter(bool              useCache,
                           int               numAssets,
                           int               nbIter,
                           int               nbSubSamples):
            MCPricesSimple(nbIter, nbSubSamples),
            isValid(false),
            changedAssets(numAssets), doingGreek(false){
            if (useCache){
                cachedHasHit = vector<BoolArraySP>(nbIter);
                for (int i = 0; i < nbIter; i++){
                    cachedHasHit[i] = BoolArraySP(new BoolArray(numAssets));
                }
                cachedNumHits = IntArraySP(new IntArray(nbIter));
                isValid = true;
            }
            // fill changedAssets with 0, 1, 2, ...
            for (int i = 0; i < numAssets; i++){
                changedAssets[i] = i;
            }
        }
        /** invoked before each greek calculation (but not before initial
            pricing run) */
        virtual void configureCache(const IntArray& changedAssets){
            if (isValid){
                this->changedAssets = changedAssets;
            }
            doingGreek = true;
        }

        virtual ~PricesRangeCounter() {}

    };

public:

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "RangeCounterSVMC::collectStateVars";
        try{
            svCollector->append(spotGen.get());             // spot level
            svCollector->append(refLevelGen.get());         // reference level
            svCollector->append(dfGen.get());               // and a DiscFactor one
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, when the future path generator is
        created and also when the past path generator is used by getEvents()) */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        static const string routine = "RangeCounterSVMC::pathGenUpdated";
        try{
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    RangeCounterSVMC(const RangeCounter*      inst,
                 const SimSeriesSP&   simSeries):
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        hasHitInHistory(nbAssets, false),
        hitsInHistory(0),
        histRefLevels(nbAssets, 0.0),
        hasHit(nbAssets, false),
        payoffPerHitCount(nbAssets+1, 0.0),
        rangeMin(nbAssets),
        rangeMax(nbAssets),
        legalRangeMin(nbAssets),
        legalRangeMax(nbAssets),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())) {
            // Precalculate interpolations to improve MC performance
            // Use the simSeries to get DateTimeArrays for each asset, to ensure that the
            // times used for the MCPath are consistent with the times used to calculate
            // the barrier levels.
            //
            // Use inst->getSampleDates() to get a per asset list of sample dates,
            // and construct the interpolations from that. This way, if the sampling
            // convention has shifted the sample date, we use the barrier for that date,
            // not the original sample date.
            // We don't have to worry about the actual model date - that should have been
            // used by the MC engine when constructing the paths. We always work in
            // sample date space.
            const DateTimeCluster& sampleDates = inst->obsMap->getSampleDates();
            for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                int iSampleSize = sampleDates[iAsset].size();
                rangeMin[iAsset].resize(iSampleSize);
                rangeMax[iAsset].resize(iSampleSize);
                legalRangeMin[iAsset].resize(iSampleSize);
                legalRangeMax[iAsset].resize(iSampleSize);

                for( int iSample=0; iSample<iSampleSize; iSample++ ) {
                    const DateTime& date = sampleDates[iAsset][iSample];
                    rangeMin[iAsset][iSample] = inst->range->interpolateMin(date);
                    rangeMax[iAsset][iSample] = inst->range->interpolateMax(date);
                    legalRangeMin[iAsset][iSample] = inst->legalRange->interpolateMin(date);
                    legalRangeMax[iAsset][iSample] = inst->legalRange->interpolateMax(date);
                }
            }

            for(int h=0; h<=nbAssets; h++) {
                // Calculate the hit penalty
                // For penaliseStayedInRange,
                //     penalty is step * number of assets that didn't hit
                // For !penaliseStayedInRange,
                //     penalty is step * number of assets that hit
                double penalty = inst->penaliseStayedInRange ?
                    (nbAssets - h) * inst->step : h * inst->step;
                payoffPerHitCount[h] = inst->notional *
                    Maths::max( inst->minRedemption,
                                inst->maxRedemption - penalty );
            }
        }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&             prices) {
        int    hits;

        PricesRangeCounter& myPrices = static_cast<PricesRangeCounter&>(prices);
        if (myPrices.cacheValid()){
            // Initialise for cached ones
            myPrices.cacheRead(pathGen->getPathIndex(),
                               hasHit,
                               hits); // this will be all hits for 'pathIdx'
        } else {
            // init for non-cached case
            hits = hitsInHistory;
            hasHit = hasHitInHistory;
        }
        const IntArray& changedAssets = myPrices.getChangedAssets();

        for(int i=0; i<changedAssets.size(); i++) {
            int iAsset = changedAssets[i];
            // update hasHit and hits for non-cached values
            if( !hasHitInHistory[iAsset] ) {
                const SVPath& path = spotSV->path(iAsset);
                // local barrier pointing to legal or risk depending on past or not
                const DoubleArray& minBarr = doingPast() ? legalRangeMin[iAsset] : rangeMin[iAsset];
                const DoubleArray& maxBarr = doingPast() ? legalRangeMax[iAsset] : rangeMax[iAsset];
                double refLevel = refLevelSV->refLevel(iAsset);

                int hitIdx;
                bool isUp;

                if( evaluateRange(iAsset, path, minBarr, maxBarr, refLevel, hitIdx, isUp) ) {
                    // this asset does hit this time
                    if (!hasHit[iAsset]) {
                        // the count does not include a hit for this asset
                        hasHit[iAsset] = true;
                        ++hits;
                    }
                } else {
                    // this asset does NOT hit this time
                    if (hasHit[iAsset]) {
                        // the count includes a hit for this asset which must be removed
                        hasHit[iAsset] = false;
                        --hits;
                    }
                }
            } else {
                hasHit[iAsset] = true;
                // hits already has a count for this
            }
        }

        // preserve values for past
        if (doingPast()) {
            // This run is historical path, so store result for MC runs
            hasHitInHistory = hasHit;
            hitsInHistory = hits;

            // In order to satisfy the BARRIER_LEVEL request we need
            // to know the ref level for each asset. This is only
            // needed for already started case.
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = refLevelSV->refLevel(iAsset);
            }
        }

        if (!doingPast() || !hasFuture()) {
            if (myPrices.cacheUpdateAllowed()){
                myPrices.cacheWrite(pathGen->getPathIndex(), hasHit, hits);
            }

            // Compute a payoff, but only when we have a "complete" situation :
            // either doingPast() and all is past, or !doingPast().
            // now scale by discount factor (notional already included)
            prices.add(payoffPerHitCount[hits] * dfSV->firstDF());

            if (doingPast() && !hasFuture()) {
                // If we're in a known state, we record known flows on
                // their known date (so no discounting).
                if (!paymentDate.empty()) {
                    knownCashFlows->addFlow(paymentDate,
                                            payoffPerHitCount[hits]);
                }
            }
        }
    }

    IMCPrices* createOrigPrices(int  nbIter,
                                        int  nbSubSamples,
                                        int  mode) {
        return new PricesRangeCounter((mode & CACHE_PRODUCT_BIT)? true: false, // use cache?
                                      nbAssets,
                                      nbIter,
                                      nbSubSamples);
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string routine = "RangeCounterSVMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }


    // any old level so that MC implied works
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "RangeCounterSVMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevelGen->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            double interpLevel  = 1.0;

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** BARRIER_LEVEL event generation.
     * Override IMCProduct::recordExtraOutput(). This is called once all the
     * MonteCarlo runs have been executed.
     * This is used in preference to implementing IHandlePaymentEvents, since
     * that already gives us PAYMENT_DATES and KNOWN_CASHFLOWS events for free
     * in IMCProduct (since this is a simple product with a single cashflow and
     * a single payment date). This way of doing things should be reviewed.
     **/
    void recordExtraOutput(CControl*     control,
                           Results*      results,
                           const IMCPrices& prices) const {
        // BARRIER_LEVEL ...
        OutputRequest*  request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        const DateTime& today = getToday();
        const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
        // Only try to satisfy this request when have some past
        if (request && (today >= lastRefDate)) {
            // This operates on a per-asset basis
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                // histRefLevels is used to report absolute barriers levels.
                // Only record events for underliers that haven't already hit
                // and we have a reference value to use.
                if (Maths::isPositive(histRefLevels[iAsset]) &&
                    !hasHitInHistory[iAsset] ) {
                    // XXX Report the ISDA adjusted levels
                    BarrierLevelArraySP levels =
                            inst->legalRange->constructBarrierLevels(
                                            today,
                                            histRefLevels[iAsset],
                                            false ); // NOT a continuous barrier
                    if (!levels->empty()) {
                        OutputRequestUtil::recordBarrierLevels(
                            control, results,
                            getMultiFactors()->assetGetTrueName(iAsset),
                            levels.get());
                    }
                }
            }
        }
    }

    // Evaluate a Path against the range levels, returning true if a hit occurs,
    // with the path index and breach direction in hitIdx and isUp.
    bool evaluateRange(const int iAsset, const SVPath& path,
        const DoubleArray& minBarr,
        const DoubleArray& maxBarr, const double refLevel,
        int& hitIdx, bool& isUp)
    {
        int    pathBegin = path.begin();
        int    pathEnd = path.end();

        if( pathEnd==pathBegin ) {
            // break early if the path is empty
            return false;
        }

        // Use the path defined per asset to accomodate different sampling
        // regimes used on each asset. This is based on the sampling dates, and
        // is converted to the actual modelling date by getModellingIndex.
        int    beginIdx = inst->obsMap->getFirstSampleIndex(iAsset, pathBegin);
        int    endIdx   = inst->obsMap->getLastSampleIndex(iAsset, pathEnd-1)+1;

        int iStep = beginIdx;
        while (iStep<endIdx) {
            // Use the barrier levels from the provided level arrays.
            // These must be consistent with the dates used to generate the path.
            double min = minBarr[iStep];
            double max = maxBarr[iStep];

            // Calculate the performance of underlying, making sure that we're using the
            // correct step for the sampling date.
            double perf = path[inst->obsMap->getModellingIndex(iAsset, iStep)] / refLevel;

            // Check barrier status
            bool isHit = Maths::isPositive( perf - max ) ||
                         Maths::isNegative( perf - min );

            if( isHit ) {
                // Populate details of the breach
                hitIdx = iStep;
                isUp = Maths::isPositive( perf - max );
                return true;
            }

            iStep++;
        }
        return false;
    }


    // Override IMCProduct::getEvents to support BarrierBreach event reporting.
    virtual void getEvents(const IMCPathGenerator*  pathGen,
                           EventResults* events,
                           const DateTime& eventDate) {
        int hits = 0;
        const DateTimeCluster& sampleCluster = inst->obsMap->getSampleDates();

        for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
            const DateTimeArray& sampleDates = sampleCluster[iAsset];
            const SVPath& path = spotSV->path(iAsset);
            // barrier breaches assessed against the legal terms only.
            const DoubleArray& minBarr = legalRangeMin[iAsset];
            const DoubleArray& maxBarr = legalRangeMax[iAsset];
            double refLevel = refLevelSV->refLevel(iAsset);

            int hitIdx;
            bool isUp;

            if( evaluateRange(iAsset, path, minBarr, maxBarr, refLevel, hitIdx, isUp) )
            {
                ++hits;

                // XXX We are issuing multiple events for the same date.
                StringArraySP assetNames(new StringArray(1));
                DoubleArraySP assetLevels(new DoubleArray(1));
                DoubleArraySP barrLevels(new DoubleArray(1));

                (*assetNames)[0] = getMultiFactors()->assetGetTrueName(iAsset);
                (*assetLevels)[0] = path[hitIdx];
                if( isUp ) {
                    (*barrLevels)[0] = maxBarr[hitIdx] * refLevel;
                } else {
                    (*barrLevels)[0] = minBarr[hitIdx] * refLevel;
                }

                // Each breach is reported seperately, since the breaches are
                // distinct events. They are reported as KNOCK_OUTs. The number
                // of unbreached underlyings is also reported.
                //
                // XXX ALERT ALERT - make sure that we use the appropriate date from
                // obsMap (has to be consistent with ISDA sample convention)
                // XXX Also be careful of movement of individual assets vs the whole. Can't
                // guarantee chronological order of output here now!!!!!!!
                events->addEvent( new BarrierBreach( sampleDates[hitIdx],
                                  "Range Counter Barrier",
                                  BarrierBreach::EUROPEAN,
                                  BarrierBreach::KNOCK_OUT,
                                  isUp, nbAssets - hits,
                                  assetNames, assetLevels, barrLevels) );
            }
        }
    }
};


//////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RangeCounter::createProduct(const MonteCarlo* model) const {
    static const string method = "RangeCounter::createProduct";

    try {
        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(assets->NbAssets()));
        // We use the modelling dates, correct?
        simSeries->addDates(obsMap->getModellingDates());

        if(model->stateVarUsed()) {
            return new RangeCounterSVMC(this, simSeries);
        } else {
            throw ModelException(method, "Non-SV Monte Carlo not supported.");
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


CClassConstSP const RangeCounter::Range::TYPE = CClass::registerClassLoadMethod(
    "RangeCounter::Range", typeid(RangeCounter::Range),
    RangeCounter::Range::load);

CClassConstSP const RangeCounter::TYPE = CClass::registerClassLoadMethod(
    "RangeCounter", typeid(RangeCounter), RangeCounter::load);

// for class loading (avoid having header file)
bool RangeCounterLoad() {
    return (RangeCounter::TYPE != 0);
}


DRLIB_END_NAMESPACE
