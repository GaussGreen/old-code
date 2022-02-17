//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Lookback.cpp
//
//   Description : Strike is max level over a set of monitoring dates
//                 The over a subsequent (necessary?) set of dates,, an average level is obtained
//                 Plain vanilla payoff = Max((AvgOut/Strike) -1, 0.0)
//                 We can floor this via lock: = Max((AvgOut/Strike) -1, lock)
//                 the lock can also be a ladder style if we set isLadder to true:
//                 lock * N < Best return < lock * N + 1
//                 and the payoff becomes: = max ((AvgOut/Strike) -1, lock * N)
//                 Final twist is there can be a cap if hasCap is set to true:
//                 payoff: Min(cap, max((AvgOut/Strike) -1, 0.0)
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/ObservationBuilder.hpp"

DRLIB_BEGIN_NAMESPACE

class Lookback : public GenericNFBase,
    virtual public IMCIntoProduct {
protected:

    /// fields ////////
    IDateBuilderSP          lookbackDates;          // all dates in the lookback period
    IDateBuilderSP          avgOutDates;            // all dates in the averaging period
    double                  lock;                   // the lock level or floor, if we don't need one, set to zero
    double                  cap;                    // maximum level for this product
    bool                    hasCap;                 // whether we have a cap
    bool                    isLadderStyleLocking;   // is the lock a ladder-style floor?
    bool                    isFlooredWithStrike;    // do we floor with strike when calculating average?

    //these are transient fields to cache the dates from teh IObservationBuilders
    DateTimeArraySP         internalLookbackDates;
    DateTimeArraySP         internalAveragingDates;

public:

    static CClassConstSP const TYPE;
    friend class LookbackSVMC;

    void Validate() {
        static const string method = "Lookback::Validate";
        try {

            // Date validation moved to Validate, to ensure that obsBuilder has been completely
            // built by GetMarket call. This means that we can't do validation until after
            // the market data has been populated.

            GenericNFBase::Validate();


            //This takes care of both checking for empty dates as well as increasing order
            DateTime::ensureStrictlyIncreasing(*internalLookbackDates,"Lookback Dates",true);
            DateTime::ensureStrictlyIncreasing(*internalAveragingDates, "Averaging Dates",true);

            // for now we do not know how to deal with getting the
            // lookback strike for assets greater than 1 it could be either
            // Max(basket in each timestep) or
            // Max ( each asset in each time step) or arrayMax (per asset)
            if ( internalLookbackDates->size()> 1 && assets->NbAssets() > 1 ) {
                throw ModelException(method,
                    "product supports only a single asset if there is more than one lookback date");
            }

            //validate to make sure averaging dates begin after last lookback date
            if(internalLookbackDates->back().isGreaterOrEqual(internalAveragingDates->front())) {
                throw ModelException(method, "averagingDates (" +
                    internalAveragingDates->front().toString() +
                    ") cannot begin before the end of lookbackDates (" +
                    internalLookbackDates->back().toString() + ")");
            }

            //Any refLevel stuff needs to be done in Validate() to ensure its after GetMarket.
            // need to make sure the ref level averaging in has finished by
            // the time we hit monitoring dates
            DateTimeArray overlap = refLevel->getFutureDates(internalLookbackDates->front());
            if (!overlap.empty()) {
                throw ModelException(method, "Lookback cannot start (" + internalLookbackDates->front().toString() +
                    ") until the reference level has been set");
            }
        }
        catch (exception& e) {
             throw ModelException(e, method);
        }


    }

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {
        static const string method = "Lookback::GetMarket";
        try {
            GenericNFBase::GetMarket(model,market);

            // Get the market for the IDateBuilder
            lookbackDates->getMarket(model, market.get());
            avgOutDates->getMarket(model, market.get());

            internalLookbackDates = lookbackDates->dates();
            internalAveragingDates = avgOutDates->dates();


        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    // validation
    void validatePop2Object(){
        static const string method = "Lookback::validatePop2Object";
        GenericNFBase::validatePop2Object();
        try {
            if (Maths::isNegative(lock)) {
                throw ModelException(method,"lock [" +
                    Format::toString(lock) +"] cannot be negative");
            }
            if (Maths::isNegative(cap)) {
                throw ModelException(method, "cap [" +
                    Format::toString(cap) +"] cannot be negative");
            }
            if(hasCap && cap == 0.0) {
                throw ModelException(method, "cap [" +
                    Format::toString(cap) + "] should be above zero if cap is used");
            }

            if(isLadderStyleLocking && lock == 0.0) {
                throw ModelException(method, "lock [" +
                    Format::toString(lock) + "]  should be above zero if LadderStyleLocking is used");
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return DateTime::merge(*internalLookbackDates, *internalAveragingDates);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
    implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    Lookback(): GenericNFBase(TYPE) {} // for reflection
    Lookback(const Lookback& rhs);     // not implemented
    Lookback& operator=(const Lookback& rhs); // not implemented

    static IObject* defaultLookback(){
        return new Lookback();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Lookback, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultLookback);
        FIELD(lookbackDates, "all dates on which levels are monitored");
        FIELD(avgOutDates, "the dates over which performance will be averaged");
        FIELD(lock, "the minimum level");
        FIELD(cap, "the maximum level");
        FIELD(hasCap, "whether we have a cap");
        FIELD(isLadderStyleLocking,"whether we are using ladder-style lock-in for the lock level");
        FIELD(isFlooredWithStrike,"do we floor with strike when calculating average?");
        FIELD_NO_DESC(internalAveragingDates);
        FIELD_NO_DESC(internalLookbackDates);
        FIELD_MAKE_TRANSIENT(internalLookbackDates);
        FIELD_MAKE_TRANSIENT(internalAveragingDates);
    }
};

/* MC product class for LookbackSV */
class LookbackSVMC : public MCProductClient,
    virtual public IMCProductLN,
    virtual public IMCProductImplied {
private:
    const Lookback* inst;         // Instrument
    int             nbAssets;     // convenient
    DoubleArray     perf;         // current perf per asset
    vector<bool>    lookbackMap;  // this will say whether we should lookback in this timestep
    vector<bool>    averagingMap; // this will say whether we should average in this timestep


    //to store state so far in case of doing past
    double          lookbackStrikeSoFar;
    int             multiplierSoFar;
    double          averageSoFar;
    int             numSoFar;

    //cache these from instrument to save walking the pointer in a tight loop
    bool            isFlooredWithStrike;
    bool            isLadder;
    bool            hasCap;
    double          notional;
    double          lock;
    double          cap;
    bool            isPaymentDateEmpty;

    // State variables and generators
    SVGenSpotSP                 spotGen;      // Generator for spot
    IRefLevel::IStateVarGenSP   refLevelGen;  // Generator for ref level
    SVGenDiscFactorSP               dfGen;        // Generator for discount factors
    SVGenSpot::IStateVarSP          spotSV;       // Spot state variable
    IRefLevel::IStateVarSP      refLevelSV;   // Ref level state variable
    SVDiscFactorSP          dfSV;         // Df state variable


public:

    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "LookbackSVMC::collectStateVars";
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
    past path generator is created, and then when the future path
    generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) {
        static const string routine = "LookbackSVMC::pathGenUpdated";
        try{
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    LookbackSVMC(const Lookback*      inst,
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
        perf(nbAssets,0.0),
        lookbackStrikeSoFar(1.0),
        averageSoFar(0.0),
        multiplierSoFar(inst->isLadderStyleLocking ? 0.0:1.0),
        numSoFar(0),
        isFlooredWithStrike(inst->isFlooredWithStrike),
        isLadder(inst->isLadderStyleLocking),
        hasCap(inst->hasCap),
        notional(inst->notional),
        lock(inst->lock),
        cap(inst->cap),
        isPaymentDateEmpty(paymentDate.empty()),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
        inst->instSettle, simSeries->getLastDate())) {
            lookbackMap = DateTime::createmapping3(simSeries->getAllDates(),*inst->internalLookbackDates);
            averagingMap = DateTime::createmapping3(simSeries->getAllDates(),*inst->internalAveragingDates);
        }

        void payoff(const IPathGenerator*  pathGen,
                    IMCPrices& prices) {

                // Same begin & end for all assets, so read from the first
                const SVPath& path = spotSV->path(0);
                double refLevel = refLevelSV->refLevel(0);
                int    beginIdx = path.begin();
                int    endIdx   = path.end();

                int iStep = beginIdx;
                double basketAverage = 0.0;
                double lookbackStrike = lookbackStrikeSoFar;
                double average = averageSoFar;
                int multiplier = multiplierSoFar;
                double num = numSoFar;

                while (iStep<endIdx) {
                    basketAverage = 0.0;
                    for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                        const SVPath& path = spotSV->path(iAsset);
                        perf[iAsset] = path[iStep] / refLevelSV->refLevel(iAsset);
                        basketAverage += perf[iAsset];
                    }
                    basketAverage /= nbAssets;

                    if (lookbackMap[iStep]) { // if it is a lookback period
                        lookbackStrike = Maths::max(lookbackStrike, basketAverage);
                    }

                    if (averagingMap[iStep]) { //if it is an averaging period
                        ++num;

                        if (isFlooredWithStrike) {
                            average = (average * (num - 1) + Maths::max(lookbackStrike, basketAverage))/num;
                        } else {
                            average = (average * (num - 1) + basketAverage)/num;
                        }

                        if (isLadder)
                        {
                            // effectively floor((avge/lookback-1.0)/cpn)
                            int floored = (average/lookbackStrike - 1.0)/lock;
                            multiplier = Maths::max(multiplier, floored );
                        }
                    }
                    ++iStep;
                }

                // preserve values for past
                if (doingPast()){
                    lookbackStrikeSoFar = lookbackStrike;
                    averageSoFar = average;
                    multiplierSoFar = multiplier;
                    numSoFar = num;
                }


                if (!doingPast() || !hasFuture()) {
                    // Compute a payoff, but only when we have a "complete" situation : either
                    // doingPast() and all is past, or !doingPast().

                    // we can do this without checking whether ladder-style locking is used
                    // because multiplier will be set to 1.0 if its not a ladder-style
                    double value = Maths::max(average/lookbackStrike - 1.0,lock * multiplier);

                    // cap the values if there is a cap
                    if (hasCap) {
                        value = Maths::min(cap,value);
                    }

                    // now scale by notional and discount factor
                    double hitValue = notional * value;

                    // If we're in a known state, we record known flows on
                    // their known dates (so no discounting).
                    if (!isPaymentDateEmpty) {
                        knownCashFlows->addFlow(paymentDate, hitValue);
                    }

                    prices.add(hitValue * dfSV->firstDF());
                }
            }

            /** Use this opportunity to do any Implied driven initialisation
            of the instrument before the main MC loop. e.g closed form barrier adjustment */
            void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

            /** Use this opportunity to do any LogNormal driven initialisation
            of the instrument before the main MC loop. e.g closed form barrier adjustment */
            void initialiseLN(const  IMCPathGenerator*  pathGen)const{
                static const string routine = "LookbackSVMC::initialiseLN";
                throw ModelException(routine, "Methodology not supported");
            }


            // any old level so that MC implied works
            CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                int                     iAsset) const {
                    static const string method = "LookbackSVMC::getVolInterp";

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

    };


//////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Lookback::createProduct(const MonteCarlo* model) const {
    static const string method = "Lookback::createProduct";

    try {
        int nbAssets = assets->NbAssets();

        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(nbAssets));

        // do we crop these dates to only be future ones??

        // do we crop these dates to only be future ones??
        simSeries->addDates(DateTime::merge(*internalAveragingDates, *internalLookbackDates));


        if(model->stateVarUsed()) {
            return new LookbackSVMC(this, simSeries);
        } else {
            throw ModelException(method, "Non-SV Monte Carlo not supported.");
        }

    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


CClassConstSP const Lookback::TYPE = CClass::registerClassLoadMethod(
    "Lookback", typeid(Lookback), Lookback::load);

// * for class loading (avoid having header file) */
bool LookbackLoad() {
    return (Lookback::TYPE != 0);
}

DRLIB_END_NAMESPACE