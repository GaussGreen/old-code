//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EnhancedStraggler.cpp
//
//   Description : A combination of a long position in a call option and a short position
//                  in a KI/KO put. The call option can lock in if it goes above a 3rd barrier.
//                  The strike for the call is set as the worst performance in a lookback period (with cap/floor)
//                  All options are on the worst performer at maturity
//
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/Events.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/DateBuilder.hpp"

DRLIB_BEGIN_NAMESPACE

class EnhancedStraggler: public GenericNFBase,
    virtual public IMCIntoProduct,
    virtual public LegalTerms::Shift,
    virtual public BarrierBreach::IEventHandler {
protected:
    /// fields ////////
    DateTime        finalLookback;      // final date in look back period
    IDateBuilderSP  monitoringDates;    // all dates after ref date to the end
    double          kiBarrier;          // the knock in barrier
    double          kiEconBarrier;      // the knock in barrier for legal terms
    double          koBarrier;          // the knock out barrier
    double          koEconBarrier;      // the knock out barrier for legal terms
    double          lockBarrier;        // the lock in barrier
    double          lockEconBarrier;    // the lock in barrier for legal terms
    double          strikeCap;          // the cap for the lookback strike
    double          strikeFloor;        // the floor for the lookback strike
    double          participation;      // the participation in the call option
    double          putStrike;          // the strike for the put option

    int             finalLookbackIdx;   // index in mon dates of last lookback

public:
    static CClassConstSP const TYPE;
    friend class EnhancedStragglerMC;
    friend class EnhancedStragglerSVMC;

    void Validate() {
        static const string method = "EnhancedStraggler::Validate";
        GenericNFBase::Validate();

        //sampleDates = DateTimeClusterConstSP(&(obsMap->getSampleDates()));

        //nbDates = 0;

        //for (int iAsset=0; iAsset< assets->NbAssets(); ++iAsset) {
        //  if ((*sampleDates)[iAsset].size()<2)
        //      throw ModelException(routine,"At least two dates required, "+
        //              Format::toString((*sampleDates)[iAsset].size())+
        //              " remains after ISDA adjustment for asset "+
        //              Format::toString(iAsset + 1) + ".");

        //}

		DateTimeArraySP monDates = monitoringDates->dates();

		// validate dates are not empty and are in order

		// validate dates are not empty - order is handled by SimSeries
		if (monDates->empty()) {
			throw ModelException(method, "No monitoring dates supplied!");
		}

		// need to make sure the ref level averaging in has finished by the time we hit monitoring dates
		DateTimeArray overlap = refLevel->getFutureDates((*monDates)[0]);
		if (!overlap.empty()) {
			throw ModelException(method, "Monitoring dates cannot start ("+ (*monDates)[0].toString()
				+") until the reference level has been set");
		}

		// check final lookback is amongst the monitoring periods
		try {
			finalLookbackIdx = finalLookback.find(*monDates);
		} catch (exception& e) {
			throw ModelException(e, method, "finalLookBackDate ("
				+ finalLookback.toString() + ") is not amongst the monitoring dates");
		}

    }

    // validation
	/** Get the asset and discount market data */
	virtual void GetMarket(const IModel*          model,
		const CMarketDataSP    market) {
		static const string method("CorrCov::GetMarket");
		try
		{
			GenericNFBase::GetMarket(model,market);
			monitoringDates->getMarket(model,market.get());
		}
		catch (exception& e) {
			throw ModelException(e,method);
		}

	}
    void validatePop2Object(){
        static const string method = "EnhancedStraggler::validatePop2Object";
        GenericNFBase::validatePop2Object();

        // check knock out barrier is above the knock in barrier
        if (!Maths::isPositive(koBarrier - kiBarrier)){
            throw ModelException(method, "KO barrier ("+
                Format::toString(koBarrier)+") should be above KI barrier ("+
                Format::toString(kiBarrier)+")");
        }
        if (!Maths::isPositive(koEconBarrier - kiEconBarrier)){
            throw ModelException(method, "Legal KO barrier ("+
                Format::toString(koEconBarrier)+") should be above legal KI barrier ("+
                Format::toString(kiEconBarrier)+")");
        }

        // check strike cap is above strike floor
        if (Maths::isNegative(strikeCap - strikeFloor)){
            throw ModelException(method, "Strike cap ("+
                Format::toString(strikeCap)+") is below strike floor ("+
                Format::toString(strikeFloor)+")");
        }

        if (!Maths::isPositive(putStrike)){
            throw ModelException(method, "put strike ("+
                Format::toString(putStrike)+") is not positive");
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
		return *(monitoringDates->dates());
    }

    // set barriers to be the economic (legal) ones
    bool sensShift(LegalTerms* shift) {
        // just replace all 3 barriers with the corresponding economic one
        kiBarrier = kiEconBarrier;
        koBarrier = koEconBarrier;
        lockBarrier = lockEconBarrier;

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
                        "Internal error - expected Monte Carlo model for EnhancedStraggler pricing");
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
    EnhancedStraggler(): GenericNFBase(TYPE) {} // for reflection
    EnhancedStraggler(const EnhancedStraggler& rhs);     // not implemented
    EnhancedStraggler& operator=(const EnhancedStraggler& rhs); // not implemented

    static IObject* defaultEnhancedStraggler(){
        return new EnhancedStraggler();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnhancedStraggler, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultEnhancedStraggler);
        FIELD(finalLookback, "last date in lookback period");
        FIELD(monitoringDates, "all dates on which levels are monitored");
        FIELD(kiBarrier, "the knock in barrier");
        FIELD(kiEconBarrier, "the knock in barrier for legal terms");
        FIELD(koBarrier, "the knock out barrier");
        FIELD(koEconBarrier, "the knock out barrier for legal terms");
        FIELD(lockBarrier, "the lock in barrier");
        FIELD(lockEconBarrier, "the lock in barrier for legal terms");
        FIELD(strikeCap, "the cap for the lookback strike");
        FIELD(strikeFloor, "the floor for the lookback strike");
        FIELD(participation, "the participation in the call option");
        FIELD(putStrike, "the strike for the put option");
        FIELD(finalLookbackIdx, "index in monitoring dates of last date in lookback period");
        FIELD_MAKE_TRANSIENT(finalLookbackIdx);
    }
    };

/* MC product class for EnhancedStraggler */
class EnhancedStragglerMC: public IMCProduct,
    virtual public IMCProductLN,
    virtual public IMCProductImplied {
private:
    const EnhancedStraggler*    inst;         // Instrument
    int                         nbAssets;     // convenient
    DoubleArray                 perf;         // current perf per asset

    // maintain state of the instrument in the past
    bool        hasLockedSoFar;
    bool        hasKISoFar;
    bool        hasKOSoFar;
    bool        allDeadSoFar;
    double      lbStrikeSoFar;
    ObservationMapSP obsMap;


public:

    EnhancedStragglerMC(const EnhancedStraggler*      inst,
        const SimSeriesSP&   simSeries):
    IMCProduct(inst->assets.get(),
        inst->valueDate,
        inst->discount.get(),
        inst->refLevel,
        simSeries,
        inst->pastValues,
        inst->instSettle.get(),
        simSeries->getLastDate()),
        inst(inst),
        nbAssets(getNumAssets()),
        perf(nbAssets, 0.0),
        hasLockedSoFar(false), hasKISoFar(false),
        hasKOSoFar(false), allDeadSoFar(false),
        lbStrikeSoFar(inst->strikeCap),
        obsMap(inst->obsMap){}

        void payoff(const IPathGenerator*  pathGen,
            IMCPrices&                prices) {

                // update with preserved values from doing past
                bool hasLocked = hasLockedSoFar;
                bool hasKI = hasKISoFar;
                bool hasKO = hasKOSoFar;
                bool allDead = allDeadSoFar;
                double lbStrike = lbStrikeSoFar;

                // indexes of where the instrument triggered the barrier
                int lockIdx = 0;
                int kiIdx = 0;
                int koIdx = 0;

                // final asset performances
                double minPerf = 0.0;

                // evaluate path against the instrument
                evaluatePath(pathGen, false,
                    hasLocked, lockIdx,
                    hasKI, kiIdx,
                    hasKO, koIdx,
                    allDead, lbStrike, minPerf);

                // preserve values for past
                if (doingPast()){
                    hasLockedSoFar = hasLocked;
                    hasKISoFar = hasKI;
                    hasKOSoFar = hasKO;
                    allDeadSoFar = allDead;
                    lbStrikeSoFar = lbStrike;
                }

                if (!doingPast() || !hasFuture()) {
                    // Compute a payoff, but only when we have a "complete" situation : either
                    // doingPast() and all is past, or !doingPast().

                    // cap and floor the lookback strike
                    lbStrike = Maths::min(Maths::max(lbStrike, inst->strikeFloor), inst->strikeCap);

                    // redemption amount
                    double value = 1.0;

                    // get the call with possible lock in
                    value += inst->participation * Maths::max(minPerf - lbStrike,
                        hasLocked ? inst->lockEconBarrier - 1.0 : 0.0);
                    if (hasKI && !hasKO) {
                        // pay the put if it has knocked in but not knocked out
                        value -= Maths::max(inst->putStrike - minPerf, 0.0);
                    }
                    // now scale by notional
                    prices.add(inst->notional * value);
                }
            }

            // Evaluate a set of asset paths against the barriers.
            // The state of these barriers (and the index of the point where they
            // triggered the barrier) is passed through the has* and *Idx parameters.
            // Also, allDead, which indicates that remaining uncrossed barrier have
            // no effect on the payoff, lbStrike, the lookback strike and minPerf,
            // the worst performer on the final step, are returned.
            // The reportAll flag makes the algorithmn run in a less optimal fashion,
            // so that barriers that have no effect on value, i.e. knock-in after
            // knock-out, are still evaluated.
            inline void evaluatePath( const IPathGenerator* pathGen, bool reportAll,
                bool& hasLocked, int& lockIdx,
                bool& hasKI, int& kiIdx,
                bool& hasKO, int& koIdx,
                bool& allDead, double& lbStrike, double& minPerf ) {
                    // local barrier pointing to legal or risk depending on past or not
                    double kiBarr = doingPast() ? inst->kiEconBarrier : inst->kiBarrier;
                    double koBarr = doingPast() ? inst->koEconBarrier : inst->koBarrier;
                    double lockBarr = doingPast() ? inst->lockEconBarrier : inst->lockBarrier;

                    int beginIdx = pathGen->begin(0); // same for all assets
                    int endIdx   = pathGen->end(0);


                    //we have to watch out for the case where one of the assets don't have a sample yet due to
                    //ISDA adjustment, this block of code caters for this, should it be inside the obsMap?
                    // if endIdx == 0 we have no past to worry about
                    if (doingPast() && endIdx != 0) {
                        // work out the last completely fixed sample set
                        int nbDatesPast = inst->obsMap->getLastSampleIndex(0/*assetidx*/, endIdx - 1);
                        for (int iAsset = 1; iAsset <nbAssets; iAsset++) {
                            nbDatesPast = Maths::min(nbDatesPast, inst->obsMap->getLastSampleIndex(iAsset, pathGen->end(iAsset)-1));
                        }
                        nbDatesPast++; // this is because we do iStep < endIdx, its a little bit unfortunate since
                        // pathGen->end() actually returns pathGen->size rather than end, we have to make the rest
                        // of the code consistent.
                        endIdx = nbDatesPast;
                    }

                    int iStep = beginIdx;

                    while ( (!allDead || reportAll) && iStep<endIdx) {
                        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                            perf[iAsset] = pathGen->Path(iAsset, 0/*iPath*/)[obsMap->getModellingIndex(iAsset,iStep)]
                            / pathGen->refLevel(iAsset, 0/*iPath*/);
                            minPerf = iAsset > 0 ? Maths::min(minPerf, perf[iAsset]) : perf[iAsset];
                        }

                        // check barrier status
                        if( !hasKI ) {
                            hasKI = !Maths::isPositive(minPerf - kiBarr);
                            kiIdx = iStep;
                        }
                        if( !hasKO ) {
                            hasKO = !Maths::isNegative(minPerf - koBarr);
                            koIdx = iStep;
                        }
                        if( !hasLocked ) {
                            hasLocked = !Maths::isNegative(minPerf - lockBarr);
                            lockIdx = iStep;
                        }

                        // do lookback if within lookback period
                        if (iStep <= inst->finalLookbackIdx) {
                            lbStrike = Maths::min(lbStrike, minPerf);
                        } else { // can only finish once lookback done
                            allDead = hasKO && hasLocked;
                        }
                        iStep++;
                    }

                    if (allDead && iStep < endIdx) {
                        // we stopped monitoring early so we need to just get final asset perfs
                        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                            perf[iAsset] = pathGen->Path(iAsset, 0/*iPath*/)[obsMap->getModellingIndex(iAsset,endIdx-1)]
                            / pathGen->refLevel(iAsset, 0/*iPath*/);
                            minPerf = iAsset > 0 ? Maths::min(minPerf, perf[iAsset]) : perf[iAsset];
                        }
                    }
                }

                /** Use this opportunity to do any LogNormal driven initialisation
                of the instrument before the main MC loop. e.g closed form barrier adjustment */
                void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

                /** Use this opportunity to do any LogNormal driven initialisation
                of the instrument before the main MC loop. e.g closed form barrier adjustment */
                void initialiseLN(const  IMCPathGenerator*  pathGen)const{
                    static const string routine = "EnhancedStragglerMC::initialiseLN";
                    throw ModelException(routine, "Methodology not supported");
                }

                // any old level so that MC implied works
                CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                    int                     iAsset) const {
                        static const string method = "EnhancedStragglerMC::getVolInterp";

                        try {
                            // one interp level per asset
                            CVolRequestLNArray reqarr(1);
                            const DateTime& today = getToday();
                            const DateTime& startDate = refLevel->getAllDates().front();
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

                    // Override IMCProduct::getEvents to support BarrierBreach event reporting.
                    virtual void getEvents(const IMCPathGenerator*  pathGen,
                        EventResults* events,
                        const DateTime& eventDate)
                    {
                        static const string method = "EnhancedStragglerMC::getEvents";
                        try {
                            // barrier state variables
                            bool hasLocked = false;
                            bool hasKI = false;
                            bool hasKO = false;
                            bool allDead = false;
                            double lbStrike = 0.0;

                            // indexes of where the instrument triggered the barrier
                            int lockIdx = 0;
                            int kiIdx = 0;
                            int koIdx = 0;

                            // final asset performances
                            double minPerf = 0.0;

                            // evaluate the historical path against the instrument
                            evaluatePath(pathGen, true,
                                hasLocked, lockIdx,
                                hasKI, kiIdx,
                                hasKO, koIdx,
                                allDead, lbStrike, minPerf);

                            // There MUST ONLY BE a trivial 1-to-1 mapping between the index
                            // and the date, so that we can use the dates for asset 0 for
                            // timestamping the events.
                            // If EHS is changed to work otherwise, then this assumption will
                            // have to be fixed
                            if( !getSimSeries()->sameDatesPerAsset() ) {
                                throw ModelException("Expected same simulation dates across all assets");
                            }
                            const DateTimeArray& simDates = getSimSeries()->getDates(0);

                            // report the barrier breach events
                            if( hasLocked ) {
                                // Lock-in occurs when ALL underlyings are ABOVE lock-in level

                                // Construct details of assets at breach
                                StringArraySP assetNames(new StringArray(nbAssets));
                                DoubleArraySP assetLevels(new DoubleArray(nbAssets));
                                DoubleArraySP barrLevels(new DoubleArray(nbAssets));

                                for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                                    (*assetNames)[iAsset] = getMultiFactors()->assetGetTrueName(iAsset);
                                    (*assetLevels)[iAsset] = pathGen->Path(iAsset, 0)[obsMap->getModellingIndex(iAsset,lockIdx)];
                                    (*barrLevels)[iAsset] = pathGen->refLevel(iAsset, 0) * inst->lockEconBarrier;
                                }

                                // Add new BarrierBreach event
                                events->addEvent( new BarrierBreach(simDates[lockIdx],
                                    "Enhanced Straggler Lock In", BarrierBreach::EUROPEAN,
                                    BarrierBreach::LOCK_IN,
                                    true, 0, assetNames, assetLevels, barrLevels) );

                            }

                            if( hasKI ) {
                                // Knock-in occurs when WORST PERFORMER is BELOW knock-in level
                                // Find worst performer
                                int worstAsset = 0;
                                double worstPerf = pathGen->Path(0, 0)[obsMap->getModellingIndex(0,kiIdx)]
                                / pathGen->refLevel(0, 0);
                                // NB Start from second asset
                                for( int iAsset=1; iAsset<nbAssets; iAsset++ ) {
                                    double perf = pathGen->Path(iAsset, 0)[obsMap->getModellingIndex(iAsset,kiIdx)]
                                    / pathGen->refLevel(iAsset, 0);
                                    if( Maths::isPositive(worstPerf-perf) ) {
                                        worstPerf = perf;
                                        worstAsset = iAsset;
                                    }
                                }
                                // Construct details of assets at breach
                                StringArraySP assetNames(new StringArray(1, getMultiFactors()->assetGetTrueName(worstAsset)));
                                DoubleArraySP assetLevels(new DoubleArray(1, pathGen->Path(worstAsset, 0)[obsMap->getModellingIndex(worstAsset,kiIdx)]));
                                DoubleArraySP barrLevels(new DoubleArray(1, pathGen->refLevel(worstAsset, 0) * inst->kiEconBarrier));

                                // Add new BarrierBreach event
                                events->addEvent( new BarrierBreach( simDates[kiIdx],
                                    "Enhanced Straggler Knock In", BarrierBreach::EUROPEAN,
                                    BarrierBreach::KNOCK_IN,
                                    false, 0, assetNames, assetLevels, barrLevels) );
                            }

                            if( hasKO ) {
                                // Knock-out occurs when ALL underlyings are ABOVE knock-out level

                                // Construct details of assets at breach
                                StringArraySP assetNames(new StringArray(nbAssets));
                                DoubleArraySP assetLevels(new DoubleArray(nbAssets));
                                DoubleArraySP barrLevels(new DoubleArray(nbAssets));

                                for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                                    (*assetNames)[iAsset] = getMultiFactors()->assetGetTrueName(iAsset);
                                    (*assetLevels)[iAsset] = pathGen->Path(iAsset, 0)[obsMap->getModellingIndex(iAsset,koIdx)];
                                    (*barrLevels)[iAsset] = pathGen->refLevel(iAsset, 0) * inst->koEconBarrier;
                                }

                                // Add new BarrierBreach event
                                events->addEvent( new BarrierBreach( simDates[koIdx],
                                    "Enhanced Straggler Knock Out", BarrierBreach::EUROPEAN,
                                    BarrierBreach::KNOCK_OUT,
                                    true, 0, assetNames, assetLevels, barrLevels) );
                            }
                        } catch (exception& e) {
                            throw ModelException(e, method);
                        }
                    }
    };

/* MC product class for EnhancedStragglerSV */
class EnhancedStragglerSVMC: public MCProductClient,
    virtual public IMCProductLN,
    virtual public IMCProductImplied {
private:
    const EnhancedStraggler*    inst;         // Instrument
    int                         nbAssets;     // convenient
    DoubleArray                 perf;         // current perf per asset

    // maintain state of the instrument in the past
    bool        hasLockedSoFar;
    bool        hasKISoFar;
    bool        hasKOSoFar;
    bool        allDeadSoFar;
    double      lbStrikeSoFar;

    // State variables and generators
    SVGenSpotSP                  spotGen;      // Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  // Generator for ref level
    SVGenDiscFactorSP            dfGen;        // Generator for discount factors
    SVGenSpot::IStateVarSP       spotSV;       // Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   // Ref level state variable
    SVDiscFactorSP dfSV;         // Df state variable

    ObservationMapSP    obsMap;
public:

    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "EnhancedStragglerSVMC::collectStateVars";
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
        static const string routine = "EnhancedStragglerSVMC::pathGenUpdated";
        try{
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    EnhancedStragglerSVMC(const EnhancedStraggler*      inst,
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
        perf(nbAssets, 0.0),
        hasLockedSoFar(false), hasKISoFar(false),
        hasKOSoFar(false), allDeadSoFar(false),
        lbStrikeSoFar(inst->strikeCap),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
        inst->instSettle, simSeries->getLastDate())),
        obsMap(inst->obsMap){
        }

        void payoff(const IPathGenerator*  pathGen,
            IMCPrices&                prices) {

                // update with preserved values from doing past
                bool hasLocked = hasLockedSoFar;
                bool hasKI = hasKISoFar;
                bool hasKO = hasKOSoFar;
                bool allDead = allDeadSoFar;
                double lbStrike = lbStrikeSoFar;

                // indexes of where the instrument triggered the barrier
                int lockIdx = 0;
                int kiIdx = 0;
                int koIdx = 0;

                // final asset performances
                double minPerf = 0.0;

                // evaluate path against the instrument
                evaluatePath(pathGen, false,
                    hasLocked, lockIdx,
                    hasKI, kiIdx,
                    hasKO, koIdx,
                    allDead, lbStrike, minPerf);

                // preserve values for past
                if (doingPast()){
                    hasLockedSoFar = hasLocked;
                    hasKISoFar = hasKI;
                    hasKOSoFar = hasKO;
                    allDeadSoFar = allDead;
                    lbStrikeSoFar = lbStrike;
                }

                if (!doingPast() || !hasFuture()) {
                    // Compute a payoff, but only when we have a "complete" situation : either
                    // doingPast() and all is past, or !doingPast().

                    // cap and floor the lookback strike
                    lbStrike = Maths::min(Maths::max(lbStrike, inst->strikeFloor), inst->strikeCap);

                    // redemption amount
                    double value = 1.0;

                    // get the call with possible lock in
                    value += inst->participation * Maths::max(minPerf - lbStrike,
                        hasLocked ? inst->lockEconBarrier - 1.0 : 0.0);
                    if (hasKI && !hasKO) {
                        // pay the put if it has knocked in but not knocked out
                        value -= Maths::max(inst->putStrike - minPerf, 0.0);
                    }
                    // now scale by notional and discount factor
                    prices.add(inst->notional * value * dfSV->firstDF());
                }
            }

            // Evaluate a set of asset paths against the barriers.
            // The state of these barriers (and the index of the point where they
            // triggered the barrier) is passed through the has* and *Idx parameters.
            // Also, allDead, which indicates that remaining uncrossed barrier have
            // no effect on the payoff, lbStrike, the lookback strike and minPerf,
            // the worst performer on the final step, are returned.
            // The reportAll flag makes the algorithmn run in a less optimal fashion,
            // so that barriers that have no effect on value, i.e. knock-in after
            // knock-out, are still evaluated.
            inline void evaluatePath( const IPathGenerator* pathGen, bool reportAll,
                bool& hasLocked, int& lockIdx,
                bool& hasKI, int& kiIdx,
                bool& hasKO, int& koIdx,
                bool& allDead, double& lbStrike, double& minPerf ) {
                    // local barrier pointing to legal or risk depending on past or not
                    double kiBarr = doingPast() ? inst->kiEconBarrier : inst->kiBarrier;
                    double koBarr = doingPast() ? inst->koEconBarrier : inst->koBarrier;
                    double lockBarr = doingPast() ? inst->lockEconBarrier : inst->lockBarrier;
                    // Same begin & end for all assets, so read from the first
                    const SVPath& path = spotSV->path(0);
                    int    beginIdx = path.begin();
                    int    endIdx   = path.end();

                    //we have to watch out for the case where one of the assets don't have a sample yet due to
                    //ISDA adjustment, this block of code caters for this, should it be inside the obsMap?
                    // if endIdx == 0 we have no past to worry about
                    if (doingPast() && endIdx != 0) {
                        // work out the last completely fixed sample set
                        int nbDatesPast = inst->obsMap->getLastSampleIndex(0/*assetidx*/, endIdx - 1);
                        for (int iAsset = 1; iAsset <nbAssets; iAsset++) {
                            nbDatesPast = Maths::min(nbDatesPast, inst->obsMap->getLastSampleIndex(iAsset,spotSV->path(iAsset).end()-1));
                        }
                        nbDatesPast++; // this is because we do iStep < endIdx, its a little bit unfortunate since
                        // pathGen->end() actually returns pathGen->size rather than end, we have to make the rest
                        // of the code consistent.
                        endIdx = nbDatesPast;
                    }

                    int iStep = beginIdx;

                    while ( (!allDead || reportAll) && iStep<endIdx) {
                        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                            const SVPath& path = spotSV->path(iAsset);
                            perf[iAsset] = path[obsMap->getModellingIndex(iAsset,iStep)] / refLevelSV->refLevel(iAsset);
                            minPerf = iAsset > 0 ? Maths::min(minPerf, perf[iAsset]) : perf[iAsset];
                        }

                        // check barrier status
                        if( !hasKI ) {
                            hasKI = !Maths::isPositive(minPerf - kiBarr);
                            kiIdx = iStep;
                        }
                        if( !hasKO ) {
                            hasKO = !Maths::isNegative(minPerf - koBarr);
                            koIdx = iStep;
                        }
                        if( !hasLocked ) {
                            hasLocked = !Maths::isNegative(minPerf - lockBarr);
                            lockIdx = iStep;
                        }

                        // do lookback if within lookback period
                        if (iStep <= inst->finalLookbackIdx) {
                            lbStrike = Maths::min(lbStrike, minPerf);
                        } else { // can only finish once lookback done
                            allDead = hasKO && hasLocked;
                        }
                        iStep++;
                    }

                    if (allDead && iStep < endIdx) {
                        // we stopped monitoring early so we need to just get final asset perfs
                        for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                            const SVPath& path = spotSV->path(iAsset);
                            perf[iAsset] = path[obsMap->getModellingIndex(iAsset,endIdx-1)] / refLevelSV->refLevel(iAsset);
                            minPerf = iAsset > 0 ? Maths::min(minPerf, perf[iAsset]) : perf[iAsset];
                        }
                    }
                }

                /** Use this opportunity to do any LogNormal driven initialisation
                of the instrument before the main MC loop. e.g closed form barrier adjustment */
                void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

                /** Use this opportunity to do any LogNormal driven initialisation
                of the instrument before the main MC loop. e.g closed form barrier adjustment */
                void initialiseLN(const  IMCPathGenerator*  pathGen)const{
                    static const string routine = "EnhancedStragglerSVMC::initialiseLN";
                    throw ModelException(routine, "Methodology not supported");
                }


                // any old level so that MC implied works
                CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                    int                     iAsset) const {
                        static const string method = "EnhancedStragglerSVMC::getVolInterp";

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

                    // Override IMCProduct::getEvents to support BarrierBreach event reporting.
                    virtual void getEvents(const IMCPathGenerator*  pathGen,
                        EventResults* events,
                        const DateTime& eventDate) {
                            static const string method = "EnhancedStragglerSVMC::getEvents";
                            try {
                                // barrier state variables
                                bool hasLocked = false;
                                bool hasKI = false;
                                bool hasKO = false;
                                bool allDead = false;
                                double lbStrike = 0.0;

                                // indexes of where the instrument triggered the barrier
                                int lockIdx = 0;
                                int kiIdx = 0;
                                int koIdx = 0;

                                // final asset performances
                                double minPerf = 0.0;

                                // evaluate the historical path against the instrument
                                evaluatePath(pathGen, true,
                                    hasLocked, lockIdx,
                                    hasKI, kiIdx,
                                    hasKO, koIdx,
                                    allDead, lbStrike, minPerf);

                                // There MUST ONLY BE a trivial 1-to-1 mapping between the index
                                // and the date, so that we can use the dates for asset 0 for
                                // timestamping the events.
                                // If EHS is changed to work otherwise, then this assumption will
                                // have to be fixed
                                if( !getSimSeries()->sameDatesPerAsset() ) {
                                    throw ModelException("Expected same simulation dates across all assets");
                                }
                                const DateTimeArray& simDates = getSimSeries()->getDates(0);

                                // report the barrier breach events
                                if( hasLocked ) {
                                    // Lock-in occurs when ALL underlyings are ABOVE lock-in level

                                    // Construct details of assets at breach
                                    StringArraySP assetNames(new StringArray(nbAssets));
                                    DoubleArraySP assetLevels(new DoubleArray(nbAssets));
                                    DoubleArraySP barrLevels(new DoubleArray(nbAssets));

                                    for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                                        const SVPath& path = spotSV->path(iAsset);
                                        (*assetNames)[iAsset] = getMultiFactors()->assetGetTrueName(iAsset);
                                        (*assetLevels)[iAsset] = path[obsMap->getModellingIndex(iAsset,lockIdx)];
                                        (*barrLevels)[iAsset] = refLevelSV->refLevel(iAsset) * inst->lockEconBarrier;
                                    }

                                    // Add new BarrierBreach event
                                    events->addEvent( new BarrierBreach( simDates[lockIdx],
                                        "Enhanced Straggler Lock In", BarrierBreach::EUROPEAN,
                                        BarrierBreach::LOCK_IN,
                                        true, 0, assetNames, assetLevels, barrLevels) );
                                }

                                if( hasKI ) {
                                    // Knock-in occurs when WORST PERFORMER is BELOW knock-in level

                                    // Find worst performer
                                    int worstAsset = 0;
                                    double worstPerf = (spotSV->path(0))[obsMap->getModellingIndex(0,kiIdx)]
                                    / refLevelSV->refLevel(0);
                                    // NB Start from second asset
                                    for( int iAsset=1; iAsset<nbAssets; iAsset++ ) {
                                        const SVPath& path = spotSV->path(iAsset);
                                        double perf = path[obsMap->getModellingIndex(iAsset,kiIdx)] / refLevelSV->refLevel(iAsset);
                                        if( Maths::isPositive(worstPerf-perf) ) {
                                            worstPerf = perf;
                                            worstAsset = iAsset;
                                        }
                                    }

                                    // Construct details of assets at breach
                                    StringArraySP assetNames(new StringArray(1, getMultiFactors()->assetGetTrueName(worstAsset)));
                                    DoubleArraySP assetLevels(new DoubleArray(1, (spotSV->path(worstAsset))[obsMap->getModellingIndex(worstAsset,kiIdx)]));
                                    DoubleArraySP barrLevels(new DoubleArray(1, refLevelSV->refLevel(worstAsset) * inst->kiEconBarrier));

                                    // Add new BarrierBreach event
                                    events->addEvent( new BarrierBreach( simDates[kiIdx],
                                        "Enhanced Straggler Knock In", BarrierBreach::EUROPEAN,
                                        BarrierBreach::KNOCK_IN,
                                        false, 0, assetNames, assetLevels, barrLevels) );
                                }

                                if( hasKO ) {
                                    // Knock-out occurs when ALL underlyings are ABOVE knock-out level

                                    // Construct details of assets at breach
                                    StringArraySP assetNames(new StringArray(nbAssets));
                                    DoubleArraySP assetLevels(new DoubleArray(nbAssets));
                                    DoubleArraySP barrLevels(new DoubleArray(nbAssets));

                                    for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                                        const SVPath& path = spotSV->path(iAsset);
                                        (*assetNames)[iAsset] = getMultiFactors()->assetGetTrueName(iAsset);
                                        (*assetLevels)[iAsset] = path[obsMap->getModellingIndex(iAsset,koIdx)];
                                        (*barrLevels)[iAsset] = refLevelSV->refLevel(iAsset) * inst->koEconBarrier;
                                    }

                                    // Add new BarrierBreach event
                                    events->addEvent( new BarrierBreach( simDates[koIdx],
                                        "Enhanced Straggler Knock Out", BarrierBreach::EUROPEAN,
                                        BarrierBreach::KNOCK_OUT,
                                        true, 0, assetNames, assetLevels, barrLevels) );
                                }
                            } catch (exception& e) {
                                throw ModelException(e, method);
                            }
                        }
    };


//////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* EnhancedStraggler::createProduct(const MonteCarlo* model) const {
    static const string method = "EnhancedStraggler::createProduct";

    try {
        int nbAssets = assets->NbAssets();

        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(nbAssets));

        // do we crop these dates to only be future ones??
        simSeries->addDates(obsMap->getModellingDates());

        if(model->stateVarUsed()) {
            return new EnhancedStragglerSVMC(this, simSeries);
        } else {
            // Otherwise, use old methodology
            return new EnhancedStragglerMC(this, simSeries);
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


CClassConstSP const EnhancedStraggler::TYPE = CClass::registerClassLoadMethod(
    "EnhancedStraggler", typeid(EnhancedStraggler), EnhancedStraggler::load);

// * for class loading (avoid having header file) */
bool EnhancedStragglerLoad() {
    return (EnhancedStraggler::TYPE != 0);
}

DRLIB_END_NAMESPACE
