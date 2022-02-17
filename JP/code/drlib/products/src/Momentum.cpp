//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Momentum.cpp
//
//   Description : Port of Momentum
//
//   Date        : Nov 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/ITimeAggregate.hpp"
#include "edginc/HandlePaymentEvents.hpp"

DRLIB_BEGIN_NAMESPACE

class Momentum: public GenericNFBase, 
                virtual public IMCIntoProduct{
protected:
    /// fields 
    //perhaps makes sense for gen time agg to contain dates itself?
    ITimeAggregateMakerSP   timeAggregate;

    // Since the floorRate (if specified, otherwise zero) is booked separately, it is deducted from the coupon here.
    double          floorRate;

    // a possible number deducted from momentum payment. Used for example in momentum option
    double          momentumStrike;

    // a multiplier to produce momentum payment
    double          participation;

    // whether coupon is floored with previous coupon or with floorRate
    bool            floorWithPreviousCoup;

    // whether it is bullish momentum;
    bool            isBullish;

    // monitoring dates and coupon payment dates
    DateTimeArray   averageOutDates;        
    
public:
    static CClassConstSP const TYPE;
    friend class MomentumMC;

    // validation
    void validatePop2Object(){
        static const string routine = "Momentum::validatePop2Object";
        GenericNFBase::validatePop2Object();
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates given!");
        }
        const DateTimeArray& periodDates = timeAggregate->getDates();
        if (periodDates.empty()) {
            throw ModelException(routine, "No couponDates given!");
        }

        // Check average dates and coupon dates "cooperate"
        // Coupon dates must be a subset of average dates
        if (!DateTime::isSubset(averageOutDates, periodDates)) {
            throw ModelException(routine, "CouponDates should be a subset of averageOutDates");
        }
        // Require some averaging before first cliquet date
        const DateTime& firstAvg = averageOutDates[0];
        const DateTime& firstCpn = periodDates[0];
        if (firstCpn < firstAvg) {
            throw ModelException(routine, "Cannot have coupon date " + firstCpn.toString() + 
                                 " before first average date " + firstAvg.toString());
        }
        // Makes no sense to have averaging after final cliquet date
        const DateTime& lastAvg = averageOutDates[averageOutDates.size()-1];
        const DateTime& lastCpn = periodDates[periodDates.size()-1];
        if (lastAvg > lastCpn) {
            throw ModelException(routine, "Cannot average on " + lastAvg.toString() + 
                                 " since after final conpon date " + lastCpn.toString());
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return averageOutDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    Momentum(): GenericNFBase(TYPE) {}    // for reflection
    Momentum(const Momentum& rhs);            // not implemented
    Momentum& operator=(const Momentum& rhs); // not implemented

    static IObject* defaultMomentum(){
        return new Momentum();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Momentum, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultMomentum);
        FIELD(timeAggregate,"How to combine values per time period");
        FIELD(floorRate, "Minimal coupon rate. Booked separately and so is deducted from coupon here");
        FIELD(momentumStrike, "Level deducted from momentum payment");
        FIELD(participation, "Participation Factor");
        FIELD(floorWithPreviousCoup, "Whether pervious coupon is needed to determine current coupon");
        FIELD(isBullish, "Whether it is a bullish momentum");
        FIELD(averageOutDates, "Average out dates for coupon determination");
    }
};

//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
class MomentumMC: public IMCProduct,
                  virtual public IMCProductLN,
                  virtual public IMCProductImplied,
                  virtual public IHandlePaymentEvents,
                  virtual public ITimeAggProductView {
private:
    const Momentum*        inst;
    int                    nbAssets;        // convenient
    int                    nbPeriods;       // convenient
    
    ITimeAggregateSP       timeAgg;         // how to turn per-period gen perfs into a single value

    // per-period values
    SimpleDoubleArray      periodPerfs;     // time dim across periods
    IntArray               periodMap;       // to track periods
    IntArray               nbAvgOutPerPeriod;  

    // per-asset values
    DoubleArray            refLevels;
    DoubleArray            priceRefs;
    DoubleArray            sum;

    // for past
    DoubleArray            sumSoFar;        //[nbAssets]
    DoubleArray            refLevelsSoFar;
    DoubleArray            priceRefsSoFar;
    SimpleDoubleArray      periodPerfsSoFar;
    int                    iPeriodSoFar;

public:

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to Momentum) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    MomentumMC(const Momentum*     inst,
               const SimSeriesSP&    simSeries):
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
        nbPeriods(inst->timeAggregate->getDates().size()),
        periodPerfs(nbPeriods,0.0),
        nbAvgOutPerPeriod(nbPeriods,0),
        refLevels(nbAssets,0.0),
        priceRefs(nbAssets,0.0),
        sum(nbAssets,0.0),
        sumSoFar(nbAssets,0.0),
        refLevelsSoFar(nbAssets,0.0),
        priceRefsSoFar(nbAssets,0.0),
        periodPerfsSoFar(nbPeriods,0.0),
        iPeriodSoFar(0){
        static const string routine = "MomentumMC::MomentumMC";

        // Tie the pieces together :
        // We need the timeAgg to give a final value
        timeAgg = ITimeAggregateSP(inst->timeAggregate->getAggregate(simSeries->getAllDates().back(),
                                                                     this,
                                                                     &periodPerfs));

        bool isTrivial;
        const DateTimeArray& allSimDates = simSeries->getAllDates();
        const DateTimeArray& periodDates = inst->timeAggregate->getDates();
        periodMap = DateTime::createMapping(allSimDates,
                                            periodDates,
                                            isTrivial);

        // we check the first averageOut date is not after first cpn date
        int iPeriod = 0; 
        for(int iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            if(inst->averageOutDates[iStep] > periodDates[iPeriod]) {
                iPeriod++;
            }
            nbAvgOutPerPeriod[iPeriod]++;
        }
    }

    // for ITimeAggProductView
    const DateTime& getValueDate() const {
        return getToday();
    }
    const YieldCurve* getYieldCurve() const {
        return discount;
    }
    DateTime settles(const DateTime& aDate) const {
        return settlement->settles(aDate, 0);
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&        prices) {
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset;
        bool   doingPast = pathGen->doingPast();

        int iPeriod = iPeriodSoFar;
        periodPerfs = periodPerfsSoFar;
        sum = sumSoFar;
        for (iAsset=0;iAsset<nbAssets;iAsset++) {
            // Special treatment for average in for first period
            if (iPeriod==0) {
                refLevels[iAsset] = pathGen->refLevel(iAsset,0);
                priceRefs[iAsset] = pathGen->refLevel(iAsset,0);
            } else {
                refLevels[iAsset] = refLevelsSoFar[iAsset];
                priceRefs[iAsset] = priceRefsSoFar[iAsset];
            }
        }

        // Form the asset basket at each Period date first
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {

            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                sum[iAsset] += pathGen->Path(iAsset,0)[iStep];
            }
 
            // do the following at a coupon payment date
            // 1. compute momentum payment
            // 2. compute Bull Factor if it is bullish
            // 3. compute coupon
            // 4. reset 
            if (periodMap[iStep]==0) {  // true iff a coupon payment date

                //1.
                double assetPerf = (sum[0]/nbAvgOutPerPeriod[iPeriod])/refLevels[0];
				double momentumPayment = (assetPerf-1.0)*Maths::sign(assetPerf-1.0);
                for(iAsset=1; iAsset<nbAssets; iAsset++) {
                    assetPerf=(sum[iAsset]/nbAvgOutPerPeriod[iPeriod])/refLevels[iAsset];
                    momentumPayment=Maths::min(momentumPayment,(assetPerf-1.0)*Maths::sign(assetPerf-1.0));
                }
                momentumPayment *= inst->participation;

                //2.
                double bullFactor = 1.0;
                if (inst->isBullish) {
                    double bullSum = 0.0;
                    for(iAsset=0; iAsset<nbAssets; iAsset++) {  
                        bullSum += pathGen->Path(iAsset,0)[iStep]/priceRefs[iAsset];
                    }
                    bullFactor=Maths::max(1.0,bullSum/nbAssets);
                }
    
                //3. 
                double value = (momentumPayment-inst->momentumStrike)*bullFactor-inst->floorRate;
                if (inst->floorWithPreviousCoup && (iPeriod>0)) {
                    periodPerfs[iPeriod] = Maths::max(periodPerfs[iPeriod-1],value);
                } else {
                    periodPerfs[iPeriod] = Maths::max(0.0,value);    
                }

                //4. 
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    refLevels[iAsset] = sum[iAsset]/nbAvgOutPerPeriod[iPeriod];
                    sum[iAsset] = 0.0;
                    if (inst->isBullish) {
                        priceRefs[iAsset] = pathGen->Path(iAsset,0)[iStep];
                    }            
                }
                iPeriod++;
            } 
        }

        if (doingPast) {
        
            periodPerfsSoFar = periodPerfs;
            iPeriodSoFar = iPeriod;
            refLevelsSoFar = refLevels;
            sumSoFar = sum;
            priceRefsSoFar = priceRefs;

        }
        if (!doingPast || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            // Perform the aggregation (across periods)
            double payoff = timeAgg->aggregate();
            prices.add(inst->notional * payoff); 
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{}

    /** Use this opportunity to do any Implied driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustments */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{}

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int             iAsset) const {
        static const string routine = "MomentumMC::getVolInterp";
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel = 1.0;    //very lazy but LN shouldn't really be going this way! XXX

        // get hold of the future strike dates
        DateTimeArray periodDates=inst->timeAggregate->getDates();
        int numLivePeriods = periodDates.size() - iPeriodSoFar;
        if (numLivePeriods<=0) {
            throw ModelException(routine, "No future periods!?");
        }
        DateTimeArray livePeriodStartDates(numLivePeriods);
        for (int iLivePeriod = 0; iLivePeriod < numLivePeriods; iLivePeriod++){
            int iPeriod = iPeriodSoFar+iLivePeriod-1;
            livePeriodStartDates[iLivePeriod] = iPeriod<0?startDate:periodDates[iPeriod];
        }

        // same strike levels per cliquet (but may need to adjust first one)
        DoubleArray  strikes(numLivePeriods, interpLevel);
        if (!fwdStarting){
            // need to set first level to absolute strike - adjusted
            // additionally for any average out samples for this cliquet
            int iStep;
            const DateTimeArray& periodDates = inst->timeAggregate->getDates();
            const DateTime& thisPeriodStart = iPeriodSoFar==0?
                startDate:periodDates[iPeriodSoFar-1];
            // find first avg date of this cliquet
            for(iStep = 0; iStep < inst->averageOutDates.size() && 
                    inst->averageOutDates[iStep] <= thisPeriodStart; iStep++) {
                ; // empty
            }
            // then walk through counting past avg dates in this cliq
            int numRemaining = nbAvgOutPerPeriod[iPeriodSoFar];
            for(; iStep < inst->averageOutDates.size() && 
                    inst->averageOutDates[iStep] <= today; iStep++) {
                numRemaining--;
            }
            if (numRemaining<=0) {
                // something wrong!
                throw ModelException(routine, "INTERNAL ERROR : numRemaining is " + Format::toString(numRemaining));
            }
            // Can't set up refLevel earlier, 'cos need PathGen. First cliq has standard ref level
            double refLevel =  iPeriodSoFar==0 ? pathGen->refLevel(iAsset, 0) : refLevelsSoFar[iAsset];
            strikes[0] = (nbAvgOutPerPeriod[iPeriodSoFar] * refLevel * interpLevel
                          - sumSoFar[iAsset])/ numRemaining;
        }

        reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
                                                           livePeriodStartDates, 
                                                           lastSimDate,
                                                           strikes));
        return reqarr;
    }


    // Satisfy IHandlePaymentEvents interface
    void recordEvents(Control* control,
                      Results* results) {
        static const string method("BoostedNFBMC::recordEvents");
        try {
            // PAYMENT_DATES is a list of all dates on which payments may occur
            // including past and potential future dates.
            // For Momentum this means asking the TimeAggregate
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request && !request->getHasFinished()) {
                OutputRequestUtil::recordPaymentDates(control,results,timeAgg->getPaymentDates());
            }
        
            // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
            // and be supplied for all past cash flows, and any future ones
            // that are determined.
            // For Momentum this means asking the TimeAggregate
            request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && !request->getHasFinished()) {
                const CashFlowArray* kfl = timeAgg->getKnownFlows();
                if (kfl && kfl->size()>0) {
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            kfl); 
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Momentum::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new MomentumMC(this, simSeries);
}


CClassConstSP const Momentum::TYPE = CClass::registerClassLoadMethod(
    "Momentum", typeid(Momentum), Momentum::load);

// * for class loading (avoid having header file) */
bool MomentumLoad() {
    return (Momentum::TYPE != 0);
}

DRLIB_END_NAMESPACE









