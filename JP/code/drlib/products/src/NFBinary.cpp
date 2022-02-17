//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : NFBinary.cpp
//
//   Description : Simple N Factor: payoff = Coupon if KI or NOT KO 
//                 
//
//   Date        : Feb 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/SVGenDiscFactor.hpp"


DRLIB_BEGIN_NAMESPACE

/** NFBinary product - Simple N Factor: payoff = Coupon if KI or not KO */
class NFBinary: public GenericNFBase, 
                virtual public IMCIntoProduct, 
                virtual public BarrierBreach::IEventHandler {
protected:
    /// fields ////////
    BarrierUnionSP barrierUnion;        //!< Currently supports BarrierSimple and BarrierSchedule
    DateTimeArray  sampleDates;         //!< Sample dates
    mutable bool   doneAmendedIsHit;    //!< Whether barrier objected has amended

public:
    static CClassConstSP const TYPE;
    friend class NFBinaryMC;
    friend class NFBinarySVMC;

    // validation
    void validatePop2Object(){
        static const string method = "NFBinary::validatePop2Object";
        GenericNFBase::validatePop2Object();
        try {
            if(!sampleDates.size()) {
                // In new type instruments, the monitoring dates are live
                // on the barrier and not on the instrument
                sampleDates = barrierUnion->getBarrier()->getBarrierSchedule()->getDateArray();
            }
            
            barrierUnion->getBarrier()->validate(assets->NbAssets());
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }
    
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return sampleDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // BarrierBreach::IEventHandler interface
    void getEvents(const BarrierBreach* breach, IModel* model, 
                   const DateTime& eventDate, EventResults* events) const;

private:
    NFBinary(): GenericNFBase(TYPE), doneAmendedIsHit(false) {} // for reflection
    NFBinary(const NFBinary& rhs);     // not implemented
    NFBinary& operator=(const NFBinary& rhs); // not implemented

    static IObject* defaultNFBinary(){
        return new NFBinary();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(NFBinary, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultNFBinary);
        FIELD(barrierUnion, "barrier Union");
        FIELD(sampleDates, "sample dates");
        FIELD(doneAmendedIsHit, "Whether isHitAmended has been set");
        FIELD_MAKE_TRANSIENT(doneAmendedIsHit);
    }
};

/* MC product class for NFBinary */
class NFBinaryMC: public IMCProduct,
                  virtual public IMCProductLN,
                  virtual public IMCProductImplied{
private:
    const NFBinary*           inst;
    BarrierSP                 barrier;   // non-const due to the adjustLN/Implied methods!...
    
    // From a past pricing; used for BARRIER_LEVEL request
    // Looking for ways to improve on this XXX
    DoubleArray               histRefLevels; // [nbAssets] 
public:

    NFBinaryMC(const NFBinary*           inst,
               const DateTimeArray&      monitorDates,
               const SimSeriesSP&        simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst),
        barrier(copy(inst->barrierUnion->getBarrier())),
        histRefLevels(getNumAssets(), 0.0){

        // create an interpolated barrier to match the sample dates
        barrier->createInterpBarrier(inst->valueDate, monitorDates);
    }

    /* payoff = Coupon  ;if KI or NOT KNO */
    /*        = 0       ;otherwise  */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        if (pathGen->doingPast()){ 
            int nbAssets = getNumAssets();
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                histRefLevels[iAsset] = pathGen->refLevel(iAsset, 0);
            }
        }
        prices.add(inst->notional * barrier->hitValue(pathGen));
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        barrier->adjustLN(pathGen,
                          this,
                          getMultiFactors(),
                          0);
    }

    /** Use this opportunity to do any Implied driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustments */
    void initialiseImplied(const  IMCPathGenerator*  pathGen)const{
        barrier->adjustImplied(pathGen,
                               this,
                               getMultiFactors(),
                               inst->discount.get(),
                               inst->instSettle.get(),
                               0);
    }

    // for the LogNormal path generator
    // XXX since there could be several interp levels needed we should
    // XXX defer construction of reqarr to a method on the barrier. TBD
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "NFBinary::getVolInterp";

        try
        {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            
            double interpLevel = barrier->getInterpLevel(pathGen, this, iAsset);

            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
            return reqarr;
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    /** For now am using this instead of implementing the IHandlePaymentEvents interface
        since I'd like to reuse the PAYMENT_DATES and KNOWN_CASHFLOWS features from
        the IMCProduct, but need to add BARRIER_LEVEL support. Should review. XXX */
    void recordExtraOutput(CControl*     control,
                           Results*      results,
                           const IMCPrices& prices) const {
        // BARRIER_LEVEL ...
        OutputRequest*  request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        const DateTime& today = getToday();
        const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
        // Only try to satisfy this request when have some past (so easy to get ref levels).
        if (request && (today >= lastRefDate)) {
            // This operates on a per-asset basis
            int nbAssets = getNumAssets();
            for(int iAsset=0; iAsset<nbAssets; iAsset++) {
                // Mostly delegate to Barrier class...
                // histRefLevels are to allow absolute barriers to be reported.
                BarrierLevelArraySP levels = barrier->reportLevels(today,
                                                                   histRefLevels[iAsset],
                                                                   iAsset);
                if (!levels->empty()) {
                    OutputRequestUtil::recordBarrierLevels(
                        control, results,
                        getMultiFactors()->assetGetTrueName(iAsset),
                        levels.get());
                }
                
            }
        }            
    }

    // go get the events based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen, 
                           EventResults* events,
                           const DateTime& eventDate) {
        int numAssets = getNumAssets();
        StringArraySP assetNames(new StringArray(numAssets));
        for (int i = 0; i < numAssets; i++) {
            (*assetNames)[i] = getMultiFactors()->assetGetTrueName(i);
        }

		// override the barrier today if necessary
		// otherwise barrier today will be the legal one
		barrier->overrideEventBarrier(eventDate);

        barrier->getEvents(pathGen, events, true, "Multi-asset barrier", 
                            assetNames.get());
    }
};


//////////////////////////////////////////////////////////////////////////

/* MC product class for NFBinarySV */
class NFBinarySVMC: public MCProductClient,
                    virtual public IMCProductLN {
private:
    const NFBinary*                   inst;             //!< Instrument
    Barrier*                          barrier;          //!< Old style barrier

    // State variables and generators
    SVGenBarrierHVStructSP            barrierStructGen; //!< Barrier generator
    SVGenBarrierHVStruct::StateVarSP  barrierStructSV;  //!< Barrier state variable
    IRefLevel::IStateVarGenSP         refLevelGen;      //!< RefLevel generator
    IRefLevel::IStateVarSP            refLevelSV;       //!< RefLevel state variable
    SVGenDiscFactorSP                    dfGen;            //!< Discount factors generator
    SVDiscFactorSP         dfSV;             //!< DF state variable

public:
    IRefLevel::IStateVarGenSP getRefLevelGen() const {
        return refLevelGen;
    }

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        static const string routine = "NFBinarySVMC::collectStateVars";
        try{
            svCollector->append(refLevelGen.get());         // reference level
            svCollector->append(barrierStructGen.get());    // barrier structure
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
        static const string routine = "NFBinarySVMC::pathGenUpdated";
        try{
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            barrierStructSV = barrierStructGen->getHitValueStructSV(
                barrierStructSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    NFBinarySVMC(const NFBinary*      inst,
                 const DateTimeArray& monitorDates,
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
        barrier(inst->barrierUnion->getBarrier()),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())) {
        
        // Convert barrier to barrier per asset generator
        SVGenBarrierHVSP barrierGen = barrier->convertBarrier(
            monitorDates, monitorDates.back(), getToday(), refLevelGen, getMultiFactors());

        // Create barrier for whole structure from barrier per asset
        barrierStructGen = SVGenBarrierHVStructSP(new SVGenBarrierHVStruct(barrierGen));
    }

    /* payoff = Coupon  ;if KI or NOT KNO */
    /*        = 0       ;otherwise  */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        
        // Get hit value and discount
        double hitValue = barrierStructSV->hitValue();

        if (doingPast() && !hasFuture()) {
            // If we're in a known state, we record known flows on
            // their known date (so no discounting).
            if (!paymentDate.empty()) {
                knownCashFlows->addFlow(paymentDate,
                                        inst->notional * hitValue);
            }
        }

        hitValue *= dfSV->firstDF();

        prices.add(inst->notional * hitValue);
    }


    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string routine = "NFBinarySVMC::initialiseLN";
        throw ModelException(routine, "Methodology not supported");
    }


    // for the LogNormal path generator
    // XXX since there could be several interp levels needed we should
    // XXX defer construction of reqarr to a method on the barrier. TBD
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string method = "NFBinarySVMC::getVolInterp";

        try {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            const DateTime& today = getToday();
            const DateTime& startDate = refLevelGen->getAllDates().front();
            const DateTime& lastSimDate = getSimSeries()->getLastDate();
            bool  fwdStarting = startDate.isGreater(today);

            // Interpolate at the last date of the barrier level
            double barrierLevel = barrierStructGen->getBarrierData()->
                getAssetBarrier(iAsset)->barrierLevels.back();
            double interpLevel  = fwdStarting ? 
                barrierLevel: barrierLevel * refLevelSV->refLevel(iAsset);

            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));

            return reqarr;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }


    /** For now am using this instead of implementing the IHandlePaymentEvents interface
        since I'd like to reuse the PAYMENT_DATES and KNOWN_CASHFLOWS features from
        the IMCProduct, but need to add BARRIER_LEVEL support. Should review. XXX */
    void recordExtraOutput(CControl*     control,
                           Results*      results,
                           const IMCPrices& prices) const {
        
        static const string method = "NFBinarySVMC::recordExtraOutput";
        
        try {
            // BARRIER_LEVEL is delegated to the state variable
            barrierStructSV->recordBarrierLevels(control, results, getMultiFactors());
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // go get the events based on this path
    virtual void getEvents(const IMCPathGenerator*  pathGen, 
                           EventResults* events,
                           const DateTime& eventDate) {
        throw ModelException("NFBinarySVMC::getEvents", 
                "Event handling not implemented for state variable NFB");
    }
};


//////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* NFBinary::createProduct(const MonteCarlo* model) const {
    static const string method = "NFBinary::createProduct";

    try {
        int nbAssets = assets->NbAssets();
        
        // Validate the barrier based on the model
        Barrier* barrier = barrierUnion->getBarrier(); // ref only
        if (!barrier) {
            throw ModelException(method, "Failed to locate barrier instance!");
        }
        barrier->validate(model, nbAssets);

        // Establish monitoring dates
        smartPtr<IMCPathConfig> pathConfig = model->getPathConfig();

        // XXX Clashes with past values validation....
        DateTimeArray monitorDates = barrier->getFutureMonitoringDates(getValueDate(),
                                                                       sampleDates,
                                                                       pathConfig);
        // Create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(nbAssets));
        simSeries->addDates(monitorDates);
    
        if(model->stateVarUsed()) {
            // State variables
            if(!doneAmendedIsHit) {
                // Amend the barrier schedule isHit values:
                // 1) Create temporary product
                // 2) Create past path generator
                // 3) Get ref level state variable
                // 4) Amend is hit flag on Barrier
                refCountPtr<NFBinarySVMC> tmpProd(new NFBinarySVMC(this, monitorDates, simSeries));
                IRefLevel::IStateVarGenSP refLevelGen = tmpProd->getRefLevelGen();
                MCPathGeneratorSP pastPathGen(pathConfig->pastPathGenerator(tmpProd.get()));
                IRefLevel::IStateVarSP oldStateVar;
                IRefLevel::IStateVarSP refLevelSV = refLevelGen->getRefLevelSV(oldStateVar, pastPathGen.get());
                DateTimeArrayArray sampleDatesPerAsset(nbAssets, monitorDates);
                barrierUnion->getBarrier()->amendIsHit(
                    tmpProd->getMultiFactors(), refLevelSV, sampleDatesPerAsset, valueDate);
                
                doneAmendedIsHit = true;
            }

            return new NFBinarySVMC(this, monitorDates, simSeries);
        } else {
            // Otherwise, use old methodology
            return new NFBinaryMC(this, monitorDates, simSeries);
        }
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

// implementation of Events interface
void NFBinary::getEvents(const BarrierBreach* breach,
                         IModel* model, 
                         const DateTime& eventDate,
                         EventResults* events) const {
    static const string method = "NFBinary::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP pastPathGenerator(
                    mc->getPathConfig()->pastPathGenerator(prod.get()));
            prod->getEvents(pastPathGenerator.get(), events, eventDate);
        } else {
            throw ModelException(method, 
                    "Internal error - expected Monte Carlo model for NFB pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const NFBinary::TYPE = CClass::registerClassLoadMethod(
    "NFBinary", typeid(NFBinary), NFBinary::load);

// * for class loading (avoid having header file) */
bool NFBinaryLoad() {
    return (NFBinary::TYPE != 0);
}

DRLIB_END_NAMESPACE
