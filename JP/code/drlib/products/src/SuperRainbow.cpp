//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SuperRainbow.cpp
//
//   Description : Test product for MC framework
//
//   Date        : May 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Performance.hpp"
#include "edginc/Format.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCProductIQuicks.hpp"

#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** SuperRainbow product - option based upon the performance of n assets
    sorted by their performance */
class SuperRainbow: public GenericNFBase, 
                    virtual public IMCIntoProduct{
protected:
    /// fields ////////
    IPerformanceSP          performance;   // gives performance
    DoubleArray             weights;       // applied to sorted performances
    bool                    isCall;
    double                  overallStrike;
    bool                    allowNonNormalWeights;

public:
    static CClassConstSP const TYPE;
    friend class SuperRainbowMC;
    friend class SuperRainbowSVMC;

    // validation
    void validatePop2Object(){
        GenericNFBase::validatePop2Object();
        if (allowNonNormalWeights) {
            // Don't require sum to 100%
            if (weights.size() != assets->NbAssets()) {
                throw ModelException("SuperRainbowMC::validatePop2Object",
                                     "Mismatch between number of assets (" + 
                                     Format::toString(assets->NbAssets()) +
                                     ") and number of weights (" + 
                                     Format::toString(weights.size()) + ")");
            }
        } else {
            // check that we've got as many weights as there are assets 
            AssetUtil::checkWeights(weights, assets->NbAssets());
        }
    }
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        throw ModelException("SuperRainbow::getSamplingDates",
                    "Internal error, shouldn't call this function since "
                    "this product has sampling dates per asset");
    }

    // most products have a single set of sampling dates so this defaults to true
    // some products e.g. SuperRainbow have a set per asset
    bool assetsShareSamplingDates() const {
        return false;
    }

    const DateTimeArrayArray* getSamplingDatesPerAsset() const {
        // a bit clunky - basically copies createProduct
        SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                    one */
        performance->recordDates(simSeries.get()); // then add dates to it

        DateTimeArrayArray* samplingDates = 
                new DateTimeArrayArray(assets->NbAssets());
        for (int i = 0; i < assets->NbAssets(); i++) {
            (*samplingDates)[i] = DateTime::merge(refLevel->getAllDates(), 
                                                  simSeries->getDates(i));
        }
        return samplingDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    SuperRainbow(): GenericNFBase(TYPE), allowNonNormalWeights(false) {} // for reflection
    SuperRainbow(const SuperRainbow& rhs); // not implemented
    SuperRainbow& operator=(const SuperRainbow& rhs); // not implemented

    static IObject* defaultSuperRainbow(){
        return new SuperRainbow();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SuperRainbow, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultSuperRainbow);
        FIELD(performance,             "Performance measurement");
        FIELD(isCall,           "is it a call option");
        FIELD(overallStrike,    "OverallStrike");
        FIELD(weights,          "Weights");
        FIELD(allowNonNormalWeights,  "True=>weights not required "
                     "to sum to 100%");
        FIELD_MAKE_OPTIONAL(allowNonNormalWeights);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
class SuperRainbowMC : public IMCProduct,
                       virtual public IMCProductLN,
                       virtual public IMCQuickGreeks,
                       virtual public IMCQuickXGamma{
private:
    const SuperRainbow*      inst; // reference to original instrument
    IPerformance::IMCPerfSP  performance; // calculates performaces for us
    DoubleArray              sortedPerf;  // working area
    SubRepriceRainbowSP      rainbowReprice; // working area (price only)
    RepriceVanillaSP         vanillaReprice; // working area (price only)
    IndexedPerfListSP        indexPerf;   // working area (price only)
public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to SuperRainbow) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
        SuperRainbowMC(const SuperRainbow*            inst,
                       const SimSeriesSP&             simSeries,
                       const IPerformance*            perf):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()), inst(inst), 
        performance(perf->getMCPerf(inst->valueDate, simSeries.get())),
        sortedPerf(getNumAssets()){
        if (performance->numPerformances() != getNumAssets()){
            throw ModelException("SuperRainbowMC",
                                 "Mismatch between number of assets and "
                                 "number of performances");
        }
    }
private:
    /** Returns the maximum delta scaling factor for each asset either
        viewed as part of the payoff (forComponents = false) or for
        each component that goes into the rainbow sum (it is the
        largest absolute derivative of the payoff/performance wrt
        each simulated asset). Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen,
                                bool                  forComponents) const{
        // first find the largest weight
        double factor = forComponents? 1.0:
            ISubReprice::Rainbow::largestAbsWeight(inst->weights);
        DoubleArray scaleFactors(inst->assets->NbAssets(), factor);
        for (int iAsset = 0; iAsset < inst->assets->NbAssets(); iAsset++) {
            scaleFactors[iAsset] *=
                performance->maxDeltaScalingFactor(futurePathGen, iAsset);
        }
        return scaleFactors;
    }
public: 
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Create IMCPrices object for first pricing call. The mode 
        parameter indicates what is needed (eg quick greeks, quick X gamma) */
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode) { // bitwise parameter
        /* Note that this method gets called for pricing the past - so it's
           worth skipping this part because if indexPerf is non null the 
           payoff is slightly slower */
        if (mode){
            // create storage space for IndexedPerfList
            indexPerf = 
                IndexedPerfListSP(IndexedPerfList::create(getNumAssets()));
            // built up our rainbow SubReprice object
            rainbowReprice = 
                SubRepriceRainbowSP(new ISubReprice::Rainbow(mode, nbIter));
            // then one to desribe our payoff
            vanillaReprice = 
                RepriceVanillaSP(new IReprice::Vanilla(rainbowReprice,
                                                       mode, nbIter, 
                                                       inst->notional));
        }
        // finally our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, vanillaReprice);
    }   

    /** Sets IMCPrices object for doing greek sens */
    virtual void setPricesForGreek(IMCPrices*               untweakedPrices,
                                   const IPathGenerator* futurePathGen,
                                   const Sensitivity*    sens) {
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        // Can we actually do quick greeks for this sensitivity?
        if (!IQuickGreeks::doQuickGreeks(sens)){
            myPrices.setMode(IQuickGreeks::NONE);
        } else {
            myPrices.setMode(IQuickGreeks::FIRST_ORDER);
        } 
        IReprice* reprice = myPrices.getReprice().get();
        reprice->setForOneSidedGreek(0); // regardless of mode
        reprice->getRepriceForCmpt()->setForOneSidedGreek(0); // ditto
    }

    /** Set IMCPrices object for doing two sided greek sens */
    virtual void setPricesForTwoSidedGreek(
        IMCPrices*                 untweakedPrices,
        const IPathGenerator*   futurePathGen,
        const ScalarShiftArray& sens){
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        // tell top IMCPrices object what's going on
        myPrices.setMode(IQuickGreeks::SECOND_ORDER);
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        IReprice* reprice = myPrices.getReprice().get();
        reprice->setForTwoSidedGreek(&params);
        // rainbow part is trivial
        reprice->getRepriceForCmpt()->setForTwoSidedGreek(0);
    }

    /** Sets IMCPrices object for doing quick x gamma. The IMCPrices
        object needs to be able to skip over paths that have zero
        x gamma. The sens array contains the shifts defining what range
        the cross gamma will be calculated over */
    virtual void setPricesForXGamma(
        IMCPrices*                  untweakedPrices,
        const IPathGenerator*    futurePathGen,
        const ScalarShiftArray&  sens) {
        static const string method("SuperRainbowMC::setPricesForXGamma");
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        myPrices.setMode(IQuickGreeks::CROSS_GAMMA);
        IReprice* reprice = myPrices.getReprice().get();
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        reprice->setForXGamma(&params);
        // then do rainbow part - need different scaleFactors though
        scaleFactors = maxDeltaFactors(futurePathGen, true);
        ISubReprice::Rainbow::XGammaParamsSet rbParams(scaleFactors,
                                                       futurePathGen,
                                                       getMultiFactors(),
                                                       sens);
        reprice->getRepriceForCmpt()->setForXGamma(&rbParams);
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        // calculate performance of each asset over the interval
        performance->calcPerf(pathGen, 0 /* iPath */, sortedPerf);
        bool recording = indexPerf.get() == 0? false: true;
        double rainbowPerf;
        if (!recording){
            // calculate rainbow performance
            rainbowPerf = 
                InstrumentUtil::rainbowPerformance(inst->weights, sortedPerf);
        } else {
            // have to manually copy performances
            indexPerf->populate(sortedPerf);
            // and then sort them
            indexPerf->sortByPerf();
            // save them
            rainbowReprice->store(*indexPerf, this, pathGen);
            // calculate rainbow performance
            rainbowPerf = indexPerf->weightedSum(inst->weights);
        }
        // apply final option formula
        rainbowPerf -= inst->overallStrike;
        if (!inst->isCall){
            rainbowPerf = -rainbowPerf;
        }
        if (recording){
            // save pre max price in vanillaReprice
            vanillaReprice->store(rainbowPerf);
        }
        // finally store true price
        prices.add(inst->notional * prices.maxWithZero(rainbowPerf));
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(
            performance->getVolInterp(this, pathGen, iAsset));
        return reqarr;
    }
};

////////////////////////////////////////////////////////////////////////////////



/* MC product class for super rainbow */
class SuperRainbowSVMC : public MCProductClient,
                         virtual public IMCProductLN,
                         virtual public IMCQuickGreeks,
                         virtual public IMCQuickXGamma{
private:
    const SuperRainbow*      inst; // reference to original instrument
    DoubleArray              sortedPerf;        // working area
    SubRepriceRainbowSP      rainbowReprice;    // working area (price only)
    RepriceVanillaSP         vanillaReprice;    // working area (price only)
    IndexedPerfListSP        indexPerf;         // working area (price only)

    // state vars
    IRefLevel::IStateVarGenSP   refLevelGen;      //!< generator for ref level
    IPerformanceSVSP            performanceGen;   //!< generator for performance state variable
    SVGenDiscFactorSP              dfGen;            //!< generator for discount factors
    
    IRefLevel::IStateVarSP      refLevelSV;       //!< ref level state variable
    IPerformanceSV::IStateVarSP performance;      //!< performance state variable
    SVDiscFactorSP   dfSV;             //!< df state variable

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get()); // and a spot one
        svCollector->append(performanceGen.get()); // and a spot one
        svCollector->append(dfGen.get()); // and a DiscFactor one
    }

                                                
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "VanillaMCSV::pathGenUpdated";

        try {
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            performance = performanceGen->getPerfStateVar(performance, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to SuperRainbow) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
        SuperRainbowSVMC(const SuperRainbow*            inst,
                       const SimSeriesSP&             simSeries,
                       const IPerformance*            perf):
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()), 
        inst(inst), 
        sortedPerf(getNumAssets()) {
        
            // Create state variable generators
            refLevelGen = IRefLevel::IStateVarGenSP(
                refLevel->createStateVarGen(getMultiFactors(), inst->valueDate));
            performanceGen = inst->performance->getStateVarGen(refLevelGen, getNumAssets());
            dfGen = SVGenDiscFactorSP(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                inst->instSettle, simSeries->getLastDate()));
    }
private:
    /** Returns the maximum delta scaling factor for each asset either
        viewed as part of the payoff (forComponents = false) or for
        each component that goes into the rainbow sum (it is the
        largest absolute derivative of the payoff/performance wrt
        each simulated asset). Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen,
                                bool                  forComponents) const{
        // first find the largest weight
        double factor = forComponents? 1.0:
            ISubReprice::Rainbow::largestAbsWeight(inst->weights);
        DoubleArray scaleFactors(inst->assets->NbAssets(), factor);
        for (int iAsset = 0; iAsset < inst->assets->NbAssets(); iAsset++) {
            scaleFactors[iAsset] *=
                performance->maxDeltaScalingFactor(iAsset);
        }
        return scaleFactors;
    }
public: 
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        static const string method = "SuperRainbowSVMC::initialiseLN";

        throw ModelException(method, "Deprecated method");
    }

    /** Create IMCPrices object for first pricing call. The mode 
        parameter indicates what is needed (eg quick greeks, quick X gamma) */
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode) { // bitwise parameter
        /* Note that this method gets called for pricing the past - so it's
           worth skipping this part because if indexPerf is non null the 
           payoff is slightly slower */
        if (mode){
            // create storage space for IndexedPerfList
            indexPerf = 
                IndexedPerfListSP(IndexedPerfList::create(getNumAssets()));
            // built up our rainbow SubReprice object
            rainbowReprice = 
                SubRepriceRainbowSP(new ISubReprice::Rainbow(mode, nbIter));
            // then one to desribe our payoff
            vanillaReprice = 
                RepriceVanillaSP(new IReprice::Vanilla(rainbowReprice,
                                                       mode, nbIter, 
                                                       inst->notional));
        }
        // finally our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, vanillaReprice);
    }   

    /** Sets IMCPrices object for doing greek sens */
    virtual void setPricesForGreek(IMCPrices*               untweakedPrices,
                                   const IPathGenerator* futurePathGen,
                                   const Sensitivity*    sens) {
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        // Can we actually do quick greeks for this sensitivity?
        if (!IQuickGreeks::doQuickGreeks(sens)){
            myPrices.setMode(IQuickGreeks::NONE);
        } else {
            myPrices.setMode(IQuickGreeks::FIRST_ORDER);
        } 
        IReprice* reprice = myPrices.getReprice().get();
        reprice->setForOneSidedGreek(0); // regardless of mode
        reprice->getRepriceForCmpt()->setForOneSidedGreek(0); // ditto
    }

    /** Set IMCPrices object for doing two sided greek sens */
    virtual void setPricesForTwoSidedGreek(
        IMCPrices*                 untweakedPrices,
        const IPathGenerator*   futurePathGen,
        const ScalarShiftArray& sens){
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        // tell top IMCPrices object what's going on
        myPrices.setMode(IQuickGreeks::SECOND_ORDER);
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        IReprice* reprice = myPrices.getReprice().get();
        reprice->setForTwoSidedGreek(&params);
        // rainbow part is trivial
        reprice->getRepriceForCmpt()->setForTwoSidedGreek(0);
    }

    /** Sets IMCPrices object for doing quick x gamma. The IMCPrices
        object needs to be able to skip over paths that have zero
        x gamma. The sens array contains the shifts defining what range
        the cross gamma will be calculated over */
    virtual void setPricesForXGamma(
        IMCPrices*                  untweakedPrices,
        const IPathGenerator*    futurePathGen,
        const ScalarShiftArray&  sens) {
        static const string method("SuperRainbowMC::setPricesForXGamma");
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        myPrices.setMode(IQuickGreeks::CROSS_GAMMA);
        IReprice* reprice = myPrices.getReprice().get();
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        reprice->setForXGamma(&params);
        // then do rainbow part - need different scaleFactors though
        scaleFactors = maxDeltaFactors(futurePathGen, true);
        ISubReprice::Rainbow::XGammaParamsSet rbParams(scaleFactors,
                                                       futurePathGen,
                                                       getMultiFactors(),
                                                       sens);
        reprice->getRepriceForCmpt()->setForXGamma(&rbParams);
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        // calculate performance of each asset over the interval
        performance->calcPerf(sortedPerf);
        bool recording = indexPerf.get() == 0? false: true;
        double rainbowPerf;
        if (!recording){
            // calculate rainbow performance
            rainbowPerf = 
                InstrumentUtil::rainbowPerformance(inst->weights, sortedPerf);
        } else {
            // have to manually copy performances
            indexPerf->populate(sortedPerf);
            // and then sort them
            indexPerf->sortByPerf();
            // save them
            rainbowReprice->store(*indexPerf, this, pathGen);
            // calculate rainbow performance
            rainbowPerf = indexPerf->weightedSum(inst->weights);
        }
        // apply final option formula
        rainbowPerf -= inst->overallStrike;
        if (!inst->isCall){
            rainbowPerf = -rainbowPerf;
        }
        
        if (doingPast() && !hasFuture()) {
            // If we're in a known state, we record known flows on
            // their known date (so no discounting).
            if (!paymentDate.empty()) {
                knownCashFlows->addFlow(paymentDate,
                                        Maths::max(inst->notional * rainbowPerf, 0.0));
            }
        }

        // Discount
        rainbowPerf *= dfSV->firstDF();

        if (recording){
            // save pre max price in vanillaReprice
            vanillaReprice->store(rainbowPerf);
        }
        // finally store true price
        prices.add(inst->notional * prices.maxWithZero(rainbowPerf));
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(
            performance->getVolInterp(this, iAsset));
        return reqarr;
    }
};


////////////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* SuperRainbow::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    performance->recordDates(simSeries.get()); // then add dates to it
    
    if (model->stateVarUsed()) {
        // State variables
        return new SuperRainbowSVMC(this, simSeries, performance.get());
    }
    // otherwise, use old methodology
    return new SuperRainbowMC(this, simSeries, performance.get());
}

CClassConstSP const SuperRainbow::TYPE = CClass::registerClassLoadMethod(
    "SuperRainbow", typeid(SuperRainbow), SuperRainbow::load);

// * for class loading (avoid having header file) */
bool SuperRainbowLoad() {
    return (SuperRainbow::TYPE != 0);
}

DRLIB_END_NAMESPACE
