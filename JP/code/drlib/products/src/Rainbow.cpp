//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Rainbow.cpp
//
//   Description : Port of DRMCRainbow - option on rainbow basket (weights
//                 are applied according to relative value, not per-asset)
//
//   Date        : May 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/MonteCarlo.hpp"
DRLIB_BEGIN_NAMESPACE

/** Rainbow product - option on a rainbow basket */
class Rainbow: public GenericNFBase, 
               virtual public IMCIntoProduct{
protected:
    /// fields ////////
    DoubleArray             weights;       // applied to sorted performances
    bool                    isCall;
    DoubleArray             strikes;       // in increasing order
    DateTimeArray           averageOutDates;

public:
    static CClassConstSP const TYPE;
    friend class RainbowMC;

    // validation
    void validatePop2Object(){
        static const string routine = "Rainbow::validatePop2Object";
        GenericNFBase::validatePop2Object();

        // check that we've got as many weights as there are assets 
        AssetUtil::checkWeights(weights, assets->NbAssets());
        
        // strikes must be increasing
        if (strikes.size()>2) {
            throw ModelException(routine,
                                 "No more than 2 strikes please : given " + 
                                 Format::toString(strikes.size()));
        } else if (strikes.size()==2) {
            if (strikes[0]>=strikes[1]) {
                throw ModelException(routine,
                                     "strikes must be increasing but : low strike " + 
                                     Format::toString(strikes[0]) + " >= high strike " +
                                     Format::toString(strikes[1]));
            }
        } else if (strikes.size()==0) {
            throw ModelException(routine,
                                 "No strikes given - require 1 or 2.");
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
    Rainbow(): GenericNFBase(TYPE) {} // for reflection
    Rainbow(const Rainbow& rhs); // not implemented
    Rainbow& operator=(const Rainbow& rhs); // not implemented

    static IObject* defaultRainbow(){
        return new Rainbow();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Rainbow, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultRainbow);
        FIELD(isCall,     "is it a call option");
        FIELD(strikes,    "strikes");
        FIELD(weights,    "Weights");
        FIELD(averageOutDates,  "averageOutDates");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
class RainbowMC: public IMCProduct,
                 virtual public IMCProductLN,
                 virtual public IMCQuickGreeks,
                 virtual public IMCQuickXGamma{
private:
    const Rainbow*           inst; // reference to original instrument
    int                      nbAssets;  // convenient
    int                      nbPath;  // num interp levels == nb strikes
    vector<DoubleArray>      sumOut;  // [iPath][iAsset] working area
    vector<DoubleArray>      sumOutSoFar; // preservation area (past)
    double                   CoP;     // +1 or -1 for call/put 
    // for quick greeks
    RepriceVanillaSP         vanillaReprice; // working area (price only)
    RepriceSpreadSP          spreadReprice;  // working area (price only)
    SubRepriceRainbowSP      rainbowReprice; // working area (price only)
    IndexedPerfListSP        indexPerf;   // working area (price only)

    enum {
        LOW_STRIKE = 0,
        HIGH_STRIKE
    };

    /** Returns the maximum delta scaling factor for each asset either
        viewed as part of the payoff (forComponents = false) or for
        each component that goes into the rainbow sum (it is the
        largest absolute derivative of the payoff/performance wrt
        each simulated asset) Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen,
                                bool                  forComponents) const{
        // first find the largest weight
        double factor = forComponents? 1.0:
            ISubReprice::Rainbow::largestAbsWeight(inst->weights);
        DoubleArray scaleFactors(nbAssets, factor);
        for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
            scaleFactors[iAsset] /= 
                futurePathGen->refLevel(iAsset, 0 /* path irrelevant*/);
        }
        return scaleFactors;
    }
public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to Rainbow) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    RainbowMC(const Rainbow*            inst,
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
        nbAssets(getNumAssets()),
        nbPath(inst->strikes.size()),
        sumOut(nbPath),
        sumOutSoFar(nbPath){
        for(int i=0; i<nbPath; i++) {
            sumOut[i] = DoubleArray(nbAssets);
            sumOutSoFar[i] = DoubleArray(nbAssets);
        }
        CoP = inst->isCall? 1 : -1;
    }

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
        IRepriceSP reprice;
        /* Note that this method gets called for pricing the past - so it's
           worth skipping this part because if indexPerf is non null the 
           payoff is slightly slower */
        if (mode){
            // create storage space for IndexedPerfList
            indexPerf = 
                IndexedPerfListSP(IndexedPerfList::create(getNumAssets()));
            // build up our rainbow SubReprice object
            rainbowReprice = 
                SubRepriceRainbowSP(new ISubReprice::Rainbow(mode, nbIter));
            if (nbPath > 1) {
                spreadReprice = RepriceSpreadSP(
                    new IReprice::Spread(rainbowReprice, mode, nbIter,
                                         inst->isCall,
                                         inst->notional, 
                                         inst->strikes[0], inst->strikes[1]));
                reprice = spreadReprice;
            } else {
                // just simple payoff
                vanillaReprice = RepriceVanillaSP(
                    new IReprice::Vanilla(rainbowReprice, mode, nbIter, 
                                          inst->notional));
                reprice = vanillaReprice;
            }
        }
        // then our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, reprice);
    }

    /** Create IMCPrices object for doing greek sens */
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
        // now sort out the vanilla/spread reprice object
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
        // now sort out the vanilla/spread reprice object
        IReprice* reprice = myPrices.getReprice().get();
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        // vanilla and spread use same params class
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
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
        // then update vanilla or spread reprice
        IReprice* reprice = myPrices.getReprice().get();
        // vanilla and spread use same params class ie TwoSidedGreekParamsSet
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen, false));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        reprice->setForXGamma(&params);
        // then do rainbow part - different scale factors though
        scaleFactors = maxDeltaFactors(futurePathGen, true);
        ISubReprice::Rainbow::XGammaParamsSet rbParams(scaleFactors,
                                                       futurePathGen,
                                                       getMultiFactors(),
                                                       sens);
        reprice->getRepriceForCmpt()->setForXGamma(&rbParams);
    }

    /* XXX With implied there's no need to distinguish the different paths/interp levels
       which here would save a sort. How can that be figured in? */
    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int    beginIdx = pathGen->begin(0); // 0 <=same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset;
 
        sumOut = sumOutSoFar;
        
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            for (int iPath=0; iPath<nbPath; iPath++) {
                for (int iStep=beginIdx; iStep<endIdx; iStep++) {
                    sumOut[iPath][iAsset] += 
                        pathGen->Path(iAsset, iPath)[iStep];
                }
            }
        }
        if (pathGen->doingPast()){ // preserve values
            sumOutSoFar = sumOut;
        }

        // Continue using sumOut as a working area - now for perfs
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            for (int iPath=0; iPath<nbPath; iPath++) {
                // calc average before passing data to rainbowReprice
                sumOut[iPath][iAsset] /= pathGen->refLevel(iAsset, iPath) *
                    inst->averageOutDates.size();
            }
        }
        double rainbowLo;
        if (!indexPerf){
            // calculate rainbow performance
            rainbowLo = InstrumentUtil::rainbowPerformance(inst->weights,
                                                           sumOut[LOW_STRIKE]);
        } else {
            // have to manually copy performances
            indexPerf->populate(sumOut[LOW_STRIKE]);
            // and then sort them
            indexPerf->sortByPerf();
            // save them
            rainbowReprice->store(*indexPerf, this, pathGen);
            // calculate rainbow performance
            rainbowLo = indexPerf->weightedSum(inst->weights);
        }
        if (spreadReprice.get()){
            // spread object stores performance
            spreadReprice->store(rainbowLo);
        }
        double payoff = CoP*(rainbowLo - inst->strikes[LOW_STRIKE]);
        if (vanillaReprice.get()){
            // must store value before max with zero (and exclude notional)
            vanillaReprice->store(payoff);
        }
        // must use Maths::max rather than prices.maxWithZero for spread
        payoff = inst->notional * (inst->strikes.size()>1?
                                   Maths::max(payoff, 0.0):
                                   prices.maxWithZero(payoff));
        if (inst->strikes.size() > 1){
            double rainbowHi = 
                InstrumentUtil::rainbowPerformance(inst->weights, 
                                                   sumOut[HIGH_STRIKE]);
            // must use Maths::max rather than prices.maxWithZero for spread
            payoff -= inst->notional *
                Maths::max(CoP*(rainbowHi-inst->strikes[HIGH_STRIKE]), 0.0);
            payoff *= CoP; // spread reversed for puts
        }
        prices.add(payoff); 
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(inst->strikes.size());
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel;
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();

        for(int i=0; i<inst->strikes.size(); i++) {
            if (fwdStarting){
                interpLevel = inst->strikes[i];
            } else {
                /* not forward starting - some samples have fixed already
                   (this includes averaging in) */
                int numDates = inst->averageOutDates.size();
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                
                interpLevel = (numDates * inst->strikes[i] * pathGen->refLevel(iAsset, i)
                               - sumOutSoFar[i][iAsset])/ numRemaining;
            }
            reqarr[i] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
        }
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Rainbow::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new RainbowMC(this, simSeries);
}

CClassConstSP const Rainbow::TYPE = CClass::registerClassLoadMethod(
    "Rainbow", typeid(Rainbow), Rainbow::load);

// * for class loading (avoid having header file) */
bool RainbowLoad() {
    return (Rainbow::TYPE != 0);
}

DRLIB_END_NAMESPACE
