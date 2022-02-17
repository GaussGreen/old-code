//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BestOrWorst.cpp
//
//   Description : Port of DRMCNFBoW - Option on Best or Worst Of N Assets
//
//   Date        : May 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/MCProductIQuicks.hpp"

DRLIB_BEGIN_NAMESPACE

/** BestOrWorst product - Option on Best or Worst Of N Assets */
class BestOrWorst: public GenericNFBase, 
                   virtual public IMCIntoProduct{
private:
    /// fields ////////
    bool                    isCall;
    bool                    isBest;
    DoubleArray             strikes;       // in increasing order
    DateTimeArray           averageOutDates;

public:
    static CClassConstSP const TYPE;
    friend class BestOrWorstMC;

    // validation
    void validatePop2Object(){
        static const string routine("BestOrWorst::validatePop2Object");
        GenericNFBase::validatePop2Object();
  
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
    BestOrWorst(): GenericNFBase(TYPE) {} // for reflection
    BestOrWorst(const BestOrWorst& rhs); // not implemented
    BestOrWorst& operator=(const BestOrWorst& rhs); // not implemented

    static IObject* defaultBestOrWorst(){
        return new BestOrWorst();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BestOrWorst, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBestOrWorst);
        FIELD(isCall,     "is it a call option");
        FIELD(isBest,     "True for Best, False for Worst.");
        FIELD(strikes,    "strikes. 2 for Call or Put Spread, have to be increasing.");
        FIELD(averageOutDates,  "averageOutDates");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super BestOrWorst */
class BestOrWorstMC: public IMCProduct,
                     virtual public IMCProductLN,
                     virtual public IMCQuickGreeks,
                     virtual public IMCQuickXGamma{
private:
    const BestOrWorst*       inst;    // reference to original instrument
    int                      nbPath;  // num interp levels == nb strikes
    vector<DoubleArray>      sumOut;  // [iPath][iAsset] working area
    vector<DoubleArray>      sumSoFar; // preservation area (past)
    double                   CoP;     // +1 or -1 for call/put 
    int                      nbAssets; // for convenience
    // for quick greeks
    RepriceVanillaSP         vanillaReprice; // working area (price only)
    RepriceSpreadSP          spreadReprice;  // working area (price only)
    SubRepriceRainbowSP      rainbowReprice; // working area (price only)

    /** Returns the maximum delta scaling factor for each asset (it is the
        largest absolute derivative of the payoff wrt each simulated asset)
        Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen) const{
        // here the largest weight, by definition, is 100%
        DoubleArray maxDeltaFactor(getNumAssets(), 1.0);
        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            maxDeltaFactor[i] /=
                futurePathGen->refLevel(i, 0 /* path irrelevant*/);
        } 
        return maxDeltaFactor;
    }
public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to BestOrWorst) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    BestOrWorstMC(const BestOrWorst*            inst,
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
        nbPath(inst->strikes.size()),
        sumOut(nbPath),
        sumSoFar(nbPath),
        nbAssets(getNumAssets()) {
        for(int i=0; i<nbPath; i++) {
            sumOut[i] = DoubleArray(nbAssets);
            sumSoFar[i] = DoubleArray(nbAssets);
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
        /* build up our rainbow SubReprice object (best of is just special
           type of rainbow) */
        if (mode){
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
        ISubRepriceSP subReprice = reprice->getRepriceForCmpt();
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
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
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
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        // vanilla and spread use same params class ie TwoSidedGreekParamsSet
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        reprice->setForXGamma(&params);
        // then do rainbow part 
        ISubReprice::Rainbow::XGammaParamsSet rbParams(scaleFactors,
                                                       futurePathGen,
                                                       getMultiFactors(),
                                                       sens);
        reprice->getRepriceForCmpt()->setForXGamma(&rbParams);
    }

    /* XXX With implied there's no need to distinguish the different
       paths/interp levels. How can that be figured in? */
    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int    iAsset;
        int    beginIdx = pathGen->begin(0); // 0 <=same for all assets
        int    endIdx   = pathGen->end(0);
        int    iPath;
        double payoff;

        // Any past samples
        sumOut = sumSoFar;
        
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            for (iPath=0; iPath<nbPath; iPath++) {
                const double* path = pathGen->Path(iAsset, iPath);
                for (int iStep=beginIdx; iStep<endIdx; iStep++) {
                    sumOut[iPath][iAsset] += path[iStep];
                }
            }
        }
        if (pathGen->doingPast()){ // preserve values from past calc
            sumSoFar = sumOut;
        }

        // Continue using sumOut as a working area - now for perfs
        for (iAsset=0; iAsset<nbAssets; iAsset++) {
            for (iPath=0; iPath<nbPath; iPath++) {
                // calc average before passing data to rainbowReprice
                sumOut[iPath][iAsset] /= 
                    pathGen->refLevel(iAsset, iPath) * 
                    inst->averageOutDates.size();
            }
        }
        
        /* Finding best (or worst) does not need a sort - just a
           single pass thru */
        double BestOrWorst[2];
        int    pickedAsset[2];
        for (iPath=0; iPath<nbPath; iPath++) {
            pickedAsset[iPath] = 0;
            //BestOrWorst[iPath] = sumOut[iPath][0];
            if (inst->isBest) {
                for (iAsset=1; iAsset<nbAssets; iAsset++) {
                    if (sumOut[iPath][iAsset] > 
                        sumOut[iPath][pickedAsset[iPath]]) {
                        pickedAsset[iPath] = iAsset;
                    } 
                }
            } else {
                for (iAsset=1; iAsset<nbAssets; iAsset++) {
                    if (sumOut[iPath][iAsset] < 
                        sumOut[iPath][pickedAsset[iPath]]) {
                        pickedAsset[iPath] = iAsset;
                    } 
                }
            }
            BestOrWorst[iPath] = sumOut[iPath][pickedAsset[iPath]];
        }
        if (rainbowReprice.get()){
            rainbowReprice->store(sumOut[0], pickedAsset[0], this, pathGen);
        }
        if (spreadReprice.get()){
            // spread object stores performance
            spreadReprice->store(BestOrWorst[0]);
        }

        payoff = CoP*(BestOrWorst[0] - inst->strikes[0]);
        if (vanillaReprice.get()){
            // must store value before max with zero
            vanillaReprice->store(payoff);
        }
        // must use Maths::max rather than prices.maxWithZero for spread
        payoff = inst->notional * (nbPath > 1? 
                                   Maths::max(payoff, 0.0):
                                   prices.maxWithZero(payoff));
        if (nbPath > 1) {
            // If more than 1 strike we have a spread to price
            // must use Maths::max rather than prices.maxWithZero for spread
            payoff -= inst->notional* 
                Maths::max(CoP*(BestOrWorst[1] - inst->strikes[1]), 0.0);
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
                               - sumSoFar[i][iAsset])/ numRemaining;
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
IMCProduct* BestOrWorst::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new BestOrWorstMC(this, simSeries);
}

CClassConstSP const BestOrWorst::TYPE = CClass::registerClassLoadMethod(
    "BestOrWorst", typeid(BestOrWorst), BestOrWorst::load);

// * for class loading (avoid having header file) */
bool BestOrWorstLoad() {
    return (BestOrWorst::TYPE != 0);
}

DRLIB_END_NAMESPACE
