//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BaskAv.cpp
//
//   Description : Port of DRMCBaskAv from EDG - option on basket with average-in
//                 and average out
//
//   Date        : April 2002
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
#include "edginc/Format.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/XCB.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCProductIQuicks.hpp"

DRLIB_BEGIN_NAMESPACE

/** BaskAv product - option based upon the performance of n assets
    sorted by their performance */
class BaskAv: public GenericNFBase, 
              virtual public IMCIntoProduct,
              public ITaxableInst::Basic {
protected:
    /// fields ////////
    string                  weightType;
    string                  smileType; // values to match XCB
    DoubleArray             weights;
    bool                    isCall;
    double                  strike;
    DateTimeArray           averageOutDates;
    bool                    isRefAvgOut;      // optional

public:
    static CClassConstSP const TYPE;
    friend class BaskAvMC;
    friend class BaskAvSVMC;
    friend class BaskAvSLMC;
    
    // validation
    void validatePop2Object(){
        static const string routine("BaskAv::validatePop2Object");
        GenericNFBase::validatePop2Object();
    
        // check that we've got as many weights as there are assets 
        // and if % weights that they sum to 100%
        if (weightType=="P") {
            AssetUtil::checkWeights(weights, assets->NbAssets());
        } else if (weightType=="U") {
            if (assets->NbAssets() != weights.size()){
                throw ModelException(routine,
                                     "Different number of assets ("+
                                     Format::toString(assets->NbAssets())+")"
                                     " to weights ("+
                                     Format::toString(weights.size())+")");
            }
        } else {
            throw ModelException(routine, "BaskAv: weightType must be U or P,"
                                 " but " + weightType + " given");
        }
        if (Maths::isNegative(strike)){
            throw ModelException(routine, "strike ("+
                                 Format::toString(strike)+") is negative");
        }
        // validate dates are not empty - order is handled by SimSeries
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates supplied!");
        }

        // This combination has potentially infinite payout
        if (isRefAvgOut && !isCall) {
            throw ModelException(routine, 
                                 "Combination of Put with isRefAvgOut True is forbidden");
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

    /** for ITaxableInst::Basic */
    const DateTime getFinalPaymentDate() const {
        return instSettle->settles(averageOutDates.back(), NULL);
    }

private:
    BaskAv(): GenericNFBase(TYPE), isRefAvgOut(false) {} // for reflection
    BaskAv(const BaskAv& rhs); // not implemented
    BaskAv& operator=(const BaskAv& rhs); // not implemented

    static IObject* defaultBaskAv(){
        return new BaskAv();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BaskAv, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBaskAv);
        FIELD(weightType,       "Unit(U) or Percentage(P) basket");
        FIELD(smileType,        "Type of smile as per XCB");
        FIELD(weights,          "Weights");
        FIELD(isCall,           "is it a call option");
        FIELD(strike,           "strike");
        FIELD(averageOutDates,  "averageOutDates");
        FIELD(isRefAvgOut,      "isRefAvgOut");
        FIELD_MAKE_OPTIONAL(isRefAvgOut);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for BaskAv - with isRefAvgOut false */
class BaskAvMC : public IMCProduct,
                 virtual public IMCProductLN,
                 virtual public IMCQuickGreeks,
                 virtual public IMCQuickXGamma{
private:
    const BaskAv*            inst;        // reference to original instrument
    int                      nbAssets;     // convenient
    DoubleArray              sumOut;       // working area - here to save alloc
    DoubleArray              sumOutSoFar;  // preservation area (past)
    bool                     isUnitBasket;
    double                   baskSumOut;
    RepriceVanillaSP         vanillaReprice; // working area (price only)

protected:
    /** Override default to allow us to switch off if we're forward 
        starting unit basket */
    virtual IMCQuickXGamma* quickXGammaSupported(){
        if (isUnitBasket){
            DateTimeArray futureRefLevelDates(inst->refLevel->
                                              getFutureDates(inst->valueDate));
            if (!futureRefLevelDates.empty()){
                // basket level is non linear in spot
                return 0;
            }
        }
        if (inst->isRefAvgOut) {
            // xgamma is not so easy
            return 0;
        }
        return this;
    }
public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    BaskAvMC(const BaskAv*            inst,
             const SimSeriesSP&       simSeries):
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
        sumOut(nbAssets),
        sumOutSoFar(nbAssets),
        baskSumOut(0.0) {
        isUnitBasket = (inst->weightType == "U");
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Create IMCPrices object for first pricing call. The storingMode 
        parameter is true if this IMCPrices object will be used for future
        tweaked pricing runs. */
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode) {
        // describe our payoff - v simple
        if (mode>0){
            if (inst->isRefAvgOut) {
                throw ModelException("BaskAvMC::createOrigPrices", 
                                     "Quick Greeks not yet supported with isRefAvgOut True.");
            }
            vanillaReprice = RepriceVanillaSP(
                new IReprice::Vanilla(mode, nbIter, inst->notional));
       }
        // then our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, vanillaReprice);
    }   

    /** Returns the maximum delta scaling factor for each asset (it is the
        largest absolute derivative of the payoff wrt each simulated asset)
        Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen) const{
        DoubleArray maxDeltaFactor(getNumAssets());
        double baskRef = 0.0;
        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            if (isUnitBasket) {
                baskRef += inst->weights[i] * 
                    futurePathGen->refLevel(i, 0/*iPath*/);
            } else {
                // basket of performances
                maxDeltaFactor[i] = fabs(inst->weights[i]) /
                    futurePathGen->refLevel(i, 0 /* path irrelevant*/);
            }
        } 
        if (isUnitBasket) {
           for (int i = 0; i < maxDeltaFactor.size(); i++) {
               maxDeltaFactor[i] = fabs(inst->weights[i]) / baskRef;
           }
        }
        return maxDeltaFactor;
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
        myPrices.getReprice()->setForOneSidedGreek(0); // regardless of mode
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
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        myPrices.getReprice()->setForTwoSidedGreek(&params);
    }

    /** Create IMCPrices object for doing quick x gamma. The IMCPrices
        object needs to be able to skip over paths that have zero
        gamma */
    virtual void setPricesForXGamma(
        IMCPrices*                 untweakedPrices,
        const IPathGenerator*   futurePathGen,
        const ScalarShiftArray& sens) {
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        myPrices.setMode(MCPricesGeneral::CROSS_GAMMA);
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        myPrices.getReprice()->setForXGamma(&params);
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("BaskAvMC::payoff");
        try {
            int    beginIdx = pathGen->begin(0); // same for all assets
            int    endIdx   = pathGen->end(0);
            double baskRef = 0;

            baskSumOut = 0;
            sumOut = sumOutSoFar;

            for (int iAsset=0; iAsset<nbAssets; iAsset++) {
                for (int iStep=beginIdx; iStep<endIdx; iStep++) {
                    sumOut[iAsset] += pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                }

                if (isUnitBasket) {
                    baskRef += inst->weights[iAsset] * pathGen->refLevel(iAsset, 0/*iPath*/);
                    baskSumOut += inst->weights[iAsset] * sumOut[iAsset];
                }
                else {
                    baskSumOut += inst->weights[iAsset] * sumOut[iAsset] / pathGen->refLevel(iAsset, 0/*iPath*/);
                }
            }
            if (isUnitBasket) {
                baskSumOut /= baskRef;
            }
            if (pathGen->doingPast()){ // preserve values
                sumOutSoFar = sumOut;
            }

            double perf = baskSumOut / inst->averageOutDates.size();
            double myPayoff;
            if (inst->isRefAvgOut) {
                // only call supported
                myPayoff = 1.0 - inst->strike / perf;
            } else {
                myPayoff = inst->isCall? 
                    (perf - inst->strike): (inst->strike - perf);
            }
            if (vanillaReprice.get()){
                vanillaReprice->store(myPayoff);
            }

            prices.add(inst->notional * prices.maxWithZero(myPayoff));
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool fwdStarting = startDate.isGreater(today);
        double interpLevel;

        if (inst->smileType == XCB::NO_SMILE) {
            reqarr[0] = CVolRequestLNSP(new ATMVolRequest());
        } else {
            if (fwdStarting){
                interpLevel = inst->strike;
            } else {
                /* not forward starting - some samples have fixed already
               (this includes averaging in) */
                int numDates = inst->averageOutDates.size();
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                // moneyness is from basket levels
                interpLevel = (numDates * inst->strike - baskSumOut)/ numRemaining;

                // here we're just using iPath = 0 since at this point
                // all the paths are the same
                // scale for each component - smileType D for now
                if (inst->smileType == XCB::FIXED_STRIKE_SMILE) {
                    interpLevel *= pathGen->refLevel(iAsset, 0);
                } else if (inst->smileType == XCB::FLOAT_STRIKE_SMILE) {
                    /*** XXX SN - not sure how correct this is! Esp during average-in! ***/
                    // Essentially we take % into component via spot, not SAS
                    interpLevel *= getMultiFactors()->assetGetSpot(iAsset);
                    // ... but
                    // the % interp is relative to baskRef, so need also to "rebase" to baskSpot
                    double assetSpot;
                    double baskRef = 0;
                    double baskSpot = 0;
                    for(int i=0;i<nbAssets; i++) {
                        assetSpot = getMultiFactors()->assetGetSpot(i);
                        if (isUnitBasket) {
                            baskRef += inst->weights[i] * pathGen->refLevel(i, 0/*iPath*/);
                            baskSpot += inst->weights[i] * assetSpot;
                        } else {
                            baskSpot += inst->weights[i] * assetSpot / pathGen->refLevel(i, 0/*iPath*/);
                        }
                    }
                    if (isUnitBasket) {
                        interpLevel *= baskRef / baskSpot;
                    } else {
                        // baskRef is 1
                        interpLevel /= baskSpot;
                    }
                }
            }
            const SimSeries* simSeries = getSimSeries();
            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     simSeries->getLastDate(),
                                                                     fwdStarting));
        }
        return reqarr;
    }
};


//////////////////////////////////////////////////////////////////////////////////

/* MC product class for BaskAv state vars - with isRefAvgOut false */
class BaskAvSVMC : public MCProductClient,
                  virtual public IMCProductLN,
                  virtual public IMCQuickGreeks,
                  virtual public IMCQuickXGamma {
private:
    const BaskAv*             inst;        // reference to original instrument
    int                       nbAssets;     // convenient
    DoubleArray               sumOut;       // working area - here to save alloc
    DoubleArray               sumOutSoFar;  // preservation area (past)
    bool                      isUnitBasket;
    double                    baskSumOut;
    RepriceVanillaSP          vanillaReprice; // working area (price only)

    SVGenSpotSP                  spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    SVGenDiscFactorSP            dfGen;        //!< Generator for discount factors
    SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVDiscFactorSP dfSV;         //!< Df state variable

protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "BaskAvSVMC::pathGenUpdated";

        try {
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

    
    /** Override default to allow us to switch off if we're forward 
        starting unit basket */
    virtual IMCQuickXGamma* quickXGammaSupported(){
        if (isUnitBasket){
            DateTimeArray futureRefLevelDates(inst->refLevel->
                                              getFutureDates(inst->valueDate));
            if (!futureRefLevelDates.empty()){
                // basket level is non linear in spot
                return 0;
            }
        }
        if (inst->isRefAvgOut) {
            // xgamma is not so easy
            return 0;
        }
        return this;
    }
public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(dfGen.get());
    }
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    BaskAvSVMC(const BaskAv*            inst,
               const SimSeriesSP&       simSeries):
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
        sumOut(nbAssets),
        sumOutSoFar(nbAssets),
        baskSumOut(0.0), spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())){
        isUnitBasket = (inst->weightType == "U");
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form
        barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Create IMCPrices object for first pricing call. The storingMode 
        parameter is true if this IMCPrices object will be used for future
        tweaked pricing runs. */
    virtual IMCPrices* createOrigPrices(int  nbIter,
                                     int  nbSubSamples,
                                     int  mode) {
        // describe our payoff - v simple
        if (mode>0){
            if (inst->isRefAvgOut) {
                throw ModelException("BaskAvSVMC::createOrigPrices", 
                                     "Quick Greeks not yet supported with "
                                     "isRefAvgOut True.");
            }
            vanillaReprice = RepriceVanillaSP(
                new IReprice::Vanilla(mode, nbIter, inst->notional));
       }
        // then our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, vanillaReprice);
    }   

    /** Returns the maximum delta scaling factor for each asset (it is the
        largest absolute derivative of the payoff wrt each simulated asset)
        Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen) const{
        DoubleArray maxDeltaFactor(getNumAssets());
        double baskRef = 0.0;
        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            if (isUnitBasket) {
                baskRef += inst->weights[i] * refLevelSV->refLevel(i);
            } else {
                // basket of performances
                maxDeltaFactor[i] = fabs(inst->weights[i]) / refLevelSV->refLevel(i);
            }
        } 
        if (isUnitBasket) {
           for (int i = 0; i < maxDeltaFactor.size(); i++) {
               maxDeltaFactor[i] = fabs(inst->weights[i]) / baskRef;
           }
        }
        return maxDeltaFactor;
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
        myPrices.getReprice()->setForOneSidedGreek(0); // regardless of mode
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
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        myPrices.getReprice()->setForTwoSidedGreek(&params);
    }

    /** Create IMCPrices object for doing quick x gamma. The IMCPrices
        object needs to be able to skip over paths that have zero
        gamma */
    virtual void setPricesForXGamma(
        IMCPrices*                 untweakedPrices,
        const IPathGenerator*   futurePathGen,
        const ScalarShiftArray& sens) {
        MCPricesGeneral& myPrices =
            dynamic_cast<MCPricesGeneral&>(*untweakedPrices);
        myPrices.setMode(MCPricesGeneral::CROSS_GAMMA);
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        myPrices.getReprice()->setForXGamma(&params);
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("BaskAvSVMC::payoff");
        try {
            double baskRef = 0.0;

            baskSumOut = 0;
            sumOut = sumOutSoFar;

            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                const SVPath& path = spotSV->path(iAsset);
                for (int iStep = path.begin(); iStep < path.end(); iStep++) {
                    sumOut[iAsset] += path[iStep];
                }
                if (isUnitBasket) {
                    baskRef += inst->weights[iAsset] * 
                        refLevelSV->refLevel(iAsset);
                    baskSumOut += inst->weights[iAsset] * sumOut[iAsset];
                } else {
                    baskSumOut += inst->weights[iAsset] * sumOut[iAsset] / 
                        refLevelSV->refLevel(iAsset);
                }
            }
            if (isUnitBasket) {
                baskSumOut /= baskRef;
            }
            if (pathGen->doingPast()){ // preserve values
                sumOutSoFar = sumOut;
            }

            double perf = baskSumOut / inst->averageOutDates.size();
            double myPayoff;
            if (inst->isRefAvgOut) {
                // only call supported
                myPayoff = 1.0 - inst->strike / perf;
            } else {
                myPayoff = inst->isCall? 
                    (perf - inst->strike): (inst->strike - perf);
            }
        
            if (doingPast() && !hasFuture()) {
                // If we're in a known state, we record known flows on
                // their known date (so no discounting).
                if (!paymentDate.empty()) {
                    knownCashFlows->addFlow(paymentDate,
                                            Maths::max(inst->notional * myPayoff, 0.0));
                }
            }

            // Discount
            myPayoff *= dfSV->firstDF();

            if (vanillaReprice.get()){
                vanillaReprice->store(myPayoff);
            }
            prices.add(inst->notional * prices.maxWithZero(myPayoff));
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool fwdStarting = startDate.isGreater(today);
        double interpLevel;

        if (inst->smileType == XCB::NO_SMILE) {
            reqarr[0] = CVolRequestLNSP(new ATMVolRequest());
        } else {
            if (fwdStarting){
                interpLevel = inst->strike;
            } else {
                /* not forward starting - some samples have fixed already
               (this includes averaging in) */
                int numDates = inst->averageOutDates.size();
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                // moneyness is from basket levels
                interpLevel = (numDates * inst->strike - baskSumOut)/
                    numRemaining;

                // here we're just using iPath = 0 since at this point
                // all the paths are the same
                // scale for each component - smileType D for now
                if (inst->smileType == XCB::FIXED_STRIKE_SMILE) {
                    interpLevel *= refLevelSV->refLevel(iAsset);
                } else if (inst->smileType == XCB::FLOAT_STRIKE_SMILE) {
                    /*** XXX SN - not sure how correct this is! Esp
                         during average-in! ***/
                    // Essentially we take % into component via spot, not SAS
                    interpLevel *= getMultiFactors()->assetGetSpot(iAsset);
                    // ... but the % interp is relative to baskRef, so
                    // need also to "rebase" to baskSpot
                    double assetSpot;
                    double baskRef = 0;
                    double baskSpot = 0;
                    for(int i=0;i<nbAssets; i++) {
                        assetSpot = getMultiFactors()->assetGetSpot(i);
                        if (isUnitBasket) {
                            baskRef += inst->weights[i]*refLevelSV->refLevel(i);
                            baskSpot += inst->weights[i] * assetSpot;
                        } else {
                            baskSpot += inst->weights[i] * assetSpot / 
                                refLevelSV->refLevel(i);
                        }
                    }
                    if (isUnitBasket) {
                        interpLevel *= baskRef / baskSpot;
                    } else {
                        // baskRef is 1
                        interpLevel /= baskSpot;
                    }
                }
            }
            const SimSeries* simSeries = getSimSeries();
            reqarr[0] = CVolRequestLNSP(
                new LinearStrikeTSVolRequest(interpLevel,
                                             startDate,
                                             simSeries->getLastDate(),
                                             fwdStarting));
        }
        return reqarr;
    }
};


//////////////////////////////////////////////////////////////////////////////////

/* stateless MC product class for BaskAv state vars - with isRefAvgOut false */
class BaskAvSLMC : public MCProductClient,
                  virtual public IMCProductLN,
                  virtual public IMCStatelessProductClient{
private:
    const BaskAv*             inst;        // reference to original instrument
    int                       nbAssets;     // convenient
    DoubleArray               sumOut;       // working area - here to save alloc
    DoubleArray               sumOutSoFar;  // preservation area (past)
    bool                      isUnitBasket;
    double                    baskSumOut;
    RepriceVanillaSP          vanillaReprice; // working area (price only)

    SVGenSpotSP                  spotGen;      //!< Generator for spot
    IRefLevel::IStateVarGenSP refLevelGen;  //!< Generator for ref level
    SVGenDiscFactorSP            dfGen;        //!< Generator for discount factors
    SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
    IRefLevel::IStateVarSP    refLevelSV;   //!< Ref level state variable
    SVDiscFactorSP dfSV;         //!< Df state variable
    int                       matIdx;       // maturity date
    DateTimeArray             pastDates;
    DateTimeArray             statelessPayOffDates;

protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "BaskAvSLMC::pathGenUpdated";

        try {
            spotSV = spotGen->getSpotSV(newPathGen);
            refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
            dfSV = dfGen->getSVDiscFactor(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    class HistoricalContext : public IHistoricalContext {
        friend class BaskAvSLMC;
        DoubleArray sumOutSoFar;
    public:
        HistoricalContext(int nbAssets):sumOutSoFar(nbAssets)
        {
            for (int i = 0; i < nbAssets; ++i)
                sumOutSoFar[i] = 0;
        }

        virtual void deepCopyTo( IHistoricalContext* destination ) const
        {
            HistoricalContext* dest = static_cast<HistoricalContext*>( destination );
            *dest = *this;
        }
    };

    typedef refCountPtr<HistoricalContext> HistoricalContextSP;
    HistoricalContextSP past; // context generated by past path generator


    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // ask for a reference level State Variable
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(dfGen.get());
    }
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    BaskAvSLMC(const BaskAv*            inst,
               const SimSeriesSP&       simSeries):
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
        sumOut(nbAssets),
        sumOutSoFar(nbAssets),
        baskSumOut(0.0), spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, simSeries->getLastDate())),
        pastDates(inst->valueDate.getPastDates(simSeries->getAllDates())),
        statelessPayOffDates(simSeries->getFutureDates(inst->valueDate)),
        past(new HistoricalContext(nbAssets)){
        isUnitBasket = (inst->weightType == "U");
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form
        barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    // IMCNewClient
    virtual DateTimeArray getPastDates()
    {
        return pastDates;
    }
    
    virtual IHistoricalContextSP createHistoricalContext()
    {
        return IHistoricalContextSP(new HistoricalContext(nbAssets));
    }

    virtual IHistoricalContextSP getInitialHistoricalContext()
    {
        return IHistoricalContextSP(past);
    }

    virtual vector<int> finalize( 
        const DateTimeArray& timeline )
    {
        // we are only interested in average out dates
        std::vector<int> keyDates;

        for (int i = 0; i < statelessPayOffDates.size(); ++i) {
            for (int j = 0; j < timeline.size(); ++j) {
                if (timeline[j] == statelessPayOffDates[i])
                    keyDates.push_back(j);
            }
        }
        matIdx = statelessPayOffDates.size() - 1;
        ASSERT(int(keyDates.size()) == statelessPayOffDates.size());
        return keyDates;
    }

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history, 
        IMCPrices& prices)
    {
        static const string routine("BaskAvSLMC::statelessPayoff");
        try {
            HistoricalContext* hContext =
                static_cast<HistoricalContext*>(history.get());

            DoubleArray spotPrices(nbAssets);
            spotSV->getSpotPrices(spotPrices);
        
            // sum up
        int iAsset;
            for (iAsset = 0; iAsset < nbAssets; ++iAsset) {
                hContext->sumOutSoFar[iAsset] += spotPrices[iAsset];
            }

            double baskRef = 0.0;
            baskSumOut = 0;
        
            bool needToPrice(false);
            if (doingPast())
                needToPrice = currentDateIdx == pastDates.size() -1;
            else
                needToPrice = currentDateIdx == matIdx;

            if (!needToPrice)
                return;

            for (iAsset = 0; iAsset < nbAssets; ++iAsset) {
                if (isUnitBasket) {
                    baskRef += inst->weights[iAsset] * 
                        refLevelSV->refLevel(iAsset);
                    baskSumOut += inst->weights[iAsset] * hContext->sumOutSoFar[iAsset];
                } else {
                    baskSumOut += inst->weights[iAsset] * hContext->sumOutSoFar[iAsset] / 
                        refLevelSV->refLevel(iAsset);
                }
            }

            if (isUnitBasket) {
                baskSumOut /= baskRef;
            }        

            double perf = baskSumOut / inst->averageOutDates.size();
            double myPayoff;
            if (inst->isRefAvgOut) {
                // only call supported
                myPayoff = 1.0 - inst->strike / perf;
            } else {
                myPayoff = inst->isCall? 
                    (perf - inst->strike): (inst->strike - perf);
            }
            // Discount
            myPayoff *= dfSV->firstDF();
            
            if (vanillaReprice.get()){
                vanillaReprice->store(myPayoff);
            }
            prices.add(inst->notional * prices.maxWithZero(myPayoff));
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) 
    {
        static const string routine("BaskAvSLMC::statelessPayoff");
        throw ModelException("payoff function called on stateless client!", routine);
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool fwdStarting = startDate.isGreater(today);
        double interpLevel;

        if (inst->smileType == XCB::NO_SMILE) {
            reqarr[0] = CVolRequestLNSP(new ATMVolRequest());
        } else {
            if (fwdStarting){
                interpLevel = inst->strike;
            } else {
                /* not forward starting - some samples have fixed already
               (this includes averaging in) */
                int numDates = inst->averageOutDates.size();
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                // moneyness is from basket levels
                interpLevel = (numDates * inst->strike - baskSumOut)/
                    numRemaining;

                // here we're just using iPath = 0 since at this point
                // all the paths are the same
                // scale for each component - smileType D for now
                if (inst->smileType == XCB::FIXED_STRIKE_SMILE) {
                    interpLevel *= refLevelSV->refLevel(iAsset);
                } else if (inst->smileType == XCB::FLOAT_STRIKE_SMILE) {
                    /*** XXX SN - not sure how correct this is! Esp
                         during average-in! ***/
                    // Essentially we take % into component via spot, not SAS
                    interpLevel *= getMultiFactors()->assetGetSpot(iAsset);
                    // ... but the % interp is relative to baskRef, so
                    // need also to "rebase" to baskSpot
                    double assetSpot;
                    double baskRef = 0;
                    double baskSpot = 0;
                    for(int i=0;i<nbAssets; i++) {
                        assetSpot = getMultiFactors()->assetGetSpot(i);
                        if (isUnitBasket) {
                            baskRef += inst->weights[i]*refLevelSV->refLevel(i);
                            baskSpot += inst->weights[i] * assetSpot;
                        } else {
                            baskSpot += inst->weights[i] * assetSpot / 
                                refLevelSV->refLevel(i);
                        }
                    }
                    if (isUnitBasket) {
                        interpLevel *= baskRef / baskSpot;
                    } else {
                        // baskRef is 1
                        interpLevel /= baskSpot;
                    }
                }
            }
            const SimSeries* simSeries = getSimSeries();
            reqarr[0] = CVolRequestLNSP(
                new LinearStrikeTSVolRequest(interpLevel,
                                             startDate,
                                             simSeries->getLastDate(),
                                             fwdStarting));
        }
        return reqarr;
    }
};

//////////////////////////////////////////////////////////////////////////////////


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* BaskAv::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    
    if(model->stateVarUsed()) {
        if (model->inStatelessMode())
            return new BaskAvSLMC(this, simSeries);
        else
            return new BaskAvSVMC(this, simSeries);
    }
    
    // otherwise, use old methodology
    return new BaskAvMC(this, simSeries);
}

CClassConstSP const BaskAv::TYPE = CClass::registerClassLoadMethod(
    "BaskAv", typeid(BaskAv), BaskAv::load);

// * for class loading (avoid having header file) */
bool BaskAvLoad() {
    return (BaskAv::TYPE != 0);
}

DRLIB_END_NAMESPACE
