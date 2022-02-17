//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Rebalanced Basket.cpp
//
//   Description : Option based upon the rebalanced n assets
//   Author      : Xiaolan zhang
//
//   Date        : April 2005It
//
//
//
//   Revision 1.2  2005/04/25 20:17:17  xZhang
//   change BaskReb to BaskRebalanced
//
//   Revision 1.1  2005/04/25 13:53:08  xZhang
//   new file for Product Rebalanced Basket
//
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/** BaskRebalanced product - option based upon the rebalanced n assets
    sorted by their performance */
class BaskRebalanced: public GenericNFBase, 
              virtual public IMCIntoProduct,
              public ITaxableInst::Basic {
protected:
    /// fields ////////
    DoubleArray             weights;
    bool                    isCall;
    double                  strike;
    DateTimeArray           rebalanceDates;

public:
    static CClassConstSP const TYPE;
    friend class BaskRebalancedMC;

    // validation
    void validatePop2Object(){
        static const string routine("BaskRebalanced::validatePop2Object");
        GenericNFBase::validatePop2Object();
  
        // check that we've got as many weights as there are assets 
        // and if % weights that they sum to 100%

        AssetUtil::checkWeights(weights, assets->NbAssets());

        if (Maths::isNegative(strike)){
            throw ModelException(routine, "strike ("+
                                 Format::toString(strike)+") is negative");
        }
        // validate dates are not empty - order is handled by SimSeries
        if (rebalanceDates.empty()) {
            throw ModelException(routine, "No rebalanceDates supplied!");
        }       
    }
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return rebalanceDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    /** for ITaxableInst::Basic */
    const DateTime getFinalPaymentDate() const {
        return instSettle->settles(rebalanceDates.back(), NULL);
    }

private:
    BaskRebalanced(): GenericNFBase(TYPE) /*, isRefAvgOut(false)*/ {} // for reflection
    BaskRebalanced(const BaskRebalanced& rhs); // not implemented
    BaskRebalanced& operator=(const BaskRebalanced& rhs); // not implemented

    static IObject* defaultBaskRebalanced(){
        return new BaskRebalanced();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BaskRebalanced, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBaskRebalanced);
        FIELD(weights,          "Weights");
        FIELD(isCall,           "is it a call option");
        FIELD(strike,           "strike");
        FIELD(rebalanceDates,  "rebalanceDates");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};


///////////////////////////////////////////////////////////////////////////
/* MC product class for BaskRebalanced - with isRefAvgOut false */
class BaskRebalancedMC : public IMCProduct,
                 virtual public IMCProductLN,
                 virtual public IMCQuickGreeks,
                 virtual public IMCQuickXGamma{
private:
    const BaskRebalanced*            inst;        // reference to original instrument
    int                      nbAssets;     // convenient
    DoubleArray              rebWeights;       // 
    DoubleArray              rebWeightsSoFar;  // preservation area (past)

    RepriceVanillaSP         vanillaReprice; // working area (price only)

protected:
    /** Override default to allow us to switch off if we're forward 
        starting unit basket */

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    BaskRebalancedMC(const BaskRebalanced*            inst,
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
        rebWeights(nbAssets),
        rebWeightsSoFar(nbAssets){
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

        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            maxDeltaFactor[i] = fabs(rebWeights[i]) ;
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
        static const string routine("BaskRebalancedMC::payoff");
        try {
            int    beginIdx = pathGen->begin(0); // same for all assets
            int    endIdx   = pathGen->end(0);
            double assetValue = 0;
            int iStep, iAsset;

            rebWeights= rebWeightsSoFar;

            for (iStep = beginIdx; iStep < endIdx; iStep++) {

                //first time  
                //if (iStep == 0  &&  pathGen->doingPast() ){
                if (iStep == 0 ){
                    assetValue = 1.0;
                }else{
                    assetValue = 0.0;
                    for (int iAsset=0; iAsset < nbAssets; iAsset++) {                       
                        assetValue += rebWeights[iAsset] * pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                    }
                }

                for (iAsset=0; iAsset < nbAssets; iAsset++) {                       
                    rebWeights[iAsset] = assetValue * inst->weights[iAsset] 
                                        / pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                }               
            }
            
            if (pathGen->doingPast()){ // preserve values
                rebWeightsSoFar = rebWeights;
            }

            double myPayoff;
            myPayoff = inst->isCall? 
                    (assetValue - inst->strike): (inst->strike - assetValue);

            if (vanillaReprice.get()){
                // save pre max price in vanillaReprice
                vanillaReprice->store(myPayoff);
            }
            prices.add(inst->notional * prices.maxWithZero(myPayoff));

        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {

        static const string method = "BaskRebalancedMC::getVolInterp";

        try
        {
            // one interp level per asset
            CVolRequestLNArray reqarr(1);
            
            double interpLevel = inst->strike * pathGen->refLevel(iAsset, 0);

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
};


//////////////////////////////////////////////////////////////////////////////////

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* BaskRebalanced::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(rebalanceDates);
    
    if(model->stateVarUsed()) {
        // State variables
        //return new BaskRebalancedSVMC(this, simSeries);
    }
    
    // otherwise, use old methodology
    return new BaskRebalancedMC(this, simSeries);
}

CClassConstSP const BaskRebalanced::TYPE = CClass::registerClassLoadMethod(
    "BaskRebalanced", typeid(BaskRebalanced), BaskRebalanced::load);

// * for class loading (avoid having header file) */
bool BaskRebalancedLoad() {
    return (BaskRebalanced::TYPE != 0);
}


DRLIB_END_NAMESPACE
