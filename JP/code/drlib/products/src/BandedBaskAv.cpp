//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BandedBaskAv.cpp
//
//   Description : Port of BandedAverage from EDG
//
//   Date        : Aug 2002
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
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/MCPrices.hpp"

DRLIB_BEGIN_NAMESPACE

class BandedBaskAv: public GenericNFBase, 
                    virtual public IMCIntoProduct {
protected:
    /// fields ////////
    string                  weightType;
    string                  smileType; // values to match XCB
    DoubleArray             weights;
    bool                    isCall;
    DoubleArray             strikes; // 1 or 2 for spread
    DateTimeArray           averageOutDates;
    double                  floor;
    bool                    isCapped;
    double                  cap;
    bool                    isBandComponents;
    bool                    isBandPerSample;

    // convenient values
    int                     nbBands; // 1 (just floor), or 2 (with cap too)
    bool                    isUnitBasket;

public:
    static CClassConstSP const TYPE;
    friend class BandedBaskAvMC;

    static const string UNIT_WEIGHTS;
    static const string PCT_WEIGHTS;

    // validation
    void validatePop2Object(){
        static const string routine("BandedBaskAv::validatePop2Object");
        GenericNFBase::validatePop2Object();

        // check that we've got as many weights as there are assets 
        // and if % weights that they sum to 100%
        if (weightType==PCT_WEIGHTS) {
            isUnitBasket = false;
            AssetUtil::checkWeights(weights, assets->NbAssets());
        } else if (weightType==UNIT_WEIGHTS) {
            isUnitBasket = true;
            if (assets->NbAssets() != weights.size()){
                throw ModelException(routine,
                                     "Different number of assets ("+
                                     Format::toString(assets->NbAssets())+")"
                                     " to weights ("+
                                     Format::toString(weights.size())+")");
            }
        } else {
            throw ModelException(routine, "weightType must be " + PCT_WEIGHTS +
                                 " or " + UNIT_WEIGHTS + ", but " + weightType+
                                 " given");
        }

        // check no more than 2 strikes, and are positive and increasing
        // Also check relative levels to floor and cap (if there) are sensible
        if (strikes.size()<1 || strikes.size()>2) {
            throw ModelException(routine, "At least 1 and no more than 2 "
                                 "strikes should be supplied (" + 
                                 Format::toString(strikes.size()) + " given)");
        }
        if (Maths::isNegative(floor)){
            throw ModelException(routine, "Floor must not be negative ("+
                                 Format::toString(floor)+" given)");
        }
        if (floor>strikes[0]) {
            throw ModelException(routine, "Floor must not be greater "
                                 "than strikes but "+
                                 Format::toString(floor)+" > " +
                                 Format::toString(strikes[0]));
        }
        // Don't really care about rough equality here
        if (strikes.size()>1 && strikes[0]>=strikes[1]) {
            throw ModelException(routine, "Strikes must be increasing but "+
                                 Format::toString(strikes[0])+" >= " +
                                 Format::toString(strikes[1]));
        }
        nbBands = 1;
        if (isCapped) {
            nbBands++;
            if (cap<strikes[strikes.size()-1]) {
                throw ModelException(routine, 
                                     "Cap must be greater than strikes but "+
                                     Format::toString(cap)+" < " +
                                     Format::toString(strikes.back()));
            }
        }

        // validate dates are not empty and are in order 
        DateTime::ensureIncreasing(averageOutDates, "averageOutDates", true);

        // Could give a meaning to this, but it's a rare case and not
        // worth the effort.
        if (isUnitBasket && isBandComponents) {
            throw ModelException(routine, "Unit baskets banding"
                                 " components is not supported!");
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
    BandedBaskAv(): GenericNFBase(TYPE) {} // for reflection
    BandedBaskAv(const BandedBaskAv& rhs); // not implemented
    BandedBaskAv& operator=(const BandedBaskAv& rhs); // not implemented

    static IObject* defaultBandedBaskAv(){
        return new BandedBaskAv();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BandedBaskAv, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultBandedBaskAv);
        FIELD(weightType,       "Unit(U) or Percentage(P) basket");
        FIELD(weights,          "Weights");
        FIELD(smileType,        "Type of smile as per XCB");
        FIELD(isCall,           "is it a call option");
        FIELD(strikes,          "strikes - 1 or 2 for spread");
        FIELD(averageOutDates,  "averageOutDates");
        FIELD(floor,  "floor");
        FIELD(isCapped,  "isCapped");
        FIELD(cap,  "cap");
        FIELD(isBandComponents,  "isBandComponents");
        FIELD(isBandPerSample,  "isBandPerSample");

        FIELD(nbBands,  "nbBands");
        FIELD_MAKE_TRANSIENT(nbBands);
        FIELD(isUnitBasket,  "isUnitBasket");
        FIELD_MAKE_TRANSIENT(isUnitBasket);

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class  */
class BandedBaskAvMC: public IMCProduct,
                      virtual public IMCProductLN,
                      /** If you want to support IQuickXGamma then you need
                          to change the use of the spread class is it
                          doesn't fully capture the optionality of the
                          banded basket which is needed for IQuickXGamma */
                      virtual public IMCQuickGreeks {
private:

    const BandedBaskAv*      inst;        // reference to original instrument
    DoubleArray              effectiveStrikes; // input strikes - floor
    int                      nbPaths; // for LN interp levels
    DoubleArray              band;
    int                      nbAvgOut; // for convenience
    DoubleMatrix             avg; // Non-banded values [iPath][iAsset]
    DoubleArray              bandedAvg;
    // These preserve values above from past to future
    DoubleMatrix             avgSoFar;
    DoubleArray              bandedAvgSoFar;
    // indexes for different bands/interps
    int                      iFloor;
    int                      iFloor2ndStrike;
    int                      iCap;

    // for quick greeks
    RepriceSpreadSP          bbSpreadReprice;  // working area (price only)
    RepriceVanillaSP         bbVanillaReprice; // working area (price only)
    RepriceSpreadSP          spreadReprice;  // working area (price only)
    RepriceVanillaSP         vanillaReprice; // working area (price only)

public:

    BandedBaskAvMC(const BandedBaskAv*            inst,
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
        effectiveStrikes(inst->strikes),
        nbPaths(inst->nbBands + (inst->strikes.size()>1?1:0)), // 1 or 2 or 3
        band(nbPaths),
        nbAvgOut(inst->averageOutDates.size()),
        avg(nbPaths,getNumAssets()),
        bandedAvg(nbPaths),
        avgSoFar(nbPaths,getNumAssets()),
        bandedAvgSoFar(nbPaths) {

        for(int i=0; i<effectiveStrikes.size(); i++) {
            effectiveStrikes[i] -= inst->floor;
        }

        // If Implied (or other model which doesn't worry about vol
        // interp) could have iFloor = iFloor2ndStrike = 0, iCap =
        // isCapped?1:0 For LN and safe but slower for all :-
        iFloor = 0;
        iFloor2ndStrike = inst->strikes.size()>1?1:0;
        iCap = inst->isCapped?iFloor2ndStrike+1:0;

        band[iFloor] = inst->floor;
        band[iFloor2ndStrike] = inst->floor; /* if extra strike =>
                                                diff interp for floor path */
        if (inst->isCapped) {
            band[iCap] = inst->cap;
        }
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
        if (mode == 0){
            // don't bother building any reprice objects
        } else {
            /* Payoff is option/spread on banded basket. Can view
               banded basket as call spread with strikes {floor, cap}.
               However viewing the banded basket as a call spread
               alters the strike(s) of the option/spread. So the final
               option/spread has {strikes} -> {strikes}-floor. Now payoff
               code actually calculates banded basket - floor. So we
               can do call spread with strikes {0, Cap-Floor} for the
               "banded basket".  If there is no cap then can view
               banded basket as call option with strike 0.
            */
            IRepriceSP bbReprice;
            if (inst->isCapped){
                // use call spread
                bbSpreadReprice = RepriceSpreadSP(
                    new IReprice::Spread(ISubRepriceSP(), mode, nbIter,
                                         true, 1.0, // irrelevant
                                         0.0, // lo strike
                                         inst->cap-inst->floor)); // hi strike
                bbReprice = bbSpreadReprice;
            } else {
                // use vanilla
                bbVanillaReprice =  RepriceVanillaSP(
                    new IReprice::Vanilla(mode, nbIter, 1.0 /* irrelevant */));
                bbReprice = bbVanillaReprice;
            }
            if (inst->strikes.size() > 1){
                // use spread
                spreadReprice = RepriceSpreadSP(
                    new IReprice::Spread(bbReprice, 
                                         mode, nbIter,
                                         inst->isCall, inst->notional,
                                         effectiveStrikes[0], // lo strike
                                         effectiveStrikes[1])); // hi strike
                reprice = spreadReprice;
            } else {
                // use vanilla
                vanillaReprice = RepriceVanillaSP(
                    new IReprice::Vanilla(bbReprice, mode,
                                          nbIter, inst->notional));
                reprice = vanillaReprice;
            }
            
        }
        // then our prices object
        return new MCPricesGeneral(nbIter, nbSubSamples, reprice);
    }   

    /** Returns the maximum delta scaling factor for each asset (it is the
        largest absolute derivative of the payoff wrt each simulated asset)
        Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen) const{
        DoubleArray maxDeltaFactor(getNumAssets(), 1.0);
        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            maxDeltaFactor[i] /=
                futurePathGen->refLevel(i, 0 /* path irrelevant*/);
            // we just scale each of the maxDeltaFactor[] by its weight. 
            maxDeltaFactor[i] *= fabs(inst->weights[i]);
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
        DoubleArray scaleFactors(maxDeltaFactors(futurePathGen));
        IReprice* reprice = myPrices.getReprice().get();
        // vanilla and spread use same params class
        IReprice::Vanilla::TwoSidedGreekParamsSet params(scaleFactors,
                                                         futurePathGen,
                                                         getMultiFactors(),
                                                         sens);
        reprice->setForTwoSidedGreek(&params);
        reprice->getRepriceForCmpt()->setForTwoSidedGreek(&params);
    }
    
    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("BandedBaskAvMC::payoff");
        int          beginIdx = pathGen->begin(0); // same for all assets
        int          endIdx   = pathGen->end(0);
        int          iPath, iAsset, iStep;
        // Any past samples
        avg = avgSoFar;
        bandedAvg = bandedAvgSoFar;

        for(iPath=0; iPath<nbPaths; iPath++) {
            if (inst->isBandComponents) {
                // Validation ensures percentage basket here
                if (inst->isBandPerSample) {
                    // Case D - ref equation 4
                    for(iAsset=0; iAsset<getNumAssets(); iAsset++) {
                        for (iStep=beginIdx; iStep<endIdx; iStep++) {
                            bandedAvg[iPath] += inst->weights[iAsset]*
                                Maths::max(pathGen->Path(iAsset, iPath)[iStep]
                                           / pathGen->refLevel(iAsset, iPath)
                                           - band[iPath], 0.0);
                        }
                    }
                } else {
                    // Case C - ref equation 3 
                    bandedAvg[iPath] = 0.0;
                    for(iAsset=0; iAsset<getNumAssets(); iAsset++) {
                        for (iStep=beginIdx; iStep<endIdx; iStep++) {
                            avg[iPath][iAsset] += 
                                pathGen->Path(iAsset, iPath)[iStep] / 
                                pathGen->refLevel(iAsset, iPath);
                        }
                        bandedAvg[iPath] += inst->weights[iAsset]*
                            Maths::max(avg[iPath][iAsset]/nbAvgOut -
                                       band[iPath],0.0);
                    }
                } 
            } else {
                // Collapse the asset dimension before banding
                double baskRef = 0.0;
                if (inst->isUnitBasket) {
                    for(iAsset=0; iAsset<getNumAssets(); iAsset++) {
                        baskRef += inst->weights[iAsset] * 
                            pathGen->refLevel(iAsset, iPath);
                    }
                }
                for (iStep=beginIdx; iStep<endIdx; iStep++) {
                    double bask = 0;
                    for(iAsset=0; iAsset<getNumAssets(); iAsset++) {
                        double comp = inst->weights[iAsset] * 
                            pathGen->Path(iAsset, iPath)[iStep];
                        if (!inst->isUnitBasket) {
                            comp /= pathGen->refLevel(iAsset, iPath);
                        }
                        bask += comp;
                    }
                    if (inst->isUnitBasket) {
                        bask /= baskRef;
                    }
                    if (inst->isBandPerSample) {
                        // Case B - ref equation 2
                        bandedAvg[iPath] += Maths::max(bask-band[iPath], 0.0);
                    } else {
                        avg[iPath][0] += bask;
                    }
                } // iStep loop 
                if (!inst->isBandPerSample) {
                    // Case A - ref equation 1
                    bandedAvg[iPath] =
                        Maths::max(avg[iPath][0]/nbAvgOut - band[iPath], 0.0);
                }
            }
        }

        if (pathGen->doingPast()) {
            avgSoFar = avg;
            bandedAvgSoFar = bandedAvg;
        }
        
        double bandedBask = (bandedAvg[iFloor] - 
                             (inst->isCapped?bandedAvg[iCap]:0)) / 
            (inst->isBandPerSample?nbAvgOut:1);
        if (bbSpreadReprice.get()){
            // the spread object takes the 'fwd' or performance
            bbSpreadReprice->store(bandedBask);
        } else if (bbVanillaReprice.get()){
            // the vanilla object takes the 'fwd' (or performance) - strike
            // here strike = 0
            bbVanillaReprice->store(bandedBask);
        }
        double fwdL = inst->isCall? (bandedBask - effectiveStrikes[0]) : 
            (effectiveStrikes[0] - bandedBask);
        if (vanillaReprice.get()){
            // the vanilla object takes the 'fwd' (or performance) - strike
            vanillaReprice->store(fwdL);
        } else if (spreadReprice.get()){
            // the spread object takes the 'fwd' or performance
            spreadReprice->store(bandedBask);
        }
        double optL = Maths::max(fwdL, 0.0);
        if (effectiveStrikes.size() == 2) {
            bandedBask = (bandedAvg[iFloor2ndStrike] - 
                          (inst->isCapped?bandedAvg[iCap]:0)) / 
                (inst->isBandPerSample?nbAvgOut:1);
            double fwdH = (inst->isCall? (bandedBask - effectiveStrikes[1]):
                           (effectiveStrikes[1] - bandedBask));
            double optH =Maths::max(fwdH, 0.0);
            // call spread or put spread
            prices.add(inst->notional * 
                       prices.maxWithZero(inst->isCall?
                                          (optL-optH): (optH-optL)));
        } else {
            prices.add(inst->notional * optL);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(nbPaths);
        
        const IRefLevel* refLevel = getRefLevel();
        const DateTime&  startDate = refLevel->getAllDates().front();
        const DateTime&  today = getToday();
        bool             fwdStarting = startDate.isGreater(today);
        int              iPath;
        double           levels[3];
        double           interpLevel;

        levels[iFloor] = inst->strikes[0];
        if (inst->strikes.size()>1) levels[iFloor2ndStrike] = inst->strikes[1];
        if (inst->isCapped) levels[iCap] = inst->cap;

        for(iPath=0; iPath<nbPaths; iPath++) {
            if (fwdStarting){
                interpLevel = levels[iPath];
            } else {
                /* not forward starting - some samples have fixed already
                   (this includes averaging in) */
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                int nbPast = nbAvgOut - numRemaining;

                if (inst->isBandPerSample)
                {
                    /* Note here we use FAvg & CAvg to adjust for the
                     * past, while !BandPerSample uses SAvg. In both
                     * cases exclusively. */
                    if (inst->isCapped && iPath==iCap) {
                        // treat the cap interp
                        interpLevel = levels[iCap];
                    } else {
                        // treat for possibly 2 interps for floor - at
                        // strike level
                        interpLevel = (nbAvgOut * levels[iPath] 
                                       - nbPast * inst->floor
                                       - bandedAvgSoFar[iPath]
                                       + (inst->isCapped?
                                          bandedAvgSoFar[iCap]:0)) / 
                            numRemaining;
                        interpLevel = Maths::max(inst->floor, interpLevel);
                    }
                }
                else
                {
                    interpLevel = (nbAvgOut * levels[iPath] - 
                                   avgSoFar[iPath][iAsset])/numRemaining;
                }
                // we have moneyness so now translate to absolute level
                // scale for each component - smileType D ...
                if (inst->smileType == XCB::FIXED_STRIKE_SMILE) {
                    interpLevel *= pathGen->refLevel(iAsset, iPath);
                } else if (inst->smileType == XCB::FLOAT_STRIKE_SMILE) {
                    /*** XXX SN - not sure how correct this is! 
                         Esp during average-in! ***/
                    // Essentially we take % into component via spot, not SAS
                    interpLevel *= getMultiFactors()->assetGetSpot(iAsset);
                    // ... but
                    // the % interp is relative to baskRef, so 
                    // need also to "rebase" to baskSpot
                    double assetSpot;
                    double baskRef = 0;
                    double baskSpot = 0;
                    for(int i=0;i<getNumAssets(); i++) {
                        assetSpot = getMultiFactors()->assetGetSpot(i);
                        if (inst->isUnitBasket) {
                            baskRef += inst->weights[i] *
                                pathGen->refLevel(i, iPath);
                            baskSpot += inst->weights[i] * assetSpot;
                        } else {
                            baskSpot += inst->weights[i] * assetSpot / 
                                pathGen->refLevel(i, iPath);
                        }
                    }
                    if (inst->isUnitBasket) {
                        interpLevel *= baskRef / baskSpot;
                    } else {
                        // baskRef is 1
                        interpLevel /= baskSpot;
                    }
                }
            }
            if (inst->smileType == XCB::NO_SMILE) {
                reqarr[iPath] = CVolRequestLNSP(new ATMVolRequest());
            } else {
                const SimSeries* simSeries = getSimSeries();
                reqarr[iPath] = CVolRequestLNSP(
                    new LinearStrikeTSVolRequest(interpLevel,
                                                 startDate,
                                                 simSeries->getLastDate(),
                                                 fwdStarting));
            }
        } // iPath
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* BandedBaskAv::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets()));//create empty one
    simSeries->addDates(averageOutDates);
    return new BandedBaskAvMC(this, simSeries);
};

CClassConstSP const BandedBaskAv::TYPE = CClass::registerClassLoadMethod(
    "BandedBaskAv", typeid(BandedBaskAv), BandedBaskAv::load);

const string BandedBaskAv::UNIT_WEIGHTS = "U";
const string BandedBaskAv::PCT_WEIGHTS = "P";

// * for class loading (avoid having header file) */
bool BandedBaskAvLoad() {
    return (BandedBaskAv::TYPE != 0);
}





DRLIB_END_NAMESPACE

