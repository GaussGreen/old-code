//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Rotator.cpp
//
//   Description : Variation on BaskAv where there's adjustment 
//                 of basket composition
//
//   Date        : July 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/**  */
class Rotator: public GenericNFBase, 
              virtual public IMCIntoProduct {
protected:
    /// fields ////////
    bool                    isCall;
    double                  strike;
    DateTimeArray           averageOutDates;
    DateTimeArray           resetDates;
    DoubleArray             initialParticipations; // len = nbSubBaskets
    BoolArray               isResetSubBasket;      // [nbSubBaskets]
    DoubleArrayArray        subBasketWeights;      // Due to IMS interface : [iAsset][jSubBasket] 
                                                   // So each asset has a DoubleArray which indicates the 
                                                   // weight it has in each sub-basket.
    bool                    isRefAvgOut;

public:
    static CClassConstSP const TYPE;
    friend class RotatorMC;

    // validation
    void validatePop2Object(){
        static const string routine("Rotator::validatePop2Object");
        GenericNFBase::validatePop2Object();

        if (Maths::isNegative(strike)){
            throw ModelException(routine, "strike ("+
                                 Format::toString(strike)+") cannot be negative");
        }
        // validate dates are not empty - order is handled by SimSeries
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates supplied!");
        }
        // Probably a mistake to try to reset after last averaging date
        if (!resetDates.empty() &&
            resetDates.back()>averageOutDates.back()) {
            throw ModelException(routine, "Final reset date " + resetDates.back().toString() +
                                 " cannot be after final averaging out date" +
                                 averageOutDates.back().toString());
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
        return DateTime::merge(averageOutDates, resetDates);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    Rotator(): GenericNFBase(TYPE), isRefAvgOut(false) {} // for reflection
    Rotator(const Rotator& rhs); // not implemented
    Rotator& operator=(const Rotator& rhs); // not implemented

    static IObject* defaultRotator(){
        return new Rotator();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Rotator, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultRotator);
        FIELD(isCall,                 "is it a call option, else put");
        FIELD(strike,                 "strike");
        FIELD(averageOutDates,        "Averaging out dates");
        FIELD(resetDates,             "Dates for reset");
        FIELD(initialParticipations,  "Starting/baseline participations");
        FIELD(isResetSubBasket,       "Per sub-basket. True=>resetting, False=>No resetting");
        FIELD(subBasketWeights,       "Array per asset for each sub-baskets definition");
        FIELD(isRefAvgOut,            "Else standard case.");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for Rotator - with isRefAvgOut false */
class RotatorMC : public IMCProduct,
                 virtual public IMCProductLN {
private:
    const Rotator*           inst;         // reference to original instrument
    int                      nbAssets;     // convenient
    int                      nbSubBaskets; // ditto
    double                   sumOutSoFar;
    DoubleArray              refLevel;     // [nbAssets] - some reset, some don't
    DoubleArray              refLevelSoFar; // as above but for past
    DoubleArray              p;            // participations[nbSubBaskets] - some reset, some don't
    DoubleArray              pSoFar;       //  as above but for past
    IntArray                 avgMap;       // easy way to flag subset of dates
    IntArray                 resetMap;     // ditto
    // these purely for reset
    IntArray                 resetSubBaskets; // [#reset sub baskets], [j] = index of jth resetting subbasket
    IntArray                 resetAssets; // [#reset assets], [k] = index of kth resetting asset
    IntArray                 noResetSubBaskets; // [#reset sub baskets], [j] = index of jth no-resetting subbasket
    IntArray                 noResetAssets; // [#reset assets], [k] = index of kth no-resetting asset
    DoubleArray              resetPartFactor; // [#reset sub baskets], [j] = p0[j-mapped] / Sum(p0)

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    RotatorMC(const Rotator*           inst,
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
        nbSubBaskets(inst->initialParticipations.size()),
        sumOutSoFar(0.0),
        refLevel(nbAssets, 0.0),
        refLevelSoFar(nbAssets, 0.0),
        p(inst->initialParticipations),
        pSoFar(inst->initialParticipations),
        resetSubBaskets(0),
        resetAssets(0),
        noResetSubBaskets(0),
        noResetAssets(0) {
        static const string routine("RotatorMC::RotatorMC");

        int  iAsset, i;
        bool isTrivial;
        resetMap = DateTime::createMapping(simSeries->getAllDates(),
                                           inst->resetDates,
                                           isTrivial);
        avgMap = DateTime::createMapping(simSeries->getAllDates(),
                                         inst->averageOutDates,
                                         isTrivial);
        
        // separate into more compute-friendly forms 
        // which allow us to loop through reset, and no reset forms independently
        // Also a natural place to check cum of weights for sub-baskets
        for(int iSubBask=0; iSubBask<nbSubBaskets; iSubBask++) {
            double sumWeights = 0.0;
            if (inst->isResetSubBasket[iSubBask]) {
                // each asset here with non-0% weight is a "resetting asset"
                resetSubBaskets.push_back(iSubBask);
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    if (!Maths::isZero(inst->subBasketWeights[iAsset][iSubBask])) {
                        resetAssets.push_back(iAsset);
                    }
                    sumWeights += inst->subBasketWeights[iAsset][iSubBask];
                }
            } else {
                 // and these are the "no resetting" varieties
                noResetSubBaskets.push_back(iSubBask);
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    if (!Maths::isZero(inst->subBasketWeights[iAsset][iSubBask])) {
                        noResetAssets.push_back(iAsset);
                    }
                    sumWeights += inst->subBasketWeights[iAsset][iSubBask];
                }
            }
            // We're implcitily percentage-type baskets
            // check weights add up to 100% (or almost anyway) - same check as in 
            // AssetUtil (if you're wondering where the magic tolerance came from)
            if (!Maths::areEqualWithinTol(sumWeights, 1.0,  0.000000001)){
                throw ModelException(routine, "Total % weights ("+
                                     Format::toString("%.12f", sumWeights)+") must "
                                     "equal 1.0");
            }
        }
        // check consistency - i.e. no assets marked for reset AND no-reset
        // This is a slight weakness in the spec of the interface (would be
        // nicer to not allow such to be possible, but that is more complicated to capture)
        for(int iResetAsset=0; iResetAsset<resetAssets.size(); iResetAsset++) {
            for(int iNoResetAsset=0; iNoResetAsset<noResetAssets.size(); iNoResetAsset++) {
                if (resetAssets[iResetAsset] == noResetAssets[iNoResetAsset]) {
                    throw ModelException(routine, "Asset #" + 
                                         Format::toString(resetAssets[iResetAsset]) +
                                         " is indicated as in both resetting and non-resetting sub-baskets");
                }
            }
        }

        // The reset requires a scaling factor for participations :-
        resetPartFactor = DoubleArray(resetSubBaskets.size(), 0.0);
        double resetPartSum = 0.0;
        for(i=0; i<resetSubBaskets.size(); i++) {
            int iSubBask = resetSubBaskets[i];
            resetPartSum += inst->initialParticipations[iSubBask];
        }
        // check resetPartSum non-zero ... be more clever about this
        // to allow all non-resetting sub-baskets. It should be safe since while
        // resetPartFactor is ill-defined, it would not be used under those 
        // circumstances
        if (resetSubBaskets.size()>0 &&
            Maths::isZero(resetPartSum)) {
            throw ModelException(routine, "Sum of resetting participations is 0!");
        }
        for(i=0; i<resetSubBaskets.size(); i++) {
            int iSubBask = resetSubBaskets[i];
            // note the resetPartfactor is indexed directly by [i]
            resetPartFactor[i] = inst->initialParticipations[iSubBask] / resetPartSum;
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static const string routine("RotatorMC::payoff");
        try {
            int    beginIdx = pathGen->begin(0); // same for all assets
            int    endIdx   = pathGen->end(0);
            int    i, j;

            double sumOut   = sumOutSoFar; 
            // refLevel (one per asset) starts being the standard but may get reset
            // We use the first SoFar value to tell us what state we're in
            if (Maths::isPositive(refLevelSoFar[0])) {
                refLevel = refLevelSoFar;
            } else {
                for(i=0;i<nbAssets;i++) {
                    refLevel[i] = pathGen->refLevel(i, 0);
                }
            }
            p = pSoFar;

            for (int iStep=beginIdx; iStep<endIdx; iStep++) {
                // each step we either get a sample of the resetting basket
                // or a reset, or both. Since any sim date is one or the other
                // and we always need the resetBask sum do that unconditionally.
                double resetBask = 0.0;
                for(i=0; i<resetSubBaskets.size(); i++) {
                    int iSubBask = resetSubBaskets[i];
                    double subBask = 0.0;
                    for (j=0; j<resetAssets.size(); j++) {
                        int iAsset = resetAssets[j];
                        double S = pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                        double wt = inst->subBasketWeights[iAsset][iSubBask];
                        subBask += p[iSubBask] * wt * S / refLevel[iAsset];
                    }           
                    resetBask += subBask;
                }
                if (avgMap[iStep]==0) {
                    double noResetBask = 0.0;
                    for(i=0; i<noResetSubBaskets.size(); i++) {
                        int iSubBask = noResetSubBaskets[i];
                        double subBask = 0.0;
                        for (j=0; j<noResetAssets.size(); j++) {
                            int iAsset = noResetAssets[j];
                            double S = pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                            double wt = inst->subBasketWeights[iAsset][iSubBask];
                            subBask += p[iSubBask] * wt * S / refLevel[iAsset];
                        }           
                        noResetBask += subBask;
                    }
                    sumOut += resetBask + noResetBask;
                }

                // handle any reset - afterwards note
                // Early validation ensures that there is a reset date before
                // any averaging date, so lastResetBask is well-defined from here.
                if (resetMap[iStep]==0) {
                    // update refLevels[] for relevant assets. Done here
                    // 'cos it's easier with the path level current for iStep.
                    for (j=0; j<resetAssets.size(); j++) {
                        int iAsset = resetAssets[j];
                        refLevel[iAsset] = pathGen->Path(iAsset, 0/*iPath*/)[iStep];
                    }
                    // Update participations based on reset-only basket level now
                    // These parts define basket at subsequent sample dates
                    for(i=0; i<resetSubBaskets.size(); i++) {
                        int iSubBask = resetSubBaskets[i];
                        // Note resetPartFactor indexed directly by [i]
                        p[iSubBask] = resetPartFactor[i] * resetBask; 
                    }
                }
            }
            if (pathGen->doingPast()){ // preserve values
                sumOutSoFar = sumOut;
                // note this is the only place refLevelSoFar is set so 
                // can use it to flag if we _were_ called for the past
                refLevelSoFar = refLevel; 
                pSoFar = p;
            }

            if (!pathGen->doingPast() || !hasFuture()) {
                // Compute a payoff, but only when we have a "complete" situation : either 
                // doingPast() and all is past, or !doingPast().
                double perf = sumOut / inst->averageOutDates.size();
                double myPayoff;
                if (inst->isRefAvgOut) {
                    // only call supported
                    myPayoff = 1.0 - inst->strike / perf;
                } else {
                    myPayoff = inst->isCall? 
                        (perf - inst->strike): (inst->strike - perf);
                }
                prices.add(inst->notional * prices.maxWithZero(myPayoff));
            }
        }
        catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    // for the LogNormal path generator. Think about details later!
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1);
        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool fwdStarting = startDate.isGreater(today);
        double interpLevel = fwdStarting? inst->strike : (inst->strike * pathGen->refLevel(iAsset, 0));

        const SimSeries* simSeries = getSimSeries();
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                 startDate,
                                                                 simSeries->getLastDate(),
                                                                 fwdStarting));
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Rotator::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    simSeries->addDates(resetDates);
    return new RotatorMC(this, simSeries);
}

CClassConstSP const Rotator::TYPE = CClass::registerClassLoadMethod(
    "Rotator", typeid(Rotator), Rotator::load);

// * for class loading (avoid having header file) */
bool RotatorLoad() {
    return (Rotator::TYPE != 0);
}

DRLIB_END_NAMESPACE
