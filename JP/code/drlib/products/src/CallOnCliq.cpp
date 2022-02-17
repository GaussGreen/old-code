//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallOnCliq.cpp
//
//   Description : 
//
//   Date        : Feb 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/Format.hpp"
#include "edginc/IAggregate.hpp"
#include "edginc/ITaxableInst.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

class CallOnCliq: public GenericNFBase, 
                  virtual public IMCIntoProduct,
                  public ITaxableInst::Basic {
private:
    /// fields ////////
    double                       globalRedemption;
    double                       globalFloor;
    double                       globalCap;
    DoubleArray                  weights;
    string                       weightType;

    // per cliquet : if only 1 entry then assume value applies to all cliquets
    DoubleArray                  floors;
    DoubleArray                  strikesPct;
    DoubleArray                  caps;
    DoubleArray                  participations;

    bool                         isRefPrevAvgOut;
    bool                         isReinvest;
    bool                         isRainbowCliquet;
    DoubleArray                  rainbowWeights;

    DateTimeArray                averageOutDates;
    DateTimeArray                cliquetDates;

public:
    static CClassConstSP const TYPE;
    friend class CallOnCliqMC;

    // validation
    void validatePop2Object(){
        static const string routine = "CallOnCliq::validatePop2Object";
        GenericNFBase::validatePop2Object();
        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates given!");
        }
        if (cliquetDates.empty()) {
            throw ModelException(routine, "No cliquetDates given!");
        }

        // Check average dates and cliquet dates "cooperate"
        // Cliquet dates must be a subset of average dates
        if (!DateTime::isSubset(averageOutDates, cliquetDates)) {
            throw ModelException(routine, "Cliquet dates should be a subset of averageOutDates");
        }
        // Require some averaging before first cliquet date
        const DateTime& firstAvg = averageOutDates[0];
        const DateTime& firstCpn = cliquetDates[0];
        if (firstCpn < firstAvg) {
            throw ModelException(routine, "Cannot have cliquet date " + firstCpn.toString() + 
                                 " before first average date " + firstAvg.toString());
        }
        // Makes no sense to have averaging after final cliquet date
        const DateTime& lastAvg = averageOutDates[averageOutDates.size()-1];
        const DateTime& lastCpn = cliquetDates[cliquetDates.size()-1];
        if (lastAvg > lastCpn) {
            throw ModelException(routine, "Cannot average on " + lastAvg.toString() + 
                                 " since after final cliquet date " + lastCpn.toString());
        }

        // XXX should be removed to a common area as per AssetUtil::checkWeights XXX
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
            throw ModelException(routine, "weightType must be U or P,"
                                 " but " + weightType + " given");
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
    CallOnCliq(): GenericNFBase(TYPE), isRefPrevAvgOut(false) {} // for reflection
    CallOnCliq(const CallOnCliq& rhs); // not implemented
    CallOnCliq& operator=(const CallOnCliq& rhs); // not implemented

    static IObject* defaultCallOnCliq(){
        return new CallOnCliq();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CallOnCliq, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultCallOnCliq);
        FIELD(globalRedemption,         "Strike of the overall option");
        FIELD(globalFloor,              "Overall option floor, relative to strike, so <=0");
        FIELD(globalCap,                "Overall option cap, relative to strike, so >=0");
        FIELD(weights,                  "Weights for asset basket");
        FIELD(weightType,               "Weight Type for asset basket: U or P");
        FIELD(floors,                   "Per-cliquet floors, length = 1 or N, relative to strike, so <=0");
        FIELD(strikesPct,               "Per-cliquet strikes, length = 1 or N");
        FIELD(caps,                     "Per-cliquet caps, length = 1 or N, relative to strike, so >=0");
        FIELD(participations,           "Per-cliquet participations, length = 1 or N");
        FIELD(isRefPrevAvgOut,          "False=>Reference for cliquet is SPOT at previous cliquet maturity");
        FIELD(isReinvest,               "True=>Returns reinvest");
        FIELD(isRainbowCliquet,         "True=>Cliquet prices are rainbowed according to rainbowWeights");
        FIELD(rainbowWeights,           "Only used if isRainbowCliquet");
        FIELD(averageOutDates,          "Sample dates for Averaging Out");
        FIELD(cliquetDates,             "Cliquet maturity dates, subset of averageOutDates");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super CallOnCliq */
class CallOnCliqMC: public IMCProduct,
                         virtual public IMCProductLN {
private:
    const CallOnCliq*        inst;      // reference to original instrument
    int                      nbAssets;  // nicer
    DoubleArray              sum;       // [nbAssets], saves alloc 
    DoubleArray              refLevel;  // [nbAssets], saves alloc 
    bool                     isUnitBasket;
    SimpleDoubleArray        cliquets;

    // historical values
    DoubleArray              sumSoFar;      // [nbAssets]
    DoubleArray              refLevelSoFar; // [nbAssets]
    SimpleDoubleArray        cliquetsSoFar;  // [nbCliquets]
    int                      iCliquetSoFar;

    // Operational aggregation and performance calcs
    IAggregateSP             timeBasket;
    IDoubleArrayModifierSP   timeBasketComponents;
    IDoubleArrayModifierSP   overallOption;
    TrivialDoubleArray       calRbw;  // just a double ... may be more natural way to phrase this?

    PerfTypeBandedPerElementMakerSP timeBasketCompMaker;
    PerfTypeSimpleBandedMakerSP     overallOptionMaker;

    IntArray                 nbAvgOutPerCliquetDate; // [nbCliquets] 
    IntArray                 cliquetMap;    // [nbAvgDates] - convenient way to track cliquet dates

public:
    
    /** equivalent to InstIntoMCProduct */
    CallOnCliqMC(const CallOnCliq*         inst,
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
        sum(nbAssets, 0.0),
        refLevel(nbAssets, 0.0),
        isUnitBasket(inst->weightType=="U"),
        cliquets(inst->cliquetDates.size(), 0.0),
        sumSoFar(nbAssets, 0.0),
        refLevelSoFar(nbAssets, 0.0),
        cliquetsSoFar(inst->cliquetDates.size(), 0.0),
        iCliquetSoFar(0),
        calRbw(0.0),
        nbAvgOutPerCliquetDate(inst->cliquetDates.size(), 0) {

        static const string routine = "CallOnCliqMC::CallOnCliqMC";
        int numCliquets = inst->cliquetDates.size();
        
        // we check the first avg date is not after first cpn date
        int iCliquet = 0; 
        for(int iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            if(inst->averageOutDates[iStep] > inst->cliquetDates[iCliquet]) {
                iCliquet++;
            }
            nbAvgOutPerCliquetDate[iCliquet]++;
        }
        bool isTrivial;
        cliquetMap = DateTime::createMapping(inst->averageOutDates,
                                            inst->cliquetDates,
                                            isTrivial);

        /////////////////
        // "Global" option -> overallOption
        // Note the "globalRedemption - 1" for reinvest case. Then it's based around 100% which is
        // more natural from a user's perspective. The "pairing" is with the ReinvestAggregators
        // which do "Product(1 + xi) - 1". Alternative is to change the behaviour of the ReinvestAggregators
        // so they just do "Product(1 + xi)", but I wonder if that's less reusable. Only time will tell!
        PerfTypeSimpleBandedMaker overallOptionMaker(inst->globalFloor,
                                                     inst->globalRedemption - (inst->isReinvest?1.0:0.0), 
                                                     inst->globalCap,
                                                     1.0); // participation
        overallOption = IDoubleArrayModifierSP(overallOptionMaker.getModifier(&calRbw));

        /////////////////
        // Cliquet "modifiers" ->timeBasketComponents
        // Should all be same length regardless
        if ((inst->floors.size() != inst->strikesPct.size())  ||
            (inst->floors.size() != inst->caps.size())  ||
            (inst->floors.size() != inst->participations.size())) {
            throw ModelException(routine, 
                                 "Must be equal: #floors=" + Format::toString(inst->floors.size()) +
                                 ", #strikes=" + Format::toString(inst->strikesPct.size()) +
                                 ", #caps=" + Format::toString(inst->caps.size()) +
                                 ", #participations=" + Format::toString(inst->participations.size()));
        }
        if (inst->floors.size() < 1) {
            throw ModelException(routine, 
                                 "Require at least one floor level!");
        }
        if (inst->floors.size() == 1) {
            // it's a simple perf modifier we need

            PerfTypeSimpleBandedMaker timeBasketCompMaker(inst->floors[0],
                                                          inst->strikesPct[0],
                                                          inst->caps[0],
                                                          inst->participations[0]);
            timeBasketComponents = IDoubleArrayModifierSP(timeBasketCompMaker.getModifier(&cliquets));
        } else {
            // per-cliquet so check length matches num cliquet dates too
            if (inst->floors.size() != numCliquets) {
                throw ModelException(routine, 
                                     "#floors=" + Format::toString(inst->floors.size()) +
                                     " must equal #cliquets=" + Format::toString(numCliquets));
            }
            PerfTypeBandedPerElementMaker timeBasketCompMaker(inst->floors,
                                                              inst->strikesPct,
                                                              inst->caps,
                                                              inst->participations);
            timeBasketComponents = IDoubleArrayModifierSP(timeBasketCompMaker.getModifier(&cliquets));
        }

        /////////////////
        // Cliquet "aggregation" -> timeBasket
        if (inst->isReinvest) {
            if (inst->isRainbowCliquet) {
                ReinvestRainbowAggregateMaker am(inst->rainbowWeights);
                timeBasket = IAggregateSP(am.getAggregate(&cliquets));
            } else {
                DoubleArray one(numCliquets, 1.0);
                ReinvestAggregateMaker am(one);
                timeBasket = IAggregateSP(am.getAggregate(&cliquets));
            }
        } else {
            if (inst->isRainbowCliquet) {
                RainbowAggregateMaker am(inst->rainbowWeights);
                timeBasket = IAggregateSP(am.getAggregate(&cliquets));
            } else {
                DoubleArray one(numCliquets, 1.0);
                BasketAggregateMaker am(one);
                timeBasket = IAggregateSP(am.getAggregate(&cliquets));
            }
        }
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) {
        int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset;

        int iCliquet = iCliquetSoFar;
        sum = sumSoFar;
        refLevel = refLevelSoFar;
        cliquets = cliquetsSoFar;
        for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
            /* Special treatment for average in for first cliq */
            if (iCliquet == 0) {
                refLevel[iAsset] = pathGen->refLevel(iAsset, 0);
            } else {
                refLevel[iAsset] = refLevelSoFar[iAsset];
            }
        }

        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            bool isCliquetDate = (cliquetMap[iStep]==0);   // true iff a cliquet date

            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
            }

            if (isCliquetDate) {
                double perf = 0.0;
                double basketRef = 0.0;

                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    double avg = sum[iAsset] / nbAvgOutPerCliquetDate[iCliquet];

                    if (isUnitBasket) {
                        perf += inst->weights[iAsset] * avg;
                        basketRef += inst->weights[iAsset] * refLevel[iAsset];
                    } else {
                        perf +=  inst->weights[iAsset] * avg / refLevel[iAsset];
                    }

                    /* Reset ref level for next cliquet*/
                    if (inst->isRefPrevAvgOut) {
                        refLevel[iAsset] = avg;
                    } else {
                        refLevel[iAsset] = pathGen->Path(iAsset, 0)[iStep];
                    }
                    sum[iAsset] = 0.0;
                }
                if (isUnitBasket) {
                    perf /= basketRef;
                }
                // then capture the value of basket return this cliquet
                cliquets[iCliquet] = perf;
                // get ready for next cliquet
                iCliquet++;
            }
        }

        // preserve values for past - before we modify/sort the cliquets.
        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            refLevelSoFar = refLevel;
            cliquetsSoFar = cliquets;
            iCliquetSoFar = iCliquet;
        }
        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff, but only when we have a "complete" situation : either 
            // doingPast() and all is past, or !doingPast().
            // Now perform the time basket - apply perf modifiers
            timeBasketComponents->apply();
            // then make a basket of them
            calRbw() = timeBasket->aggregate();
            overallOption->apply();
            prices.add(inst->notional * calRbw()); 
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string routine = "CallOnCliqMC::getVolInterp";
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel = inst->strikesPct[0]; // very lazy but LN shouldn't really be going this way! XXX

        // get hold of the future strike dates
        int numLiveCliqs = inst->cliquetDates.size() - iCliquetSoFar;
        if (numLiveCliqs<=0) {
            throw ModelException(routine, "No future cliquets!?");
        }
        DateTimeArray liveCliqStartDates(numLiveCliqs);
        for (int iLiveCliq = 0; iLiveCliq < numLiveCliqs; iLiveCliq++){
            int iCliquet = iCliquetSoFar+iLiveCliq-1;
            liveCliqStartDates[iLiveCliq] = iCliquet<0?startDate:inst->cliquetDates[iCliquet];
        }

        // same strike levels per cliquet (but may need to adjust first one)
        DoubleArray  strikes(numLiveCliqs, interpLevel);
        if (!fwdStarting){
            // need to set first level to absolute strike - adjusted
            // additionally for any average out samples for this cliquet
            int iStep;
            const DateTime& thisCliquetStart = iCliquetSoFar==0?
                startDate:inst->cliquetDates[iCliquetSoFar-1];
            // find first avg date of this cliquet
            for(iStep = 0; iStep < inst->averageOutDates.size() && 
                    inst->averageOutDates[iStep] <= thisCliquetStart; iStep++) {
                ; // empty
            }
            // then walk through counting past avg dates in this cliq
            int numRemaining = nbAvgOutPerCliquetDate[iCliquetSoFar];
            for(; iStep < inst->averageOutDates.size() && 
                    inst->averageOutDates[iStep] <= today; iStep++) {
                numRemaining--;
            }
            if (numRemaining<=0) {
                // something wrong!
                throw ModelException(routine, "INTERNAL ERROR : numRemaining is " + Format::toString(numRemaining));
            }
            // Can't set up refLevel earlier, 'cos need PathGen. First cliq has standard ref level
            double refLevel =  iCliquetSoFar==0 ? pathGen->refLevel(iAsset, 0) : refLevelSoFar[iAsset];
            strikes[0] = (nbAvgOutPerCliquetDate[iCliquetSoFar] * refLevel * interpLevel
                          - sumSoFar[iAsset])/ numRemaining;
        }
        reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
                                                           liveCliqStartDates, 
                                                           lastSimDate,
                                                           strikes));
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CallOnCliq::createProduct(const MonteCarlo* model) const {

    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new CallOnCliqMC(this, simSeries);
}

CClassConstSP const CallOnCliq::TYPE = CClass::registerClassLoadMethod(
    "CallOnCliq", typeid(CallOnCliq), CallOnCliq::load);

// * for class loading (avoid having header file) */
bool CallOnCliqLoad() {
    return (CallOnCliq::TYPE != 0);
}

DRLIB_END_NAMESPACE





