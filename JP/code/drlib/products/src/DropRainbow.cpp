//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DropRainbow.cpp
//
//   Description : Port of DropRainbow without pick feature. 
//                 At selected dates can drop a number of the best/worst
//                 performing assets, and a rainbow basket is formed from those
//                 not dropped. Option on that. Extends slightly by allowing
//                 drop from middle of asset perfs
//
//   Date        : Oct 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

class DropRainbow: public GenericNFBase, 
                   virtual public IMCIntoProduct {
protected:
    /// fields ////////
    DoubleArray             weights;            // [nbAssets-nbDroppedAssets]
    bool                    isCall;
    double                  strike;
    DateTimeArray           averageOutDates;
    IntArray                averageFlags;       /* 0: AvgOut, 1: drop date*/
// This provides the "spot perf" idea too, since if avgFromStart==false and there 
// are only drop dates it is automatically a spot perf.
    bool                    avgFromStart; 
    int                     nbDropBest;
    int                     nbDropMid;
    int                     nbDropWorst;
    bool                    isMidBest;
    int                     nbDropDates;        // transient

    enum FlagType {
        AVGDATE = 0,
        DROPDATE };

public:
    static CClassConstSP const TYPE;
    friend class DropRainbowMC;

    // validation
    void validatePop2Object(){
        static const string routine("DropRainbow::validatePop2Object");
        GenericNFBase::validatePop2Object();
        Maths::checkNonNegative(strike, "strike");

        // validate dates are not empty and are in order 
        DateTime::ensureIncreasing(averageOutDates, "averageOutDates", true);

        // Count number of drop dates and require that last date must be drop
        if (averageOutDates.size() != averageFlags.size()) {
            throw ModelException(routine, "Must have equal number of dates ("+
                                 Format::toString(averageOutDates.size())+
                                 ") and flags ("+
                                 Format::toString(averageFlags.size())+
                                 ")");
        }
        int iDate;
        // Note we allow 0 drop dates
        nbDropDates = 0;
        for (iDate = 0; iDate < averageOutDates.size(); iDate++) {
            if (averageFlags[iDate] != DROPDATE &&
                averageFlags[iDate] != AVGDATE) {
                throw ModelException(routine,
                                     "Flag number " + 
                                     Format::toString(iDate+1) + 
                                     " has value " + 
                                     Format::toString(averageFlags[iDate]) +
                                     ". Must be either " +
                                     Format::toString(AVGDATE) +
                                     " or " + Format::toString(DROPDATE));

            }
            if (averageFlags[iDate] == DROPDATE) {
                nbDropDates++;
            }
        }

        Maths::checkNonNegative(nbDropBest, "nbDropBest");
        Maths::checkNonNegative(nbDropMid, "nbDropMid");
        Maths::checkNonNegative(nbDropWorst, "nbDropWorst");

        int nbDropped = (nbDropBest+nbDropWorst+nbDropMid)*nbDropDates;
        if (nbDropped > assets->NbAssets()) {
            throw ModelException(routine, "(nbDropBest[" + Format::toString(nbDropBest) + 
                                 "] + nbDropWorst[" + Format::toString(nbDropWorst) +
                                 "] + nbDropMid[" + Format::toString(nbDropMid) +
                                 "]) * nbDropDates[" + Format::toString(nbDropDates) +
                                 "] = " + Format::toString(nbDropped) +
                                 " which must not exceed number of assets (" +
                                 Format::toString(assets->NbAssets()) + ")");
        }
        int nbRemaining = assets->NbAssets() - nbDropped;
        // check num weights 
        if (weights.size() != nbRemaining) {
            throw ModelException(routine, "Require #weights (" + Format::toString(weights.size()) +
                                 ") equal to #non-dropped assets (" + Format::toString(nbRemaining) +
                                 ")");
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
    DropRainbow(): GenericNFBase(TYPE), nbDropDates(0) {} // for reflection
    DropRainbow(const DropRainbow& rhs); // not implemented
    DropRainbow& operator=(const DropRainbow& rhs); // not implemented

    static IObject* defaultDropRainbow(){
        return new DropRainbow();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DropRainbow, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultDropRainbow);
        FIELD(weights,          "Weights");
        FIELD(isCall,           "is it a call option");
        FIELD(strike,           "strike");
        FIELD(averageOutDates,  "averageOutDates");
        FIELD(averageFlags,     "average flags: 0 - avg, 1 - avg&drop");
        FIELD(avgFromStart,     "false=> avg from previous drop date");
        FIELD(nbDropBest,       "nbDropBest");
        FIELD(nbDropMid,        "nbDropMid");
        FIELD(nbDropWorst,      "nbDropWorst");
        FIELD(isMidBest,        "isMidBest");
        FIELD(nbDropDates,      "nbDropDates");
        FIELD_MAKE_TRANSIENT(nbDropDates);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
class DropRainbowMC : public IMCProduct,
                      virtual public IMCProductLN {
private:
    int                      nbAssets;            // convenient
    const DropRainbow*       inst;                // reference to original instrument
    DoubleArray              sum;                 // per asset - class member saves alloc in payoff
    IntArray                 nbAvgOutPerDropDate; // [nbDropDates+1] the final being nbAvgDates since last drop
    int                      nbDroppedPerDate;
    int                      midOffset;           // convenient index offset to locate mid-droppers

    struct IndexedPerfs {
        double perf;
        int    iAsset;
    };
    struct IndexedPerfs      swap;

    // These preserve values from past to future
    DoubleArray              sumSoFar;
    int                      iDropDateSoFar;
    int                      nbDroppedSoFar;

    vector<IndexedPerfs>     perfs;
    class OrderPerfs {
    public:
        int operator() (const IndexedPerfs& p1, const IndexedPerfs& p2) {
            return p1.perf > p2.perf;
        }
    };

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    DropRainbowMC(const DropRainbow*       inst,
                  const SimSeriesSP&       simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        nbAssets(getNumAssets()),
        inst(inst),
        sum(nbAssets),
        nbAvgOutPerDropDate(inst->nbDropDates+1, 0),
        nbDroppedPerDate(inst->nbDropBest + inst->nbDropMid + inst->nbDropWorst),
        midOffset(inst->nbDropMid-2+(inst->isMidBest?0:1)),
        sumSoFar(nbAssets),
        iDropDateSoFar(0),
        nbDroppedSoFar(0),
        perfs(nbAssets){

        static const string routine("DropRainbowMC::DropRainbowMC");
        int iAsset;

        // Start the indexing for the perfs 
        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            perfs[iAsset].iAsset = iAsset;
        }

        int iStep;
        int iDrop = 0;
        for(iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            nbAvgOutPerDropDate[iDrop]++;
            if(inst->averageFlags[iStep] == DropRainbow::DROPDATE) {
                iDrop++;
            }
        }
        // cumulative number of AvgOut Dates
        if(inst->avgFromStart) {
            for(iDrop = 1; iDrop < inst->nbDropDates+1; iDrop++) {
                nbAvgOutPerDropDate[iDrop] += nbAvgOutPerDropDate[iDrop-1];
            }
        } else if (nbAvgOutPerDropDate[inst->nbDropDates]==0) {
            // Special case when final sample date is also a drop date
            // so we re-record nb samples for averaging
            nbAvgOutPerDropDate[inst->nbDropDates] = nbAvgOutPerDropDate[inst->nbDropDates-1];
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

        int    iAsset, iPerf;
        int    beginIdx = pathGen->begin(0); // same for all assets
        int    endIdx   = pathGen->end(0);

        int    iDropDate = iDropDateSoFar;
        int    nbDropped = nbDroppedSoFar;
        sum = sumSoFar;
        
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            // Only need to worry about those not yet dropped
            for (iPerf=0; iPerf<nbAssets-nbDropped; iPerf++) {

                iAsset = perfs[iPerf].iAsset;
                sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                
                if (inst->averageFlags[iStep] == DropRainbow::DROPDATE) {
                    perfs[iPerf].perf = sum[iAsset] / pathGen->refLevel(iAsset, 0) / 
                        nbAvgOutPerDropDate[iDropDate];
                    if (!inst->avgFromStart) {
                        sum[iAsset] = 0.0;
                    }
                }
            }

            // Hereafter all is perf related, not asset. We're concerned with best/worst and not
            // which asset it actually is.
            if (inst->averageFlags[iStep] == DropRainbow::DROPDATE) {
                // The dropped items will be shunted to the end of the array, and 
                // successive sorting will be done for the front not-yet-dropped ones.
                // Sorting is best to worst.
                sort(perfs.begin(), perfs.end()-nbDropped, OrderPerfs());
                // "Rotate" the dropped assets into appropriate places at the end of the array
                // The worst are already in position.
                // Deal with mid-droppers :-
                int  count;
                int  nbLeft = nbAssets - nbDropped; // we start this drop date with this many perfs
                // will move backwards with both these indexes (more natural since  
                // we're shifting entries to the end of the array)
                int  iWrite = nbLeft-1-inst->nbDropWorst; 
                int  iRead = (nbLeft+midOffset)/2; // integer arithmetic
                for(count=0;count<inst->nbDropMid;count++, iRead--, iWrite--) {
                    swap = perfs[iWrite];
                    perfs[iWrite] = perfs[iRead];
                    perfs[iRead] = swap;
                }
                // Finally drop any best performers :-
                // iWrite continues from loop above
                iRead = inst->nbDropBest - 1;
                for(count=0;count<inst->nbDropBest;count++, iRead--, iWrite--) {
                    swap = perfs[iWrite];
                    perfs[iWrite] = perfs[iRead];
                    perfs[iRead] = swap;
                }
                // This is a deterministic Dropping - we Drop the same number each Drop date
                nbDropped += nbDroppedPerDate;
                iDropDate++;
            }
        }

        if (pathGen->doingPast()){ 
            // preserve values
            sumSoFar = sum;
            iDropDateSoFar = iDropDate;
            nbDroppedSoFar = nbDropped;
        }
        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff. 
            // Do this only when we have a "complete" situation : either doingPast() and all is past, or 
            // !doingPast().
            // The perfs for the final basket are of the remaining assets perfs. Compute these if
            // the final sample date is not a drop date.
            for (iPerf=0; inst->averageFlags[endIdx-1] != DropRainbow::DROPDATE && 
                     iPerf<nbAssets-nbDropped; iPerf++) {
                iAsset = perfs[iPerf].iAsset;
                perfs[iPerf].perf = sum[iAsset] / pathGen->refLevel(iAsset, 0) / 
                    nbAvgOutPerDropDate[inst->nbDropDates];
            }
            // The non-dropped perfs are all collected at the start of perfs[].
            // We need one final sort for the rainbow of the basket. 
            sort(perfs.begin(), perfs.end()-nbDropped, OrderPerfs());
            double basket = 0.;
            for (iPerf=0; iPerf<inst->weights.size(); iPerf++) {
                basket += inst->weights[iPerf] * perfs[iPerf].perf;
            }
            double payoff = inst->isCall? (basket - inst->strike) : (inst->strike - basket);
            prices.add(payoff>0.0 ? inst->notional * payoff : 0.0);
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
        if (fwdStarting){
            interpLevel = inst->strike;
        } else {
            /* not forward starting - some samples have fixed already
               (this includes averaging in) */
            if (inst->avgFromStart) {
                int numDates = inst->averageOutDates.size();
                int numRemaining = 
                    today.numFutureDates(inst->averageOutDates);
                double soFar = sumSoFar[iAsset] / pathGen->refLevel(iAsset, 0);
                // moneyness is from basket levels
                interpLevel = (numDates * inst->strike - soFar)/ numRemaining;
            } else {
                interpLevel = inst->strike;
            }
            // here we're just using iPath = 0 since at this point
            // all the paths are the same
            interpLevel *= pathGen->refLevel(iAsset, 0);
        }
        const SimSeries* simSeries = getSimSeries();
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                 startDate,
                                                                 simSeries->getLastDate(),
                                                                 fwdStarting));
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* DropRainbow::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new DropRainbowMC(this, simSeries);
}

CClassConstSP const DropRainbow::TYPE = CClass::registerClassLoadMethod(
    "DropRainbow", typeid(DropRainbow), DropRainbow::load);

// * for class loading (avoid having header file) */
bool DropRainbowLoad() {
    return (DropRainbow::TYPE != 0);
}

DRLIB_END_NAMESPACE
