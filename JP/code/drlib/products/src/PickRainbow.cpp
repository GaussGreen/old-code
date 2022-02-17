//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PickRainbow.cpp
//
//   Description : Port of DropRainbow and PickRainbowFloor models from EDG
//                 At selected dates can pick or drop a number of the best/worst
//                 performing assets, and a rainbow basket is formed from those
//                 picked/dropped. Option on that.
//
//   Date        : Sep 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Check.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

class PickRainbow: public GenericNFBase, 
                   virtual public IMCIntoProduct {
protected:
    /// fields ////////
    DoubleArray             weights;            // [nbPickedAssets]
    bool                    isCall;
    double                  strike;
    DateTimeArray           averageOutDates;
    DoubleArray             floors;             // [nbPickDates], allowed optional for ease of testing
    bool                    isCapped;
    DoubleArray             caps;               // [nbPickDates], allowed optional for ease of testing
    IntArray                averageFlags;       /* 0: AvgOut, 1: pick date*/
// This provides the "spot perf" idea too, since if avgFromStart==false and there 
// are only pick dates it is automatically a spot perf.
    bool                    avgFromStart; 
    int                     nbPickBest;
    int                     nbPickWorst;
    int                     nbPickDates;        // transient

    enum FlagType {
        AVGDATE = 0,
        PICKDATE };

public:
    static CClassConstSP const TYPE;
    friend class PickRainbowMC;

    // validation
    void validatePop2Object(){
        static const string routine("PickRainbow::validatePop2Object");
        GenericNFBase::validatePop2Object();
        Maths::checkNonNegative(strike, "strike");

        // validate dates are not empty and are in order 
        DateTime::ensureIncreasing(averageOutDates, "averageOutDates", true);

        // Count number of pick dates and require that last date must be pick
        if (averageOutDates.size() != averageFlags.size()) {
            throw ModelException(routine, "Must have equal number of dates ("+
                                 Format::toString(averageOutDates.size())+
                                 ") and flags ("+
                                 Format::toString(averageFlags.size())+
                                 ")");
        }
        int iDate;
        nbPickDates = 0;
        for (iDate = 0; iDate < averageOutDates.size(); iDate++) {
            if (averageFlags[iDate] != PICKDATE &&
                averageFlags[iDate] != AVGDATE) {
                throw ModelException(routine,
                                     "Flag number " + 
                                     Format::toString(iDate+1) + 
                                     " has value " + 
                                     Format::toString(averageFlags[iDate]) +
                                     ". Must be either " +
                                     Format::toString(AVGDATE) +
                                     " or " + Format::toString(PICKDATE));

            }
            if (averageFlags[iDate] == PICKDATE) {
                nbPickDates++;
            }
        }

        if (nbPickDates<1) {
            throw ModelException(routine, "Must have at least 1 pick date!");
        }

        if (averageFlags[averageOutDates.size()-1] != PICKDATE) {
            throw ModelException(routine, "Last flag must always be a Pick Date");
        }

        // for ease of using ported regression tests allow floors to be optional
        // and catch here - create a 0% floor
        if (floors.size() == 0) {
            floors = DoubleArray(nbPickDates, 0.0);
        }

        if (floors.size() != nbPickDates) {
            throw ModelException(routine, "#floors (" + Format::toString(floors.size()) +
                                 ") must equal nbPickDates (" + Format::toString(nbPickDates) +")");
        }

        // caps per pick date; flagged via isCapped
        if (isCapped) {
            if (caps.size() != nbPickDates) {
                throw ModelException(routine, "#caps (" + Format::toString(caps.size()) +
                                     ") must equal nbPickDates (" + Format::toString(nbPickDates) +")");
            }
            for(int i=0; i<caps.size(); i++) {
                if (floors[i]>caps[i]) {
                    throw ModelException(routine, "At pick date #" + Format::toString(i+1) + " Floor (" +
                                         Format::toString(floors[i]) + ") > Cap (" + Format::toString(caps[i]) + ")");
                }
            }
        }

        Maths::checkNonNegative(nbPickBest, "nbPickBest");
        Maths::checkNonNegative(nbPickWorst, "nbPickBest");

        int nbPicked = (nbPickBest+nbPickWorst)*nbPickDates;
        if (nbPicked > assets->NbAssets()) {
            throw ModelException(routine, "(nbPickBest[" + Format::toString(nbPickBest) + 
                                 "] + nbPickWorst[" + Format::toString(nbPickWorst) +
                                 "]) * nbPickDates[" + Format::toString(nbPickDates) +
                                 "] = " + Format::toString(nbPicked) +
                                 " which must not exceed number of assets (" +
                                 Format::toString(assets->NbAssets()) + ")");
        }
        // check num weights and that they sum to 100%
        Check::percWeights(weights, nbPicked,
                           "total picked assets");
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
    PickRainbow(): GenericNFBase(TYPE), isCapped(false), nbPickDates(0) {} // for reflection
    PickRainbow(const PickRainbow& rhs); // not implemented
    PickRainbow& operator=(const PickRainbow& rhs); // not implemented

    static IObject* defaultPickRainbow(){
        return new PickRainbow();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PickRainbow, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultPickRainbow);
        FIELD(weights,          "Weights");
        FIELD(isCall,           "is it a call option");
        FIELD(strike,           "strike");
        FIELD(averageOutDates,  "averageOutDates");
        FIELD(averageFlags,     "average flags: 0 - avg, 1 - avg&pick");
        FIELD(avgFromStart,     "false=> avg from previous pick date");
        FIELD(nbPickBest,       "nbPickBest");
        FIELD(nbPickWorst,      "nbPickWorst");
        FIELD(floors,           "floors");
        FIELD_MAKE_OPTIONAL(floors);
        FIELD(isCapped,         "isCapped");
        FIELD_MAKE_OPTIONAL(isCapped);
        FIELD(caps,             "caps");
        FIELD_MAKE_OPTIONAL(caps);
        FIELD(nbPickDates,      "nbPickDates");
        FIELD_MAKE_TRANSIENT(nbPickDates);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
class PickRainbowMC : public IMCProduct,
                      virtual public IMCProductLN {
private:
    int                      nbAssets;            // convenient
    const PickRainbow*       inst;                // reference to original instrument
    DoubleArray              sum;                 // per asset - class member saves alloc in payoff
    IntArray                 nbAvgOutPerPickDate;
    int                      nbPickedTotal;       // convenient

    // These preserve values from past to future
    DoubleArray              sumSoFar;
    int                      iPickDateSoFar;
    int                      nbPickedBestSoFar;
    int                      nbPickedWorstSoFar;

    struct IndexedPerfs {
        double perf;
        int    iAsset;
    };
    vector<IndexedPerfs>     perfs;
    class OrderPerfs {
    public:
        int operator() (const IndexedPerfs& p1, const IndexedPerfs& p2) {
            return p1.perf > p2.perf;
        }
    };
    DoubleArray              myFloors;  // rearranged to index easier in payoff
    DoubleArray              myCaps;  // rearranged to index easier in payoff
    DoubleArray              finalPerfs;
    class OrderDoubles {
    public:
        int operator() (double d1, double d2) {
            return d1 > d2;
        }
    };

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    PickRainbowMC(const PickRainbow*       inst,
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
        nbAvgOutPerPickDate(inst->nbPickDates, 0),
        nbPickedTotal(inst->nbPickDates*(inst->nbPickBest+inst->nbPickWorst)),
        sumSoFar(nbAssets),
        iPickDateSoFar(0),
        nbPickedBestSoFar(0),
        nbPickedWorstSoFar(0),
        perfs(nbAssets),
        myFloors(nbPickedTotal),
        myCaps(nbPickedTotal),
        finalPerfs(nbPickedTotal){

        static const string routine("PickRainbowMC::PickRainbowMC");
        int iAsset;

        // Start the indexing for the perfs 
        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            perfs[iAsset].iAsset = iAsset;
        }

        int iStep;
        int iPick = 0;
        for(iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            nbAvgOutPerPickDate[iPick]++;
            if(inst->averageFlags[iStep] == PickRainbow::PICKDATE) {
                iPick++;
            }
        }
        // cumulative number of AvgOut Dates
        if(inst->avgFromStart) {
            for(iPick = 1; iPick < inst->nbPickDates; iPick++) {
                nbAvgOutPerPickDate[iPick] += nbAvgOutPerPickDate[iPick-1];
            }
        }

        // To ease use in payoff, map floors into an array (myFloors) that can be indexed by  
        // iPerf (that is, according to the position in the array of all picked asset perfs).
        // inst->floors is indexed by pick date, and the array of picked asset perfs is
        // organised with the best placed at the beginning counting up (in blocks of nbPickBest 
        // per pick date), and the worst at the end counting down (in blocks of nbPickWorst
        // per pick date). 
        // The best picked are the entries for iPerf = 0, 1, etc
        int iPerf;
        for(iPerf=0; iPerf<inst->nbPickDates*inst->nbPickBest; iPerf++) {
            // Each pick date nbPickBest are picked, so the mapping is easy
            iPick = iPerf / inst->nbPickBest;  // integer arithmetic
            myFloors[iPerf] = inst->floors[iPick];
            if (inst->isCapped) {
                myCaps[iPerf] = inst->caps[iPick];
            }
        }
        // The worst are entries from nbPickedTotal-1, nbPickedTotal-2, etc
        for(;  // continues from last loop
            iPerf<nbPickedTotal; 
            iPerf++) {
            // Counting down means the "picked" order is the reverse of 
            // the "perf" order, hence use of "-iPerf"
            iPick = (nbPickedTotal-1-iPerf) / inst->nbPickWorst;  // integer arithmetic
            myFloors[iPerf] = inst->floors[iPick];
            if (inst->isCapped) {
                myCaps[iPerf] = inst->caps[iPick];
            }
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

        int    iPickDate = iPickDateSoFar;
        int    nbPickedBest = nbPickedBestSoFar;
        int    nbPickedWorst = nbPickedWorstSoFar;
        sum = sumSoFar;
        
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            for (iPerf=nbPickedBest; iPerf<nbAssets-nbPickedWorst; iPerf++) {

                iAsset = perfs[iPerf].iAsset;
                sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                
                if (inst->averageFlags[iStep] == PickRainbow::PICKDATE) {
                    perfs[iPerf].perf = sum[iAsset] / pathGen->refLevel(iAsset, 0) / 
                        nbAvgOutPerPickDate[iPickDate];
                    if (!inst->avgFromStart) {
                        sum[iAsset] = 0.0;
                    }
                }
            }

            // Hereafter all is perf related, not asset. We're concerned with best/worst and not
            // which asset it actually is.
            if (inst->averageFlags[iStep] == PickRainbow::PICKDATE) {
                // "Picking" is achieved by shrinking the array of perfs from the appropriate end
                // Sorted best to worst so picking best moves the start of the array
                sort(perfs.begin()+nbPickedBest, perfs.end()-nbPickedWorst, OrderPerfs());
                // This is a deterministic picking - we pick the same number each pick date
                nbPickedBest+=inst->nbPickBest;
                nbPickedWorst+=inst->nbPickWorst;
                iPickDate++;
            }
        }

        if (pathGen->doingPast()){ 
            // preserve values
            sumSoFar = sum;
            iPickDateSoFar = iPickDate;
            nbPickedBestSoFar = nbPickedBest;
            nbPickedWorstSoFar = nbPickedWorst;
        }
        if (!pathGen->doingPast() || !hasFuture()) {
            // Compute a payoff. 
            // Do this only when we have a "complete" situation : either doingPast() and all is past, or 
            // !doingPast().
            // We need to bring the picked perfs into a contiguous block
            // At this stage we don't need the index information, so can just sort 
            // doubles. This is probably faster. Also, there is an implicit assumption
            // in the way we use perfs[] that the list of iAssets is complete for the next iteration
            // (we don't care about the actual order, just that every iAsset appears exactly once)
            for(iPerf=0; iPerf<nbPickedBest; iPerf++) {
                finalPerfs[iPerf] = Maths::max(myFloors[iPerf], perfs[iPerf].perf);
                if (inst->isCapped) {
                    finalPerfs[iPerf] = Maths::min(myCaps[iPerf], finalPerfs[iPerf]);
                }
            }
            int iPickedWorst = nbAssets-nbPickedWorst;
            for(iPerf=nbPickedBest; iPerf<nbPickedTotal; iPerf++, iPickedWorst++) {
                finalPerfs[iPerf] = Maths::max(myFloors[iPerf], perfs[iPickedWorst].perf);
                if (inst->isCapped) {
                    finalPerfs[iPerf] = Maths::min(myCaps[iPerf], finalPerfs[iPerf]);
                }
            }
            // Then we need one final sort for the rainbow of the basket. 
            sort(finalPerfs.begin(), finalPerfs.end(), OrderDoubles());
            double basket = 0.;
            for (iPerf=0; iPerf<nbPickedTotal; iPerf++) {
                basket += inst->weights[iPerf] * finalPerfs[iPerf];
            }
            double payoff = inst->isCall? (basket - inst->strike) : (inst->strike - basket);
            prices.add(payoff>0.0 ? inst->notional * payoff : 0.0);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        /* 1. Not sure how easy this adjustment is in EDR since we'd need to intercept 
           path generation, which properly means creating a new pathGen, or adding
           a new method. For now perhaps rely on using Implied, so all we need to
           really do is get the date info correct.  

           Past pickings adjustment applies to the whole rainbow basket, not per asset.
           This requires a guess at the weights of any final basket. We reuse the payoff 
           function after projecting each asset forward at spot level in a dummy path. 
           This is very similar to the PastFromProduct idea, but we take 
           from it the weights of the final rainbow basket. 

           2. I shall do the adjust for past averaging samples, just to ease 
           migration of tests. But the averaging from previous pick date will interp
           without any adjustment.
        */
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
IMCProduct* PickRainbow::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new PickRainbowMC(this, simSeries);
}

CClassConstSP const PickRainbow::TYPE = CClass::registerClassLoadMethod(
    "PickRainbow", typeid(PickRainbow), PickRainbow::load);

// * for class loading (avoid having header file) */
bool PickRainbowLoad() {
    return (PickRainbow::TYPE != 0);
}

DRLIB_END_NAMESPACE
