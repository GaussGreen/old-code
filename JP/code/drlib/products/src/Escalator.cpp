//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Escalator.cpp
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
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include <algorithm>


DRLIB_BEGIN_NAMESPACE

// typedef array<DoubleArraySP, DoubleArray> DoubleArrayArray;

// CClassConstSP const DoubleArrayArray::TYPE = CClass::registerClassLoadMethod(
//    "DoubleArrayArray", typeid(DoubleArrayArray), load);

// typedef smartPtr<DoubleArrayArray> DoubleArrayArraySP;
// typedef smartConstPtr<DoubleArrayArray> DoubleArrayArrayConstSP;


/** Escalator product - option based upon sum of coupons depending on 
    the performance of n assets */
class Escalator: public GenericNFBase, 
                 virtual public IMCIntoProduct{
protected:
    /// fields ////////
    
    DateTimeArray           averageOutDates;
    DoubleMatrix            coupons;
    DoubleArray             couponsNoBreach;

    IntArray                averageFlags;       /* 0: AvgOut, 1: coupon date*/
    
    double                  barrier;
    double                  highKSpread;
    double                  lowKSpread;
    bool                    digInPast;

    double                  overallStrike;
    bool                    avgFromStart;
    bool                    dynamicSpreads;

    bool                    useGlobalFloor;
    int                     nbDropBest;
    int                     nbDropWorst;

    static const int         COUPONDATE;
    static const int         AVGDATE;


public:
    static CClassConstSP const TYPE;
    friend class EscalatorMC;

    // validation
    void validatePop2Object(){
        static const string method = "Escalator::validatePop2Object";
        GenericNFBase::validatePop2Object();
        try {
            if (averageFlags.size() != averageOutDates.size()){
                throw ModelException(method, "Must have equal number of"
                                     "dates ("+
                                     Format::toString(averageOutDates.size())+
                                     ") and flags ("+
                                     Format::toString(averageFlags.size())+
                                     ")");
            }
                
            int iDate, nbCouponDates = 0;
            for(iDate = 0; iDate < averageOutDates.size(); iDate++) {
                if(averageFlags[iDate] != Escalator::COUPONDATE &&
                   averageFlags[iDate] != Escalator::AVGDATE) {
                    throw ModelException(method,
                                         "Flag number " + 
                                         Format::toString(iDate+1) + 
                                         " has value " + 
                                         Format::toString(averageFlags[iDate]) +
                                         ". Must be either " +
                                         Format::toString(Escalator::AVGDATE) +
                                         " or " + Format::toString(Escalator::COUPONDATE));

                }
                if(averageFlags[iDate] == COUPONDATE) {
                    nbCouponDates += 1;
                }
            }
            
            if(averageFlags[averageOutDates.size()-1] != COUPONDATE) {
                throw ModelException("Last flag must always be a CouponDate");
            }
            
            if(nbCouponDates != couponsNoBreach.size()) {
                throw ModelException("Number of coupon dates is " +
                                     Format::toString(couponsNoBreach.size()) +
                                     " and number of average flags with value " +
                                     Format::toString(Escalator::COUPONDATE) + 
                                     " is " +
                                     Format::toString(nbCouponDates) +
                                     ". Must be equal.");
            }

            if(couponsNoBreach.size() != coupons.numRows()) {
                throw ModelException("Size of coupons-no-breach" + 
                                     Format::toString(couponsNoBreach.size()) +
                                     " must equal number of rows in coupons matrix " +
                                     Format::toString(coupons.numRows()));
            }

            if(assets->NbAssets() != coupons.numCols()) {
                throw ModelException("Number of assets " +
                                     Format::toString(assets->NbAssets()) +
                                     " must be equal to number of columns in coupons matrix " +
                                     Format::toString(coupons.numCols()));
            }
            
            DateTime::ensureIncreasing(averageOutDates, "averageOutDates", true);
            
            if(Maths::isNegative(barrier)) {
                throw ModelException("Barrier must be non-negative. Here "+
                                     Format::toString(barrier));
            }
            
            if(Maths::isNegative(lowKSpread)) {
                throw ModelException("Low spread must be non-negative. Here "+
                                     Format::toString(lowKSpread));
            }

            if(Maths::isNegative(highKSpread)) {
                throw ModelException("High spread must be non-negative. Here "+
                                     Format::toString(highKSpread));
            }

            if(Maths::isZero(lowKSpread) && Maths::isZero(highKSpread)) {
                throw ModelException("Both high and low spreads are zero. At least one must be positive.");
            }

            if(!dynamicSpreads && Maths::isNegative(barrier - lowKSpread)) {
                throw ModelException("Using normal spreads. Barrier is " +
                                     Format::toString(barrier) +
                                     " and low spread is " +
                                     Format::toString(lowKSpread) +
                                     ". Their difference is negative");
            }

            if (nbCouponDates*(nbDropBest+nbDropWorst)>assets->NbAssets()) {
                throw ModelException("(nbDropBest+nbDropWorst) * nbCouponDates (" 
                                     + Format::toString(nbCouponDates*(nbDropBest+nbDropWorst)) + 
                                     ") must be <= numAssets (" 
                                     + Format::toString(assets->NbAssets()) + ")");
            }

        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
        
        
        return;
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
    Escalator(): GenericNFBase(TYPE), overallStrike(0.0), avgFromStart(false), 
        dynamicSpreads(false), useGlobalFloor(false), 
        nbDropBest(0), nbDropWorst(0) {}         // for reflection
    Escalator(const Escalator& rhs);            // not implemented
    Escalator& operator=(const Escalator& rhs); // not implemented

    static IObject* defaultEscalator(){
        return new Escalator();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Escalator, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultEscalator);
        FIELD(averageOutDates,  "the average out dates");
        FIELD(coupons,          "coupons");
        FIELD(couponsNoBreach,  "couponsNoBreach");
        FIELD(averageFlags,     "average flags");
        FIELD(barrier,          "barrier");
        FIELD(highKSpread,      "highKSpread");
        FIELD(lowKSpread,       "lowKSpread");
        FIELD(digInPast,        "digInPast");
        FIELD(overallStrike,    "overallStrike");
        FIELD_MAKE_OPTIONAL(overallStrike);
        FIELD(avgFromStart,     "avgFromStart");
        FIELD_MAKE_OPTIONAL(avgFromStart);
        FIELD(dynamicSpreads,   "dynamicSpreads");
        FIELD_MAKE_OPTIONAL(dynamicSpreads);
        FIELD(useGlobalFloor,   "use global floor");
        FIELD_MAKE_OPTIONAL(useGlobalFloor);
        FIELD(nbDropBest,   "nbDropBest (default 0)");
        FIELD_MAKE_OPTIONAL(nbDropBest);
        FIELD(nbDropWorst,   "nbDropWorst (default 0)");
        FIELD_MAKE_OPTIONAL(nbDropWorst);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
class EscalatorMC : public IMCProduct,
                    virtual public IMCProductLN {
private:
    int                      numAssets;
    const Escalator*         inst; // reference to original instrument

    DoubleMatrix             incremCoupons;
    int                      nbCouponDates;

    DoubleArray              sumOut;
    class IndexedPerfs {
    public:
        double perf;
        int    iAsset;
    };
    vector<IndexedPerfs>            perfs;
    class orderPerfs {
    public:
        int operator() (const IndexedPerfs& p1, const IndexedPerfs& p2) {
            return p1.perf > p2.perf;
        }
    };
    int                      nbDroppedBest;
    int                      nbDroppedBestSoFar;
    int                      nbDroppedWorst;
    int                      nbDroppedWorstSoFar;
    DoubleArray              callSpreads;
    IntArray                 nbAvgOutPerCouponDate;

    int                      currentCouponDate;         // past
    double                   historicCouponPayment;     // past
    DoubleArray              historicSumOut;            // past

    double                   loBarrier;
    double                   hiBarrier;

    static const int         COUPONDATE;
    static const int         AVGDATE;
    static const double      TINY_SPREAD;

public:
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to SuperRainbow) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    EscalatorMC(const Escalator*               inst,
                const SimSeriesSP&             simSeries):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        numAssets(getNumAssets()),
        inst(inst), 
        nbCouponDates(inst->couponsNoBreach.size()),
        sumOut(numAssets), 
        perfs(numAssets), 
        nbDroppedBest(0),
        nbDroppedBestSoFar(0),
        nbDroppedWorst(0),
        nbDroppedWorstSoFar(0),
        callSpreads(numAssets),
        nbAvgOutPerCouponDate(nbCouponDates),
        currentCouponDate(0),
        historicCouponPayment(0.0),
        historicSumOut(numAssets),
        loBarrier(inst->barrier - inst->lowKSpread), 
        hiBarrier(inst->barrier + inst->highKSpread) {

        int iAsset;

        // Start the indexing for the perfs 
        for(iAsset=0; iAsset<numAssets; iAsset++) {
            perfs[iAsset].iAsset = iAsset;
        }

        incremCoupons = DoubleMatrix(numAssets + 1, nbCouponDates);
        
        // incremental coupons
        int iCouponDate;
        for(iCouponDate = 0; iCouponDate < nbCouponDates; iCouponDate++) {
            incremCoupons[0][iCouponDate] = inst->couponsNoBreach[iCouponDate];
        }
        
        double maxCouponIncrement = 0.0;
        for(iAsset = 0; iAsset < numAssets; iAsset++) {
            for(iCouponDate = 0; iCouponDate < nbCouponDates; iCouponDate++) {
                // incremental coupons
                if(iAsset == 0) {
                    incremCoupons[iAsset+1][iCouponDate] = inst->coupons[iAsset][iCouponDate] - 
                                                         inst->couponsNoBreach[iCouponDate];
                } else {
                    incremCoupons[iAsset+1][iCouponDate] = inst->coupons[iAsset][iCouponDate] - 
                                                           inst->coupons[iAsset-1][iCouponDate];
                }
                
                if(Maths::isNegative(incremCoupons[iAsset+1][iCouponDate])) {
                    throw ModelException("Coupon numbers " +
                                         Format::toString(iAsset) +
                                         " and " +
                                         Format::toString(iAsset+1) +
                                         " for coupon date number " +
                                         Format::toString(iCouponDate+1) +
                                         " are decreasing. Must be non-decreasing.");
                }

                if(Maths::isPositive(incremCoupons[iAsset+1][iCouponDate] - maxCouponIncrement)) {
                    maxCouponIncrement = incremCoupons[iAsset+1][iCouponDate];
                }
            }
        }

        if(inst->dynamicSpreads && 
           Maths::isNegative(inst->barrier - inst->lowKSpread * maxCouponIncrement)) {
            throw ModelException("Dynamic spreads. Barrier level of " +
                                 Format::toString(inst->barrier) + 
                                 " minus low spread of " +
                                 Format::toString(inst->lowKSpread) +
                                 " times maximum coupon increment " +
                                 Format::toString(maxCouponIncrement) +
                                 " is negative. Must be positive.");
        }
        
        int iStep;
        iCouponDate = 0;
        for(iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            nbAvgOutPerCouponDate[iCouponDate] += 1;
            if(inst->averageFlags[iStep] == COUPONDATE) {
                iCouponDate += 1;
            }
        }

        // cumulative number of AvgOut Dates
        if(inst->avgFromStart) {
            for(iCouponDate = 1; iCouponDate < nbCouponDates; iCouponDate++) {
                nbAvgOutPerCouponDate[iCouponDate] += nbAvgOutPerCouponDate[iCouponDate-1];
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

        int    beginIdx = pathGen->begin(0); // 0 <=same for all assets
        int    endIdx   = pathGen->end(0);
        int    iAsset, iPerf;

        sumOut = historicSumOut;
        nbDroppedBest = nbDroppedBestSoFar;
        nbDroppedWorst = nbDroppedWorstSoFar;
        double totalCouponPayment = historicCouponPayment;
        int couponDate = currentCouponDate;

        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            for (iPerf=nbDroppedBest; iPerf<numAssets-nbDroppedWorst; iPerf++) {
                iAsset = perfs[iPerf].iAsset;

                sumOut[iAsset] += pathGen->Path(iAsset, 0)[iStep];
                
                if(inst->averageFlags[iStep] == COUPONDATE) {
                    perfs[iPerf].perf = sumOut[iAsset] / pathGen->refLevel(iAsset, 0) / 
                        nbAvgOutPerCouponDate[couponDate];
                    if(!inst->avgFromStart) {
                        sumOut[iAsset] = 0.0;
                    }
                }
            }

            // Hereafter all is perf related, not asset. We're concerned with best/worst and not
            // which asset it actually is.
            // Also support all being dropped, when we contribute nothing to payoff
            if(inst->averageFlags[iStep] == COUPONDATE &&
               nbDroppedBest<numAssets-nbDroppedWorst) {
                // "Dropping" is achieved by shrinking the array of perfs from the appropriate end
                // Sorted best to worst so dropping best moves the start of the array
                sort(perfs.begin()+nbDroppedBest, perfs.end()-nbDroppedWorst, orderPerfs());
                for (iPerf=nbDroppedBest; iPerf < numAssets-nbDroppedWorst; iPerf++) {
                    if(pathGen->doingPast() && inst->digInPast) {
                        if(Maths::isNegative(perfs[iPerf].perf - inst->barrier)) {
                            callSpreads[iPerf] = 0.0;
                        } else {
                            callSpreads[iPerf] = 1.0;
                        }
                    } else {
                        double loStrike, hiStrike;
                        int incremCouponIndex;
                        if(inst->dynamicSpreads) {
                            incremCouponIndex = iPerf + 1;
                            loStrike = inst->barrier - inst->lowKSpread * 
                                       incremCoupons[incremCouponIndex][couponDate] * 100.0;

                            hiStrike = inst->barrier + inst->highKSpread * 
                                       incremCoupons[incremCouponIndex][couponDate] * 100.0;

                            if(Maths::equals(loStrike, hiStrike)) {
                                loStrike *= (1.0 - TINY_SPREAD);
                                hiStrike *= (1.0 + TINY_SPREAD);
                            }
                        } else {
                            loStrike = loBarrier;
                            hiStrike = hiBarrier;
                        }
                    
                        callSpreads[iPerf] = (Maths::max(perfs[iPerf].perf - loStrike, 0.0) - 
                                               Maths::max(perfs[iPerf].perf - hiStrike, 0.0)) / 
                                              (hiStrike - loStrike);
                    }
                }

                // Here need a 0-based index to reference incremCoupons, hence the shifting of iPerf
                totalCouponPayment += incremCoupons[0][couponDate];
                for (iPerf=0; iPerf < numAssets-nbDroppedBest-nbDroppedWorst; iPerf++) {
                    totalCouponPayment += incremCoupons[iPerf+1][couponDate] * callSpreads[iPerf+nbDroppedBest];
                }
                // Anything to drop? Check if the best perf so far contributed - if it 
                // did then drop as many as indicated.
                if (Maths::isPositive(callSpreads[nbDroppedBest])) {
                    // yes!
                    nbDroppedBest+=inst->nbDropBest;
                    nbDroppedWorst+=inst->nbDropWorst;
                }
                
                couponDate += 1;
            }
        }

        if (pathGen->doingPast()){ // preserve values
            historicSumOut = sumOut;
            currentCouponDate = couponDate;
            historicCouponPayment = totalCouponPayment;
            nbDroppedBestSoFar = nbDroppedBest;
            nbDroppedWorstSoFar = nbDroppedWorst;
        }
        
        double payoff;
        if(inst->useGlobalFloor) {
            payoff = Maths::max(totalCouponPayment - inst->overallStrike, 0.0);
        } else {
            payoff = totalCouponPayment;
        }

        prices.add(inst->notional * payoff);
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here

        const IRefLevel* refLevel = getRefLevel();
        const DateTime& startDate = refLevel->getAllDates().front();
        const DateTime& today = getToday();
        bool fwdStarting = startDate.isGreater(today);
        
        // just to be consistent with LN MC in EDG
        double interpLevel = inst->barrier + inst->highKSpread;
        if (!fwdStarting){
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
IMCProduct* Escalator::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new EscalatorMC(this, simSeries);
}

CClassConstSP const Escalator::TYPE = CClass::registerClassLoadMethod(
    "Escalator", typeid(Escalator), Escalator::load);

// * for class loading (avoid having header file) */
bool EscalatorLoad() {
    return (Escalator::TYPE != 0);
}

const int Escalator::COUPONDATE = 1;
const int Escalator::AVGDATE = 0;

const int EscalatorMC::COUPONDATE = Escalator::COUPONDATE;
const int EscalatorMC::AVGDATE = Escalator::AVGDATE;
const double EscalatorMC::TINY_SPREAD = 0.001;

DRLIB_END_NAMESPACE
