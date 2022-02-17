//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CalendarRainbow.cpp
//
//   Description : Port of EDG CalendarRainbow
//
//   Date        : Nov 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/Format.hpp"
#include "edginc/IAggregate.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

/** CalendarRainbow product - option on doubly rainbowed perfs - across time and asset */
class CalendarRainbow: public GenericNFBase, 
                       virtual public IMCIntoProduct{
private: // was protected - any reason why?
    /// fields ////////
    IAggregateMakerSP            timeBasket;
    IDoubleArrayModifierMakerSP  timeBasketComponents;
    IAggregateMakerSP            assetBasket;
    IDoubleArrayModifierMakerSP  assetBasketComponents;
    IDoubleArrayModifierMakerSP  overallOption;
    bool                         avgFromStart;
    DateTimeArray                averageOutDates;
    DateTimeArray                couponDates;
    bool                         isCliquetStyle;
    bool                         isRefPrevAvgOut;

// perhaps a flag for which of time or asset dimension gets rainbowed first?

public:
    static CClassConstSP const TYPE;
    friend class CalendarRainbowMC;

    // validation
    void validatePop2Object(){
        static const string routine = "CalendarRainbow::validatePop2Object";
        GenericNFBase::validatePop2Object();

        if (averageOutDates.empty()) {
            throw ModelException(routine, "No averageOutDates given!");
        }
        if (couponDates.empty()) {
            throw ModelException(routine, "No couponDates given!");
        }

        // Check average dates and coupon dates "cooperate"
        // Coupon dates must be a subset of average dates
        if (!DateTime::isSubset(averageOutDates, couponDates)) {
            throw ModelException(routine, "Coupon dates should be a subset of averageOutDates");
        }
        // Require some averaging before first coupon date
        const DateTime& firstAvg = averageOutDates[0];
        const DateTime& firstCpn = couponDates[0];
        if (firstCpn < firstAvg) {
            throw ModelException(routine, "Cannot have coupon date " + firstCpn.toString() + 
                                 " before first average date " + firstAvg.toString());
        }
        // Makes no sense to have averaging after final coupon date
        const DateTime& lastAvg = averageOutDates[averageOutDates.size()-1];
        const DateTime& lastCpn = couponDates[couponDates.size()-1];
        if (lastAvg > lastCpn) {
            throw ModelException(routine, "Cannot average on " + lastAvg.toString() + 
                                 " since after final coupon date " + lastCpn.toString());
        }

        // Don't see how this can make sense so forbid for now
        if (isCliquetStyle) {
            if (avgFromStart) {
                throw ModelException(routine, "Cannot avgFromStart if isCliquetStyle");
            }
        } else {
            if (isRefPrevAvgOut) {
                // should/can we check anything here?
                // The paramater is meaningless if not cliquet...
            }
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
    CalendarRainbow(): GenericNFBase(TYPE), isRefPrevAvgOut(false) {} // for reflection
    CalendarRainbow(const CalendarRainbow& rhs); // not implemented
    CalendarRainbow& operator=(const CalendarRainbow& rhs); // not implemented

    static IObject* defaultCalendarRainbow(){
        return new CalendarRainbow();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CalendarRainbow, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultCalendarRainbow);
        FIELD(timeBasket,             "timeBasket");
        FIELD(timeBasketComponents,   "timeBasketComponents");
        FIELD(assetBasket,            "assetBasket");
        FIELD(assetBasketComponents,  "assetBasketComponents");
        FIELD(overallOption,          "overallOption");
        FIELD(avgFromStart,           "avgFromStart");
        FIELD(averageOutDates,        "averageOutDates");
        FIELD(couponDates,           "couponDates");
        FIELD(isCliquetStyle,         "isCliquetStyle");
        FIELD(isRefPrevAvgOut,         "isRefPrevAvgOut");
        FIELD_MAKE_OPTIONAL(isRefPrevAvgOut);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super CalendarRainbow */
class CalendarRainbowMC: public IMCProduct,
                         virtual public IMCProductLN {
private:
    const CalendarRainbow*   inst;      // reference to original instrument
    int                      nbAssets;  // nicer
    DoubleArray              sum;       // [nbAssets], saves alloc 
    DoubleArray              refLevel;  // [nbAssets], saves alloc 
    SimpleDoubleArray        assetComps;
    SimpleDoubleArray        coupons;

    // historical values
    DoubleArray              sumSoFar;      // [nbAssets]
    DoubleArray              refLevelSoFar; // [nbAssets]
    SimpleDoubleArray        couponsSoFar;  // [nbCoupons]
    int                      iCouponSoFar;

    // Operational aggregation and performance calcs
    IAggregateSP             timeBasket;
    IDoubleArrayModifierSP   timeBasketComponents;
    IAggregateSP             assetBasket;
    IDoubleArrayModifierSP   assetBasketComponents;
    IDoubleArrayModifierSP   overallOption;
    TrivialDoubleArray       calRbw;  // just a double ... may be more natural way to phrase this?

    IntArray                 nbAvgOutPerCouponDate; // [nbCoupons] 
    IntArray                 couponMap;    // [nbAvgDates] - convenient way to track coupon dates

public:
    
    /** equivalent to InstIntoMCProduct */
    CalendarRainbowMC(const CalendarRainbow*         inst,
                      const SimSeriesSP&             simSeries):
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
        assetComps(nbAssets, 0.0),
        coupons(inst->couponDates.size(), 0.0),
        sumSoFar(nbAssets, 0.0),
        refLevelSoFar(nbAssets, 0.0),
        couponsSoFar(inst->couponDates.size(), 0.0),
        iCouponSoFar(0),
        calRbw(0.0),
        nbAvgOutPerCouponDate(inst->couponDates.size(), 0) {
        
        // we check the first avg date is not after first cpn date
        int iCoupon = 0; 
        for(int iStep = 0; iStep < inst->averageOutDates.size(); iStep++) {    
            // number of AvgOut per period
            if(inst->averageOutDates[iStep] > inst->couponDates[iCoupon]) {
                iCoupon++;
            }
            nbAvgOutPerCouponDate[iCoupon]++;
        }
        // cumulative number of AvgOut Dates
        if(inst->avgFromStart) {
            for(iCoupon = 1; iCoupon < inst->couponDates.size(); iCoupon++) {
                nbAvgOutPerCouponDate[iCoupon] += nbAvgOutPerCouponDate[iCoupon-1];
            }
        } 
        bool isTrivial;
        couponMap = DateTime::createMapping(inst->averageOutDates,
                                            inst->couponDates,
                                            isTrivial);

        timeBasketComponents = IDoubleArrayModifierSP(inst->timeBasketComponents->getModifier(&coupons));
        timeBasket = IAggregateSP(inst->timeBasket->getAggregate(&coupons));
        assetBasketComponents = IDoubleArrayModifierSP(inst->assetBasketComponents->getModifier(&assetComps));
        assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
        overallOption = IDoubleArrayModifierSP(inst->overallOption->getModifier(&calRbw));
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

        int iCoupon = iCouponSoFar;
        sum = sumSoFar;
        refLevel = refLevelSoFar;
        coupons = couponsSoFar;
        for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
            /* Special treatment for average in for first cliq */
            if (!inst->isCliquetStyle || iCoupon == 0) {
                refLevel[iAsset] = pathGen->refLevel(iAsset, 0);
            } else {
                refLevel[iAsset] = refLevelSoFar[iAsset];
            }
        }

        // Form the asset basket at each coupon date first
        for (int iStep=beginIdx; iStep<endIdx; iStep++) {
            bool isCouponDate = (couponMap[iStep]==0);   // true iff a coupon date

            for(iAsset=0; iAsset<nbAssets; iAsset++) {
                sum[iAsset] += pathGen->Path(iAsset, 0)[iStep];

                if (isCouponDate) { 
                    double avg = sum[iAsset] / nbAvgOutPerCouponDate[iCoupon];
                    assetComps[iAsset] = avg / refLevel[iAsset];
                    if (inst->isCliquetStyle)
                    {
                        /* Reset ref level for next coupon - done per asset note*/
                        if (inst->isRefPrevAvgOut) {
                            refLevel[iAsset] = avg;
                        } else {
                            refLevel[iAsset] = pathGen->Path(iAsset, 0)[iStep];
                        }
                    }
                    if (!inst->avgFromStart) {
                        sum[iAsset] = 0.0;
                    }
                }
            }

            if (isCouponDate) {
                // apply any perf modifiers
                assetBasketComponents->apply();
                // then form into an appropriate basket
                coupons[iCoupon] = assetBasket->aggregate();
                // get ready for next coupon
                iCoupon++;
            }
        }

        // preserve values for past - before we modify/sort the coupons.
        if (pathGen->doingPast()){ 
            sumSoFar = sum;
            refLevelSoFar = refLevel;
            couponsSoFar = coupons;
            iCouponSoFar = iCoupon;
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
        static const string routine = "CalendarRainbowMC::getVolInterp";
        CVolRequestLNArray reqarr(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        const DateTime&    today = getToday();
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();
        bool               fwdStarting = startDate.isGreater(today);
        double             interpLevel = inst->assetBasketComponents->getInterpLevel(iAsset);

        if (inst->isCliquetStyle) {
            // get hold of the future strike dates
            int numLiveCliqs = inst->couponDates.size() - iCouponSoFar;
            if (numLiveCliqs<=0) {
                throw ModelException(routine, "No future coupons!?");
            }
            DateTimeArray liveCliqStartDates(numLiveCliqs);
            for (int iCliquet = 0; iCliquet < numLiveCliqs; iCliquet++){
                int iCoupon = iCouponSoFar+iCliquet-1;
                liveCliqStartDates[iCliquet] = iCoupon<0?startDate:inst->couponDates[iCoupon];
            }

            // same strike levels per cliquet (but may need to adjust first one)
            DoubleArray  strikes(numLiveCliqs, interpLevel);
            if (!fwdStarting){
                // need to set first level to absolute strike - adjusted
                // additionally for any average out samples for this cliquet
                int iStep;
                const DateTime& thisCouponStart = iCouponSoFar==0?
                    startDate:inst->couponDates[iCouponSoFar-1];
                // find first avg date of this coupon
                for(iStep = 0; iStep < inst->averageOutDates.size() && 
                        inst->averageOutDates[iStep] <= thisCouponStart; iStep++) {
                    ; // empty
                }
                // then walk through counting past avg dates in this cliq
                int numRemaining = nbAvgOutPerCouponDate[iCouponSoFar];
                for(; iStep < inst->averageOutDates.size() && 
                        inst->averageOutDates[iStep] <= today; iStep++) {
                    numRemaining--;
                }
                if (numRemaining<=0) {
                    // something wrong!
                    throw ModelException(routine, "INTERNAL ERROR : numRemaining is " + Format::toString(numRemaining));
                }
                // Can't set up refLevel earlier, 'cos need PathGen. First cliq has standard ref level
                double refLevel =  iCouponSoFar==0 ? pathGen->refLevel(iAsset, 0) : refLevelSoFar[iAsset];
                strikes[0] = (nbAvgOutPerCouponDate[iCouponSoFar] * refLevel * interpLevel
                              - sumSoFar[iAsset])/ numRemaining;
            }
            reqarr[0] =  CVolRequestLNSP(new CliquetVolRequest(fwdStarting, 
                                                               liveCliqStartDates, 
                                                               lastSimDate,
                                                               strikes));
        } else {
            // per asset 
            if (!fwdStarting){
                // some samples have fixed already (this includes averaging in)
                if (inst->avgFromStart) {
                    int numDates = inst->averageOutDates.size();
                    int numRemaining = 
                        today.numFutureDates(inst->averageOutDates);
                    
                    interpLevel = (numDates * interpLevel * pathGen->refLevel(iAsset, 0)
                                   - sumSoFar[iAsset])/ numRemaining;
                } else {
                    interpLevel *= pathGen->refLevel(iAsset, 0);
                }
            }
            reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                     startDate,
                                                                     lastSimDate,
                                                                     fwdStarting));
        }
        return reqarr;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* CalendarRainbow::createProduct(const MonteCarlo* model) const {

    // XXX Cliquet style resetting is not supported with implied - how enforce that?

    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(averageOutDates);
    return new CalendarRainbowMC(this, simSeries);
}

CClassConstSP const CalendarRainbow::TYPE = CClass::registerClassLoadMethod(
    "CalendarRainbow", typeid(CalendarRainbow), CalendarRainbow::load);

// * for class loading (avoid having header file) */
bool CalendarRainbowLoad() {
    return (CalendarRainbow::TYPE != 0);
}

DRLIB_END_NAMESPACE





