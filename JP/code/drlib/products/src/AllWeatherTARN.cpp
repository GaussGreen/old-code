//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AllWeatherTARN.cpp
//
//   Description : As TargetRedemptionNote but with the following additional features:

//                 - a mandatory lower barrier which, if breached by all assets, will result
//                   in a one-off x% coupon and early redemption.
//
//                 - an optional upper barrier which, if breached by all assets, will result
//                   in a special perf-dependent coupon but no early redemption.
//                   The same aggregation method is used as for the no-breach coupon but a
//                   different generalised performance may be specified.
//                   The upper barrier is persistent - it applies on subsequent coupon dates 
//                   even if already breached.
//                   
//
//                 The all-weather coupons (both up and down) are subject to the targetLevel
//
//   Date        : Mar 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/TargetRedemptionNote.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/

/** AllWeatherTARN product - as TARN but with upper and lower barriers.
                             If all asset performances breach these barriers
                             there is a special coupon and possible early redemption */
class AllWeatherTARN: public TargetRedemptionNote,
                      virtual public LegalTerms::Shift,
                      virtual public BarrierBreach::IEventHandler,
                      virtual public IMCIntoProduct{
private:
    /// fields ////////

    // couponPerf used in interface only. Assigned to 'performance' in validatePop2Object()
    IDoubleArrayModifierMakerSP  couponPerf;   

    // ALL-WEATHER fields
    // All assets must breach the up/downBarrier for the special all-weather event to be triggered
    DoubleArray                  AWBarriers;        // [0=down (mand), 1=up (opt)] all-weather barriers
    DoubleArray                  ecoAWBarriers;     // [0=down (opt), 1=up (opt)] economic all-weather barriers

    double                       downBarrier;       // transient
    double                       ecoDownBarrier;    // transient, defaulted to downBarrier
    bool                         isUpBarrier;       // transient
    double                       upBarrier;         // transient
    double                       ecoUpBarrier;      // transient, defaulted to upBarrier

    // Special all-weather coupon may be flat or dependent on asset-performance
    IDoubleArrayModifierMakerSP  upPerf;
    double                       downCoupon;
    bool                         payDoubleCoupon;   // whether to pay both the regular coupon and the downCoupon
                                                    // when a lower barrier breach triggers early redemption
public:
    static CClassConstSP const TYPE;
    friend class AllWeatherTARNMC;

    // validation
    void validatePop2Object(){
    static const string method = "AllWeatherTARN::validatePop2Object";
    try{

        // Assign couponPerf to performance
        performance = couponPerf;

        // parent
        TargetRedemptionNote::validatePop2Object();

        // Just validate all-weather fields
        // Barriers
        if ((AWBarriers.size() > 2) || AWBarriers.size() <= 0) {
            throw ModelException(method, "1 (down) or 2 (down and up) all-weather barriers must be provided");
        }

        if (ecoAWBarriers.size() > AWBarriers.size()) {
            throw ModelException(method, "There cannot be more economic barriers (" +
                                         Format::toString(ecoAWBarriers.size()) +
                                         ") than all-weather barriers (" +
                                         Format::toString(AWBarriers.size()) + ")");
        }

        // Set internal barrier fields
        downBarrier = AWBarriers[0];
        if (2==AWBarriers.size()) {
            isUpBarrier = true;
            upBarrier = AWBarriers[1];
        }

        // Default economic barriers if necessary
        for (int i=0; i<AWBarriers.size(); i++) {
            if (!(i<ecoAWBarriers.size())) {
                ecoAWBarriers.push_back(AWBarriers[i]);
            }
        }

        ecoDownBarrier = ecoAWBarriers[0];
        if (isUpBarrier) {
            ecoUpBarrier = ecoAWBarriers[1];
        }

        // Down-barrier
        if (downBarrier>1.0 || downBarrier<0.0) {
            throw ModelException(method, "Down barrier (" + Format::toString(downBarrier) +
                                         ") must be between 0.0 and 1.0");
        }

        // Economic down-barrier
        if (ecoDownBarrier>1.0 || ecoDownBarrier<0.0) {
            throw ModelException(method, "Economic down barrier (" + Format::toString(downBarrier) +
                                         ") must be between 0.0 and 1.0");
        }

        if (isUpBarrier) {
            if (upBarrier<1.0) {
                throw ModelException(method, "Up barrier (" + Format::toString(upBarrier) +
                                     ") if supplied cannot be less than 1.0");
            }

            if (ecoUpBarrier<1.0) {
                throw ModelException(method, "Economic up barrier (" + Format::toString(upBarrier) +
                                     ") if supplied cannot be less than 1.0");
            }

            if (!upPerf) {
                throw ModelException(method, "Up-coupon performance must be provided if an up barrier is provided");
            }

        }

        // If cliquet-style it is not clear whether the reference level for computation of the 
        // sunny coupon should be the initial spot price or the reset spot level. Disable for now.
        // As it stands the initial spot price would be used for both all-weather barrier testing and
        // computation of sunny coupon.
        if (isCliquetStyle) {
            throw ModelException(method, "Coupons may not be cliquet-style in the All-Weather TARN");
        }

    } catch (exception& e){
        throw ModelException(e, method);
        }
    }

    /** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift) {
        // Set the barriers for pricing equal to the economic barriers
        upBarrier = ecoUpBarrier;
        downBarrier = ecoDownBarrier;
        targetLevel = ecoTargetLevel;

        return false;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
    implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

    // Implementation of the BarrierBreach::IEventHandler interface
    virtual void getEvents(const BarrierBreach*     breach,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const;

private:
    AllWeatherTARN(): TargetRedemptionNote(TYPE),
                      downBarrier(2.0), ecoDownBarrier(2.0), isUpBarrier(false),
                      upBarrier(-1.0), ecoUpBarrier(-1.0), upPerf(0) {} // for reflection

    AllWeatherTARN(const AllWeatherTARN& rhs); // not implemented
    AllWeatherTARN& operator=(const AllWeatherTARN& rhs); // not implemented

    static IObject* defaultAllWeatherTARN(){
        return new AllWeatherTARN();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AllWeatherTARN, clazz);
        SUPERCLASS(TargetRedemptionNote);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(BarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultAllWeatherTARN);

        // All-weather Fields
        // couponPerf is used in interface only, assigned to 'performance' in validatePop2Object
        FIELD(couponPerf, "performance to apply to the aggregated asset basket");   
        FIELD(AWBarriers, "Lower and upper all-weather barriers");
        FIELD(ecoAWBarriers, "Economic lower and upper all-weather barriers");

        FIELD(downBarrier, "Lower barrier");
        FIELD_MAKE_TRANSIENT(downBarrier);
        FIELD(ecoDownBarrier, "Lower economic barrier");
        FIELD_MAKE_TRANSIENT(ecoDownBarrier);
        FIELD(isUpBarrier, "Whether there is an upper barrier");
        FIELD_MAKE_TRANSIENT(isUpBarrier);
        FIELD(upBarrier, "Upper barrier");
        FIELD_MAKE_TRANSIENT(upBarrier);
        FIELD(ecoUpBarrier, "Upper economic barrier");
        FIELD_MAKE_TRANSIENT(ecoUpBarrier);

        FIELD(upPerf, "Coupon payable if all assets breach upper barrier");
        FIELD_MAKE_OPTIONAL(upPerf);
        FIELD(downCoupon, "Coupon payable if all assets breach lower barrier");
        FIELD(payDoubleCoupon, "Whether the regular coupon is also paid upon lower barrier breach");

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super TargetRedemptionNote */
class AllWeatherTARNMC: public TargetRedemptionNoteMC {
private:
    // ALL-WEATHER fields
    // All assets must breach the up/downBarrier for the special all-weather event to be triggered
    double                 downBarrier;
    double                 ecoDownBarrier;
    bool                   isUpBarrier;
    double                 upBarrier;
    double                 ecoUpBarrier;

    // On upper breach by all assets how to aggregate the performances
    IAggregateSP           upAggregate;

    // Simple performance to apply to aggregate to get up-coupon
    IDoubleArrayModifierSP upPerf;

    double                 downCoupon;
    SimpleDoubleArray      upCoupon;

    SimpleDoubleArray      assetSimple;             // simple asset performances from start to coupon date: Si(T)/Si(0)

    bool                   payDoubleCoupon;         // Whether the regular coupon is also paid upon lower barrier breach

    // Needed for BARRIER_LEVEL output request
    DateTimeArray          barrierDate;            // Barrier dates
    DoubleArray            startLevel;             // Initial ref levels for each asset

    // BarrierBreach event details
    struct BarrierBreachDetails {
        int iCpnIdx;
        bool isSunny;
        DoubleArray levels;
        DoubleArray refLevels;
    };
    vector<struct BarrierBreachDetails> cachedBreaches;

    // Modify coupons according to all-weather conditions
    // Returns the coupon idx corresponding to the redemption date
    // (target is not yet considered)
    virtual int couponsOverride(const IPathGenerator*  pathGen,
                                DoubleArray& newCoupons,
                                bool& earlyRedeemed ) {

        int    endIdx   = pathGen->end(0);
        int redeemCpn = inst->couponDates.size() - 1;
        double myDownBarrier = downBarrier;
        double myUpBarrier = upBarrier;
        DoubleArray levels(nbAssets);

        if (pathGen->doingPast()) {
            // Override with economic barriers
            myDownBarrier = ecoDownBarrier;
            myUpBarrier = ecoUpBarrier;
        }

        for(int iAsset = 0; iAsset<nbAssets; iAsset++) {
            // startLevel is later needed for barrier reporting
            startLevel[iAsset] = pathGen->refLevel(iAsset, 0);
        }

        // Clear any cached barrier breach information.
        cachedBreaches.clear();

        bool downBreach, upBreach;
        int iCoupon = 0;  
        for (int iStep=0; iStep<endIdx; iStep++) {
      
            if (couponMap[iStep]==0) {
                // A coupon date

                if (earlyRedeemed) {
                    // All assets have already breached the lower barrier
                    newCoupons[iCoupon] = 0.0;
                }
                else {

                    downBreach = true;
                    upBreach = true;
                    for(int iAsset = 0; iAsset<nbAssets; iAsset++) {
                        levels[iAsset] = pathGen->Path(iAsset, 0)[iStep];

                        assetSimple[iAsset] = levels[iAsset] / startLevel[iAsset];

                        if (assetSimple[iAsset] > myDownBarrier) {
                            downBreach = false;
                        }

                        if (!isUpBarrier || (assetSimple[iAsset] < myUpBarrier)) {
                            upBreach = false;
                        }

                        if (!downBreach && !upBreach) {
                            // There is no all-weather barrier breach - exit loop early
                            break;
                        }
                    }

                    // Apply all-weather conditions:
                    // 1. If all assets have breached the upper barrier, pay a performance-dependent coupon
                    // 2. If all assets have breached the lower barrier, pay a one off coupon and early redeem the note
                    // 3. Otherwise apply the standard TARN no-breach coupon

                    if (upBreach) {
                        upCoupon[0] = upAggregate->aggregate();
                        upPerf->apply();
                        newCoupons[iCoupon] = upCoupon[0];
                    }
                    else if (downBreach) {
                        // Note: may 'early' redeem on last coupon date
                        earlyRedeemed = true;
                        redeemCpn = iCoupon;
                        if (payDoubleCoupon) {
                            newCoupons[iCoupon] = downCoupon + coupons[iCoupon];
                        } else {
                            newCoupons[iCoupon] = downCoupon;
                        }
                    } else {
                        newCoupons[iCoupon] = coupons[iCoupon];
                    }

                    if (pathGen->doingPast() && (upBreach || downBreach)) {
                        // Only record breach details on past, not simulation runs
                        struct BarrierBreachDetails details;

                        details.iCpnIdx = iCoupon;
                        details.isSunny = upBreach;
                        details.levels = levels;
                        details.refLevels = startLevel;

                        cachedBreaches.push_back(details);
                    }
                }
                iCoupon++;
            }
        }
        return redeemCpn;
    }

    // Store barrier dates (if applicable)
    virtual void barrierDates(const int earlyRedeemCpn) {
        static const string method("AllWeatherTARNMC::barrierDates");
        try {
            // need to report barrier levels over a date range
            DateTime toDate = BarrierLevel::barrierWindow(Today);

            // All weather breaches can only occur on coupon dates
            // earlyRedeemCpn is < inst->couponDates.size()
            for (int i=0; i<=earlyRedeemCpn && (inst->couponDates[i]<=toDate); i++) {
                if (inst->couponDates[i] >= Today) {
                    barrierDate.push_back(inst->couponDates[i]);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // Record barriers (if applicable)
    virtual void recordBarriers(Control* control,
                        Results* results) const {
        static const string method("AllWeatherTARNMC::recordBarriers");
        try {
            OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
            if (request && !request->getHasFinished())
            {
                BarrierLevelArray levels(0);
                for (int iAsset=0; iAsset<nbAssets; iAsset++) {

                    for (int iDate=0; iDate<barrierDate.size(); iDate++) {
                        // Down barrier is mandatory
                        BarrierLevel lowerBarrier(false,                    // isUp
                                                  barrierDate[iDate],       // date
                                                  ecoDownBarrier * startLevel[iAsset], // level
                                                  false);                   // isContinuous
                        levels.push_back(lowerBarrier);

                        // Up barrier is optional
                        if (isUpBarrier) {
                            BarrierLevel upperBarrier(true,                 // isUp
                                                      barrierDate[iDate],   // date
                                                      ecoUpBarrier * startLevel[iAsset], // level
                                                      false);               // isContinuous
                            levels.push_back(upperBarrier);
                        }
                    }

                    if (!levels.empty()) {

                        OutputRequestUtil::recordBarrierLevels(control,
                                                               results,
                                                               getMultiFactors()->assetGetTrueName(iAsset),
                                                               &levels);
                    }
                    levels.clear();
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

public:

    /** equivalent to InstIntoMCProduct */
    AllWeatherTARNMC(const AllWeatherTARN* inst,
                     const SimSeriesSP& simSeries):
            TargetRedemptionNoteMC(inst, simSeries),
            downBarrier(inst->downBarrier),
            ecoDownBarrier(inst->ecoDownBarrier),
            isUpBarrier(inst->isUpBarrier),
            upBarrier(inst->upBarrier),
            ecoUpBarrier(inst->ecoUpBarrier),
            downCoupon(inst->downCoupon),
            upCoupon(1, 0.0),       // performance will be applied for each up coupon in turn so just length 1
            assetSimple(nbAssets, 0.0),
            payDoubleCoupon(inst->payDoubleCoupon),
            barrierDate(0),
            startLevel(nbAssets)
    {

        // Use normal no-breach coupon aggregation method against assetSimple performances in
        // computing up-coupon
        // Attach upPerf to up coupon performance
        upAggregate = IAggregateSP(inst->assetBasket->getAggregate(&assetSimple));
        upPerf = IDoubleArrayModifierSP(inst->upPerf->getModifier(&upCoupon));
    }

    // Override TargetRedemptionNote::retrieveEvents.
    // If we have knocked out due to a rainy day condition, then there cannot be a target redemption
    // after it, even if the rainy day coupon causes the target to be hit (and is capped by it).
    virtual void retrieveEvents(EventResults* events) const {
        bool hasRainy = false;
        int iRainyIdx;
        for (unsigned int i=0; i<cachedBreaches.size(); i++) {
            const struct BarrierBreachDetails& breach = cachedBreaches[i];
            if( !breach.isSunny ) {
                hasRainy = true;
                iRainyIdx = breach.iCpnIdx;
            }
        }

        if( hasTargetRedeemed &&
            !(hasRainy && iRainyIdx<=cachedRedemptionDetails.couponIdx) )
        {
            // Emit TargetRedemption event
            const DateTime& redemptionDate = inst->couponDates[cachedRedemptionDetails.couponIdx];
            
            string floatingLegType = TargetRedemption::NOT_APPLICABLE;
            // Handle floating leg details
            if( inst->hasFloater ) {
                if( inst->isKO ) {
                    floatingLegType = TargetRedemption::KNOCK_OUT;
                } else {
                    floatingLegType = TargetRedemption::KNOCK_IN;
                }
            }

            events->addEvent(new TargetRedemption(redemptionDate,
                                                  cachedRedemptionDetails.finalCoupon, 
                                                  cachedRedemptionDetails.totalCoupon,
                                                  cachedRedemptionDetails.target, 
                                                  cachedRedemptionDetails.bonus,
                                                  cachedRedemptionDetails.redemption,
                                                  floatingLegType));
        }
    }

    void retrieveBarrierBreachEvents(EventResults* events) const {
        for (unsigned int i=0; i<cachedBreaches.size(); i++) {
            const struct BarrierBreachDetails& breach = cachedBreaches[i];
            if( !(hasTargetRedeemed && breach.iCpnIdx > cachedRedemptionDetails.couponIdx ) ) {
                // Ignore all sunny-side/rainy-side coupons that occur after early redemption
                const DateTime& eDate = inst->couponDates[breach.iCpnIdx];
                
                string barrDesc = breach.isSunny ? "All Weather TARN - Sunny side barrier" :
                                                   "All Weather TARN - Rainy side barrier";
                string barrType = breach.isSunny ? BarrierBreach::NOT_APPLICABLE : // Special coupon
                                                   BarrierBreach::KNOCK_OUT; // Early redeems
                string monType = BarrierBreach::EUROPEAN;
                int hitsLeft = 0; // Not used

                StringArraySP assetNames(new StringArray(nbAssets));
                DoubleArraySP assetLevels(new DoubleArray(nbAssets));
                DoubleArraySP barrLevels(new DoubleArray(nbAssets));
                double barr = breach.isSunny ? ecoUpBarrier :
                                               ecoDownBarrier;
                for( int iAsset=0; iAsset<nbAssets; iAsset++ ) {
                    (*assetNames)[iAsset] = getMultiFactors()->assetGetTrueName(iAsset);
                    (*assetLevels)[iAsset] = breach.levels[iAsset];
                    (*barrLevels)[iAsset] = breach.refLevels[iAsset] * barr;
                }

                events->addEvent( new BarrierBreach( eDate, barrDesc, monType, barrType,
                    breach.isSunny, hitsLeft, assetNames, assetLevels, barrLevels ) );
            }
        }
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* AllWeatherTARN::createProduct(const MonteCarlo* model) const {

    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
    one */
    simSeries->addDates(averageOutDates);
    return new AllWeatherTARNMC(this, simSeries);
}

    // Implementation of the BarrierBreach::IEventHandler interface
void AllWeatherTARN::getEvents(const BarrierBreach*     breach,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const
{
    static const string method = "AllWeatherTARN::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP past = prod->runPast(mc);
            
            // Little bit of nastiness. Down cast the IMCProduct to the AllWeatherTARNMC,
            // so that we can invoke the barrier breach retrieval method
            AllWeatherTARNMC *awt = dynamic_cast<AllWeatherTARNMC*>(prod.get());
            if (awt) {
                awt->retrieveBarrierBreachEvents(events);
            } else {
                throw ModelException( method,
                    "Internal error - expected AllWeatherTARNMC to be product for BarrierBreach events");
            }
        } else {
            throw ModelException(method, 
                "Internal error - expected Monte Carlo model for AllWeatherTARN pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
};

CClassConstSP const AllWeatherTARN::TYPE = CClass::registerClassLoadMethod(
    "AllWeatherTARN", typeid(AllWeatherTARN), AllWeatherTARN::load);

// * for class loading (avoid having header file) */
bool AllWeatherTARNLoad() {
    return (AllWeatherTARN::TYPE != 0);
}

DRLIB_END_NAMESPACE
