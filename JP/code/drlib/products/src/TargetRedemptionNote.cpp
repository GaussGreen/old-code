//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TargetRedemptionNote.cpp
//
//   Description : 
//
//   Date        : Nov 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/TargetRedemptionNote.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
TargetRedemptionNote::TargetRedemptionNote(CClassConstSP clazz): GenericNFBase(clazz), 
                            assetPerf(0), performance(0), isCliquetStyle(false), isRefPrevAvgOut(false),  
                            isMakeWhole(false), ecoTargetLevel(0.0), hasFloater(false), floater(0), 
                            floorWithPreviousCoup(false), receiveOvershoot(false), isKO(true), 
                            optionAtMat(false), onePayment(false) 
{} 

// validation
void TargetRedemptionNote::validatePop2Object(){
    static const string routine = "TargetRedemptionNote::validatePop2Object";
    GenericNFBase::validatePop2Object();   

    if (Maths::isZero(targetLevel)) {
        // Could also interpret a zero target as no target, but this is not transparent
        // Better to validate against and for user to enter an unachieveable target if there is no target
        throw ModelException(routine, "No target level supplied!");
    }

    if (Maths::isZero(ecoTargetLevel)) {
        // No economic target level has been specified
        ecoTargetLevel = targetLevel;
    }

    if (averageOutDates.empty()) {
        throw ModelException(routine, "No average out dates given!");
    }
    if (couponDates.empty()) {
        throw ModelException(routine, "No coupon dates given!");
    }
        
    // Check average dates and coupon dates "cooperate"
    // Coupon dates must be a subset of average dates
    if (!DateTime::isSubset(averageOutDates, couponDates)) {
        throw ModelException(routine, "Coupon dates should be a subset of averageOutDates");
    }

    // Require some averaging before first coupon date
    const DateTime& firstAvg = averageOutDates[0];
    const DateTime& firstCpn = couponDates[0];

    // Does this make sense?  couponDates is already a subset of averageOutDates.
    // Did you mean == here?
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
    
    // the following validations should be removed as they are tested (cliquet, avg)
    // or implemented (libor) properly
    
    if (bonusCoupons.size() != 0 && bonusCoupons.size() != couponDates.size()) {
        throw ModelException(routine, "Number of bonus coupons = " + Format::toString(bonusCoupons.size())
            + " does not match number of coupon dates = " + Format::toString(couponDates.size()) + ".");
        
    }
    
    // validate that hasFloater has floater
    if (hasFloater && !floater.get()){
        throw ModelException(routine,
            "hasFloater = true, but floater is not found");
    }

    // validate that the settlement is not physical
    if (instSettle->isPhysical()) {
        throw ModelException(routine, 
                             "physical settlement not supported");
    }

    if (!assetPerf) {
        // No asset-level performance is specified so default to a zero-strike forward 
        // (equivalent to not applying a performance)
        assetPerf = IDoubleArrayModifierMakerSP(new PerfTypeSimpleMaker(
                            "F",        // perfType 
                            0.0,        // strike
                            1.0));      // participation
        assetPerf->validatePop2Object();
    }
}

/** Get the asset and discount market data */
void TargetRedemptionNote::GetMarket(const IModel*          model, 
                                     const CMarketDataSP    market) {
    // parent
    GenericNFBase::GetMarket(model, market);
    // relevant elements of self
    if (hasFloater) {
        floater->getMarket(model, market.get());
    }
}

bool TargetRedemptionNote::sensShift(Theta* shift) {
    static const string method = "TargetRedemptionNote::sensShift";
    try  {
        const DateTime& newDate = shift->rollDate(valueDate);
        // and fixings
        if (hasFloater) {
            floater->setFixingforThetaShift(valueDate,
                discount.get(),
                newDate);
        }
        
        // roll the parent (updates value date etc)
        GenericNFBase::sensShift(shift);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }    
    return true;
    // When 'couponCurve' is live can consider letting it take over.
}

/** Satisfy LegalTerms::Shift interface */
bool TargetRedemptionNote::sensShift(LegalTerms* shift) {
    // Set the barriers for pricing equal to the economic barriers
    targetLevel = ecoTargetLevel;

    return false;
}

// TargetRedemption::IEventHandler interface
void TargetRedemptionNote::getEvents(const TargetRedemption* target, IModel* model, 
                                     const DateTime& eventDate, EventResults* events) const {
    static const string method = "TargetRedemptionNote::getEvents";

    try {
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            auto_ptr<IMCProduct> prod(createProduct(mc));
            MCPathGeneratorSP past = prod->runPast(mc);
            prod->retrieveEvents(events);
        } else {
            throw ModelException(method, 
                    "Internal error - expected Monte Carlo model for TargetRedemptionNote pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
 }

// returns the array of all sampling dates the instrument will ever need
// excluding (possibly) the ref level dates (handled in GenericNFBase)
const DateTimeArray TargetRedemptionNote::samplingDates() const {
    return averageOutDates;
}

class TargetRedemptionNoteHelper {
public:
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TargetRedemptionNote, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LegalTerms::Shift);
        IMPLEMENTS(TargetRedemption::IEventHandler);
        EMPTY_SHELL_METHOD(defaultTargetRedemptionNote);
        FIELD(assetPerf, "performance to apply to raw asset performances before aggregation");
        FIELD_MAKE_OPTIONAL(assetPerf);
        FIELD(performance, "performance to apply to the aggregated asset basket");
        FIELD_MAKE_OPTIONAL(performance); // Optional for AWT compatibility
        FIELD(assetBasket,            "assetBasket");
        FIELD(avgFromStart,           "avgFromStart");
        FIELD(averageOutDates,        "averageOutDates");
        FIELD(couponDates,           "couponDates");
        FIELD(frontCoupons,    "guaranteed coupons recieved before the coupons based on the performance.");
        FIELD_MAKE_OPTIONAL(frontCoupons); // to be removed
        FIELD(targetLevel,    "level which serves as a barrier for the cum coupons paid");
        FIELD(ecoTargetLevel,    "Economic target level");
        FIELD_MAKE_OPTIONAL(ecoTargetLevel);
        FIELD(isMakeWhole,     "True if no early redemption means that targetLevel - cum coupon payments are made at maturity.");
        FIELD(isCliquetStyle,         "Whether the note is cliquet style");
        FIELD_MAKE_OPTIONAL(isCliquetStyle);
        FIELD(isRefPrevAvgOut,         "true if the reference level is reset to previous average level");
        FIELD_MAKE_OPTIONAL(isRefPrevAvgOut);
        FIELD(hasFloater, "true if contract is swap style, false if note style");
        FIELD(floater, "libor leg");
        FIELD_MAKE_OPTIONAL(floater);
        FIELD(floorWithPreviousCoup, "true: all coupons after the first one are floored with previous coupon.");
        FIELD_MAKE_OPTIONAL(floorWithPreviousCoup);
        FIELD(receiveOvershoot,         "true: do not cap cum coupons at the target level when breached.");
        FIELD_MAKE_OPTIONAL(receiveOvershoot);
        FIELD(isKO,         "true/false: libor knocks out/in on early redemption");
        FIELD_MAKE_OPTIONAL(isKO);
        FIELD(optionAtMat,         "true: if note does not early redeem, investor receives additional option at maturity.");
        FIELD_MAKE_OPTIONAL(optionAtMat);
        FIELD(matPerformance,   "performance at maturity to apply to the aggregated asset basket in optionAtMat case.");
        FIELD_MAKE_OPTIONAL(matPerformance);
        FIELD(bonusCoupons, "bonus coupons schedule. The investor receives one bonus coupon when the target is breached.");
        FIELD_MAKE_OPTIONAL(bonusCoupons);
        FIELD(onePayment, "if true, all the coupons are paid on the redemption date");
        FIELD_MAKE_OPTIONAL(onePayment);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultTargetRedemptionNote(){
        return new TargetRedemptionNote();
    }
};

/** equivalent to InstIntoMCProduct */
TargetRedemptionNoteMC::TargetRedemptionNoteMC(const TargetRedemptionNote*         inst,
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
            coupons(inst->couponDates.size(), 0.0),
            sum(nbAssets, 0.0),
            refLevel(nbAssets, 0.0),
            assetComps(nbAssets, 0.0),
            sumSoFar(nbAssets, 0.0),
            refLevelSoFar(nbAssets, 0.0),
            couponsSoFar(inst->couponDates.size(), 0.0),
            iCouponSoFar(0),
            nbAvgOutPerCouponDate(inst->couponDates.size(), 0),
            iFirstUnKnownFloater(0),
            hasTargetRedeemed(false)
{
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
    couponMap = DateTime::createMapping(simSeries->getAllDates(),
        inst->couponDates,
        isTrivial);
    
    // set up couponsFV: forward value factors for rainbow coupons
    couponsFV = DoubleArray(inst->couponDates.size());
    DateTime matSettlementDate = settlement->settles(inst->couponDates.back(),0);
    for(iCoupon = 0; iCoupon < inst->couponDates.size(); iCoupon++) {
        
        DateTime couponPayDate = settlement->settles(inst->couponDates[iCoupon],0); 
        // While we're at it, may as well store pay dates if each coupon is paid 
        // if one payment only the last coupon date is stored
        if ((!inst->onePayment) || ((iCoupon+1) == inst->couponDates.size())) {
            instPaymentDates.push_back(couponPayDate);
        }

        if (couponPayDate <= inst->valueDate)
        {
            couponsFV[iCoupon] = 0.0;
        }
        else
        {
            couponsFV[iCoupon] = 1.0 / discount->pv(couponPayDate, matSettlementDate);
        }
    }
    
    // compute liborParticipation: amount of libor holder will receive if early redemption occurs at coupon i.
    if (inst->hasFloater)
    {
        iFirstUnKnownFloater = 0; // init to 0 in case inst is not a KO

        liborParticipation = DoubleArray(inst->couponDates.size());
        liborParticipation = *computeLiborParticipation(matSettlementDate).get();
        
        // add libor payment dates to instPaymentDates
        PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
        DateTimeArraySP payDates = fltpay->paymentDates();
        instPaymentDates.insert(instPaymentDates.end(), payDates->begin(),payDates->end());

        // prepare the KNONW_CASHFLOW.  All cashflow in floating leg, 
        // which are scheduled before the first coupon date,
        // are known to be paid, so add to KNOWN_CASHFLOW.
        // equity linked coupon are not calculated, because it should be calculated in cacheKnownCFS.
        if (inst->isKO)
        {
            PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
    
            CashFlowArrayConstSP floatcf(fltpay->knownCashflows(inst->valueDate,
                                                                0,
                                                                false,
                                                                inst->valueDate,
                                                                0.0));
    
            DateTime firstCpnDate = inst->couponDates[0];                
            int i = 0;
            while (i < floatcf->size() && (*floatcf)[i].date <= firstCpnDate){
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow((*floatcf)[i].date, (*floatcf)[i].amount * inst->notional));
                i ++;
            }
            iFirstUnKnownFloater = i;                
        }
    }
    else
    {
        liborParticipation = DoubleArray(0);
    }
    
    // now attach assetComps to asset-level performance and aggregation method, 
    // and coupons to coupon-level performance
    assetPerf = IDoubleArrayModifierSP(inst->assetPerf->getModifier(&assetComps));
    assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
    performance = IDoubleArrayModifierSP(inst->performance->getModifier(&coupons));
}

// control variate is done here
void TargetRedemptionNoteMC::recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const {
    
    if (!inst->isKO && inst->hasFloater)
    {
        double liborPrice = inst->notional*
                    inst->floater->getPV(inst->valueDate, inst->discount.get());
        
        results->storePrice(liborPrice + results->retrievePrice(),
            results->getCcyName());
        
    }
}

void TargetRedemptionNoteMC::payoff(const IPathGenerator*  pathGen,
                                    IMCPrices&                prices) {
    int    beginIdx = pathGen->begin(0); // 0 <- same for all assets
    int    endIdx   = pathGen->end(0);
    int    iAsset;

    // update with preserved values from doing past
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

            // Apply asset-level performances to assetComps[]
            assetPerf->apply();
        
            // then compute aggregate (e.g. worst performance)
            coupons[iCoupon] = assetBasket->aggregate();
            
            // get ready for next coupon
            iCoupon++;
        }
    }

    // preserve values for past - before we modify/sort the coupons.
    if (doingPast()){ 
        sumSoFar = sum;
        refLevelSoFar = refLevel;
        couponsSoFar = coupons;         
        iCouponSoFar = iCoupon;
    }

    // Value of final coupon. Used in option calculations.
    double aggAtMat = coupons[coupons.size()-1];   // store this before we apply the performance
    
    // apply perf modifier to the coupon values.
    performance->apply();

    if (doingPast()){ 
        cacheKnownCFs(aggAtMat, pathGen);

        // change paymentCashFlows if has breach
    }

    if (!doingPast() || !hasFuture()) {
        // Compute a payoff, but only when we have a "complete" situation : either 
        // doingPast() and all is past, or !doingPast().

        double value = computeTargetRedemptionPayoff(aggAtMat, pathGen);
    
        prices.add(inst->notional*(value)); 

    }
}

// Modify coupons according to any special (e.g. all-weather) conditions
// Returns the coupon idx corresponding to the redemption date (ignoring target)
int TargetRedemptionNoteMC::couponsOverride(const IPathGenerator*  pathGen,
                                            DoubleArray& newCoupons,
                                            bool& earlyRedeemed) {

    // No overriding on the vanilla TARN
    earlyRedeemed = false;
    int endIdx   = pathGen->end(0);
    int iCoupon=0;
    for (int iStep=0; iStep<endIdx; iStep++) {

        if (couponMap[iStep]==0) {
        // A coupon date

            newCoupons[iCoupon] = coupons[iCoupon];
            iCoupon ++;
        }
    }

    // Ignoring target there is no early redemption
    return inst->couponDates.size() - 1;
}
    

// Returns a double array of same size as the couponDates array.
// The i'th value represents the maturity settlement value of the unpaid funding payments
// participated in, given that early redemption occurs at coupon i.
// Note: if early redemption occurs at coupon i, then the last funding payment participated in (paid or unpaid)
// is the first one on or after the coupon date when early redemption occurs.
DoubleArraySP TargetRedemptionNoteMC::computeLiborParticipation(const DateTime& matSettlementDate)
{
    CashFlowArrayConstSP cf(inst->floater->getCashFlowArray(inst->valueDate, inst->discount.get()));
    DoubleArraySP libParticipation(new DoubleArray(inst->couponDates.size()));
    int iCoup = 0;
    for (iCoup = 0; iCoup < inst->couponDates.size(); iCoup++)
    {
        double cum = 0.0;
        int iLib = 0;
        bool done = false;
        while (!done && iLib < cf->size())
        {
            DateTime thisLiborDate = (*cf)[iLib].date;
            if (thisLiborDate > inst->valueDate)
            {   // only include libor payments that haven't been payed yet
                cum += (*cf)[iLib].amount / 
                    inst->discount->pv(thisLiborDate,matSettlementDate);
            }
            if (thisLiborDate.getDate() >= inst->couponDates[iCoup].getDate())
            {
                done = true;
            }
            
            iLib ++;
        }
        (*libParticipation)[iCoup] = cum * (inst->isKO? 1.0 : -1.0);
    }
    return libParticipation;
}

double TargetRedemptionNoteMC::computeTargetRedemptionPayoff(
        double aggregateAtMaturity,
        const IPathGenerator*  pathGen)
{
    double maturityValue = 0.0; // we need forward value at maturity as MC engine pv's.
    DoubleArray newCoupons(coupons.size());
    bool earlyRedeemed = false;
    struct RedemptionDetails details;
    
    // compute newCoupons and couponIdxAtBreach
    int earlyRedeemCpn = 
        processCoupons(aggregateAtMaturity, newCoupons, earlyRedeemed, details, pathGen);

    // now add to maturityValue with couponsFV.
    for (int i = 0; i <= earlyRedeemCpn; i++)
    {
        maturityValue += newCoupons[i] * couponsFV[i];      
    }
    
    // for floater case, need to add the forward valued libor cashflows, which is stored in liborParticipation.
    if (inst->hasFloater)
    {
        maturityValue += liborParticipation[earlyRedeemCpn];
    }
    
    return maturityValue;       
}

void TargetRedemptionNoteMC::cacheKnownCFs(
        double aggregateAtMaturity,
        const IPathGenerator*  pathGen)
{
    DoubleArray newCoupons(coupons.size());
    bool earlyRedeemed = false;
    struct RedemptionDetails details;
    
    // compute newCoupons and couponIdxAtBreach
    int earlyRedeemCpn = 
        processCoupons(aggregateAtMaturity, newCoupons, earlyRedeemed, details, pathGen);

    // if has early redeemed then we know exactly the payment dates
    // so the payment dates are calculated again
    if (earlyRedeemed) {
        instPaymentDates.clear();
    }

    // add the settlement date for each coupon date to the knownCFs
    int iCoup = 0; 
    int nextCoup = 0;   // used below for fees kcfs.
    
    // Add known cashflows up to and including last coupon before today
    // Add payment dates up to and including couponIdxAtBreach (which maybe 
    // a particular future coupon date if flooring on previous coupons for example)
    while(iCoup <= earlyRedeemCpn)
    {
        if (!inst->onePayment) {
            if (inst->couponDates[iCoup] <= inst->valueDate) {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow(settlement->settles(inst->couponDates[iCoup],0), 
                                                   inst->notional*newCoupons[iCoup]));
            }
            // instPaymentDates is empty if the instrument has early redeemed
            // so we add again all the coupon date
            if (earlyRedeemed) {
                instPaymentDates.push_back(settlement->settles(inst->couponDates[iCoup],0));
            }
        }
        else if (earlyRedeemed && iCoup==earlyRedeemCpn) {
            // if one Payment, we need to had the breach date and cashflow if before today
            if (inst->couponDates[iCoup] <= inst->valueDate) {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow(settlement->settles(inst->couponDates[iCoup],0), 
                                                   inst->notional*newCoupons[iCoup]));
            }
            instPaymentDates.push_back(settlement->settles(inst->couponDates[iCoup],0));
        } 

        // Obtain the index of the earlier of next coupon after valueDate and next coupon after early redemption
        if (inst->couponDates[iCoup] <= inst->valueDate) {
            nextCoup ++;
        }
        
        iCoup++;
    }

    // todo: in case of one coupon left with make whole, !receiveOvershoot, it is actually known.

    // Store barrier dates (if applicable)
    barrierDates(earlyRedeemCpn);
    
    // grab any floaters
    if (inst->hasFloater)
    {
        // compute lastGuaranteedCouponDate, to be used in capturing the fees
        DateTime lastGuaranteedCouponDate(0,0);
        
        // if early redeemed: nextCoup equals the coupon index after early-redeem coupon
        // if not early redeemed: nextCoup equals 1 + the index of the last coupon < value date
        if (earlyRedeemed || nextCoup  == inst->couponDates.size())
        {   // don't participate in any coupons after breached
            lastGuaranteedCouponDate = inst->couponDates[nextCoup-1]; 
        }
        else
        {   // if not breached, get at least one more coupon
            lastGuaranteedCouponDate = inst->couponDates[nextCoup];
        }
        
        PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
        
        CashFlowArrayConstSP floatcf(fltpay->knownCashflows(inst->valueDate,
            0,
            false,
            inst->valueDate,
            0.0));
        
        // for isKO, libor payments are only guaranteed until the first libor payment on or after the coupon
        // corresponding to lastGuaranteedCouponDate.  for !isKO, the payments are known after this point, only if
        // we have early redeemed.
        bool feeAlive = true;
        for (int i = iFirstUnKnownFloater; i < floatcf->size(); i++) {
            if ((feeAlive && inst->isKO) || (!feeAlive && !inst->isKO && earlyRedeemed))
            {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow((*floatcf)[i].date, (*floatcf)[i].amount * inst->notional));

                // instPaymentDates is empty if the instrument has early redeemed
                // so we add again all the coupon date                    
                if (earlyRedeemed) {
                    int j = (i == iFirstUnKnownFloater)? 0 : i;
                    while (j <= i){
                        DateTimeArraySP payDates(new DateTimeArray(1,(*floatcf)[j].date));
                        instPaymentDates.insert(instPaymentDates.end(), payDates->begin(),payDates->end());
                        j++;
                    }
                }
            }
            if ((*floatcf)[i].date >= lastGuaranteedCouponDate)
            {
                feeAlive = false; // game over
            }
        }              
    }

    // Flag whether early redemption has occured and cache details if it has
    hasTargetRedeemed = earlyRedeemed;
    if( hasTargetRedeemed ) {
        cachedRedemptionDetails = details;
    }
}

// isReporting indicates whether coupons are being processed for reporting (KNOWN_CASHFLOWS and BARRIER_LEVEL)
// or risk/pricing purposes
int TargetRedemptionNoteMC::processCoupons(
        double aggregateAtMaturity,
        DoubleArray &newCoupons, /* M */
        bool &earlyRedeemed,
        struct RedemptionDetails& details,
        const IPathGenerator*  pathGen) {        /* M */

    // Modify coupons for any special (e.g. all-weather conditions)
    int earlyRedeemCpn = couponsOverride(pathGen, newCoupons, earlyRedeemed);

    double myTarget = inst->targetLevel;
    
    if (doingPast()) {
        // Override with economic target for reporting purposes
        myTarget = inst->ecoTargetLevel;
    }

    applyTarget(aggregateAtMaturity, earlyRedeemCpn, earlyRedeemed, details, newCoupons, myTarget);

    return earlyRedeemCpn;

}

// This function sets 'newCoupons' to the actual cash amounts of the payoff, based on the performances in 'coupons'.
// It also sets 'earlyRedeemCpn' to the index of the coupon at the point of knockout.  If there is no knockout, then
// couponIdxAtBreach is set to the coupon.size() - 1.
void TargetRedemptionNoteMC::applyTarget(
        double aggregateAtMaturity,
        int &earlyRedeemCpn, // earlier of maturity and early redemption
        bool &earlyRedeemed,        /* M */
        struct RedemptionDetails& details,
        DoubleArray& newCoupons,
        double myTarget)   
{
    bool breached = false;
    int couponIdxAtBreach = newCoupons.size()-1;
    
    double cum = 0.0;
    double previousCoup = 0.0;  // store this for floorWithPreviousCoup
    double sumCoupon = 0.0;     // store this for onePayment
    int i;

    // Initialise redemption details
    details.couponIdx = 0;
    details.bonus = 0.0;
    details.finalCoupon = 0.0;
    details.totalCoupon = 0.0;
    details.target = 0.0;
    details.redemption = 0.0;

    // Compute coupons up to earlyRedeemCpn
    for(i=0; i<=earlyRedeemCpn; i++) {

        if (!breached) {
            double thisCoup = newCoupons[i];

            if (inst->floorWithPreviousCoup && i > 0 && thisCoup < previousCoup)
            {
                thisCoup = previousCoup;
            }
            cum += thisCoup; // todo: spot at start !!!

            // Allow small tolerance so the "==" is not at the whim of numerical accuracy
            if (!Maths::isNegative(cum - myTarget)) { 
                breached = true;
                couponIdxAtBreach=i;
                if (!inst->receiveOvershoot)
                {
                    thisCoup -= cum - myTarget;
                }
            }

            previousCoup = thisCoup;
            newCoupons[i] = thisCoup;

            // Keep track of the total coupon paid (adjusted for overshoot)
            sumCoupon += thisCoup;
            // if one payment sum the coupon
            if (inst->onePayment) {
                newCoupons[i] = 0.;
            }
        }
        else {
            newCoupons[i] = 0.0;
        }
    }

    if (earlyRedeemCpn < newCoupons.size()-1) {
        // Instrument has redeemed before last coupon date, set subsequent coupons to 0
        for (i=earlyRedeemCpn+1; i<newCoupons.size(); i++) {
            newCoupons[i] = 0.0;
        }
    }

    // Redemption occurs at earlier of target breach and lower all-weather risk barrier breach
    earlyRedeemed = earlyRedeemed || breached;
    earlyRedeemCpn = Maths::min(couponIdxAtBreach, earlyRedeemCpn);

    // need to sort out the cash bit        
    // now add Make Whole amount, if necessary.  
    if (inst->isMakeWhole && !earlyRedeemed)
    {
        newCoupons.back() += (myTarget - cum);
    }
    
    // if onePayment, pay the sum of coupon at the last date
    if (inst->onePayment) {
        newCoupons[earlyRedeemCpn] += sumCoupon;
    }

    // now add any bonus coupons
    if (breached && inst->bonusCoupons.size() != 0) {
        newCoupons[earlyRedeemCpn] += inst->bonusCoupons[earlyRedeemCpn];
        sumCoupon += inst->bonusCoupons[earlyRedeemCpn];

        // store details of bonus in redemption information
        details.bonus =  inst->bonusCoupons[earlyRedeemCpn];
    }

    if (earlyRedeemed) {
        // Record details of redemption
        details.couponIdx = earlyRedeemCpn;
        details.finalCoupon = newCoupons[earlyRedeemCpn];
        details.totalCoupon = sumCoupon;
        details.target = myTarget;
    }

    // now add principal, at point of early redeem, otherwise at maturity. 
    // Need last FV factor to properly handle settlement.  (todo: verify this!)
    if (!inst->hasFloater)
    {
        newCoupons[earlyRedeemCpn] += 1.0;
    }
    
    if (earlyRedeemed) {
        // Update the complete redemption amount
        details.redemption = newCoupons[earlyRedeemCpn];
    }

    // add option at maturity
    TrivialDoubleArray lastPerformance = aggregateAtMaturity;
    if (inst->optionAtMat && !earlyRedeemed)
    {
        IDoubleArrayModifierSP performance(inst->matPerformance->getModifier(&lastPerformance));
        performance->apply();
        newCoupons.back() += lastPerformance(); 
    }     
}

// Satisfy IHandlePaymentEvents interface
// for Tarn specific record of events (pay dates, cashflows etc)
void TargetRedemptionNoteMC::recordEvents(Control* control,
    Results* results) {
    static const string method("TargetRedemptionNoteMC::recordEvents");
    try {
        
        // PAYMENT_DATES is a list of all dates on which payments may occur
        // including past and potential future dates.
        // For TargetRedemptionNote this means all funding payment dates
        // as well as the settle date for each coupon date.
        // instPaymentDates were populated in the constructor for TargetRedemptionNoteMC
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            OutputRequestUtil::recordPaymentDates(control,results,&instPaymentDates);
        }
        
        // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
        // and be supplied for all past cash flows, and any future ones
        // that are determined.
        // For TargetRedemptionNote this means any past paid fees from the funding component, as well as 
        // past coupons.  Moreover, if the instrument is currently alive, then we are guaranteed the funding payment,
        // which may be known, corresponding to the next coupon.  Finally, there is a 'pathological' case where one coupon
        // remains and whose amount is known in the 'make whole', 'don't receive overshoot' situation.
        //
        // Note: the function cacheKnownCFs has already computed knownCFs.
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished())
        {                   
            knownCFs.recordKnownCashFlows(control, results);            
        }

        // Record BARRIER_LEVELs (if any)
        recordBarriers(control, results);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Store barriers (if applicable)
void TargetRedemptionNoteMC::barrierDates(const int couponIdxAtBreach)
{}

// Record barriers (if applicable)
void TargetRedemptionNoteMC::recordBarriers(Control* control,
                                            Results* results) const
{}

// for the LogNormal path generator
CVolRequestLNArray TargetRedemptionNoteMC::getVolInterp(const IMCPathGenerator* pathGen,
    int                     iAsset) const {
    static const string routine = "TargetRedemptionNoteMC::getVolInterp";
    CVolRequestLNArray reqarr(1);
    const DateTime&    startDate = getRefLevel()->getAllDates().front();
    const DateTime&    today = getToday();
    const DateTime&    lastSimDate = getSimSeries()->getLastDate();
    bool               fwdStarting = startDate.isGreater(today);
    double             interpLevel = pathGen->refLevel(iAsset, 0); //todo review
    
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

// Override MCProduct::retrieveEvents.
// Mechanism to support TargetRedemption event reporting.
void TargetRedemptionNoteMC::retrieveEvents(EventResults* events) const {
    if( hasTargetRedeemed )
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

// -----------------------------------------------------------------------------------
// State var version methods 

/** equivalent to InstIntoMCProduct */
TargetRedemptionNoteSVMC::TargetRedemptionNoteSVMC(
    const TargetRedemptionNote*         inst,
    const SimSeriesSP&             simSeries) :
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
        coupons(inst->couponDates.size(), 0.0),
        sum(nbAssets, 0.0),
        refLevel(nbAssets, 0.0),
        assetComps(nbAssets, 0.0),
        sumSoFar(nbAssets, 0.0),
        refLevelSoFar(nbAssets, 0.0),
        couponsSoFar(inst->couponDates.size(), 0.0),
        iCouponSoFar(0),
        nbAvgOutPerCouponDate(inst->couponDates.size(), 0),
        iFirstUnKnownFloater(0),
        hasTargetRedeemed(false),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                    getToday())),
        matDfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                                inst->instSettle, simSeries->getLastDate()))
{

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
    couponMap = DateTime::createMapping(simSeries->getAllDates(),
        inst->couponDates,
        isTrivial);
    
    // set up couponsFV: forward value factors for rainbow coupons
    DateTime matSettlementDate = settlement->settles(inst->couponDates.back(),0);
    for(iCoupon = 0; iCoupon < inst->couponDates.size(); iCoupon++) {
        
        DateTime couponPayDate = settlement->settles(inst->couponDates[iCoupon],0); 
        // While we're at it, may as well store pay dates if each coupon is paid 
        // if one payment only the last coupon date is stored
        if ((!inst->onePayment) || ((iCoupon+1) == inst->couponDates.size())) {
            instPaymentDates.push_back(couponPayDate);
        }
    }
    couponDfGen = SVGenDiscFactorSP(new SVGenDiscFactor(inst->valueDate,
                                                  inst->discount.getSP(),
                                                  inst->instSettle, 
                                                  inst->couponDates));
    
    // compute liborParticipation: amount of libor holder will recieve if early redemption occurs at coupon i.
    if (inst->hasFloater)
    {
        // Build the SV Generator
        liborLegGen = LiborLeg::LiborLegSVGenSP(inst->floater->createLiborLegSVGen(inst->discount.getSP()));

        iFirstUnKnownFloater = 0; // init to 0 in case inst is not a KO

        // add libor payment dates to instPaymentDates
        PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
        DateTimeArraySP payDates = fltpay->paymentDates();
        instPaymentDates.insert(instPaymentDates.end(), payDates->begin(),payDates->end());

        // prepare the KNONW_CASHFLOW.  All cashflow in floating leg, 
        // which are scheduled before the first coupon date,
        // are known to be paid, so add to KNOWN_CASHFLOW.
        // equity linked coupon are not calculated, because it should be calculated in cacheKnownCFS.
        if (inst->isKO)
        {
            PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
    
            CashFlowArrayConstSP floatcf(fltpay->knownCashflows(inst->valueDate,
                                                                0,
                                                                false,
                                                                inst->valueDate,
                                                                0.0));
    
            DateTime firstCpnDate = inst->couponDates[0];                
            int i = 0;
            while (i < floatcf->size() && (*floatcf)[i].date <= firstCpnDate){
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow((*floatcf)[i].date, (*floatcf)[i].amount * inst->notional));
                i ++;
            }
            iFirstUnKnownFloater = i;                
        }
    }
    
    // now attach assetComps to asset-level performance and aggregation method, 
    // and coupons to coupon-level performance
    assetPerf = IDoubleArrayModifierSP(inst->assetPerf->getModifier(&assetComps));
    assetBasket = IAggregateSP(inst->assetBasket->getAggregate(&assetComps));
    performance = IDoubleArrayModifierSP(inst->performance->getModifier(&coupons));
}

void TargetRedemptionNoteSVMC::recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const {
    // not used but required 
}

void TargetRedemptionNoteSVMC::payoff(const IPathGenerator*  pathGen,
                                      IMCPrices&                prices) {
    int    iAsset;

    // update with preserved values from doing past
    int iCoupon = iCouponSoFar;
    sum = sumSoFar;
    refLevel = refLevelSoFar;
    coupons = couponsSoFar;

    for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
        /* Special treatment for average in for first cliq */
        if (!inst->isCliquetStyle || iCoupon == 0) {
            refLevel[iAsset] = refLevelSV->refLevel(iAsset);
        } else {
            refLevel[iAsset] = refLevelSoFar[iAsset];
        }
    }

    // Same begin & end for all assets, so read from the first
    const SVPath& path = spotSV->path(0/*iAsset*/);
    int    beginIdx = path.begin();
    int    endIdx   = path.end();

    // Form the asset basket at each coupon date first
    for (int iStep=beginIdx; iStep<endIdx; iStep++) {
        bool isCouponDate = (couponMap[iStep]==0);   // true iff a coupon date
    
        for(iAsset=0; iAsset<nbAssets; iAsset++) {
            const SVPath& path = spotSV->path(iAsset);
            sum[iAsset] += path[iStep];
        
            if (isCouponDate) { 
                double avg = sum[iAsset] / nbAvgOutPerCouponDate[iCoupon];
                assetComps[iAsset] = avg / refLevel[iAsset];
                if (inst->isCliquetStyle)
                {
                    /* Reset ref level for next coupon - done per asset note*/
                    if (inst->isRefPrevAvgOut) {
                        refLevel[iAsset] = avg;
                    } else {
                        refLevel[iAsset] = path[iStep];
                    }
                }
                if (!inst->avgFromStart) {
                    sum[iAsset] = 0.0;
                }
            }
        }
    
        if (isCouponDate) {

            // Apply asset-level performances to assetComps[]
            assetPerf->apply();
        
            // then compute aggregate (e.g. worst performance)
            coupons[iCoupon] = assetBasket->aggregate();
            
            // get ready for next coupon
            iCoupon++;
        }
    }

    // preserve values for past - before we modify/sort the coupons.
    if (doingPast()){ 
        sumSoFar = sum;
        refLevelSoFar = refLevel;
        couponsSoFar = coupons;         
        iCouponSoFar = iCoupon;
    }

    // Value of final coupon. Used in option calculations.
    double aggAtMat = coupons[coupons.size()-1];   // store this before we apply the performance

    // apply perf modifier to coupon values.
    performance->apply();

    if (doingPast()){ 
        cacheKnownCFs(aggAtMat);

        // change paymentCashFlows if has breach
    }

    if (!doingPast() || !hasFuture()) {
        // Compute a payoff, but only when we have a "complete" situation : either 
        // doingPast() and all is past, or !doingPast().

        double value = computeTargetRedemptionPayoff(aggAtMat);

        // discounting done by each part so no need here
        // value *= matDfSV->firstDF();

        if (!inst->isKO && inst->hasFloater)
        {
            value += liborLegSV->getPV(inst->valueDate);
        }
        prices.add(inst->notional*value); 
    }
}
// Stateless payoff version

class TRNHistoricalContext : public IHistoricalContext {
    friend class TargetRedemptionNoteSVMC;
 
    DoubleArray sumSoFar;
    DoubleArray refLevelSoFar;
    SimpleDoubleArray couponsSoFar;         
    int iCouponSoFar;
    IHistoricalContextSP refLevelContext;

public:
    TRNHistoricalContext(int nbAssets, IHistoricalContextSP refLevelHC): 
      sumSoFar(nbAssets,0.0),refLevelSoFar(nbAssets,0.0), couponsSoFar(nbAssets,0.0),
      refLevelContext(refLevelHC)
    {
    }

    virtual void deepCopyTo( IHistoricalContext* destination ) const 
    {
        TRNHistoricalContext* dest = static_cast<TRNHistoricalContext*>( destination );
        dest->sumSoFar = sumSoFar;
        dest->refLevelSoFar = refLevelSoFar;
        dest->couponsSoFar = couponsSoFar;
        dest->iCouponSoFar = iCouponSoFar;
        refLevelContext->deepCopyTo(dest->refLevelContext.get());
    }

};


IHistoricalContextSP TargetRedemptionNoteSVMC::createHistoricalContext()    
{
    return IHistoricalContextSP(new TRNHistoricalContext(
                        nbAssets, 
                        this->refLevelSV->createHistoricalContext()));
}

IHistoricalContextSP TargetRedemptionNoteSVMC::getInitialHistoricalContext()    
{
    return IHistoricalContextSP(new TRNHistoricalContext(
                        nbAssets, 
                        this->refLevelSV->getInitialHistoricalContext()));
}



DateTimeArray TargetRedemptionNoteSVMC::getPastDates()
{
    return DateTimeArray();
}

vector<int> TargetRedemptionNoteSVMC::finalize( const DateTimeArray& simDates )
{
    return DateTime::createMapping2( 
        simDates, 
        inst->averageOutDates, 
        inst->valueDate, 
        lastPastDateIdx);
}

void TargetRedemptionNoteSVMC::statelessPayOff(
    int currentDateIdx,
    IHistoricalContextSP history,
    IMCPrices& prices ) 
{
    
    TRNHistoricalContext* hContext = 
        static_cast<TRNHistoricalContext*>(history.get()); 

    // update with preserved values from doing past
    int iCoupon = hContext->iCouponSoFar;
    sum = hContext->sumSoFar;
    refLevel = hContext->refLevelSoFar;
    coupons = hContext->couponsSoFar;
    
    int iAsset;
    for (iAsset = 0; iAsset < nbAssets; iAsset ++) {
        // Special treatment for average in for first cliq 
        if (!inst->isCliquetStyle || iCoupon == 0) {
          refLevel[iAsset] = refLevelSV->refLevel2(iAsset, hContext->refLevelContext);
        }
    }

    // Form the asset basket at each coupon date first
    bool isCouponDate = (couponMap[currentDateIdx]==0);   // true iff a coupon date

    for(iAsset=0; iAsset<nbAssets; iAsset++) {
        double spot = spotSV->getSpotPrice(iAsset);
        sum[iAsset] += spot;

        if (isCouponDate) { 
            double avg = sum[iAsset] / nbAvgOutPerCouponDate[iCoupon];
            assetComps[iAsset] = avg / refLevel[iAsset];
            if (inst->isCliquetStyle)
            {
                // Reset ref level for next coupon - done per asset note
                if (inst->isRefPrevAvgOut) {
                    refLevel[iAsset] = avg;
                } else {
                    refLevel[iAsset] = spot;
                }
            }
            if (!inst->avgFromStart) {
                sum[iAsset] = 0.0;
            }
        }
    }

    if (isCouponDate) {

        // Apply asset-level performances to assetComps[]
        assetPerf->apply();

        // then compute aggregate (e.g. worst performance)
        coupons[iCoupon] = assetBasket->aggregate();

        // get ready for next coupon
        iCoupon++;
    }

    // preserve values for past - before we modify/sort the coupons.
    hContext->sumSoFar = sum;
    hContext->refLevelSoFar = refLevel;
    hContext->couponsSoFar = coupons;         
    hContext->iCouponSoFar = iCoupon;

    double aggAtMat = coupons[coupons.size()-1];   // store this before we apply the performance

    // apply perf modifier to coupons
    performance->apply();

    if ( doingPast() && lastPastDateIdx == currentDateIdx ) {
        cacheKnownCFs(aggAtMat);

    }

    // change paymentCashFlows if has breach

    if ( currentDateIdx == inst->averageOutDates.size() - 1 ) {
    //if (!doingPast() || !hasFuture()) {
        // Compute a payoff, but only when we have a "complete" situation : either 
        // doingPast() and all is past, or !doingPast().

        double value = computeTargetRedemptionPayoff(aggAtMat);

        // discounting done by each part so no need here
        // value *= matDfSV->firstDF();

        if (!inst->isKO && inst->hasFloater)
        {
            value += liborLegSV->getPV(inst->valueDate);
        }
        prices.add(inst->notional*value); 
 
    }
}

// Modify coupons according to any special (e.g. all-weather) conditions
// Returns the coupon idx corresponding to the redemption date (ignoring target)
int TargetRedemptionNoteSVMC::couponsOverride(DoubleArray& newCoupons,
                                              bool& earlyRedeemed) {

    // No overriding on the vanilla TARN
    earlyRedeemed = false;
    const SVPath& path = spotSV->path(0/*iAsset*/);
    int    endIdx = path.end();

    int iCoupon=0;
    for (int iStep=0; iStep<endIdx; iStep++) {

        if (couponMap[iStep]==0) {
        // A coupon date

            newCoupons[iCoupon] = coupons[iCoupon];
            iCoupon ++;
        }
    }

    // Ignoring target there is no early redemption
    return inst->couponDates.size() - 1;
}
    

// Returns a double array of same size as the couponDates array.
// The i'th value represents the maturity settlement value of the unpaid funding payments
// participated in, given that early redemption occurs at coupon i.
// Note: if early redemption occurs at coupon i, then the last funding payment participated in (paid or unpaid)
// is the first one on or after the coupon date when early redemption occurs.
DoubleArraySP TargetRedemptionNoteSVMC::computeLiborParticipation(const DateTime& matSettlementDate)
{
    throw ModelException("TargetRedemptionNoteSVMC::computeLiborParticipation",
                         "Meaningless with stochastic rates");
}

// SV version does it's PV here so not needed in payoff()
double TargetRedemptionNoteSVMC::computeTargetRedemptionPayoff(
        double aggregateAtMaturity)
{
    double presentValue = 0.0; 
    DoubleArray newCoupons(coupons.size());
    bool earlyRedeemed = false;
    struct RedemptionDetails details;

    // compute newCoupons and earlyRedeemCpn
    int earlyRedeemCpn = 
        processCoupons(aggregateAtMaturity, newCoupons, earlyRedeemed, details);
    
    // now add to presentValue with PV from coupon dates.
    for (int i = 0; i <= earlyRedeemCpn; i++)
    {
        presentValue += newCoupons[i] * couponDfSV->path()[i];      
    }
    
    // for floater case, need to add the libor cashflows
    if (inst->hasFloater) {
        const CashFlowArray* cf = liborLegSV->getCashFlowArray(inst->valueDate);
        // Since we need to compute a single one of the liborParticipations here
        // that being the [couponIdxAtBreach] we can in principle know where
        // the "cum" summing should start from to save recomputing each time.
        // For now I'm being lazy...
        double cum;
        for (int iCoup = 0; iCoup <= earlyRedeemCpn; iCoup++) {
            cum = 0.0;
            int iLib = 0;
            bool done = false;
            while (!done && iLib < cf->size()) {
                const DateTime& thisLiborDate = (*cf)[iLib].date;
                if (thisLiborDate > inst->valueDate)
                {   // only include libor payments that haven't been paid yet
                    cum += (*cf)[iLib].amount * liborLegSV->payDatePV(iLib);
                }
                if (thisLiborDate.getDate() >= inst->couponDates[iCoup].getDate()) {
                    done = true;
                }
                iLib ++;
            }
        }
        double liborPart = cum * (inst->isKO? 1.0 : -1.0);

        presentValue += liborPart;
    }
    return presentValue;       
}

void TargetRedemptionNoteSVMC::cacheKnownCFs(
        double aggregateAtMaturity)
{
    DoubleArray newCoupons(coupons.size());
    bool earlyRedeemed = false;
    struct RedemptionDetails details;

    // compute newCoupons and earlyRedeemCpn
    int earlyRedeemCpn = 
        processCoupons(aggregateAtMaturity, newCoupons, earlyRedeemed, details);

    // if has early redeemed then we know exactly the payment dates
    // so the payment dates are calculated again
    if (earlyRedeemed) {
        instPaymentDates.clear();
    }

    // add the settlement date for each coupon date to the knownCFs
    int iCoup = 0; 
    int nextCoup = 0;   // used below for fees kcfs.
    
    // Add known cashflows up to and including last coupon before today
    // Add payment dates up to and including couponIdxAtBreach (which maybe 
    // a particular future coupon date if flooring on previous coupons for example)
    while(iCoup <= earlyRedeemCpn)
    {
        if (!inst->onePayment) {
            if (inst->couponDates[iCoup] <= inst->valueDate) {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow(settlement->settles(inst->couponDates[iCoup],0), 
                                                   inst->notional*newCoupons[iCoup]));
            }
            // instPaymentDates is empty if the instrument has early redeemed
            // so we add again all the coupon date
            if (earlyRedeemed) {
                instPaymentDates.push_back(settlement->settles(inst->couponDates[iCoup],0));
            }
        }
        else if (earlyRedeemed && iCoup==earlyRedeemCpn) {
            // if one Payment, we need to had the breach date and cashflow if before today
            if (inst->couponDates[iCoup] <= inst->valueDate) {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow(settlement->settles(inst->couponDates[iCoup],0), 
                                                   inst->notional*newCoupons[iCoup]));
            }
            instPaymentDates.push_back(settlement->settles(inst->couponDates[iCoup],0));
        } 

        // Obtain the index of the earlier of next coupon after valueDate and next coupon after early redemption
        if (inst->couponDates[iCoup] <= inst->valueDate) {
            nextCoup ++;
        }
        
        iCoup++;
    }

    // todo: in case of one coupon left with make whole, !receiveOvershoot, it is actually known.

    // Store barrier dates (if applicable)
    barrierDates(earlyRedeemCpn);
    
    // grab any floaters
    if (inst->hasFloater)
    {
        // compute lastGuaranteedCouponDate, to be used in capturing the fees
        DateTime lastGuaranteedCouponDate(0,0);
        
        // if early redeemed: nextCoup equals the coupon index after early-redeem coupon
        // if not early redeemed: nextCoup equals 1 + the index of the last coupon < value date
        if (earlyRedeemed || nextCoup  == inst->couponDates.size())
        {   // don't participate in any coupons after breached
            lastGuaranteedCouponDate = inst->couponDates[nextCoup-1]; 
        }
        else
        {   // if not breached, get at least one more coupon
            lastGuaranteedCouponDate = inst->couponDates[nextCoup];
        }
        
        PayStreamSP fltpay(inst->floater->makePayStream(inst->discount.get()));
        
        CashFlowArrayConstSP floatcf(fltpay->knownCashflows(inst->valueDate,
            0,
            false,
            inst->valueDate,
            0.0));
        
        // for isKO, libor payments are only guaranteed until the first libor payment on or after the coupon
        // corresponding to lastGuaranteedCouponDate.  for !isKO, the payments are known after this point, only if
        // we have early redeemed.
        bool feeAlive = true;
        for (int i = iFirstUnKnownFloater; i < floatcf->size(); i++) {
            if ((feeAlive && inst->isKO) || (!feeAlive && !inst->isKO && earlyRedeemed))
            {
                knownCFs.addKnownCashFlow(inst->discount->getCcy(), 
                                          CashFlow((*floatcf)[i].date, (*floatcf)[i].amount * inst->notional));

                // instPaymentDates is empty if the instrument has early redeemed
                // so we add again all the coupon date                    
                if (earlyRedeemed) {
                    int j = (i == iFirstUnKnownFloater)? 0 : i;
                    while (j <= i){
                        DateTimeArraySP payDates(new DateTimeArray(1,(*floatcf)[j].date));
                        instPaymentDates.insert(instPaymentDates.end(), payDates->begin(),payDates->end());
                        j++;
                    }
                }
            }
            if ((*floatcf)[i].date >= lastGuaranteedCouponDate)
            {
                feeAlive = false; // game over
            }
        }              
    }

    // Flag whether early redemption has occured and cache details if it has
    hasTargetRedeemed = earlyRedeemed;
    if( hasTargetRedeemed ) {
        cachedRedemptionDetails = details;
    }
}

// isReporting indicates whether coupons are being processed for reporting (KNOWN_CASHFLOWS and BARRIER_LEVEL)
// or risk/pricing purposes
int TargetRedemptionNoteSVMC::processCoupons(
        double aggregateAtMaturity,
        DoubleArray &newCoupons, /* M */
        bool &earlyRedeemed,
        struct RedemptionDetails& details) {

    // Modify coupons for any special (e.g. all-weather conditions)
    int earlyRedeemCpn = couponsOverride(newCoupons, earlyRedeemed);

    double myTarget = inst->targetLevel;
    
    if (doingPast()) {
        // Override with economic target for reporting purposes
        myTarget = inst->ecoTargetLevel;
    }

    applyTarget(aggregateAtMaturity, earlyRedeemCpn, earlyRedeemed, details, newCoupons, myTarget);

    return earlyRedeemCpn;
}

// This function sets 'newCoupons' to the actual cash amounts of the payoff, based on the performances in 'coupons'.
// It also sets 'earlyRedeemCpn' to the index of the coupon at the point of knockout.  If there is no knockout, then
// couponIdxAtBreach is set to the coupon.size() - 1.
void TargetRedemptionNoteSVMC::applyTarget(
        double aggregateAtMaturity,
        int &earlyRedeemCpn, // earlier of maturity and early redemption
        bool &earlyRedeemed,        /* M */
        struct RedemptionDetails& details,
        DoubleArray& newCoupons,
        double myTarget)   
{
    bool breached = false;
    int couponIdxAtBreach = newCoupons.size()-1;
    
    double cum = 0.0;
    double previousCoup = 0.0;  // store this for floorWithPreviousCoup
    double sumCoupon = 0.0;     // store this for onePayment
    int i;

    // Initialise redemption details
    details.couponIdx = 0;
    details.bonus = 0.0;
    details.finalCoupon = 0.0;
    details.totalCoupon = 0.0;
    details.target = 0.0;
    details.redemption = 0.0;

    // Compute coupons up to earlyRedeemCpn
    for(i=0; i<=earlyRedeemCpn; i++) {

        if (!breached) {
            double thisCoup = newCoupons[i];

            if (inst->floorWithPreviousCoup && i > 0 && thisCoup < previousCoup)
            {
                thisCoup = previousCoup;
            }
            cum += thisCoup; // todo: spot at start !!!

            // Allow small tolerance so the "==" is not at the whim of numerical accuracy
            if (!Maths::isNegative(cum - myTarget)) { 
                breached = true;
                couponIdxAtBreach=i;
                if (!inst->receiveOvershoot)
                {
                    thisCoup -= cum - myTarget;
                }
            }

            previousCoup = thisCoup;
            newCoupons[i] = thisCoup;
        
            // Keep track of the total coupon paid (adjusted for overshoot)
            sumCoupon += thisCoup;
            // if one payment sum the coupon
            if (inst->onePayment) {
                newCoupons[i] = 0.;
            }
        }
        else {
            newCoupons[i] = 0.0;
        }
    }

    if (earlyRedeemCpn < newCoupons.size()-1) {
        // Instrument has redeemed before last coupon date, set subsequent coupons to 0
        for (i=earlyRedeemCpn+1; i<newCoupons.size(); i++) {
            newCoupons[i] = 0.0;
        }
    }

    // Redemption occurs at earlier of target breach and lower all-weather risk barrier breach
    earlyRedeemed = earlyRedeemed || breached;
    earlyRedeemCpn = Maths::min(couponIdxAtBreach, earlyRedeemCpn);
    
    // need to sort out the cash bit        
    // now add Make Whole amount, if necessary.  
    if (inst->isMakeWhole && !earlyRedeemed)
    {
        newCoupons.back() += (myTarget - cum);
    }
    
    // if onePayment, pay the sum of coupon at the last date
    if (inst->onePayment) {
        newCoupons[earlyRedeemCpn] += sumCoupon;
    }

    // now add any bonus coupons
    if (breached && inst->bonusCoupons.size() != 0) {
        newCoupons[earlyRedeemCpn] += inst->bonusCoupons[earlyRedeemCpn];
        sumCoupon += inst->bonusCoupons[earlyRedeemCpn];

        // store details of bonus in redemption information
        details.bonus =  inst->bonusCoupons[earlyRedeemCpn];
    }
    
    if (earlyRedeemed) {
        // Record details of redemption
        details.couponIdx = earlyRedeemCpn;
        details.finalCoupon = newCoupons[earlyRedeemCpn];
        details.totalCoupon = sumCoupon;
        details.target = myTarget;
    }


    // now add principal, at point of early redeem, otherwise at maturity. 
    // Need last FV factor to properly handle settlement.  (todo: verify this!)
    if (!inst->hasFloater)
    {
        newCoupons[earlyRedeemCpn] += 1.0;
    }
    
    if (earlyRedeemed) {
        // Update the complete redemption amount
        details.redemption = newCoupons[earlyRedeemCpn];
    }

    // add option at maturity
    TrivialDoubleArray lastPerformance = aggregateAtMaturity;
    if (inst->optionAtMat && !earlyRedeemed)
    {
        IDoubleArrayModifierSP performance(inst->matPerformance->getModifier(&lastPerformance));
        performance->apply();
        newCoupons.back() += lastPerformance(); 
    }     
}

// Satisfy IHandlePaymentEvents interface
// for Tarn specific record of events (pay dates, cashflows etc)
void TargetRedemptionNoteSVMC::recordEvents(Control* control,
    Results* results) {
    static const string method("TargetRedemptionNoteSVMC::recordEvents");
    try {
        
        // PAYMENT_DATES is a list of all dates on which payments may occur
        // including past and potential future dates.
        // For TargetRedemptionNote this means all funding payment dates
        // as well as the settle date for each coupon date.
        // instPaymentDates were populated in the constructor for TargetRedemptionNoteSVMC
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            OutputRequestUtil::recordPaymentDates(control,results,&instPaymentDates);
        }
        
        // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
        // and be supplied for all past cash flows, and any future ones
        // that are determined.
        // For TargetRedemptionNote this means any past paid fees from the funding component, as well as 
        // past coupons.  Moreover, if the instrument is currently alive, then we are guaranteed the funding payment,
        // which may be known, corresponding to the next coupon.  Finally, there is a 'pathological' case where one coupon
        // remains and whose amount is known in the 'make whole', 'don't receive overshoot' situation.
        //
        // Note: the function cacheKnownCFs has already computed knownCFs.
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished())
        {                   
            knownCFs.recordKnownCashFlows(control, results);            
        }

        // Record BARRIER_LEVELs (if any)
        recordBarriers(control, results);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// Store barriers (if applicable)
void TargetRedemptionNoteSVMC::barrierDates(const int couponIdxAtBreach)
{}

// Record barriers (if applicable)
void TargetRedemptionNoteSVMC::recordBarriers(Control* control,
                                            Results* results) const
{}
    
// for the LogNormal path generator
CVolRequestLNArray TargetRedemptionNoteSVMC::getVolInterp(const IMCPathGenerator* pathGen,
    int                     iAsset) const {
    static const string routine = "TargetRedemptionNoteSVMC::getVolInterp";
    CVolRequestLNArray reqarr(1);
    const DateTime&    startDate = getRefLevel()->getAllDates().front();
    const DateTime&    today = getToday();
    const DateTime&    lastSimDate = getSimSeries()->getLastDate();
    bool               fwdStarting = startDate.isGreater(today);
    double             interpLevel = refLevelSV->refLevel(iAsset);
    
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
            double refLevel =  iCouponSoFar==0 ? refLevelSV->refLevel(iAsset) : refLevelSoFar[iAsset];
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
                
                interpLevel = (numDates * interpLevel * refLevelSV->refLevel(iAsset)
                    - sumSoFar[iAsset])/ numRemaining;
            } else {
                interpLevel *= refLevelSV->refLevel(iAsset);
            }
        }
        reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
            startDate,
            lastSimDate,
            fwdStarting));
    }
    return reqarr;
}


void TargetRedemptionNoteSVMC::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    static const string routine = "TargetRedemptionNoteSVMC::pathGenUpdated";
    
    try {
        spotSV = spotGen->getSpotSV(newPathGen);
        refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
        matDfSV = matDfGen->getSVDiscFactor(newPathGen);
        couponDfSV = couponDfGen->getSVDiscFactor(newPathGen);
        if (inst->hasFloater) {
            liborLegSV = liborLegGen->getLiborLegSV(liborLegSV, newPathGen);
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void TargetRedemptionNoteSVMC::collectStateVars(IStateVariableCollectorSP svCollector) const{
    // ask for a reference level State Variable
    svCollector->append(refLevelGen.get());
    svCollector->append(spotGen.get());
    svCollector->append(matDfGen.get());
    svCollector->append(couponDfGen.get());
    if (inst->hasFloater) {
        svCollector->append(liborLegGen.get());
    }
}

// Override MCProduct::retrieveEvents.
// Mechanism to support TargetRedemption event reporting.
void TargetRedemptionNoteSVMC::retrieveEvents(EventResults* events) const {
    if( hasTargetRedeemed )
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

//--------------------------------------------------------------------------

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* TargetRedemptionNote::createProduct(const MonteCarlo* model) const {
    
    // XXX Cliquet style resetting is not supported with implied - how enforce that?
    
    // we create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
    one */
    simSeries->addDates(averageOutDates);

    if(model->stateVarUsed()) {
        // State variables
        return new TargetRedemptionNoteSVMC(this, simSeries);
    }
    
    // otherwise, use old methodology
    return new TargetRedemptionNoteMC(this, simSeries);
}

CClassConstSP const TargetRedemptionNote::TYPE = 
    CClass::registerClassLoadMethod("TargetRedemptionNote", 
                                    typeid(TargetRedemptionNote), 
                                    TargetRedemptionNoteHelper::load);
bool  TargetRedemptionNoteLoad() {
    return (TargetRedemptionNote::TYPE != 0);
}


DRLIB_END_NAMESPACE
