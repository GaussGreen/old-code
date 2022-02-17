//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIRunTime.cpp
//
//   Description : Run time classes for SPI and Rainbow SPI
//                 Handles the algorithm in a central place
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIRunTime.hpp"
#include "edginc/SyntheticPortfolioInsurance.hpp"

DRLIB_BEGIN_NAMESPACE

SPICouponsRT::SPICouponsRT(ICouponsSPISP coupons) : coupons(coupons), B(0.0) {}

double SPICouponsRT::getCoupon(SPIRunTime* dynBask, double BF, int iStep) {
    // default  implementation. 
    double coupon = 0.0;
    B = dynBask->B;

    if (coupons->isCouponStep(iStep)) {
        coupon = coupons->getCoupon(B, BF, iStep);
    }
    return coupon;
}

double SPICouponsRT::getCumSE() { 
    return 0.0;
}

double SPICouponsRT::getCumTE() { 
    return 0.0;
}

double SPICouponsRT::getCumB() { 
    return B;
}

SPICouponsRT::~SPICouponsRT() {};

/*****************************************************************************/
SPICouponsBasketRT::SPICouponsBasketRT(ICouponsSPISP coupons) : SPICouponsRT(coupons) {}

double SPICouponsBasketRT::getCoupon(SPIRunTime* dynBask, double BF, int iStep) {
    double coupon = 0.0;
    B = dynBask->B;

    if (coupons->isCouponStep(iStep) &&
        coupons->isThresholdBreached(B, iStep)) {
        // only if basket level above threshold
        coupon = coupons->getCoupon(B, BF, iStep);
    }
    return coupon;
}

SPICouponsBasketRT::~SPICouponsBasketRT() {};

/*****************************************************************************/
// coupons contingent on TE level
SPICouponsExposureRT::SPICouponsExposureRT(ICouponsSPISP coupons) : SPICouponsRT(coupons), SE(0.0), TE(0.0) {}

double SPICouponsExposureRT::getCoupon(SPIRunTime* dynBask, double BF, int iStep) {
    double coupon = 0.0;
    SE = TE = 0.0; // clear old results - for reporting purposes
    B = dynBask->B;

    if (coupons->isCouponStep(iStep)) {
        dynBask->sustLoss(iStep, BF);
        SE = dynBask->SE; // copy so we can report cum value
        TE = dynBask->algo->targetExp(SE);
        if (coupons->isThresholdBreached(TE, iStep)) {
            // only if TE above threshold
            coupon = coupons->getCoupon(B, BF, iStep);
        }
    }
    return coupon;
}

double SPICouponsExposureRT::getCumSE() { 
    return SE;
}

double SPICouponsExposureRT::getCumTE() { 
    return TE;
}

 SPICouponsExposureRT::~SPICouponsExposureRT() {};

 /****************************************************************************/
 SPIFeesRT::SPIFeesRT(SPIRunTime* dynBasket) : 
    fees(dynBasket->fees), dynBask(dynBasket),
    firstStep(dynBasket->iStepFirstRebal), 
    lastStep(dynBasket->iStepLastRebal) {}

SPIFeesRT::~SPIFeesRT() {};

double SPIFeesRT::calculateFees(int iStep, double& feeToPay, double pBF, double pBL) {
    double currentFee = 0.0;
    if (iStep > firstStep && iStep <= lastStep) {
        // get the non-contingent fees
        currentFee = fees->getFeeAmount(iStep, dynBask->nZ, dynBask->Z(iStep), 
                                dynBask->nE, dynBask->E);
        feeToPay += fees->getPaidFeeAmount(iStep, dynBask->nZ, dynBask->Z(iStep), 
                                dynBask->nE, dynBask->E);
        // now get the contingent ones (if any)
        // note: SE is calculated on basket before any fees deducted if there are ctgt fees
        currentFee += getCtgtFee(iStep, pBF, pBL);
        feeToPay += getCtgtPaidFee(iStep);
    }
    return currentFee;
}

double SPIFeesRT::getCtgtFee(int iStep, double pBF, double pBL) {
    // default  implementation. 
    return 0.0;
}

double SPIFeesRT::getCtgtPaidFee(int iStep) {
    // default  implementation. 
    return 0.0;
}

/*****************************************************************************/

// fees contingent on SE level
SPIFeesExposureRT::SPIFeesExposureRT(SPIRunTime* dynBasket,
                    ISPIBondFloorSP bondFloor, ILockInSPI* lockIn) 
    : SPIFeesRT(dynBasket), bFloor(bondFloor), lock(lockIn) {}

SPIFeesExposureRT::~SPIFeesExposureRT() {};

double SPIFeesExposureRT::getCtgtFee(int iStep, double pBF, double pBL) {
    double fee = 0.0;
    thresholdBreached = false;

    // let's calculate the current SE
    lock->apply(pBL, dynBask->B, pBF, iStep);
    double BF = bFloor->getLevel(pBL, iStep);
    dynBask->sustLoss(iStep, BF);
    double SE = dynBask->SE; // copy so we can report cum value

    if (fees->isThresholdBreached(SE)) {
        // only if TE above threshold
        thresholdBreached = true;
        fee += fees->getContingentFeeAmount(iStep, dynBask->nZ, dynBask->Z(iStep), 
                                dynBask->nE, dynBask->E);
    }
    return fee;
}

/******
// note to save time we rely on getFee having been called first to set the threshold flag
*******/
double SPIFeesExposureRT::getCtgtPaidFee(int iStep) {
    double fee = 0.0;
    if (thresholdBreached) { // calculated in getFee
        fee += fees->getContingentPaidFeeAmount(iStep, dynBask->nZ, dynBask->Z(iStep), 
                                dynBask->nE, dynBask->E);
    }
    return fee;
}

/*****************************************************************************/
// Run time basket

// start values from the "SoFar" 
// I think I need a BSoFar to handle after mat lag-adj-only use of value()
void SPIRunTime::init(int  startIdx,
                      bool doingPast) {
    // The lock-in now depends on BF, so we need to initialise outside the loop rather
    // than relying on a value being defined inside. If there is no previous step (so first 
    // rebal is first step) we just use that again.
    BF = bondFloor->getLevel(BL, startIdx>0?startIdx-1:startIdx); 
    B  = BSoFar;
    BL = BLSoFar;
    nZ = nZSoFar;
    nE = nESoFar;
    LC = LCSoFar;
    LB = LBSoFar;
    ELag = ELagSoFar;
    nEDiffLag = nEDiffLagSoFar;
    isRebalLag = isRebalLagSoFar;
    A = ASoFar;
    RC = RCSoFar;
    cutoffStep = 0;
    sumB = sumBSoFar;
    sumCoupons = sumCouponsSoFar;
    isCutoff = isCutoffSoFar;
    accumulatedFeeToPay = accumulatedFeeToPaySoFar;
    currentFeeAmount = 0.0;
    for(int iBucket=0; iBucket<gapRiskByBucket.size(); iBucket++) {
        gapRiskByBucket[iBucket] = 0.0;
    }
}

// preserve values from past into "SoFar"
void SPIRunTime::setSoFar() {
    BSoFar  = B;
    BLSoFar  = BL;
    nZSoFar = nZ;
    nESoFar = nE;
    LCSoFar = LC;
    LBSoFar = LB;
    ELagSoFar = ELag;
    nEDiffLagSoFar = nEDiffLag;
    isRebalLagSoFar = isRebalLag;
    ASoFar = A;
    RCSoFar = RC;
    sumBSoFar = sumB;
    sumCouponsSoFar = sumCoupons;
    isCutoffSoFar = isCutoff;
    accumulatedFeeToPaySoFar = accumulatedFeeToPay;
}

SPIRunTime::SPIRunTime(const SyntheticPortfolioInsurance*  inst,
                       ILockInSPI*              lockIn,
                       InstrumentSettlementConstSP  settlement,
                       const DateTime*          instPaymentDate,
                       const DateTimeArray&     feePayDates):
    settlement(settlement),
    dynBask(inst->basket.get()),
    algo(dynBask->algorithm->getAlgorithmSPI()),
    bond(dynBask->bond->getBondSPI()),
    fees(dynBask->feesSPI->getFeesSPI()),
    coupons(dynBask->couponsSPI->getCouponsSPI()),
    loanCost(dynBask->loanCost->getLoanCostSPI()),    
    cutoff(inst->cutoff->getCutoffSPI()),
    lockIn(lockIn),
    numAlgAssets(algo->getNumRiskyAssets()),
    disc(inst->discount.get()),
    numSteps(dynBask->getRebalanceDates().size()),
    B(0.), BF(0.), BL(0.), SC(0.), UE(0.), UC(0.), SE(0.), TE(0.),
    bondPrices(numSteps), ZSpecial(0.),
    haveZSpecial(-1),
    nZ(0.),
    nE(numAlgAssets, 0.),
    A(numAlgAssets, 0.),
    LC(0.),
    LB(0.),
    RC(numAlgAssets, 0.),
    RCFactor(numSteps),
    cutoffStep(0),
    isCutoff(false),
    sumB(0.), couponAmt(0.),
    payoff(0.), sumCoupons(0.),
    accumulatedFeeToPay(0.0),
    paidFeeToPay(0.),
    currentFeeAmount(0.0),
    isRebal(false),
    debugOn(inst->debugOn),
    notional(inst->notional),
    matCashFlow(inst->matCashFlow),
    E(numAlgAssets, 0.),
    Einit(numAlgAssets, 0.),
    iStepLastRebal(0),
    loanCostAccrueFactorArray(numSteps),
    feeTimeFactorArray(numSteps),
    instKnownCashFlows(new SPIKnownFlows(inst->valueDate, disc, instPaymentDate)),
    BSoFar(dynBask->initialBasketLevel),
    BLSoFar(lockIn->getInitialLockIn()),
    nZSoFar(0.0),
    LCSoFar(0.),
    LBSoFar(0.),
    sumBSoFar(0.0),
    sumCouponsSoFar(0.0),
    isCutoffSoFar(false),
    accumulatedFeeToPaySoFar(0.0),
    terminatedEarly(false), terminalStep(-1),
    terminalCoupon(0.0), terminalTotalCoupon(0.0),
    targetCoupon(0.0), bonusCoupon(0.0), terminalTotalRedemption(0.0) {
    
    static const string routine = "SPIRunTime::SPIRunTime";
    try{

        iFirstFutureStep = inst->valueDate.findUpper(dynBask->getRebalanceDates());

        int             i, iStep;
        int             numDates = inst->basket->getRebalanceDates().size();
        int             iLastRebal = numDates - inst->basket->maxExecLagDays - 1;
        
        // To minimise memory use we form a single IntArray which flags the meaning of each date
        // XXX This is becoming somewhat irrelevant since so many other variations which have arrays
        // XXX along sample dates, so should encapsulate this. 
        sampleFlags = IntArray(inst->basket->getRebalanceDates().size());
        bool isTrivial;
        IntArray     feeNotifMap = DateTime::createMapping(inst->basket->getRebalanceDates(),
                                                            feePayDates,
                                                            isTrivial);
        
        IntArray     avgMap = DateTime::createMapping(inst->basket->getRebalanceDates(),
                                                        inst->averageOutDates,
                                                        isTrivial);

        // Establish 'sampleFlags'...
        /////////////////////////////
        // the last maxLagDays entries are just samples to allow the lag adjustment to be computed
        // and have no rebalance, fee payment, or lock-in
        for(iStep=0; iStep<inst->basket->maxPubLagDays; iStep++) {
            sampleFlags[iStep] = SAMPLE_ONLY;
        }
        for(iStep=inst->basket->maxPubLagDays; iStep<numDates - inst->basket->maxExecLagDays; iStep++) {
            sampleFlags[iStep] = REBALANCE;
            if (feeNotifMap[iStep]==0) {
                sampleFlags[iStep] |= FEE_NOTIFICATION;
            }
            if (avgMap[iStep]==0) {
                sampleFlags[iStep] |= AVERAGE;
            }
        }
        sampleFlags[iLastRebal] |= FINAL_REBAL;
        for(iStep=iLastRebal+1; iStep<numDates; iStep++) {
            sampleFlags[iStep] = SAMPLE_ONLY;
        }
        sampleFlags[numDates-1] |= AVERAGE; // checked this above
        
        int lagSize = Maths::max(dynBask->maxExecLagDays, 1); // always have at least 1 - no lag handled carefully
        ELag = DoubleMatrix(numAlgAssets, lagSize);
        nEDiffLag = DoubleMatrix(numAlgAssets, lagSize);
        isRebalLag = BoolArray(lagSize);
        ELagSoFar = DoubleMatrix(numAlgAssets, lagSize);
        nEDiffLagSoFar = DoubleMatrix(numAlgAssets, lagSize);
        isRebalLagSoFar = BoolArray(lagSize);
        nESoFar = DoubleArray(numAlgAssets, 0.0);
        ASoFar = DoubleArray(numAlgAssets, 0.0);
        RCSoFar = DoubleArray(numAlgAssets, 0.0);
        
        // Check algo is consistent with dynamic basket - for now just numAlgAssets
        if (dynBask->numPubLagDays.size() != numAlgAssets) {
            throw ModelException(routine, "Must have a publication lag supplied per asset, but " + 
                                    Format::toString(dynBask->numPubLagDays.size()) + 
                                    " pub lags given, and have " +
                                    Format::toString(numAlgAssets) + 
                                    " assets");
        }
        
        /* FIRST REBAL is a few days after the first sample date (to allow publication lag to work in) 
            Take the largest of all the asset's pub lags */
        iStepFirstRebal = dynBask->maxPubLagDays;
        
        /* FINAL_REBAL is a few days before the final sample date (to allow execution lag to play itself out)
            tf[i] is most recent t[i] such that t[i] is fee payment date : for bondFee calc this does not include current
            date; for bond floor it does. */
        iStepLastRebal = numSteps - dynBask->maxExecLagDays - 1;
        // sanity check
        if (iStepFirstRebal >= numSteps) {
            throw ModelException(routine, 
                                    "maxPubLagDays (" +
                                    Format::toString(iStepFirstRebal) + 
                                    " exceeds the number of rebalance dates (" +
                                    Format::toString(numSteps) + "!");
        }
        if (iStepLastRebal < 0) {
            throw ModelException(routine, 
                                    "maxExecLagDays (" +
                                    Format::toString(dynBask->maxExecLagDays) + 
                                    " exceeds the number of rebalance dates (" +
                                    Format::toString(numSteps) + "!");
        }
        if (!(sampleFlags[iStepLastRebal] & FINAL_REBAL)) {
            throw ModelException(routine, "Internal error");
        }
        const DateTime& finalRebalDate = dynBask->getLastRebalDate();
        // Give bond some useful info : today and maturity date. Last seems a bit nonsensical...
        bond->init(inst->valueDate, finalRebalDate, &dynBask->getRebalanceDates());
        bond->getBondPrices(disc, bondPrices);
        
        // Overwrite bond price today with that got from the special curve attached to the bond
        // "today" means the latest determined sample - but that will come from past bond values
        // To make meaningful use of the EURIBOR curve need to use it for the next future date
        // This is the value that will be fed back into the model (by Pyramd) as the
        // closing for the next daily rebalance date 
        // only bother doing this if the trade is alive and started
        if (inst->valueDate <= dynBask->getFinalDate()) {
            const DateTime& nextDailyRebalDate = dynBask->getNextRebalanceDate(inst->valueDate);
            ZSpecial = bond->getBondPriceToday(nextDailyRebalDate);
            for(i=0; i<numSteps;i++) {
                if (dynBask->getRebalanceDates()[i].getDate()==inst->valueDate.getDate()) {
                    bondPrices[i] = ZSpecial;
                    haveZSpecial = i;
                    break;
                } else if (dynBask->getRebalanceDates()[i].getDate()>inst->valueDate.getDate()) {
                    break;
                }
            }
        }
 
        // XXX This is nasty. There should be a natural place to init algo but we
        // XXX have a dependency on bondPrices now so has to be not earlier than here
        algo->init(&bondPrices); 
        loanCost->init(inst->valueDate, bond->getYC(), inst->dayCountBasis,
                       dynBask->getRebalanceDates());
        fees->init(numAlgAssets, dynBask->getRebalanceDates(), iStepFirstRebal, 
                    inst->dayCountBasis, loanCost.get(), disc);
        coupons->init(dynBask->getRebalanceDates());

        for(i=0; i<iStepFirstRebal; i++) {
            loanCostAccrueFactorArray[i] = 0.0;
        }
        for(i=iStepFirstRebal; i<numSteps;i++) {
            // capture the (Rji+SL)*(ti-tr)/Basis up-front
            loanCostAccrueFactorArray[i] = (i==0)? 0. : 
                loanCost->getAccrueFactor(disc, dynBask->getRebalanceDates()[i-1], 
                                            dynBask->getRebalanceDates()[i], 
                                            inst->dayCountBasis);
        }

        // now set up rebal cost past/future
        for (i = 0; i < numSteps; i++) {
            if (inst->valueDate.getDate() < dynBask->getRebalanceDates()[i].getDate()) {
                RCFactor[i] = dynBask->futureRebalCostRate;
            } else {
                RCFactor[i] = dynBask->pastRebalCostRate;
            }
        }

        // now make the run time coupon helper
        string type = coupons->getThresholdType(); 
        if (type == SPI_THRESHOLD_TYPE_BASKET) {
            couponsRT = SPICouponsRTSP(new SPICouponsBasketRT(coupons));
        }
        else if (type == SPI_THRESHOLD_TYPE_TE) {
            couponsRT = SPICouponsRTSP(new SPICouponsExposureRT(coupons));
        }
        else {
            couponsRT = SPICouponsRTSP(new SPICouponsRT(coupons));
        }

        // feeDF. The actual payment takes place a few days after the "paymentdate" 
        // Strictly only need feeDF on dates flagged as "fee payment dates" but
        // this is simpler (sorry)
        const DateTime& matSettleDate = settlement->settles(inst->basket->getRebalanceDates().back(), 
                                                            0); // asset optional
        double matPV = disc->pv(matSettleDate);
        int numRebalDates = inst->basket->getRebalanceDates().size();
        feeDFArray = DoubleArray(numRebalDates);
        couponDFArray = DoubleArray(numRebalDates, 0.0);
        for(iStep=0; iStep<numRebalDates; iStep++) {
            const DateTime feePayDate = fees->getFeePayDate(inst->basket->getRebalanceDates()[iStep]);
            feeDFArray[iStep] = 1.0 / disc->pv(feePayDate, matSettleDate);
            couponDFArray[iStep] = coupons->getCouponPVFactor(iStep, disc) / matPV;
        }

        // Internal class to handle interdependencies for bond floor calculation
        bondFloor = ISPIBondFloorSP(SPIBondFloor::make(this,
                                                    inst->dayCountBasis,
                                                    inst->valueDate));

        if (fees->hasContingentFees()) {
            feesRT = SPIFeesRTSP(new SPIFeesExposureRT(this, bondFloor, lockIn));
        } else {
            feesRT = SPIFeesRTSP(new SPIFeesRT(this));
        }

        postCutoff = inst->cutoff->getCutoffSPI()->getPostCutoff(iLastRebal);

        /** Gap Risk **/
        // Translate gapRiskBucketOffsets(offsets from valueDate) into gapRiskBucketDates
        // Always put in first date today. This means I can trust all values will be captured
        // and by definition we have 0 gap risk before today.
        StringArray offsets = inst->gapRiskBucketOffsets;
        if (offsets.size()<1) {
            // We provide a default list if none is provided
            string defaultOffsets[9] = {"1Y","2Y","3Y","4Y","5Y","7Y","10Y","20Y","30Y"};
            int num = sizeof(defaultOffsets)/sizeof(string);
            offsets = StringArray(num);
            for(i=0; i<num; i++) {
                offsets[i] = defaultOffsets[i];
            }
        }
        gapRiskBucketDates = DateTimeArray(1+offsets.size());
        gapRiskBucketDates[0] = inst->valueDate;
        for(i=0; i<offsets.size(); i++) {
            MaturityPeriod m(offsets[i]);
            gapRiskBucketDates[i+1] = m.toDate(inst->valueDate);
        }
        DateTime::ensureIncreasing(gapRiskBucketDates,
                                    "Gap risk offsets must produce increasing date list", 
                                    true /*failIfEmpty*/);
        
        // Only needed if requested but we don't know that yet so always allocate
        gapRiskBucketIdx = IntArray(numDates);
        gapRiskByBucket = DoubleArray(gapRiskBucketDates.size());
        // Contribution in bucket i means gap risk aggregated between dates BucketDate[i] and BucketDate[i+1]
        // Last bucket has holds all future contributions
        int lastBucket = gapRiskBucketDates.size()-1;
        int iBucket = 0;
        for(iStep=0; iStep<numDates; iStep++) {
            // as soon as we're putting values into the last bucket we stay doing so
            if (iBucket < lastBucket &&
                inst->basket->getRebalanceDates()[iStep] > gapRiskBucketDates[iBucket+1]) {
                iBucket++;
            }
            gapRiskBucketIdx[iStep] = iBucket;
        }

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

// single algorithm step - returns true if we terminate early at this step
bool SPIRunTime::singleSPIStep(bool doingPast, int iStep, 
                               int endIdx, double& futurePayoff) {
    static const string routine("SPIRunTime::singleSPIStep");

    double UL = 0.0, SL = 0.0;
    int flag = sampleFlags[iStep];
    bool isFirstRebalStep = (iStep == iStepFirstRebal);
    bool isRebalancing = (iStep > iStepFirstRebal);
    // reset
    payoff = 0.0;
    paidFeeToPay = 0.0;

    // Compute current dynamic basket value. This has not yet paid any coupons.
    B = value(iStep);

    // now remove fees
    currentFeeAmount = feesRT->calculateFees(iStep,  accumulatedFeeToPay, 
                                             BF, BL);
    B -= currentFeeAmount;

    // Any lock-in. Note this uses the values of BL & BF from the previous step.
    lockIn->apply(BL, B, BF, iStep);

    // Compute current bond floor level
    BF = bondFloor->getLevel(BL, iStep);

    // Account for coupons. Note the circularity here - Basket(ex coupons) depends on BF; which depends
    // on BL, which depends on Basket(cum coupon). The BF may itself also be affected by coupons (any min coupon>0)
    // but that is FUTURE coupons only. May be cleaner to have a method on runTimeBask which
    // adjusts "B" by paying the coupon.
    bool isCouponTargetMet = false;
    couponAmt = couponsRT->getCoupon(this, BF, iStep);
    if (!Maths::isZero(couponAmt)) {
        // Target coupon sum and possible early redepmtion
        double earlyRedemptionAmt = 0.0;
        // note this function will update the sumCoupons and possibly
        // amend the actual coupon amount if the cap is applied
        isCouponTargetMet = coupons->isTargetMet(iStep, couponAmt, sumCoupons,
                                                 earlyRedemptionAmt);
        // now subtract the finalised coupon from the basket
        B -= couponAmt;

        // if we've breached we get the maturity cash flow early
        earlyRedemptionAmt += matCashFlow;
        if (doingPast) {
            double knownPayment = notional * 
                (couponAmt+ (isCouponTargetMet?earlyRedemptionAmt:0.));
            // This also takes care of notified payments which have not yet been paid.
            instKnownCashFlows->addFlow(coupons->getCouponPayDate(iStep),
                                        knownPayment);
            if (isCouponTargetMet) {
                // store some stuff for events
                terminalStep = iStep;
                terminalCoupon = couponAmt;
                terminalTotalCoupon = sumCoupons;
                targetCoupon = coupons->getTarget();
                bonusCoupon = earlyRedemptionAmt - matCashFlow;
                // Same value as knownPayment (above)
                terminalTotalRedemption = couponAmt + earlyRedemptionAmt;

                if (!Maths::isZero( accumulatedFeeToPay)) {
                    const DateTime& payDate = 
                        fees->getFeePayDate(dynBask->getRebalanceDates()[iStep]);
                    instKnownCashFlows->addFlow(payDate,
                                                notional *  accumulatedFeeToPay);
                }
            }
        } else {
            futurePayoff += couponAmt * couponDF(iStep);
        }
        // If we've hit the target we can just finish now
        if (isCouponTargetMet) {
            if(!doingPast) {
                futurePayoff += couponDF(iStep) * earlyRedemptionAmt // paid at coupon pay date
                            + feeDF(iStep) *  accumulatedFeeToPay; // paid immediately
            }
            else {
                terminatedEarly = true; // ensures we don't try to simulate future
            }
            // make sure debug info is correct - stops here
            if (debugOn || doingPast) {
                // approximate since not all paid on same day (see above)
                payoff = earlyRedemptionAmt; // coupon handled separately
                paidFeeToPay = accumulatedFeeToPay;
            }
            return true;
        }
    }

    // any fees to take out at this step?
    // Once a fee is "notified" it is no longer included in the SPI but
    // rather becomes a simple cash flow pending its payment. Note this can
    // get a bit complicated if there are overlapping notification & payment
    // dates. We need to keep track of which fees have been notified and not
    // paid.
    // We model accumulation of fees between notification dates, hence the reset here.
    if (flag & FEE_NOTIFICATION) {
        if (doingPast) {
            // The idea is that when a flow is past we know what it is (so can enter it
            // as a KNOWN_CASHFLOW) and it then no longer counts in the option payoff.
            // Performance not really an issue since this is done for past only
            // XXX Would be nicer to reuse instPaymentDates...
            const DateTime& payDate = fees->getFeePayDate(dynBask->getRebalanceDates()[iStep]);
            instKnownCashFlows->addFlow(payDate,
                                        notional *  accumulatedFeeToPay);
        } else {
            futurePayoff +=  accumulatedFeeToPay * feeDF(iStep);
        }
        // Start accumulating for next fee payment date
        paidFeeToPay =  accumulatedFeeToPay; // so client report is complete on fee notif dates
        accumulatedFeeToPay = 0.0;
    }

    // Compute Unbalanced position values
    UE = isCutoff ? postCutoff->getPostCutoffExposure(B, BF, iStep) : 
            isFirstRebalStep ? 0.0 : unbalExposure(iStep);
    UC = isCutoff ? postCutoff->getPostCutoffCrash() : 
            isFirstRebalStep ? 0.0 : algo->unbalCrash(nE, E, iStep);
    UL = UE * UC;

    SL = sustLoss(iStep, BF);
    TE = algo->targetExp(SE);
    
    // calculate the average after coupon has been deducted
    if (flag & AVERAGE) {
        sumB += B;
    }

    // Decide if cutoff (wait until we're rebalancing first, else no exposure) or rebalance. 
    bool isCutoffNow = isRebalancing && cutoff->isCutoff(UE, 
                                                         TE, B, BF);
    if (isCutoffNow && !cutoffStep) {
        cutoffStep = iStep; // record first cutoff
    }
    isRebal = (isFirstRebalStep || 
                    (flag & FINAL_REBAL) ||
                    ((flag & REBALANCE) &&
                        !isCutoff &&
                        (algo->doRebalance(UL, SL, iStep) || isCutoffNow)));
    
    if (isCutoffNow) {  // the cutoff decision needed for the target exposure calc is now or before
        isCutoff = true;
    }

    // Update dynamic basket (rebalance)
    double reweightExposure;
    double reweightCrash;
    if (isRebal) { // rebalance according to sustainable target exposure 
        reweightExposure = isCutoff? postCutoff->getPostCutoffExposure(B, BF, iStep) : TE;
        reweightCrash = isCutoff? postCutoff->getPostCutoffCrash() : SC;
    } else {       // rebalance according to unbalanced exposure
        reweightExposure = UE;
        reweightCrash = UC;
    }

    rebalanceAndPrepareForNextStep(iStep, isRebal, 
                                                reweightCrash, reweightExposure);

    // Since only need to calc this on a pricing run it would 
    // be cleaner and faster to condition the code out, but
    // not sure how I can decide. XXX
    // NB Done AFTER nE has been rebalanced for this step, and no need to
    // measure on the final rebal date (hence the strict "<")
    if (iStep<iStepLastRebal) {
        double stepGapRisk = algo->gapRiskAtStep(B, BF, 
                                                 equityComponent(),
                                                 iStep);
        // If skipping dates need to scale
        dynBask->scaleGapRisk(iStep, stepGapRisk);
        gapRiskByBucket[gapRiskBucketIdx[iStep]] += stepGapRisk;
    }

    return false;
}

bool SPIRunTime::getBondClosing(double& closing) {
    closing = ZSpecial;
    return haveZSpecial>=0;
}

// the i is step index
double SPIRunTime::value(int i) {
    if (i==iStepFirstRebal) {
        return (B=dynBask->initialBasketLevel);
    }
    int iAsset;
    if (i<=iStepLastRebal) { // before and on mat
        // LC : not sure about the spec here... 
        // might be nicer to treat LC at [i], rather than refer to [i-1]?
        LC = 0.0;
        if (sampleFlags[i] & REBALANCE) {
            // accrue more cost each rebalance date
            // spec says accrual is on LB at previous rebal date, not at prev sample date
            // but then the Rji are a matrix of rates. resolve what is wanted! XXX 
            LC += LB * loanCostAccrueFactor(i);
        }
        
        B = nZ * Z(i) - LC - LB;
        for(iAsset=0; iAsset<numAlgAssets; iAsset++) {
            B += nE[iAsset] * E[iAsset];
        }
    } 
    
    // trading lag adjusts - note that post FINAL_REBAL this routine is purely an adjustment, and
    // there is no assignment to B. 
    for(iAsset=0; iAsset<numAlgAssets; iAsset++) {
        B -= lagAdj(i, iAsset);
        B -= rebalCost(i, iAsset);
    }
    
    return B;
}

// needed for gap risk only!
double SPIRunTime::equityComponent() {
    double eq = 0.0;
    for(int iAsset=0; iAsset<numAlgAssets; iAsset++) {
        eq += nE[iAsset] * E[iAsset];
    }
    return eq;
}

double SPIRunTime::rebalCost(int iStep, int iAsset) {
    RC[iAsset] = 0.; 
    // turn off cost if skipping - not if lag is zero??
    if (dynBask->doRebalCost(iStep)) {
        int iLag = dynBask->maxExecLagDays > 0 ? 
            (iStep - dynBask->numExecLagDays[iAsset]) % dynBask->maxExecLagDays : 0;
        // While this is <0 we are in the initial steps where we are still gathering
        // info about levels for lag adj. We do NOT perform any rebal cost here. Note 
        // the condition is applied PER ASSET, so one asset may be being lag adjusted
        // before another. One exception to this is if we've overridden the rebal costs
        // for the initial/final days - we then do no adjustment
        if (iStep == dynBask->numExecLagDays[iAsset] && dynBask->overrideInitialRebalCost) {
            iLag = -1;
        }
        if (iStep == iStepLastRebal + dynBask->numExecLagDays[iAsset] 
                && dynBask->overrideFinalRebalCost) {
            iLag = -1;
        }
        if (iLag>=0 && isRebalLag[iLag]) {
            // we want (nE(i-g) - nE(i-g-1)) * E(i)
            RC[iAsset] = RCFactor[iStep] * fabs(nEDiffLag[iAsset][iLag] * E[iAsset]);
        }
    }
    return RC[iAsset];
}

// this could be const if it wasn't for the debug array A[]!
double SPIRunTime::lagAdj(int iStep, int iAsset) {
    A[iAsset] = 0.; 
    // no lag means no lag adjustment; also not if skipping
    if (dynBask->doLagAdj(iStep)) {
        int iLag = (iStep - dynBask->numExecLagDays[iAsset]) % dynBask->maxExecLagDays;
        // While this is <0 we are in the initial steps where we are still gathering
        // info about levels for lag adj. We do NOT perform any lag adj here. Note 
        // the condition is applied PER ASSET, so one asset may be being lag adjusted
        // before another. One exception to this is if we've overridden the lag adjustment
        // for the initial day we then do no adjustment
        if (iStep == dynBask->numExecLagDays[iAsset] && dynBask->overrideInitialLag) {
            iLag = -1;
        }
        if (iLag>=0 && isRebalLag[iLag]) {
            A[iAsset] = nEDiffLag[iAsset][iLag] * (E[iAsset] - ELag[iAsset][iLag]);
        }
    }
    return A[iAsset];
}

double SPIRunTime::unbalExposure(int i) {
    double z = nZ * Z(i);
    double e = 0.;
    for(int iAsset=0; iAsset<numAlgAssets; iAsset++) {
        e += nE[iAsset] * E[iAsset];
    }
    if (Maths::isZero(z+e)) {
        // can happen after last rebal date
        return 0.0;
    }
    // XXX The next 2 clauses (shame to have ifs) try to cope with a big crash situation
    // XXX where B->0 and possibly beyond 0
    if (Maths::isZero(e)) {
        return 0.0;
    }
    if (!Maths::isPositive(B)) {
        // XXX Could have negative B too, which makes nonsense UE. This limit is more
        // XXX sensible : shame it is arbitrary... Would be nice to think that the sequence
        // XXX of rebal/cutoff events would mean this line is not actually executed!
        return VERY_BIG;
    }
    double UE = e/(z+e) * (B+LB)/B;
    return UE;
}

double SPIRunTime::sustLoss(int iStep, double BF) {
    // Cope with extreme case when B may have gone -ve due to extreme equity crash and high fees to pay
    double buffer = Maths::isPositive(B)? 1.-BF/B : 0.0;
    SC = algo->sustCrash(buffer, iStep);
    SE = Maths::max(buffer/SC, 0.);
    return SE * SC;
}

void SPIRunTime::rebalanceAndPrepareForNextStep(int iStep, bool isRebal, 
                                    double crash, double exposure) {
    // recording HERE for use later (so no "-g"). Also cope with no lag
    int iLag = dynBask->maxExecLagDays>0 ? (iStep % dynBask->maxExecLagDays) : 0; 

    // only defined for 1 or 2 assets
    const DoubleArray& pE = algo->rebalWeights(crash);

    bool   afterFinalRebal = (iStep>=iStepLastRebal); //includes final rebal date
    nZ = afterFinalRebal ? 0. : Maths::max(0., 1. - exposure)*B/Z(iStep);
    for(int iAsset=0; iAsset<numAlgAssets; iAsset++) {

        if (isRebal) {
            ELag[iAsset][iLag] = E[iAsset];
            nEDiffLag[iAsset][iLag] = -nE[iAsset];
        }
        nE[iAsset] = afterFinalRebal ? 0. : pE[iAsset]*exposure*B/E[iAsset];

        // complete the job - so they really are diffs
        if (isRebal) {
            nEDiffLag[iAsset][iLag] += nE[iAsset];
        }
    }
    // after mat we have validated that there are no rebal flags
    isRebalLag[iLag] = isRebal;

    // update class vars to know we've moved a step on
    // note that here since we are updating LB at the very end of the current sample date, its
    // value throughout the above is actually for the previous sample date
    LB = B * Maths::max(exposure - 1., 0.);
}

// Vol interp = Cash(i)/( H(i) * Z(i) ) with H(i) = nE1(i)*E1(0)+nE2(i)*E2(0)
// and Cash(i) = BF/nZ.Z/F/LC etc
double SPIRunTime::volInterp(int endIdx,      // endIdx of the past path gen
                            double BFSoFar,
                            double strike) {

    double interp;
    if (endIdx<=iStepFirstRebal) { // no past yet
        interp = Maths::isZero(strike)? 1.0 : strike; // XXX DO ME!
    } else {
        int    i = endIdx-1; // last past step

        // compute Cash ...
        double myLC = 0.0;
        if (sampleFlags[i] & REBALANCE) {
            // accrue more cost each rebalance date
            // spec says accrual is on LB at previous rebal date, not at prev sample date
            // but then the Rji are a matrix of rates. resolve what is wanted! XXX 
            myLC += LBSoFar * loanCostAccrueFactor(i);
        }
        // fee - note ignoring equity fees for this purpose
        DoubleArray noEquityAlloc(numAlgAssets, 0.0);
        double F = fees->getFeeAmount(i, nZ, Z(i), noEquityAlloc, Einit) + myLC;

        double pvFactor = disc->pv(dynBask->getRebalanceDates().back());
        // Poor excuse time ... should properly use "strike * pvFactor" but that changes
        // lots of tests which have daft past bond prices. Those levels mean the equivalent discounting
        // factor is stupid and there's a nonsensical mismatch : Z[i] is much < pvFactor. BUT 
        // the definition of bondFloorFactor means bff should be > pvFactor (because of fees).
        // So, for reasonable data this change (addition of Max()) should make no difference.
        // Rather than fuss with loads of test cases I've decided to use Z[i] instead of 
        // pvFactor. Then no diffs and the "error" is small. Maybe fix properly later... XXX
        double myStrike = Maths::max(BFSoFar, strike*Z(i) /*pvFactor*/);
        double cash = myStrike - nZSoFar * Z(i) + F + LBSoFar;
        
        // compute normalising factor
        double H = 0.0;
        for(int iAsset=0; iAsset<numAlgAssets; iAsset++) {
            H += nE[iAsset] * Einit[iAsset];
        }
        if (Maths::isZero(H)) {
            // if no equity alloc we don't really care ...
            interp = VERY_BIG;
        } else {
            interp = cash / (H * pvFactor);
        }
    }
    return interp;
}

// Isolate IR-dependent elements behind a function
double SPIRunTime::Z(int iStep) const {
    return bondPrices[iStep];
}

double SPIRunTime::loanCostAccrueFactor(int iStep) {
    return loanCostAccrueFactorArray[iStep];
}

double SPIRunTime::feeTimeFactor(int iStep) {
    throw ModelException("SPIRunTime::feeTimeFactor", "got here");
    return feeTimeFactorArray[iStep];
}

double SPIRunTime::feeDF(int iStep) {
    return feeDFArray[iStep];
}

double SPIRunTime::couponDF(int iStep) {
    return couponDFArray[iStep];
}

void SPIRunTime::refreshDiscountFactorsQuotients(DoubleArray&      DFQuot) const {
    throw ModelException("SPIRunTime::refreshDiscountFactorsQuotients",
                         "Internal error - method should be over-ridden!");
}


/** --------------------------------------------------------------- */
/*  StateVar compliant form */
/*  --------------------------------------------------------------- **/

SPIRunTimeSV::SPIRunTimeSV(const SyntheticPortfolioInsurance* inst,
                           ILockInSPI*                        lockIn,
                           InstrumentSettlementConstSP        settlement,
                           const DateTimeArray&               feePayDates):
    SPIRunTime(inst, lockIn, settlement, 
               0,  // use default instPaymentDate == today, so no PV adj 
               feePayDates) {
    
    static const string routine = "SPIRunTimeSV::SPIRunTimeSV";
    try{
        DateTimeArray couponPayDates = coupons->getPaymentDates();
        if (couponPayDates.size()>0) {
            dfGenCoupon = SVGenDiscFactorSP(new SVGenDiscFactor(inst->valueDate,
                                                          inst->discount.getSP(),
                                                          couponPayDates));
        }
        // may need a feeDF on a coupon date if there's early termination
        DateTimeArray couponDates = coupons->getEssentialDates();
        DateTimeArray feeDFDates = DateTime::merge(couponDates, feePayDates);
        DateTimeArray feeDFPayDates = DateTimeArray(feeDFDates.size());
        for (int j = 0; j < feeDFDates.size(); ++j) {
            feeDFPayDates[j] = fees->getFeePayDate(feeDFDates[j]);
        }
        if (feeDFDates.size()>0) {
            dfGenFee = SVGenDiscFactorSP(new SVGenDiscFactor(inst->valueDate,
                                                       inst->discount.getSP(),
                                                       feeDFPayDates));
        }
        feeDatesMap.clear(); // let's not leak
        feeDatesMap = DateTime::getProjection(inst->basket->getRebalanceDates(),
                                              feeDFDates);

        // These may be sparse
        const DateTimeArray& rebalDates = dynBask->getRebalanceDates();
        if (rebalDates.size()<1) {
            throw ModelException(__FUNCTION__,
                                 "No rebalance dates!");
        }
        dfGenRebalDates.clear();
        dfGenRebalDates = vector<SVGenExpectedDiscFactorSP>(rebalDates.size()-1);
        for (unsigned int i=0; i<dfGenRebalDates.size(); i++) {
            dfGenRebalDates[i] = SVGenExpectedDiscFactorSP(
                new SVGenExpectedDiscFactor(rebalDates[i], // calcDate
                                            rebalDates[i], // pvDate - SRM requires this equal to calcDate
                                            inst->discount.getSP(),
                                            DateTimeArray(1, rebalDates[i+1]),
                                            false /*computeLog*/));
        }
        dfSVRebalDates.clear();
        dfSVRebalDates = vector<SVExpectedDiscFactorSP>(dfGenRebalDates.size());
        // bond is a different set of expected DFs - from mat to each sim date
        DateTimeArray bondMat(1, dynBask->bond->getBondSPI()->getEssentialDates().back());
        dfGenBond.clear();
        for (int i=0; 
             i<rebalDates.size() && rebalDates[i].getDate()<bondMat[0].getDate(); i++) {
            dfGenBond.push_back(SVGenExpectedDiscFactorSP(
                                    new SVGenExpectedDiscFactor(rebalDates[i], // calcDate
                                                                rebalDates[i], // pvDate - SRM requires this equal to calcDate
                                                                inst->discount.getSP(),
                                                                bondMat,
                                                                false /*computeLog*/)));
        }
        dfSVBond.clear();
        dfSVBond = vector<SVExpectedDiscFactorSP>(dfGenBond.size());
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

void SPIRunTimeSV::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    static const string routine = "SPIRunTimeSV::pathGenUpdated";

    try {
        if (dfGenCoupon.get()) {
            dfSVCoupon = dfGenCoupon->getSVDiscFactor(newPathGen);
        }
        if (dfGenFee.get()) {
            dfSVFee = dfGenFee->getSVDiscFactor(newPathGen);
        }
        for(unsigned int i=0; i<dfSVRebalDates.size(); i++) {
            dfSVRebalDates[i] = dfGenRebalDates[i]->getSVExpectedDiscFactor(newPathGen);
        }
        for(unsigned int i=0; i<dfSVBond.size(); i++) {
            dfSVBond[i] = dfGenBond[i]->getSVExpectedDiscFactor(newPathGen);
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
void SPIRunTimeSV::collectStateVars(IStateVariableCollectorSP svCollector) const{
    if (dfGenCoupon.get()) {
        svCollector->append(dfGenCoupon.get());
    }
    if (dfGenFee.get()) {
        svCollector->append(dfGenFee.get());
    }
    for(unsigned int i=0; i<dfSVRebalDates.size(); i++) {
        svCollector->append(dfGenRebalDates[i].get());
    }
    for(unsigned int i=0; i<dfSVBond.size(); i++) {
        svCollector->append(dfGenBond[i].get());
    }
}

// Called at start of each payoff() call
void SPIRunTimeSV::init(int  startIdx,
                        bool doingPast) {
    // non-SV stuff
    SPIRunTime::init(startIdx, doingPast);

    if (!doingPast) {
        // then need to refresh bond prices etc
        bond->getFutureBondPrices(dfSVBond, iFirstFutureStep,
                                  bondPrices);
        if (haveZSpecial>=0) {
            // questionable but for now use the deterministic one....?
            bondPrices[haveZSpecial] = ZSpecial;
        }
        algo->init(&bondPrices); 
        loanCost->refresh(&dfSVRebalDates);
        if (dfGenFee.get()) {
            fees->refresh(iFirstFutureStep);
        }
        // only future steps
        for(int i=iFirstFutureStep+1; i<numSteps;i++) {
            loanCostAccrueFactorArray[i] = 
                loanCost->getFutureAccrueFactor(i);
        }
        bondFloor->refresh(iFirstFutureStep);
    }
}

double SPIRunTimeSV::Z(int iStep) const {
    return bondPrices[iStep];
}

double SPIRunTimeSV::loanCostAccrueFactor(int iStep) {
    return loanCostAccrueFactorArray[iStep];
}

double SPIRunTimeSV::feeTimeFactor(int iStep) {
    throw ModelException("SPIRunTimeSV::feeTimeFactor",
                         "not implemented");
}

double SPIRunTimeSV::feeDF(int iStep) {
    if (dfSVFee.get() &&
        feeDatesMap[iStep]>=0) {
        return dfSVFee->path()[feeDatesMap[iStep]];
    }
    // perhaps throw exception?
    return 0.0;
}

double SPIRunTimeSV::couponDF(int iStep) {
    if (dfSVCoupon.get()) {
        return coupons->getCouponPVFactor(iStep,
                                          dfSVCoupon);
    }
    // should be irrelevant
    return 0.0;
}

void SPIRunTimeSV::refreshDiscountFactorsQuotients(DoubleArray&      DFQuot) const {
    for (int i = 0; i < numSteps - 1 ; i++) {
        DFQuot[i] = dfSVRebalDates[i]->firstDF(); // expected DF from [i+1] to [i]
    }
}


DRLIB_END_NAMESPACE
