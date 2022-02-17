//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRibTurbo.cpp
//
//   Description : Rib Turbo component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KRibTurbo.hpp"

DRLIB_BEGIN_NAMESPACE

class KRibTurboTree : public FDProduct {
public:
    /******************** variables ********************/
    KRibTurboConstSP  inst;
    FDProductArray    legFXprods;
    // support for different floating index per coupon per leg
    // indexProds[i][j] where i = leg number, j = coupon number
    vector<FDProductArray > legIndexProds;
    CouponSchedDates& sched;
    KRib2TreeSP ribProd;

    /******************** methods ********************/
    KRibTurboTree(const KRibTurboConstSP &inst, FDModel* model);

    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int &step, FDProduct::UpdateType updateType);
    virtual const TreeSlice& getValue(int step, DateTime eventDate) const;

    virtual string getOutputName() const { return inst->outputName; }
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results); 

protected:
    void addCoupon(TreeSlice &accuSlice, int step, int cpnIdx) const;
    void getFullLegValue(int step, DateTime eventDate, TreeSlice &fullLegSliceL) const;
    void getForwardValue(int step, DateTime eventDate, TreeSlice &fwdMainSliceL) const;
    void calculateFullCoupon(int step, int cpnIdx, DateTime resetDate,
                             TreeSlice& couponSliceL) const;

private:
    TreeSliceSP mainSlice;
    vector<TreeSliceSP> couponSlice;  // to store each coupon until added to turbo

    // slices used in getValue/stub calculations - never DEVd
    TreeSliceSP getValueSlice;
    TreeSliceSP stubCpnSlice;
    
    // principal indicators to simplify calc
    vector<vector<bool> >payCpnPrincipal;
    vector<bool> payInitCpnPrincipal;

    // domestic discounting curve name 
    string discountCurveName;

    // calculate this once and store as member variable rather than
    // recalculating numerous times.  min full accrual date is 
    // min(fullAccrDatesArray) where fullAccrDates are 
    // min(couponReset, couponAccStart)
    DateTime minFullAccrDate;

};

/********************************** KRibTurboTree *******************************/

KRibTurboTree::KRibTurboTree(const KRibTurboConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst), sched(*inst->sched) {
    try {
        // Create the underlying reset FDProducts - if one IProdCreator is provided
        // then we assume same index product for all reset dates.  If separate IProdCreators
        // provided, we have to create separate FDProducts

        // only support one leg if instrument payment/discount curve is not the same as
        // model discount curve, and leg must be same currency/discount curve as instrument
        discountCurveName = inst->getActualDiscount()->getName();
        int nbLegs = inst->legs->size();
        if (discountCurveName != inst->discount->getName()) {
            if (nbLegs != 1)
                throw ModelException("Product tree pricer only supports one turbo leg when "
                    "instrument discount curve " + inst->discount->getName() + " is different "
                    "from the model discount/pricing curve " + discountCurveName + 
                    ". Current instrument supplied with " + Format::toString(nbLegs) + " legs");
            int legIdx;
            for (legIdx=0; legIdx < nbLegs; ++legIdx) {
                KRibTurbo::IndexLeg &leg = *(*inst->legs)[legIdx];
                if (leg.discount->getName() != inst->discount->getName())
                    throw ModelException("Product tree pricer only supports the turbo leg "
                        "defined with the same discount curve as the instrument discount "
                        "curve (" + inst->discount->getName() + ") when the instrument discount "
                        "curve is not the same as the model discount/pricing curve " +
                        discountCurveName + ". Leg discount curve supplied = " + 
                        leg.discount->getName());
            }
        }

        int nbCoupons = sched.nbCoupons;
        legFXprods.resize(nbLegs);  // only 1 FX per leg required
        legIndexProds.resize(nbLegs);
        int legIdx;
        for (legIdx=0; legIdx < nbLegs; ++legIdx) {

            KRibTurbo::IndexLeg &leg = *(*inst->legs)[legIdx];
            // FX IndexSpecs
            legFXprods[legIdx] = model->createProduct(leg.fx);
            legFXprods[legIdx]->addModelResetDates(sched.reset, sched.pay);
            // Fixing indices
            legIndexProds[legIdx].resize(nbCoupons);
            if (leg.indexes.size() == 1) {
                FDProductSP prod = model->createProduct(leg.indexes[0]);
                for (int i = 0; i < nbCoupons; i++)
                    legIndexProds[legIdx][i] = prod; 
            }
            else if (leg.indexes.size()) {
                ASSERT(leg.indexes.size() == nbCoupons);
                for (int i = 0; i < nbCoupons; i++)
                    legIndexProds[legIdx][i] = model->createProduct(leg.indexes[i]); 
            }
            // register modelResetDates
            for (int cpnIdx = 0; cpnIdx < nbCoupons; cpnIdx++) {
                if (!Maths::isZero((*leg.weights)[cpnIdx])) {
                    legIndexProds[legIdx][cpnIdx]->addModelResetDate(
                        sched.reset[cpnIdx], sched.accEnd[cpnIdx]);
                }
            }
        }

        // create the RIB product if needed
        if (inst->rib.get()) {
            FDProductSP prod = model->createProduct(inst->rib);
            ribProd = DYNAMIC_POINTER_CAST<KRib2Tree>(prod);
        }

        // account for most general case where the reset dates may not be
        // strictly increasing, so the first coupon reset date may not be the 
        // first reset date of the turbo
        minFullAccrDate = min(sched.accStart[0], sched.reset[0]);
        for (int i = 1; i < nbCoupons; i++) {
            DateTime fullAccrDate = min(sched.accStart[i], sched.reset[i]);
            if (fullAccrDate < minFullAccrDate)
                minFullAccrDate = fullAccrDate;
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

// register all dates, and the model remove's those that are <= valueDate
void KRibTurboTree::init(Control*) const {
    try {
        DateTime today = model->getToday();
        DateTime valueDate = model->getValueDate();

        model->addCritDates(sched.reset);
        model->addCritDates(sched.pay);
        // add accrual start dates to simplify calc logic in cases where 
        // coupon resetDate > accStartDate, so there is a guaranteed critical
        // date to add the coupon to the turbo
        model->addCritDates(sched.accStart);

        // register discounting required for coupon payments - for any payments with reset
        // dates <= today, register zero start as today.  Discounting slice is expected 
        // to be available between and including the relevant zeroUseDate and pay dates
        int i, j;
        int nbCoupons = sched.nbCoupons;
        for (i = 0; i < nbCoupons; i++) {
            if (valueDate < sched.pay[i]) {
                DateTime zeroUseDate = max(sched.reset[i], today);
                model->registerZero(zeroUseDate, sched.pay[i], inst->discount->getName());
            }
        }
        // also account for special case of any delayed payments of fully accrued coupons
        // that we calculate on the today date (ie. getValue(0) for full turbo value 
        // (never used in the fwdTurboValue calculation)
        for (i = 0; i < nbCoupons; i++) {
            if (sched.accEnd[i] <= valueDate && valueDate < sched.pay[i] )
                model->registerZero(today, sched.pay[i], inst->discount->getName());
        }

        // register discounters with model required for principal payment discounting
        // between the principal payment dates and the associated coupon reset dates 
        // (ie. the dates where the principal payments are actually added to the turbo swap)
        // for those coupons that have a principal payment
        int nbLegs = inst->legs->size();
        for (i = 0; i < nbCoupons; i++) {
            for (j = 0; j < nbLegs; j++) {
                DateTime zeroUseDate;
                if (valueDate < sched.pay[i]) {
                    zeroUseDate = max(sched.reset[i], today);
                    KRibTurbo::IndexLeg *leg = (*inst->legs)[j].get();
                    double principalPayment = leg->principalPayments[i];
                    if (!Maths::isZero(principalPayment)) {
                        model->registerZero(zeroUseDate, sched.pay[i], leg->discount->getName());
                    }
                }
            }
        }

        // no need to pre-register discounting for initial principal payments as 
        // they are applied/calculated with no discounting period, so the tree
        // getZero function returns 1 (domestic) or FX (foreign) when useDate=matDate
    }
    catch (exception& e) { 
        throw makeException(e, __FUNCTION__);
    }
}

void KRibTurboTree::initProd(void) {
    try {
        int i, j;
        int nbCoupons = sched.nbCoupons;
        DateTime today = model->getToday();
        DateTime valueDate = model->getValueDate();

        mainSlice = model->createSlice(discountCurveName);
        mainSlice->name = inst->outputName + "_mainSlice";
        *mainSlice = 0.0;
        startDEV(mainSlice);  
        
        // create one slice per coupon payment then clear the slices in the update 
        // function when they are no longer required.  The default createSlice function does not
        // allocate any size for the slice (this is automatically allocated on assignent in 
        // the calc), so there is no memory overhead in creating them all up front here
        couponSlice.resize(nbCoupons);

        for (i = 0; i < nbCoupons; i++) {
            couponSlice[i] = model->createSlice(discountCurveName);
            couponSlice[i]->name = inst->outputName + "_couponSlice_" + 
                                   Format::toString(i);
        }

        // getValue/stub slices - never DEVd
        getValueSlice = model->createSlice();
        getValueSlice->name = inst->outputName + "_getValueSlice";
        stubCpnSlice = model->createSlice();
        stubCpnSlice->name = inst->outputName + "_stubCpnSlicee";

        // variables to simplify update/calc function
        int nbLegs = inst->legs->size();
        payCpnPrincipal.resize(nbLegs);
        for (i = 0; i < nbLegs; i++)
            payCpnPrincipal[i].resize(nbCoupons, false);

        payInitCpnPrincipal.resize(nbLegs, false);

        for (i = 0; i < nbLegs; i++) {
            for (j = 0; j < nbCoupons; j++) {
                DateTime zeroUseDate;
                if (valueDate < sched.pay[j]) {
                    KRibTurbo::IndexLeg *leg = (*inst->legs)[i].get();
                    double principalPayment = leg->principalPayments[j];
                    if (!Maths::isZero(principalPayment))
                        payCpnPrincipal[i][j] = true;
                }
            }
        }
        
        for (i = 0; i < nbLegs; i++) {
            const KRibTurbo::IndexLeg *leg = (*inst->legs)[i].get();
            if (leg->payInitialPrincipal) {
                if (!Maths::isZero((*leg->notionals)[0]))
                    payInitCpnPrincipal[i] = true;
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__); 
    }
}


// apart from treeStep == 0, getValue returns the value of the forward turbo (as defined 
// by the accrual period boundary dates, to allow generic offseting of reset/payment dates).  
// A forward turbo coupon is presumed non active (ie. dropped) if the accrual end date 
// is <= (currentDate or valueDate).
// At treeStep = 0, this function returns the full leg value, NOT the value of the forward
// turbo.  The full leg value includes any future payments from non active coupons (which 
// would be dropped in a forward turbo calculation) and assumes full payment of any part
// accrued coupons (thus does not take account of leg's stub convention). 
// In summary, treeStep == 0 is a special case of the getValue function
// The pricing framework guarantees that pricing events like call/KO exercise that require a
// forward turbo value are not included on <= valueDate, so getValue(treeStep == 0) is
// guaranteed to be asking for the full turbo value, not the forward turbo value.
const TreeSlice& KRibTurboTree::getValue(int step, DateTime eventDate) const {
   try {
        DateTime currentDate = model->getDate(step);

        if (eventDate != currentDate) {
            throw ModelException("Cannot be valued at a date != currentDate");
        }

        // check trivial case
        if (currentDate <= minFullAccrDate) {
            // all coupons have been added to mainSlice by update function at this
            // stage if all coupons have fully accrued and reset
            *getValueSlice = *mainSlice;
            return *getValueSlice;
        }

        if (step == 0)
            getFullLegValue(step, eventDate, *getValueSlice);
        else 
            getForwardValue(step, eventDate, *getValueSlice);

        return *getValueSlice;
    }
    catch (exception& e) { 
        throw makeException(e, __FUNCTION__); 
    }
}

void KRibTurboTree::getFullLegValue(int step, DateTime eventDate, TreeSlice &fullLegSliceL) const {
   try {
        DateTime valueDate = model->getValueDate();
        DateTime currentDate = model->getDate(step);

        // calculate the full leg value, including any future payments
        // no matter if those payments are from non active coupons

        fullLegSliceL = *mainSlice;

        int nbCoupons = sched.nbCoupons;
        for (int cpnIdx = 0; cpnIdx < nbCoupons; cpnIdx++) {

            DateTime accStartDate = sched.accStart[cpnIdx];
            DateTime accEndDate = sched.accEnd[cpnIdx];
            DateTime resetDate = sched.reset[cpnIdx];
            DateTime payDate = sched.pay[cpnIdx];
            DateTime fullAccrDate = min(resetDate, accStartDate);

            // only include cashflows after the value date
            if (payDate <= valueDate)
                continue;

            // check if coupon has been added by update function already
            if (currentDate <= fullAccrDate && valueDate < accEndDate )
                continue;

            addCoupon(fullLegSliceL, step, cpnIdx);
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}


void KRibTurboTree::getForwardValue(int step, DateTime eventDate, TreeSlice &fwdMainSliceL) const { 
    try {
/*
        // if being called to value itself as a past value/fixing ie. before the tree
        // start date, defer to the instrument itself which can do the fully
        // deterministic model independent calculation
        // ??? speculative at this stage - no guarantee logic is sensible
        DateTime modelStartDate = model->model->getDate(0)
        if (eventDate < modelStartDate) {
            CashflowInfo cfi;
            double pastValue = inst->getValue(eventDate, cfi);
            mainSliceL = pastValue;
            return mainSlice;
        }
*/
        DateTime valueDate = model->getValueDate();
        DateTime currentDate = model->getDate(step);
        DateTime accLastDate = sched.accEnd.back();

        // simple case first: forward turbo has no value if on or after last
        // accrual date as no coupons are active in terms of the forward turbo,
        // even though there may be payments after this date
        if (accLastDate <= currentDate) {
            fwdMainSliceL = 0.0;
            return;
        }

        // take local copy of fully accrued and reset leg so far
        fwdMainSliceL = *mainSlice;

        // Calculate stub coupon and add any fully accrued coupons that have
        // not yet been added in the update function
        int nbCoupons = sched.nbCoupons;
        for (int cpnIdx = 0; cpnIdx < nbCoupons; cpnIdx++) {

            DateTime accStartDate = sched.accStart[cpnIdx];
            DateTime accEndDate = sched.accEnd[cpnIdx];
            DateTime resetDate = sched.reset[cpnIdx];
            DateTime payDate = sched.pay[cpnIdx];
            DateTime fullAccrDate = min(resetDate, accStartDate);

            // check if current coupon is active.  
            // If not active, any payments associated with this coupon are not 
            // considered part of the fwdturbo value even if they are after currentDate
            if (accEndDate <= valueDate || accEndDate <= currentDate)
                continue;

            // check if coupon has been added by update function already
            if (currentDate <= fullAccrDate && valueDate < accEndDate)
                continue;

            if (currentDate <= accStartDate) {  
                // coupon has fully accrued
                addCoupon(fwdMainSliceL, step, cpnIdx);
                continue;
            }

            // Coupon has not fully accrued.
            // This covers cases where a stub is required as the current date is
            // within the accrual period of the current active coupon
            if (inst->stubType == StubType::STUB_NOT_ALLOWED) {
                throw ModelException("Value requested on "+currentDate.toString()
                    +" implies a using a stub for coupon #"
                    +Format::toString(cpnIdx)+", but stubType == STUB_NOT_ALLOWED");
            }

            if (inst->stubType == StubType::STUB_NONE) {
                // STUB_NONE stub convention pays the full coupon with the cashflow date being
                // the coupon payment date
                addCoupon(fwdMainSliceL, step, cpnIdx);
                continue;
            }

            // for other stub types, will have to implement a calculatePartialCoupon function
            // as the day count conventions for each indexLeg may be different so can't just 
            // scale the full coupon value
            // ??? need to sort out what exercise principal payments mean in these cases
            throw ModelException("stubTypes except STUB_NONE are not yet supported - "
                                "getValue date = " + currentDate.toString() + 
                                ", coupon index = " + Format::toString(cpnIdx) + 
                                ", coupon accrual startDate = " + accStartDate.toString() + 
                                ", coupon accrual endDate = " + accEndDate.toString());
        }
        // ??? when eventDate != currentDate, will have to discount this initial principal
        // payment between event and currentDate - todo.  Put additional check here to force 
        // addition of code to cover this case and thus ensure this case won't be forgotten about
        if (currentDate != eventDate)
            throw ModelException("Adding initial Prinicpal payment at exercise logic is only "
                "implemented for the case where the eventDate " + eventDate.toString() + 
                " equals the current getValue date " + currentDate.toString() +
                ".  Additional code required to request a forward starting forward turbo value");

        // finally any initialPrincipalPayments for each leg to the forward value on the 
        // current (exercise) date if not already included in the forward value (ie. exercise date < turbo
        // start date).  Initial prinicpal assumed to pay on the exercise date, so the
        // payment amount added directly to forward value
        int nbLegs = inst->legs->size();
        DateTime legStartDate = sched.accStart[0];
        if (currentDate > legStartDate) {
            int i;
            for (i = 0; i < nbLegs; i++) {
                const KRibTurbo::IndexLeg &leg = *(*inst->legs)[i];
                // only case to worry about here is if currentDate > accStartDate[0], as logic
                // above has already added intial principal payment of currentDate <= accStartDate[0]
                if (payInitCpnPrincipal[i] && leg.payInitialPrincipalAtExercise) {
                    TreeSliceSP stubZeroSlice;
                    // use the zero to elegantly account for currency of principal
                    model->getZero(currentDate, currentDate, leg.discount->getName(), stubZeroSlice);
                    double principalPayment = (*leg.notionals)[0];
                    fwdMainSliceL -= principalPayment * *stubZeroSlice;
                }
            }
        }
    }
    catch (exception& e) { 
        throw makeException(e, __FUNCTION__); 
    }
}

// calculates individual undiscounted turbo coupon, including cap/floor modification
void KRibTurboTree::calculateFullCoupon(int step, int cpnIdx, DateTime resetDate,
                                        TreeSlice& couponSliceL) const {
   try {
        DateTime currentDate = model->getDate(step);
        int i;

        couponSliceL = 0.0;
        // calculate coupon payment as sum of the turbo indexlegs
        int nbLegs = inst->legs->size();
        for (i = 0; i < nbLegs; i++) {
            const KRibTurbo::IndexLeg &leg = *(*inst->legs)[i];
            double spread   = (*leg.spreads)[cpnIdx];
            double dcf      = leg.dcfs[cpnIdx];
            double notional = (*leg.notionals)[cpnIdx];
            double weight   = (*leg.weights)[cpnIdx];

            const TreeSlice &fx = legFXprods[i]->getValue(step, resetDate);

            if (Maths::isZero(weight)) {
                if (leg.rateType == RateType::SIMPLE) {
                    couponSliceL += (notional * dcf * spread) * fx;
                }
                else {
                    couponSliceL += notional * (::pow (1. + spread, dcf) - 1.) * fx;
                }
            }
            else {
                const TreeSlice &indexSlice =   
                    legIndexProds[i][cpnIdx]->getValue(step, resetDate);

                if (leg.rateType == RateType::SIMPLE) {
                    couponSliceL += (notional * dcf * weight) * (indexSlice + spread/weight) * 
                                    fx;
                } 
                else {
                    couponSliceL += notional * (pow (1. + (weight * indexSlice + spread), dcf) - 1.) * 
                                    fx;
                }
            }
        }
        // apply any provided cap and floor to coupon.  Cap and floor are applied on 
        // coupon payment value rather than rates as the turbo boths adds indexLegs of different
        // currencies together in the coupon calculation, and some of the rates may be compound
        // This is OK as we apply the cap/floor before discounting the coupon payment
        double cap;
        double floor;
        double capFloorNotional;
        double capFloorDcf;
        double pmt;

        if (inst->capFloorNotionals.get()) {
            capFloorNotional = (*inst->capFloorNotionals)[cpnIdx];
            capFloorDcf = inst->capFloorDcfs[cpnIdx];

            if (inst->capRates.get()) {  // cap requested
                cap = (*inst->capRates)[cpnIdx];
                pmt = cap * capFloorNotional * capFloorDcf;
                // if notional -ve, cap amount becomes floor amount
                if (capFloorNotional < 0.0)
                    couponSliceL = smax(couponSliceL, pmt);
                else
                    couponSliceL = smin(couponSliceL, pmt);
            }
            if (inst->floorRates.get()) {  // floor requested
                floor = (*inst->floorRates)[cpnIdx];
                pmt = floor * capFloorNotional * capFloorDcf;
                if (capFloorNotional < 0.0)
                    couponSliceL = smin(couponSliceL, pmt);
                else
                    couponSliceL = smax(couponSliceL, pmt);
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/*  ??? something to think about
void KRibTurboTree::calculateFullBinaryCoupon(int step, int cpnIdx, DateTime resetDate,
                                              TreeSlice& couponSliceL) const {

}*/


void KRibTurboTree::addCoupon(TreeSlice &accuSlice, int step, int cpnIdx) const {
    try {
        DateTime valueDate = model->getValueDate();
        DateTime currentDate = model->getDate(step);
        DateTime accStartDate = sched.accStart[cpnIdx];
        DateTime resetDate = sched.reset[cpnIdx];
        DateTime payDate = sched.pay[cpnIdx];

        if (currentDate <= resetDate) {
            // The value of the coupon is available 
            // (computed by update, even if it is not added to mainSlice yet)
            if (!ribProd) {
                accuSlice += *couponSlice[cpnIdx];
            }
            else {
                TreeSliceSP ribCoupon;
                ribProd->getRibCoupon(step, cpnIdx, ribCoupon);
                accuSlice += (*ribCoupon);
            }
        }
        else {
            // reset date is before current date.
            // To calculate this, we need to use either a state variable or to
            // approximate the unknown reset rate by using the deterministic ratio method.
            // The model decides which method it is going to do - we simply call
            // the getValue function of the underlying with the eventDate = resetDate
            TreeSlice &stubCpnSliceL = *stubCpnSlice;
            calculateFullCoupon(step, cpnIdx, resetDate, stubCpnSliceL);

            TreeSliceSP stubZeroSlice;
            if (!ribProd) {
                model->getZero(currentDate, payDate, inst->discount->getName(), stubZeroSlice);
            }
            else {
                ribProd->getRibCoupon(step, cpnIdx, stubZeroSlice);
            }
            accuSlice += (stubCpnSliceL * (*stubZeroSlice));
            stubCpnSliceL.clear();
        }

        // add any principal payments
        int i;
        int nbLegs = inst->legs->size();
        for (i = 0; i < nbLegs; i++) {
            if (payCpnPrincipal[i][cpnIdx]) {
                TreeSliceSP stubZeroSlice;
                const KRibTurbo::IndexLeg &leg = *(*inst->legs)[i];
                model->getZero(currentDate, payDate, leg.discount->getName(), stubZeroSlice);
                double principalPayment = leg.principalPayments[cpnIdx];
                accuSlice += principalPayment * (*stubZeroSlice);
            }
        }

        // also have to account for initial principal payment
        for (i = 0; i < nbLegs; i++) {
            if (cpnIdx == 0 &&
                payInitCpnPrincipal[i] &&
                valueDate < accStartDate &&
                currentDate <= accStartDate) {

                TreeSliceSP stubZeroSlice;
                const KRibTurbo::IndexLeg &leg = *(*inst->legs)[i];
                model->getZero(currentDate, currentDate, leg.discount->getName(), stubZeroSlice);
            
                accuSlice -= (*leg.notionals)[0] * (*stubZeroSlice);
            }
        }
    }
    catch (exception& e) { 
        throw makeException(e, __FUNCTION__); 
    }
}


// Function calculates and adds any fully accrued coupons to the mainSlice on
// the full accrual date, if that date is > valueDate.
// which is defined as min(accStartDate, resetDate) 
void KRibTurboTree::update(int &step, FDProduct::UpdateType updateType) {
    try {
        DateTime currentDate = model->getDate(step);
        DateTime valueDate = model->getValueDate();
        int cpnIdx;
        int nbCoupons = sched.nbCoupons;
        int nbLegs = inst->legs->size();
        // note: product currency can be different from model pricing currency
        string discountCurve = inst->discount->getName();  

        // active coupons are calculated on the reset date - code supports the case
        // of multiple coupon resets on the same reset date
        for (cpnIdx = 0; cpnIdx < nbCoupons; cpnIdx++) {

            DateTime resetDate = sched.reset[cpnIdx];
            DateTime accStartDate = sched.accStart[cpnIdx];
            DateTime accEndDate = sched.accEnd[cpnIdx];
            DateTime payDate = sched.pay[cpnIdx];
            TreeSlice& couponSliceL = *couponSlice[cpnIdx];
            DateTime fullAccrDate = min(resetDate, accStartDate);

            // "set couponSlice[]": calculate full coupon
            if (currentDate == resetDate) {
                TreeSliceSP zeroSlice;
                // if resetDate == payDate, getZero returns discount slice of const value 1.0
                model->getZero(currentDate, payDate, discountCurve, zeroSlice);
                calculateFullCoupon(step, cpnIdx, currentDate, couponSliceL);
                if (ribProd.get()) {
                    ribProd->couponHasReset(couponSlice[cpnIdx], cpnIdx);
                }
                couponSliceL *= (*zeroSlice);  // discount payment
                // must DEV coupon slice between resetDate and fullAccrDate
                startDEV(couponSlice[cpnIdx]);
            }

            if (accEndDate <= valueDate)  {
                // coupon not active - has fully accrued, but may or may not have paid yet.
                // Adding fully accrued but future payment coupons is handled in getValue
                continue;
            }

            if (payDate <= valueDate) {
                // ??? catch this obscure case for now - should probably support in future.
                // The payment date is before the valueDate and is thus dropped, but the 
                // coupon is active as the accrual end date is >= valueDate.  This may occur in
                // a modified following situation on the payment date, but error for now if it 
                // crosses the valueDate.  Will still price OK if doesn't cross the valueDate
                // so should either support fully or not at all
                throw ModelException("Coupon index no. " + Format::toString(cpnIdx) + " pays (" +
                                     payDate.toString() + ") <= the value date " +
                                     valueDate.toString() + " but the accrual end Date " + 
                                     accEndDate.toString() + " is after the value date, which "
                                     "means the coupon is active.  This coupon configuration "
                                     "is not supported at this stage if it crosses the valueDate");
            }

            if (currentDate != fullAccrDate)
                continue;

            // set mainSlice += principal, coupon
            addCoupon(*mainSlice, step, cpnIdx);

            // [begin] delete slices that are no more needed

            // couponSlice
            couponSlice[cpnIdx]->clear();
            stopDEV(couponSlice[cpnIdx]);

            // principalPaymentSlice
  //          principalPaymentSlice[cpnIdx]->clear();
  //          stopDEV(principalPaymentSlice[cpnIdx]);

            // initialPrincPaymentSlice
 //           initialPrincPaymentSlice->clear();
 //           stopDEV(initialPrincPaymentSlice);

            // rib coupon
            if (ribProd.get())
                ribProd->deleteCoupon(cpnIdx);

            // [end] delete slices that are no more needed
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KRibTurboTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        string ccy;
        if (inst->discount.get()) 
            ccy = inst->discount->getCcy();

        recordSliceToOutputName(ctrl, results, model, inst->isTopComponent(), 
                                inst->outputName, ccy, getValue(0,model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/*****************************************************************************************/

static void checkArraySize(const string &name, int arraySize, int nbCoupons) {
    if (arraySize != nbCoupons) {
        throw ModelException(name + " array size (" + Format::toString(arraySize) +
                             ") should equal number of coupons (" + 
                             Format::toString(nbCoupons) + ")");
    }
}

/************************************ KRibTurbo::IndexLeg *********************************/

void KRibTurbo::IndexLeg::validatePop2Object(void) {
    try {
        int nbCoupons = notionals->size();

        if (nbCoupons<1) 
            throw ModelException("At least one coupon must be specified in instrument");

        if (!weights.get())
            weights.reset(new DoubleArray(nbCoupons, 0.0));  // create dummy spread value 0.0
        checkArraySize("weights", weights->size(), nbCoupons);

        if (!spreads.get()) 
            spreads.reset(new DoubleArray(nbCoupons, 0.0));  // create dummy spread value 0.0
        checkArraySize("spreads", spreads->size(), nbCoupons);

        int nbIndexes = indexes.size();
        if (nbIndexes == 0) {
            for (int i=0; i<nbCoupons; ++i) {
                if (!Maths::isZero((*weights)[i])) 
                    throw ModelException("No index is supplied but the weight (for coupon " +
                                         Format::toString(i) + ") is " +
                                         Format::toString((*weights)[i])+" !=0");
            }
        } 
        else if (nbIndexes != 1 && nbIndexes != nbCoupons) 
            throw ModelException("If index supplied, there must be either 1 only (assumed to apply "
                                 "to all reset dates), or 1 index per reset date.  Supplied " +
                                 Format::toString(nbIndexes) + " indexes but leg has " +
                                 Format::toString(nbCoupons) + " reset dates/coupons");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void KRibTurbo::IndexLeg::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(KRibTurbo::IndexLeg, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(notionals, "Notional schedule for each coupon");
    FIELD(spreads, "Fixed spread for each coupon (as decimal)"); 
    FIELD_MAKE_OPTIONAL(spreads);
    FIELD(weights, "Weight for each coupon to apply to index");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(discount, "");
    FIELD(payInitialPrincipal, "Pay initial principal at accrual start date of first coupon");
    FIELD_MAKE_OPTIONAL(payInitialPrincipal);
    FIELD(payPrincipal, "Pay amortization payments and final notional as principal payment");
    FIELD_MAKE_OPTIONAL(payPrincipal);
    FIELD(payInitialPrincipalAtExercise, "At exercise/getValue requested date of leg, include "
          "any initial principal payments if enabled by the payInitialPrincipal flag as part of "
          "the forward value of the leg.  For example, if payInitialPrincipal is true and "
          "payPrincipalAtExercise is true, the initial prinicpal payment will be included in "
          "the forward leg value and will be assumed to pay on the exercise date (ie. moved "
          "forward from the leg start date), even if the exercise/getValue date is > start date "
          "of the leg.");
    FIELD_MAKE_OPTIONAL(payInitialPrincipalAtExercise);
            
    FIELD(rateType, "Rate type");
    FIELD(dcc, "day count convention for leg");
    FIELD_MAKE_OPTIONAL(dcc);
    FIELD(dcfs, ""); 
    FIELD_MAKE_OPTIONAL(dcfs);
    FIELD(indexes, "reset index/underlying component - supply either 1 to apply for all reset dates"
          " or 1 component per coupon"); 
    FIELD_MAKE_OPTIONAL(indexes);
    FIELD(fx, "IndexSpecFX to be used to convert this leg to the turbo payment ccy"); 
    FIELD_MAKE_OPTIONAL(fx);

    // transient fields to store and copy calculated values
    FIELD(principalPayments,""); 
    FIELD_MAKE_TRANSIENT(principalPayments);
}


CClassConstSP const KRibTurbo::IndexLeg::TYPE = CClass::registerClassLoadMethod(
    "KRibTurbo::IndexLeg", typeid(KRibTurbo::IndexLeg), KRibTurbo::IndexLeg::load);


typedef KRibTurbo::IndexLegArray KRibTurboIndexLegArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("KRibTurbo::IndexLegArray", KRibTurboIndexLegArray);


/************************************ KRibTurbo *********************************/

FDProductSP KRibTurbo::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KRibTurboTree(KRibTurboConstSP(this), model));
}

// called via KComponent in the getMarket function
void KRibTurbo::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        int nbCoupons = sched->nbCoupons;
        int i;

        sched->setup(model, market);

        if (sched->reset.size() == 0)
            throw ModelException("Product (outputName = " + outputName +
                                 " requires reset date even if all legs are fixed rate");

        if (capFloorNotionals.get()) {
            checkArraySize("capFloorNotionals", capFloorNotionals->size(), nbCoupons);
            if (floorRates.get()) 
                checkArraySize("floorRates", floorRates->size(), nbCoupons);
            if (capRates.get()) 
                checkArraySize("capRates", capRates->size(), nbCoupons);
            if (capFloorDcc.get()) {
                if (capFloorDcfs.size())
                    throw ModelException("capFloorDcc and capFloorDcfs must not be provided at the same time");
                capFloorDcfs.resize(nbCoupons);
                sched->calcDcfs(*capFloorDcc, capFloorDcfs);
            }
            else {
                if (!capFloorDcfs.size())
                    throw ModelException("capFloorDcc or capFloorDcfs required if capFloorNotionals provided");
                if (capFloorDcfs.size()!=nbCoupons)
                    throw ModelException("capFloorDcfs must have nbCoupons ("+Format::toString(nbCoupons)+") entries");
            }
            if (capRates.get() && floorRates.get()) {
                for (i = 0; i < capRates->size(); i++) {
                    if ((*capRates)[i] <= (*floorRates)[i])
                        throw ModelException("Supplied capRate " + Format::toString((*capRates)[i]) + 
                            " for coupon index " + Format::toString(i) + " must be greater than " 
                            "supplied floor rate " + Format::toString((*floorRates)[i]));
                }
            }   
        } 
        else {
            if (floorRates.get())
                throw ModelException("If capFloorNotionals not provided, must not provide floorRates");
            if (capRates.get())
                throw ModelException("If capFloorNotionals not provided, must not provide capRates");
            if (capFloorDcc.get())
                throw ModelException("If capFloorNotionals not provided, must not provide capFloorDcc");
            if (capFloorDcfs.size())
                throw ModelException("If capFloorNotionals not provided, must not provide capFloorDcfs");
        }

        if (legs->size()<1) 
            throw ModelException("Must provide at least one leg in turbo for pricing");

        // legs checking
        int nbLegs = legs->size();
        for (int legI=0; legI<nbLegs; ++legI) {
            try { 
                IndexLeg *leg = (*legs)[legI].get();
                if (!leg)
                    throw ModelException("This leg is not a valid object - NULL pointer");

                { // check array sizes
                    int nbNotionals = leg->notionals->size();
                    if (nbNotionals != nbCoupons) {
                        throw ModelException("Number of coupon dates "
                        + Format::toString(nbCoupons)
                        + " does not equal number of notionals: " 
                        + Format::toString(nbNotionals));
                    }
                    // spreads and weights already checked in validatePop2Object
                }

                // calculate the day count fractions if needed
                if (leg->dcc.get()) {
                    if (leg->dcfs.size())
                        throw ModelException("dcc and dcfs must not be provided at the same time");
                    sched->calcDcfs(*leg->dcc, leg->dcfs);
                }
                else {
                    if (!leg->dcfs.size())
                        throw ModelException("dcc or dcfs required");
                    if (leg->dcfs.size()!=nbCoupons)
                        throw ModelException("dcfs must have nbCoupons ("+Format::toString(nbCoupons)+") entries");
                }

                // default currency is domestic
                if (leg->discount.getName().empty()) {
                    leg->discount.setName(discount.getName());
                }

                // creates fx KIndexSpecs based on the currencies of the legs
                if (!leg->fx) {
                    leg->fx.reset(new IndexSpecFX(leg->discount->getName(), discount->getName()));
                }
                leg->fx->addResetDates(sched->reset);
                leg->fx->setup(model, market);
                if (leg->fx->sameCurve) {
                    // cannot be wrong because it means that fx has been constructed implicitly
                }
                else {
                    // check that the correct factor was provided
                    string fxRiskName = leg->fx->factor->getRiskCcy()->getName();
                    string fxBaseName = leg->fx->factor->getBaseCcy()->getName();
                    if (fxRiskName != leg->discount->getName()) {
                        throw ModelException("fx->factor \""+leg->fx->factor->getName()
                            +"\" has risk curve \""
                            +fxRiskName
                            +"\" when it should be "
                            +leg->discount->getName());
                    }
                    if (fxBaseName != discount->getName()) {
                        throw ModelException("fx->factor \""+leg->fx->factor->getName()
                            +"\" has base curve \""
                            +fxBaseName
                            +" when it should be "
                            +discount->getName());
                    }
                }
                // send reset dates to the underlying indexes
                if (leg->indexes.size() == 1) {
                    DateTimeArray fixingDates;
                    FlexDates::filterOnWeight(sched->reset, *leg->weights, fixingDates);
                    leg->indexes[0]->addResetDates(fixingDates);
                    leg->indexes[0]->setup(model, market);
                }
                else if (leg->indexes.size() == nbCoupons) {

                    DateTimeArray fixingDates(1);
                    for (i = 0; i < leg->indexes.size(); i++) {
                        if (Maths::isZero((*leg->weights)[i]))
                            continue;
                        fixingDates[0] = sched->reset[i];
                        leg->indexes[i]->addResetDates(fixingDates);
                        leg->indexes[i]->setup(model, market);
                    }
                }

                // precalculate principal payments to simplify the update function
                leg->principalPayments.resize(nbCoupons, 0.0);

                if (leg->payPrincipal) {
                    for (i = 0; i < nbCoupons-1; i++) {
                        // for each coupon
                        double principalPmt = (*leg->notionals)[i] - (*leg->notionals)[i+1];
                        leg->principalPayments[i] = principalPmt;
                    }
                    // add final payment
                    double lastPmt = (*leg->notionals).back();
                    leg->principalPayments[nbCoupons-1] = lastPmt;
                }
            }
            catch (exception& e) {
                throw ModelException(e, "In leg["+Format::toString(legI)+"]");
            }
        }
        // RIB
        if (rib.get()) {
            rib->discountCurveName = getActualDiscount()->getName();
            rib->accStart = sched->accStart;
            rib->accEnd = sched->accEnd;
            rib->resetDates = sched->reset;
            rib->payDates = sched->pay;
            rib->setup(model, market);
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

DateTime KRibTurbo::getLastDate(void) const {
    return sched->pay.back();
}

void KRibTurbo::reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const {
    try {
        int nbLegs = legs->size();
        int nbCoupons = sched->nbCoupons;

        // principal (initial & coupon) payments
        for (int cpnIdx=-1; cpnIdx < nbCoupons; ++cpnIdx) { // -1 for initial principal payment
            CashflowInfoSP cfi(new CashflowInfo);
            cfi->componentName = outputName;
            cfi->cfType = CashflowInfo::CfType::PRINCIPAL;
            cfi->date = (cpnIdx==-1 ? sched->accStart[0]
                                    : sched->pay[cpnIdx]);
            for (int legIdx=0; legIdx < nbLegs; ++legIdx) {
                IndexLeg &leg = *(*legs)[legIdx];

                if (cpnIdx==-1 && !leg.payInitialPrincipal)
                    continue;

                if (cpnIdx!=-1 && !leg.payPrincipal) 
                    continue;

                if (!amountsNeeded)
                    continue;

                double notional = (cpnIdx==-1 ? -(*leg.notionals)[0]
                                              : leg.principalPayments[cpnIdx]);
                double fx = leg.fx->getValue(cfi->date, *cfi);

                cfi->amount += notional * fx;
                cfi->push("leg["+Format::toString(legIdx)+"]", notional);
            }
            if (!Maths::isZero(cfi->amount))
                cashflowInfos.push_back(cfi);
        }

        // coupon payments
        for (int cpnIdx = 0; cpnIdx < nbCoupons; cpnIdx++) {
            CashflowInfoSP cpnCfi(new CashflowInfo);
            cashflowInfos.push_back(cpnCfi);

            cpnCfi->componentName = outputName;
            cpnCfi->cfType = CashflowInfo::CfType::COUPON;
            cpnCfi->date = sched->pay[cpnIdx];

            if (!amountsNeeded)
                continue;

            double ribFraction = (!rib ? 1.
                                       : rib->getRibFraction(cpnIdx, *cpnCfi));

            if (cpnCfi->amountType == CashflowInfo::AmountType::UNKNOWN)
                continue;

            cpnCfi->push("accStart", sched->accStart[cpnIdx]);
            cpnCfi->push("accEnd", sched->accEnd[cpnIdx]);

            for (int legIdx=0; legIdx < nbLegs; ++legIdx) {
                IndexLeg &leg = *(*legs)[legIdx];
                CashflowInfo cfiLoc;

                double fx       = leg.fx->getValue(cpnCfi->date, cfiLoc);
                double notional = (*leg.notionals)[cpnIdx];
                double dcf      = leg.dcfs[cpnIdx];
                double weight   = (*leg.weights)[cpnIdx];
                double spread   = (*leg.spreads)[cpnIdx];
                double reset    = 0.;

                if (!Maths::isZero(weight)) {
                    int indexesIdx  = (leg.indexes.size()==1 ? 0 : cpnIdx);
                    reset = leg.indexes[indexesIdx]->getValue(sched->reset[cpnIdx], cfiLoc);
                }

                double tmp;
                if (leg.rateType == RateType::SIMPLE)
                     tmp = (reset * weight + spread)  * dcf;
                else tmp = ::pow(1 + (reset * weight + spread), dcf) - 1.;

                cfiLoc.push("dcf", dcf);
                cfiLoc.push("dcfFrac", sched->calcDcfFrac(cpnIdx, *leg.dcc));
                cfiLoc.push("notional", notional);
                cfiLoc.push("weight", weight);
                cfiLoc.push("spread", spread);
                cfiLoc.push("rateType", RateTypeBoxedEnum::toString(leg.rateType));

                cpnCfi->merge(cfiLoc, "leg["+Format::toString(legIdx)+"]");
                cpnCfi->amount += notional  * ribFraction * fx * tmp;
            }
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KRibTurbo::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(KRibTurbo, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(legs, ""); 

    FIELD(sched, "");
    FIELD(rib, "");
    FIELD_MAKE_OPTIONAL(rib);

    FIELD(capFloorNotionals, ""); 
    FIELD_MAKE_OPTIONAL(capFloorNotionals);
    FIELD(floorRates, ""); 
    FIELD_MAKE_OPTIONAL(floorRates);
    FIELD(capRates, ""); 
    FIELD_MAKE_OPTIONAL(capRates);
    FIELD(capFloorDcc, ""); 
    FIELD_MAKE_OPTIONAL(capFloorDcc);
    FIELD(capFloorDcfs, ""); 
    FIELD_MAKE_OPTIONAL(capFloorDcfs);
    FIELD(stubType, "Stub convention when exercised/knocked out "
                 "etc during accrual period");
 }

CClassConstSP const KRibTurbo::TYPE = CClass::registerClassLoadMethod(
    "KRibTurbo", typeid(KRibTurbo), KRibTurbo::load);

DEFINE_TEMPLATE_TYPE(KRibTurboArray);

bool KRibTurboLoad() { 
    return KRibTurbo::TYPE != 0; 
}

DRLIB_END_NAMESPACE
