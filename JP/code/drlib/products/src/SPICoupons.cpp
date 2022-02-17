//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPICoupons.cpp
//
//   Description : Coupons interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPICoupons.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/SVPath.hpp"

DRLIB_BEGIN_NAMESPACE

/*****************************************************************************/
// this is the external interface for Coupons so that users can bolt in 
// any Coupons type they want - note this includes the SPICouponsWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ICouponsSPI as soon as possible
void ICouponsSPIInterface::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ICouponsSPIInterface, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}
CClassConstSP const ICouponsSPIInterface::TYPE = CClass::registerInterfaceLoadMethod(
    "ICouponsSPIInterface", typeid(ICouponsSPIInterface), ICouponsSPIInterface::load);

/*****************************************************************************/

ICouponsSPI::~ICouponsSPI() {}

SPICouponsNone::SPICouponsNone(): SPIInterfaceIMS(TYPE) {} // for reflection

ICouponsSPISP SPICouponsNone::getCouponsSPI() {
    static const string routine = "SPICouponsNone::getCouponsSPI";
    // essentially just returns itself
    ICouponsSPISP theCoupons = ICouponsSPISP(this, NullDeleter()); // FIXME can be a problem
    return theCoupons;
}

bool SPICouponsNone::doesNothing() const {
    return true;
}

void SPICouponsNone::init(const DateTimeArray& rebalDates) {}

double SPICouponsNone::getCoupon(double B,
                         double BF,
                         int iStep) const {
    return 0.0;
}

double SPICouponsNone::getCouponPVFactor(int               iStep,
                                 const YieldCurve* disc) const {
    return 0.0;
}

double SPICouponsNone::getCouponPVFactor(int               iStep,
                                         SVDiscFactorSP dfSV) const {
    return 0.0;
}

const DateTime& SPICouponsNone::getCouponPayDate(int iStep) const {
    throw ModelException("SPICouponsNone::getCouponPayDate",
                         "Internal error!");
}

bool SPICouponsNone::isCouponDate(const DateTime& myDate) const {
    return false;
}

bool SPICouponsNone::isCouponStep(int iStep) const {
    return false;
}

bool SPICouponsNone::isThresholdBreached(double level, int iStep) const {
    return false;
}

string SPICouponsNone::getThresholdType() const {
    return SPI_THRESHOLD_TYPE_NONE;
}

bool SPICouponsNone::guaranteedCoupons() const { return false; }

double SPICouponsNone::getGuaranteedCoupon(int iStep) const {
    return 0.0;
}

bool SPICouponsNone::onlyFixedCoupons() const { return true; }

bool SPICouponsNone::hasEarlyRedemption() const {
    return false;
}

double SPICouponsNone::getTarget() const {
    return 0.0;
}

bool SPICouponsNone::isTargetMet(int     iStep,
                                 double&  couponAmt,
                                 double&  sumCouponsSoFar,
                                 double&  bonusCoupon) const {
    return false;
}

DateTimeArray SPICouponsNone::getPaymentDates() const {
    DateTimeArray dummyDates(0);
    return dummyDates;
}

DateTimeArray SPICouponsNone::getEssentialDates() const {
    DateTimeArray dummyDates(0);
    return dummyDates;
}

IObject* SPICouponsNone::defaultSPICouponsNone(){
    return new SPICouponsNone();
}

// post getMarket validation (handles settlement of coupons)
void SPICouponsNone::Validate(const InstrumentSettlement* instSettle) {
    // do nothing
}


/** Invoked when Class is 'loaded' */
void SPICouponsNone::load(CClassSP& clazz){
    REGISTER(SPICouponsNone, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    IMPLEMENTS(ICouponsSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPICouponsNone);
    FIELD(dummyString, "dummyString");
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(dummyString);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPICouponsNone::TYPE = CClass::registerClassLoadMethod(
    "SPICouponsNone", typeid(SPICouponsNone), SPICouponsNone::load);

/*****************************************************************************/

SPICouponsStd::SPICouponsStd(): SPIInterfaceIMS(TYPE), minCoupon(0), fixedCoupon(0),
					thresholdType(SPI_THRESHOLD_TYPE_BASKET), couponStrike(1.0),
					struckAtPrevious(false), struckEx(false),
                    allowEarlyRedemption(false), targetLevel(0.0), bonusCoupons(0),
                    hasCouponCap(false), couponCap(0.0), currentStrike(1.0) {} // for reflection

// validation
void SPICouponsStd::validatePop2Object(){
    static const string routine = "SPICouponsStd::validatePop2Object";

    int len = couponDates.size();
    if (len<1) {
        isOK = false;
        err = "Require at least one coupon date!";
        return;
    }
    // min coupons only appear in the Super SPI interface so could be empty
    // if we came in via SPI (set them all to zero)
    if (minCoupon.empty()) {
        minCoupon = DoubleArray(len, 0.0);
    }
    if (minCoupon.size() != len) {
        isOK = false;
        err = "Min Coupon array length (=" + Format::toString(minCoupon.size()) +
            ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
        return;
    }
    if (maxCoupon.size() != len) {
        isOK = false;
        err = "Max Coupon array length (=" + Format::toString(maxCoupon.size()) +
            ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
        return;
    }
    //fixed coupons were added later so could be empty (set them all to zero)
    if (fixedCoupon.empty()) {
        fixedCoupon = DoubleArray(len, 0.0);
    }
    else if (fixedCoupon.size() != len) {
        isOK = false;
        err = "Fixed Coupon array length (=" + Format::toString(minCoupon.size()) +
            ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
        return;
    }
    if (participations.size() != len) {
        isOK = false;
        err = "Participations array length (=" + Format::toString(participations.size()) +
            ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
        return;
    }
    if (threshold.size() != len) {
        isOK = false;
        err = "Threshold array length (=" + Format::toString(threshold.size()) +
            ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
        return;
    }
    if (!(thresholdType==SPI_THRESHOLD_TYPE_BASKET
        || thresholdType==SPI_THRESHOLD_TYPE_TE)) {
        isOK = false;
        err = "Unrecognised coupon threshold type " + thresholdType +
               ". Expected " + SPI_THRESHOLD_TYPE_BASKET + ", " + SPI_THRESHOLD_TYPE_TE;
    }
    if (Maths::isNegative(couponStrike)) {
        isOK = false;
        err = "Coupon strike (= " + Format::toString(couponStrike) +
            ") must be non-negative";
        return;
    }
    if (allowEarlyRedemption) {
        // min coupons only appear in the Super SPI interface so could be empty
        // if we came in via SPI (set them all to zero)
        if (bonusCoupons.empty()) {
            bonusCoupons = DoubleArray(len, 0.0);
        }
        else if (bonusCoupons.size() != len) {
            isOK = false;
            err = "Bonus coupon array length (=" + Format::toString(bonusCoupons.size()) +
                ") must equal #Coupon Dates (=" + Format::toString(len) + ")";
            return;
        }
        // Meaningless to have targetLevel > Sum(maxCoupon) - never reach
        // or targetLevel < Sum(minCoupon) - always reach
        double sumMaxCoupon = 0.0;
        double sumMinCoupon = 0.0;
        for(int j=0; j<maxCoupon.size(); j++) {
            sumMaxCoupon += maxCoupon[j];
            sumMinCoupon += minCoupon[j];
        }
        if (Maths::isPositive(targetLevel - sumMaxCoupon)) {
            isOK = false;
            err = "Target level (" +
                Format::toString(targetLevel) +
                ") must not exceed the sum of max coupons (" +
                Format::toString(sumMaxCoupon) + ")";
            return;
        }
        if (!Maths::isPositive(targetLevel - sumMinCoupon)) {
            isOK = false;
            err = "Target level (" +
                Format::toString(targetLevel) +
                ") must be greater than the sum of min coupons (" +
                Format::toString(sumMinCoupon) + ")";
            return;
        }
        if (hasCouponCap && Maths::isNegative(couponCap - targetLevel)) {
            isOK = false;
            err = "Total coupon cap (" +
                Format::toString(couponCap) +
                ") cannot be less than the target level (" +
                Format::toString(targetLevel) + ")";
            return;
        }
    }
    for(int i=0; i<minCoupon.size(); i++) {
        // what would negative floor mean?
        if (Maths::isNegative(minCoupon[i])) {
            isOK = false;
            err = "Min coupon #" + Format::toString(i+1) +
                "=" + Format::toString(minCoupon[i]) +
                " is negative!";
        }
        if (Maths::isPositive(minCoupon[i]-maxCoupon[i])) {
            isOK = false;
            err = "Min coupon #" + Format::toString(i+1) +
                "=" + Format::toString(minCoupon[i]) +
                " >= Max coupon " + Format::toString(maxCoupon[i]);
        }
    }
}

// post getMarket validation (handles settlement of coupons)
void SPICouponsStd::Validate(const InstrumentSettlement* instSettle) {
    // sort out the settlement

    //first get the holidays from the settlement object
    HolidayCollectorSP holidayVisitor = HolidayCollectorSP(
            new HolidayCollector());
    instSettle->accept(holidayVisitor.get());
    HolidayConstSP settleHols = !holidayVisitor->getHoliday() ?
                                HolidayConstSP(Holiday::weekendsOnly()) :
                                holidayVisitor->getHoliday();

    // now build a settlement object based on these holidays
    CashSettlePeriod cpnSettle(settlePeriod, settleHols.get());

    // then adjust the coupon dates
    couponPayDates = DateTimeArray(couponDates.size());
    for(int j=0; j<couponPayDates.size(); j++) {
        couponPayDates[j] = cpnSettle.settles(couponDates[j],0/*asset*/);
    }
}

ICouponsSPISP SPICouponsStd::getCouponsSPI() {
    static const string routine = "SPICouponsStd::getCouponsSPI";
    // essentially just returns itself
    ICouponsSPISP theCoupons = ICouponsSPISP(this, NullDeleter()); // FIXME can be a problem
    return theCoupons;
}

bool SPICouponsStd::doesNothing() const {
    return false;
}

void SPICouponsStd::init(const DateTimeArray& rebalDates) {
    // create a map which gives iCoupon from iStep
    iCouponMap = IntArray(rebalDates.size(), -1); // -1 flags not a coupon date
    int iCoupon = 0;
    for (int iDate=0; iDate<rebalDates.size() && iCoupon<couponDates.size(); iDate++) {
        if (couponDates[iCoupon]==rebalDates[iDate]) {
            iCouponMap[iDate] = iCoupon;
            iCoupon++;
        }
    }
    if (iCoupon<couponDates.size()) {
        throw ModelException("SPICouponsStd::init",
                             "Supplied coupon date ("+
                             couponDates[iCoupon].toString()+")"
                             " not found in list of all dates");
    }
}

bool SPICouponsStd::guaranteedCoupons() const {
    // any coupons with 0 threshold and >0 minCoupon?
    for(int i=0; i<minCoupon.size(); i++) {
        if (Maths::isPositive(minCoupon[i]) &&
            Maths::isZero(threshold[i])) {
            return true;
        }
    }
    return false;
}

double SPICouponsStd::getGuaranteedCoupon(int iStep) const {
    int iCoupon = iCouponMap[iStep];
    if (iCoupon>=0 && Maths::isZero(threshold[iCoupon])) {
        return minCoupon[iCoupon];
    }
    return 0.0;
}

bool SPICouponsStd::onlyFixedCoupons() const {
    for(int i=0; i<minCoupon.size(); i++) {
        if (!Maths::equals(minCoupon[i], maxCoupon[i])){
            return false;
        }
    }
    return true;
}

// note only called if iStep is a coupon step
double SPICouponsStd::getCoupon(double B,
                         double BF,
                         int iStep) const {
    double coup = 0.0;
    int iCoupon = iCouponMap[iStep];
    if (iCoupon == 0) {
        // reset current strike if first coupon
        currentStrike = couponStrike;
    }
    coup = Maths::max(minCoupon[iCoupon],
                      Maths::min(B-BF,
                                 Maths::min(maxCoupon[iCoupon],
                                            fixedCoupon[iCoupon] +
                                                Maths::max(participations[iCoupon] * (B - currentStrike), 0.0))));
    // if we need to change the strike do it ready for next coupon
    if (struckAtPrevious) {
        currentStrike = struckEx ? B - coup : B;
    }
    return coup;
}

const DateTime& SPICouponsStd::getCouponPayDate(int iStep) const {
    int iCoupon = iCouponMap[iStep];
    if (iCoupon>=0) {
        return couponPayDates[iCoupon];
    }
    throw ModelException("SPICouponsStd::getCouponPayDate",
                         "Internal error : No coupon date at step " + Format::toString(iStep));
}

double SPICouponsStd::getCouponPVFactor(int               iStep,
                                 const YieldCurve* disc) const {
    int iCoupon = iCouponMap[iStep];
    if (iCoupon>=0) {
        return disc->pv(couponPayDates[iCoupon]);
    }
    return 0.0;
}

double SPICouponsStd::getCouponPVFactor(int               iStep,
                                        SVDiscFactorSP dfSV) const {
    int iCoupon = iCouponMap[iStep];
    if (iCoupon>=0) {
        return dfSV->getDF(iCoupon);
    }
    return 0.0;
}

bool SPICouponsStd::isCouponDate(const DateTime& myDate) const {
    try {
        myDate.find(couponDates);
        return true;
    } catch (exception&) {
        return false;
    }
}

bool SPICouponsStd::hasEarlyRedemption() const {
    return allowEarlyRedemption;
}

double SPICouponsStd::getTarget() const {
    return targetLevel;
}

bool SPICouponsStd::isCouponStep(int iStep) const {
    return iCouponMap[iStep] >= 0;
}

bool SPICouponsStd::isThresholdBreached(double level, int iStep) const {
    // only called for a coupon step
    return Maths::isPositive(level-threshold[iCouponMap[iStep]]);
}

string SPICouponsStd::getThresholdType() const {
    for (int i = 0; i < threshold.size(); i++) {
        if (Maths::isPositive(threshold[i])) {
            return thresholdType;
        }
    }
    return SPI_THRESHOLD_TYPE_NONE;
}

// checks whether the target is met. NB: if it is we also check the overall cap
// We update the sum of coupons (and the coupon if the cap is triggered)
bool SPICouponsStd::isTargetMet(int     iStep,
                                double&  couponAmt,
                                double&  sumCouponsSoFar,
                                double&  bonusCoupon) const {
    if (!allowEarlyRedemption) {
        return false;
    }
    bool isMet = false;
    sumCouponsSoFar += couponAmt;
    // this check may be unnecessary if called always when it's known to be a coupon paying step
    int iCoupon = iCouponMap[iStep];
    if (iCoupon>=0 &&
        !Maths::isNegative(sumCouponsSoFar - targetLevel)) {
        isMet = true;
        bonusCoupon = bonusCoupons[iCoupon];
        // now apply cap if necessary
        if (hasCouponCap && Maths::isPositive(sumCouponsSoFar - couponCap)) {
            couponAmt -= (sumCouponsSoFar - couponCap);
            sumCouponsSoFar = couponCap;
        }
    }
    return isMet;
}

DateTimeArray SPICouponsStd::getPaymentDates() const {
    return couponPayDates;
}

DateTimeArray SPICouponsStd::getEssentialDates() const {
    return couponDates;
}

SPICouponsStd::SPICouponsStd(CClassConstSP clazz): SPIInterfaceIMS(clazz), minCoupon(0), fixedCoupon(0),
                     thresholdType(SPI_THRESHOLD_TYPE_BASKET),
                     couponStrike(1.0), struckAtPrevious(false), struckEx(false),
                     allowEarlyRedemption(false), targetLevel(0.0), bonusCoupons(0),
                     hasCouponCap(false), couponCap(0.0) {}

IObject* SPICouponsStd::defaultSPICouponsStd(){
    return new SPICouponsStd();
}

/** Invoked when Class is 'loaded' */
void SPICouponsStd::load(CClassSP& clazz){
    REGISTER(SPICouponsStd, clazz);
    SUPERCLASS(SPIInterfaceIMS);
    EMPTY_SHELL_METHOD(defaultSPICouponsStd);
    IMPLEMENTS(ICouponsSPIInterface);
    FIELD(couponDates,"Notification dates of coupons");
    FIELD(minCoupon,"Floor - one per coupon date");
    FIELD(maxCoupon,"Cap - one per coupon date");
    FIELD(fixedCoupon,"Fixed piece of the coupon - one per coupon date");
    FIELD(participations,"One per coupon date");
    FIELD(threshold,"Coupon only paid if thresholdType exceeds this level. One per coupon date");
    FIELD(thresholdType,"Either Basket or Target Exposure");
    FIELD(couponStrike,"Coupon strike");
    FIELD(struckAtPrevious,"Coupon struck at basket level at previous coupon date");
    FIELD(struckEx,"If struckAtPrevious is level ex coupon");
    FIELD(settlePeriod,"#Days settlement for coupons");
    FIELD(allowEarlyRedemption,"true=>switches on target level and bonus coupons");
    FIELD(targetLevel,"If allowEarlyRedemption, level for sum of coupons which triggers early redemption");
    FIELD(bonusCoupons,"If allowEarlyRedemption, the amount redeemed");
    FIELD(hasCouponCap, "Is there an overall cap for the total sum of coupons");
    FIELD(couponCap, "The overall cap for the total total sum of coupons");
    FIELD_NO_DESC(couponPayDates);
    FIELD_MAKE_TRANSIENT(couponPayDates);
// All must be optional to allow the IMS interface to be selective
    FIELD_MAKE_OPTIONAL(couponDates);
    FIELD_MAKE_OPTIONAL(minCoupon);
    FIELD_MAKE_OPTIONAL(maxCoupon);
    FIELD_MAKE_OPTIONAL(fixedCoupon);
    FIELD_MAKE_OPTIONAL(participations);
    FIELD_MAKE_OPTIONAL(threshold);
    FIELD_MAKE_OPTIONAL(thresholdType);
    FIELD_MAKE_OPTIONAL(couponStrike);
    FIELD_MAKE_OPTIONAL(struckAtPrevious);
    FIELD_MAKE_OPTIONAL(struckEx);
    FIELD_MAKE_OPTIONAL(settlePeriod);
    FIELD_MAKE_OPTIONAL(allowEarlyRedemption);
    FIELD_MAKE_OPTIONAL(targetLevel);
    FIELD_MAKE_OPTIONAL(bonusCoupons);
    FIELD_MAKE_OPTIONAL(hasCouponCap);
    FIELD_MAKE_OPTIONAL(couponCap);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPICouponsStd::TYPE = CClass::registerClassLoadMethod(
    "SPICouponsStd", typeid(SPICouponsStd), SPICouponsStd::load);

/***********************************************************************************/

SPICouponsStdSuper::SPICouponsStdSuper(): SPICouponsStd(TYPE) {} // for reflection

IObject* SPICouponsStdSuper::defaultSPICouponsStdSuper(){
    return new SPICouponsStdSuper();
}

/** Invoked when Class is 'loaded' */
void SPICouponsStdSuper::load(CClassSP& clazz){
    REGISTER(SPICouponsStdSuper, clazz);
    SUPERCLASS(SPICouponsStd);
    EMPTY_SHELL_METHOD(defaultSPICouponsStdSuper);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPICouponsStdSuper::TYPE = CClass::registerClassLoadMethod(
    "SPICouponsStdSuper", typeid(SPICouponsStdSuper), SPICouponsStdSuper::load);

/***********************************************************************************/
SPICouponsWrapper::SPICouponsWrapper(CClassConstSP clazz): CObject(clazz){}

ICouponsSPISP SPICouponsWrapper::getCouponsSPI() {
    static const string routine = "SPICouponsWrapper::getCouponsSPI";
    ICouponsSPISP theCoupons;

    if (SPICouponsType.empty() ||
        SPICouponsType==SPI_COUPON_TYPE_NONE) {
        if (!couponsNone.get()) {
            couponsNone = SPICouponsNoneSP(new SPICouponsNone());
        }
        theCoupons = ICouponsSPISP(couponsNone.get(), NullDeleter()); //FIXME
    } else if (SPICouponsType==SPI_COUPON_TYPE_STD) {
        if (!couponsStd.get()) {
            throw ModelException(routine, "Expected couponsStd but none supplied!");
        }
        theCoupons = ICouponsSPISP(couponsStd.get(), NullDeleter()); // FIXME
    } else {
        throw ModelException(routine, "Unrecognised SPICouponsType " + SPICouponsType +
                             ". Expected " + SPI_COUPON_TYPE_STD + ", " +
                             SPI_COUPON_TYPE_NONE);
    }

	// need to cast for validity check
	SPIInterfaceIMS* actualCoupons = dynamic_cast<SPIInterfaceIMS*>(theCoupons.get());
	if (!actualCoupons) {
        throw ModelException(routine, "Internal error in SPI coupon validation");
	}
    if (!actualCoupons->isValid()) {
        throw ModelException(routine, "Invalid SPICoupons (type " + SPICouponsType +
                             ") : " + actualCoupons->errString());
    }
    return theCoupons;
}

bool SPICouponsWrapper::doesNothing() const {
    return (SPICouponsType.empty() || SPICouponsType==SPI_COUPON_TYPE_NONE);
}


/** Invoked when Class is 'loaded' */
void SPICouponsWrapper::load(CClassSP& clazz){
    REGISTER(SPICouponsWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICouponsSPIInterface);
    EMPTY_SHELL_METHOD(defaultSPICouponsWrapper);
    FIELD(SPICouponsType, "SPICouponsType");
    FIELD(couponsNone,  "Coupons - None");
    FIELD_MAKE_OPTIONAL(couponsNone);
    FIELD(couponsStd,  "Coupons Standard");
    FIELD_MAKE_OPTIONAL(couponsStd);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// for reflection
SPICouponsWrapper::SPICouponsWrapper(): CObject(TYPE){}

IObject* SPICouponsWrapper::defaultSPICouponsWrapper(){
    return new SPICouponsWrapper();
}
CClassConstSP const SPICouponsWrapper::TYPE = CClass::registerClassLoadMethod(
    "SPICouponsWrapper", typeid(SPICouponsWrapper), SPICouponsWrapper::load);

/****************************************************************************/

// empty derived class which is just an SPICouponsWrapper
// this enables us to split the interface at the IMS level
// into SPI and Super SPI. We restrict access to the full set of
// params for the former at the Aladdin/IMS level
SPICouponsWrapperSuper::SPICouponsWrapperSuper(): SPICouponsWrapper(TYPE) {} // for reflection

IObject* SPICouponsWrapperSuper::defaultSPICouponsWrapperSuper(){
    return new SPICouponsWrapperSuper();
}

/** Invoked when Class is 'loaded' */
void SPICouponsWrapperSuper::load(CClassSP& clazz){
    REGISTER(SPICouponsWrapperSuper, clazz);
    SUPERCLASS(SPICouponsWrapper);
    EMPTY_SHELL_METHOD(defaultSPICouponsWrapperSuper);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPICouponsWrapperSuper::TYPE = CClass::registerClassLoadMethod(
    "SPICouponsWrapperSuper", typeid(SPICouponsWrapperSuper), SPICouponsWrapperSuper::load);

DRLIB_END_NAMESPACE
