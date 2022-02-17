//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPICoupons.hpp
//
//   Description : Coupons interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_CPN_HPP
#define EDR_SPI_CPN_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

#define SPI_THRESHOLD_TYPE_BASKET   "Basket" // the default for backwards compatibility
#define SPI_THRESHOLD_TYPE_TE       "Target Exposure"
#define SPI_THRESHOLD_TYPE_NONE     "No Threshold" // for internal use - denotes zero threshold

/*****************************************************************************/
// this is the actual interface of what a coupon object 
// needs to do internally
class PRODUCTS_DLL ICouponsSPI {
public:

    virtual void init(const DateTimeArray& rebalDates) = 0;

    virtual double getCoupon(double B,
                             double BF,
                             int iStep) const = 0;

    virtual const DateTime& getCouponPayDate(int iStep) const = 0;
    // Present value from coupon pay date. Having it as part of the interface
    // means easier to keep conditions out of main code.
    virtual double getCouponPVFactor(int               iStep,
                                     const YieldCurve* disc) const = 0;

    virtual double getCouponPVFactor(int               iStep,
                                     SVDiscFactorSP dfSV) const = 0;

    virtual bool isCouponDate(const DateTime& myDate) const = 0;

    virtual bool isCouponStep(int iStep) const = 0;

    virtual bool isThresholdBreached(double level, int iStep) const = 0;

    virtual string getThresholdType() const = 0;

    // indicates whether bond floor is to be adjusted
    virtual bool guaranteedCoupons() const = 0;
    // provides the min coupon - if not a coupon date will simply be 0
    virtual double getGuaranteedCoupon(int iStep) const = 0;

    // indicates (for Rainbow SPI whether there are only guaranteed coupons)
    virtual bool onlyFixedCoupons() const = 0;

    virtual bool hasEarlyRedemption() const = 0;

    virtual double getTarget() const = 0;

    // Could possibly do all this inside Coupon class but that would need
    // methods to reset the sum and be aware of past etc. Not now ... XXX
    // If this returns true then bonusCoupon will be set appropriately.
    virtual bool isTargetMet(int     iStep,
                             double&  couponAmt,
                             double&  sumCouponsSoFar,
                             double&  bonusCoupon) const = 0;

    // not sure why this yet... KNOWN_CASHFLOWS at some stage??
    virtual DateTimeArray getPaymentDates() const = 0;

    // for building of sim date list
    virtual DateTimeArray getEssentialDates() const = 0;

    // post getMarket validation (handles settlement of coupons)
    virtual void Validate(const InstrumentSettlement* instSettle) = 0;

    virtual ~ICouponsSPI();
};
DECLARE_REF_COUNT(ICouponsSPI);

/*****************************************************************************/
// this is the external interface for abstraction so that users can bolt in 
// any Coupons type they want - note this includes the SPICouponsWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ICouponsSPI as soon as possible
class PRODUCTS_DLL ICouponsSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual ICouponsSPISP getCouponsSPI() = 0;

    virtual bool doesNothing() const = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<ICouponsSPIInterface> ICouponsSPIInterfaceSP;

/*****************************************************************************/
/** A placeholder to allow coupons to be optional
 */
class PRODUCTS_DLL SPICouponsNone : public SPIInterfaceIMS,
                      virtual public ICouponsSPI,
                      virtual public ICouponsSPIInterface {
public:
    static CClassConstSP const TYPE;
    SPICouponsNone(); // for reflection

    virtual ICouponsSPISP getCouponsSPI();

    virtual bool doesNothing() const;

    virtual void init(const DateTimeArray& rebalDates);

    virtual double getCoupon(double B,
                             double BF,
                             int iStep) const;

    virtual double getCouponPVFactor(int               iStep,
                                     const YieldCurve* disc) const;

    virtual double getCouponPVFactor(int               iStep,
                                     SVDiscFactorSP dfSV) const;

    virtual const DateTime& getCouponPayDate(int iStep) const;

    virtual bool isCouponDate(const DateTime& myDate) const;

    virtual bool isCouponStep(int iStep) const;

    virtual bool isThresholdBreached(double level, int iStep) const;

    virtual string getThresholdType() const;

    virtual bool guaranteedCoupons() const;

    virtual double getGuaranteedCoupon(int iStep) const;

    virtual bool onlyFixedCoupons() const;

    virtual bool hasEarlyRedemption() const;

    virtual double getTarget() const;

    virtual bool isTargetMet(int     iStep,
                             double&  couponAmt,
                             double&  sumCouponsSoFar,
                             double&  bonusCoupon) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual DateTimeArray getEssentialDates() const;

    // post getMarket validation (handles settlement of coupons)
    virtual void Validate(const InstrumentSettlement* instSettle);

private:
    // I think IMS requires something
    string               dummyString;

    SPICouponsNone(const SPICouponsNone& rhs); // not implemented
    SPICouponsNone& operator=(const SPICouponsNone& rhs); // not implemented

    static IObject* defaultSPICouponsNone();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPICouponsNone> SPICouponsNoneSP;

class PRODUCTS_DLL SPICouponsStd : public SPIInterfaceIMS,
                      virtual public ICouponsSPI,
                      virtual public ICouponsSPIInterface {
private:
    // all arrays have this length
    DateTimeArray couponDates;       // "notification" - paid a bit later
    DoubleArray   minCoupon;         // for each coupon
    DoubleArray   maxCoupon;
    DoubleArray   fixedCoupon;      // a fixed piece of the coupon
    DoubleArray   participations;
    DoubleArray   threshold;
    string        thresholdType;    // coupons contingent on basket or TE above some level
    double        couponStrike;
    bool          struckAtPrevious;  // can override strike with basket level at previous coupon
    bool          struckEx;          // if struckAtPrevious true is basket level ex (i.e. post) coupon
    int           settlePeriod;      // will become a CashSettlePeriod with same hols as inst

    // Target level for sum of coupons so far, and if met an early redemption with principal plus bonus
    bool          allowEarlyRedemption;  // true=> have target and bonus; else ignore these
    double        targetLevel;
    DoubleArray   bonusCoupons;
    bool          hasCouponCap;
    double        couponCap; //cap for the total value of all coupons

    mutable double currentStrike; // strike for current coupon. Could be couponStrike or basket value $unregistered
                                 // at last coupon date (cum or ex)
    // transient
    IntArray      iCouponMap; // from iStep -> iCoupon $unregistered
    DateTimeArray couponPayDates;

public:
    static CClassConstSP const TYPE;
    SPICouponsStd();// for reflection

    // validation
    void validatePop2Object();

    virtual ICouponsSPISP getCouponsSPI();

    virtual bool doesNothing() const;

    void init(const DateTimeArray& rebalDates);

    virtual bool guaranteedCoupons() const;

    virtual double getGuaranteedCoupon(int iStep) const;

    virtual bool onlyFixedCoupons() const;

    // note only called if iStep is a coupon step
    virtual double getCoupon(double B,
                             double BF,
                             int iStep) const;

    virtual const DateTime& getCouponPayDate(int iStep) const;

    virtual double getCouponPVFactor(int               iStep,
                                     const YieldCurve* disc) const;

    virtual double getCouponPVFactor(int               iStep,
                                     SVDiscFactorSP dfSV) const;

    virtual bool isCouponDate(const DateTime& myDate) const;

    virtual bool hasEarlyRedemption() const;

    virtual double getTarget() const;

    virtual bool isCouponStep(int iStep) const;

    virtual bool isThresholdBreached(double level, int iStep) const;

    virtual string getThresholdType() const;

    virtual bool isTargetMet(int     iStep,
                             double&  couponAmt,
                             double&  sumCouponsSoFar,
                             double&  bonusCoupon) const;

    virtual DateTimeArray getPaymentDates() const;

    virtual DateTimeArray getEssentialDates() const;

    // post getMarket validation (handles settlement of coupons)
    virtual void Validate(const InstrumentSettlement* instSettle);

protected:

    SPICouponsStd(CClassConstSP clazz);

private:
    SPICouponsStd(const SPICouponsStd& rhs); // not implemented
    SPICouponsStd& operator=(const SPICouponsStd& rhs); // not implemented

    static IObject* defaultSPICouponsStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPICouponsStd> SPICouponsStdSP;

/***********************************************************************************/

// empty derived class which is just an SPICouponsStd
// this enables us to split the interface at the IMS level
// into SPI and Super SPI. We restrict access to the full set of
// params for the former at the Aladdin/IMS level
class PRODUCTS_DLL SPICouponsStdSuper : public SPICouponsStd {
public:
    static CClassConstSP const TYPE;
    SPICouponsStdSuper(); // for reflection

private:
    SPICouponsStdSuper(const SPICouponsStdSuper& rhs); // not implemented
    SPICouponsStdSuper& operator=(const SPICouponsStdSuper& rhs); // not implemented

    static IObject* defaultSPICouponsStdSuper();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

/***********************************************************************************/

#define SPI_COUPON_TYPE_NONE   "None"
#define SPI_COUPON_TYPE_STD    "Standard"

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPICouponsWrapper : public CObject,
                      virtual public ICouponsSPIInterface  {
public:
    static CClassConstSP const TYPE;

    string                SPICouponsType;
    SPICouponsNoneSP      couponsNone;
    SPICouponsStdSP       couponsStd;

protected:
    SPICouponsWrapper(CClassConstSP clazz);

public:
    virtual ICouponsSPISP getCouponsSPI();

    virtual bool doesNothing() const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    // for reflection
    SPICouponsWrapper();

    static IObject* defaultSPICouponsWrapper();
};
typedef smartPtr<SPICouponsWrapper> SPICouponsWrapperSP;

/****************************************************************************/

// empty derived class which is just an SPICouponsWrapper
// this enables us to split the interface at the IMS level
// into SPI and Super SPI. We restrict access to the full set of
// params for the former at the Aladdin/IMS level
class PRODUCTS_DLL SPICouponsWrapperSuper : public SPICouponsWrapper {
public:
    static CClassConstSP const TYPE;
    SPICouponsWrapperSuper(); // for reflection

private:
    SPICouponsWrapperSuper(const SPICouponsWrapperSuper& rhs); // not implemented
    SPICouponsWrapperSuper& operator=(const SPICouponsWrapperSuper& rhs); // not implemented

    static IObject* defaultSPICouponsWrapperSuper();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
