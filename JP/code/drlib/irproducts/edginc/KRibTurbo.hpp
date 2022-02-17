//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRibTurbo.hpp
//
//   Description : combined RIB and floating/turbo component
//
//----------------------------------------------------------------------------

#ifndef KRIBTURBO_HPP
#define KRIBTURBO_HPP

#include "edginc/config.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/CouponSched.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/KRib2.hpp"

DRLIB_BEGIN_NAMESPACE


class KRibTurbo : public KComponent,
                  virtual public FDModel::IIntoProduct {

public:
    // A bit opaque calling this IndexLeg, but the generalized turbo is 
    // a special case because it is more than just and index (it controls 
    // principal payments, contains notional schedules for each coupon etc) 
    // but it is not a collection of swap legs as the FX is applied at the 
    // index reset date, NOT the payment date which means it can NOT be used
    // to model the sum of various currency payment swap legs that happen to 
    // share the same schedule
    class IndexLeg : public CObject {
    public:
        static CClassConstSP const TYPE;
        friend class KRibTurboTree;
        friend class KRibTurbo;
        friend class SimpleRibMQProd;
        friend class SimpleRib;

        /************ exported fields ************/
    public:
        YieldCurveWrapper  discount;/**< currency of the leg */

    protected:
        DoubleArraySP  notionals;    /**< Notional/coupon (in currency) */
        DoubleArraySP  spreads;      /**< Fixed spread/coupon */
        DoubleArraySP  weights;      /**< Weight/coupon for floating index */
        
        bool payInitialPrincipal;    /**< Pay initial principal on leg start date */
        bool payPrincipal;           /**< Pay any interim and final principals */
        /**< Pay any principals set to true in above 2 payPrinicpal fields if value of
             leg requested during life of the turbo */
        bool payInitialPrincipalAtExercise;

        RateType::Enum rateType;       /**< Simple/Continuous rate */
        DayCountConventionSP dcc;    /**< day count convention for leg, used if dcfs not provided */
        DoubleArray  dcfs;           /**< Day count fraction/coupon, used if dcc not provided */
        CModel::IProdCreatorArray indexes;   /**< reset object */
        IndexSpecFXSP fx;            /**<  convert from leg currency to pricing currency */

        /******** transient ********/
        DoubleArray principalPayments;  /*< calculated coupon principal payments */

        /******************** methods *******************/
    public:
        virtual void validatePop2Object(void);

    protected:
        IndexLeg(CClassConstSP const &type) :
            CObject(type), payInitialPrincipal(false), payPrincipal(false), 
            payInitialPrincipalAtExercise(false) {}

    private:
        static IObject* defaultConstructor(void) { return new IndexLeg(TYPE); }
        static void load(CClassSP& clazz);
    };
    DECLARE(IndexLeg);

    /****************************** variables ********************************/
    static CClassConstSP const TYPE;
    friend class KRibTurboTree;
    friend class SimpleRibMQProd;

    IndexLegArraySP legs;
    CouponSchedDatesSP sched;
    KRib2SP rib;

protected:    
    // payment cap/floor notional and dcc fields - may be same as one of the 
    // turbo legs or different, so defined here at the swap/payment level
    DoubleArraySP floorRates;
    DoubleArraySP capRates;
    DoubleArraySP capFloorNotionals;
    DayCountConventionSP capFloorDcc; // optional. If not provided, capFloorDcfs must be provided
    DoubleArray capFloorDcfs; // optional. If not provided, capFloorDcc must be provided

    StubType::Enum stubType;  /**< stub type if leg exercised/knocked out */
    
public:

    /* KComponent:: */
    virtual DateTime getLastDate(void) const;


protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* KComponent:: */
    virtual void reportCashFlows(CashflowInfoArray &cashflowInfos, 
        bool amountsNeeded ) const;

    KRibTurbo(CClassConstSP const &type) : KComponent(type) {};

private:
    static void load(CClassSP& clazz); 
    static IObject* defaultConstructor(void) { return new KRibTurbo(TYPE); }
};
DECLARE(KRibTurbo);

DRLIB_END_NAMESPACE
#endif
