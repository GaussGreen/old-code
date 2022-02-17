//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KRibFloatLeg.hpp
//
//   Description : experimental combined RIB and floating leg component
//
//----------------------------------------------------------------------------
#ifndef _KRIBFLOATLEG_HPP
#define _KRIBFLOATLEG_HPP

#include "edginc/CouponSched.hpp"
#include "edginc/KRib.hpp"

DRLIB_BEGIN_NAMESPACE

class KRibFloatLeg : public KComponent,
                  virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KRibFloatLegTree;

    /****************** exported fields ************/
public: // can be modified by shell instruments
    // for each coupon - begin
    CouponSchedDatesSP sched;
    bool payInitialPrincipal; // pay initial notional value on accrual start date
    bool payPrincipal;        // pay the difference of the notional at each coupon
    IProdCreatorSP index;     // underlying index for floater

    // KRib is used as an interface component, and does not create a RIB FDProduct
    KRibSP   rib;

protected:
    DoubleArraySP   notionals;  // notional
    DoubleArraySP   weights;    // coupon is weights[]*index + spreads[]
    DoubleArraySP   spreads;    
    DoubleArraySP   dcfs;       // If dcfs are provided, dcc is ignored
    DoubleArraySP   capRates;    // cap rate as decimal
    DoubleArraySP   floorRates;  // floor rate as decimal
    // for each coupon - end

    RateType::Enum rateType;    // simple or continuous compounding rate
    DayCountConventionSP dcc; // day count convention

    /****************** transient fields ************/
private:
    // principal payment fields calculated in the setup function
    DateTimeArray principalDates;
    DoubleArray principalPayments;

    /************************ methods **********************/
public:
    /* KComponent:: */
    virtual DateTime getLastDate(void) const;
    /* KComponent:: */
    virtual void getAccEndDates(DateTimeArray &accEndDates) const;
    /* CObject:: */
    virtual void validatePop2Object(void);

protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* KComponent:: */
    void reportEvents(const KnownCashflows*, IModel* model,
        const DateTime& eDate, EventResults* events) const;

    KRibFloatLeg(CClassConstSP const &type) 
        : KComponent(type), payInitialPrincipal(false), payPrincipal(true),
        rateType(RateType::SIMPLE) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KRibFloatLeg(TYPE); }
};
typedef smartPtr<KRibFloatLeg> KRibFloatLegSP;
typedef smartConstPtr<KRibFloatLeg> KRibFloatLegConstSP;

DRLIB_END_NAMESPACE

#endif


