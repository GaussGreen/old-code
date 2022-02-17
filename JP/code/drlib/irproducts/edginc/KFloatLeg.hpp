//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KFloatLeg.hpp
//
//   Description : floater component
//
//----------------------------------------------------------------------------
#ifndef _KFLOATLEG_HPP
#define _KFLOATLEG_HPP

#include "edginc/KComponent.hpp"
#include "edginc/CouponSched.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MonteCarlo.hpp"

DRLIB_BEGIN_NAMESPACE

class KFloatLeg : public KComponent,
                  public virtual FDModel::IIntoProduct,
                  public virtual IMCIntoProduct
{
public:
    static CClassConstSP const TYPE;
    friend class KFloatLegTree;

    /****************** exported fields ************/
public: // can be modified by shell instruments
    // for each coupon - begin
    CouponSchedDatesSP sched;
    bool payInitialPrincipal; // pay initial notional value on accrual start date
    bool payPrincipal;        // pay the difference of the notional at each coupon
    IProdCreatorSP index;     // underlying index for floater

protected:
    DoubleArraySP   notionals;  // notional
    DoubleArraySP   weights;    // coupon is weights[]*index + spreads[]
    DoubleArraySP   spreads;    
    DoubleArraySP   dcfs;       // If dcfs are provided, dcc is ignored
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

    KFloatLeg(
            const string &discount, 
            const string &outputName, 
            CouponSchedDatesSP sched,
            DoubleArraySP   notionals,
            DoubleArraySP   weights,
            DoubleArraySP   spreads,
            const RateType::Enum& rateType,
            bool payInitialPrincipal,
            bool payPrincipal,
            DayCountConventionSP dcc,
            IProdCreatorSP index)
        :   KComponent(discount, outputName, TYPE), 
            sched(sched),
            payInitialPrincipal(payInitialPrincipal),
            payPrincipal(payPrincipal),
            index(index),
            notionals(notionals),
            weights(weights),
            spreads(spreads),
            rateType(rateType),
            dcc(dcc)
        { validatePop2Object(); }

protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /*IMCIntoProduct interface method */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* KComponent:: */
    virtual void reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const;

    KFloatLeg(CClassConstSP const &type) 
        : KComponent(type), payInitialPrincipal(false), payPrincipal(true),
        rateType(RateType::SIMPLE) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KFloatLeg(TYPE); }
};
DECLARE(KFloatLeg);

DRLIB_END_NAMESPACE

#endif
