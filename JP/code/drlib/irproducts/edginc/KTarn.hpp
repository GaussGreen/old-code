//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KTarn.hpp
//
//   Description : Tarn component
//
//----------------------------------------------------------------------------
#ifndef _KTarn_HPP
#define _KTarn_HPP

#include "edginc/KComponent.hpp"

DRLIB_BEGIN_NAMESPACE

// tarn component
class KTarn : public KComponent,
              virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;    

    void setup(const IModel* model, const MarketData* market);

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    KTarn(CClassConstSP const &type = TYPE) 
    :   KComponent(type), targetLevel(0), knockInOut(KnockType::KNOCK_OUT), 
        finalPayment(true), payoffSmoothing(false), pastRealizedAmount(0.), 
        excessCoupon(false), nbStateVar(1), stateVarNbStdDeviations(0.) 
    {}

    struct KnockType {
        enum Enum {KNOCK_IN, KNOCK_OUT};
    };
    typedef BoxedEnum<KnockType::Enum> KnockTypeBoxedEnum;

private:
    friend class KTarnTree;

    /****************** exported fields ************/
    double          targetLevel;
    KnockType::Enum knockInOut;
    bool            finalPayment;
    bool            payoffSmoothing;
    double          pastRealizedAmount;  

    bool            excessCoupon; // if this is turned on, a full complex coupon is paid upon knock-out. Otherwise, the last complex payment equals redemption level minus sum of past coupons.

    int		        nbStateVar;   // number of states to price this product
    double          stateVarNbStdDeviations;


    IProdCreatorSP  complexLeg; // the component that generates the monitored cashflows 
    IProdCreatorSP  fundingLeg; // the component that generates the monitored cashflows 

    /****************** transient fields ************/
    DateTimeArray   fundAccEndDates;
    DateTimeArray   cplxAccEndDates;
    DateTimeArray   cplxPay;
    /****************** methods ************/
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KTarn(); }
};
typedef smartPtr<KTarn> KTarnSP;
typedef smartConstPtr<KTarn> KTarnConstSP;

DRLIB_END_NAMESPACE

#endif


