//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KCashflow.hpp
//
//   Description : Array of cashflow payments on specific date
//               : Payments are assumed to all be in same currency
//----------------------------------------------------------------------------

#ifndef KCASHFLOW_HPP
#define KCASHFLOW_HPP

#include "edginc/config.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/IndexSpecFX.hpp"

DRLIB_BEGIN_NAMESPACE

class KCashflow : public KComponent,
                  virtual public FDModel::IIntoProduct 
{
    friend class KCashflowTree;

    /************* reflection ************/
public:
    static CClassConstSP const TYPE;
    KCashflow(CClassConstSP const &type) : KComponent(type), cfType(CashflowInfo::CfType::UNSET) {}
private:
    static void load(CClassSP& clazz); 
    static IObject* defaultConstructor(void) { return new KCashflow(TYPE); }

    /******** exported variables *********/
    DateTimeArray dates;
    DoubleArray amounts;
    CashflowInfo::CfType::Enum cfType;

    /******** transient variables ********/
    IndexSpecFXSP fx; /**< convert from leg currency to pricing currency */

    /************ functions **************/
public:
    /* KComponent:: */
    virtual DateTime getLastDate(void) const;

protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* KComponent:: */
    virtual void reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const;
};

typedef smartPtr<KCashflow> KCashflowSP;
typedef smartConstPtr<KCashflow> KCashflowConstSP;
typedef array<KCashflowSP, KCashflow> KCashflowArray;

DRLIB_END_NAMESPACE
#endif
