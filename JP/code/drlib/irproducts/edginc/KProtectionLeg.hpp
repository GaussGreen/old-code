//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : KProtectionLeg.hpp
//
//   Description : KRiskyComponent that pays out losses between specified
//                 attachment points (i.e. CDS tranche protection) for a
//                 specified period. Forward settlement is simple (i.e.
//                 not like a forward tranche!) and may be conditional
//                 or unconditional. The recovery rate may be determined
//                 by the market, or may be digital.
//
//   Author      : Charles Morcom
//
//   Date        : May 15, 2006
//
//----------------------------------------------------------------------------
#ifndef QR_KPROTECTIONLEG_HPP
#define QR_KPROTECTIONLEG_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/KRiskyComponent.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(KProtectionLeg)

/**A KProtectionLeg pays out contingent on default. If settlementIsUnconditional, losses before
  * the start date are paid at the start date. All other losses are paid as they are incurred. */
class KProtectionLeg : public KRiskyComponent,
                  virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const  TYPE;
    /* DATA MEMBERS */
    double      notional;
    DateTime    protectionStartDate;
    DateTime    protectionEndDate;
    bool        isMarketRecovery;
    double      digitalRecoveryRate;
    double      kLo;
    double      kHi;
    bool        settlementIsUnconditional;
    friend class KProtectionLegTree;


public:
    /* KComponent:: */
    virtual DateTime getLastDate() const;
    /* CObject:: */
    virtual void validatePop2Object();

    KProtectionLeg(
            const string& discount,
            const string& cdsParSpreads,
            const string& outputName, 
            double notional,
            const DateTime& protectionStartDate,
            const DateTime& protectionEndDate,
            bool settlementIsConditional = true,
            bool marketRecovery = true,
            double digitalRecovery = 0.0,
            double lowerAttachment = 0.0,
            double upperDetachment = 1.0)
        :   KRiskyComponent(discount, cdsParSpreads, outputName, TYPE),
        notional(notional),
        protectionStartDate(protectionStartDate),
        protectionEndDate(protectionEndDate),
        isMarketRecovery(marketRecovery),
        digitalRecoveryRate(digitalRecovery),
        kLo(lowerAttachment),
        kHi(upperDetachment)
    { validatePop2Object(); }

protected:
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* KComponent:: */
    void reportEvents(const KnownCashflows*, IModel* model,
        const DateTime& eDate, EventResults* events) const;

    KProtectionLeg(CClassConstSP const &type) 
        : KRiskyComponent("","","",type), notional(0),
        isMarketRecovery(true), digitalRecoveryRate(0.0), kLo(0.0), kHi(1.0), 
        settlementIsUnconditional(true) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KProtectionLeg(TYPE); }
};

DRLIB_END_NAMESPACE

#endif


