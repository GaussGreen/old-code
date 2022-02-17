//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashFlowStream.hpp
//
//   Description : CashFlow instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 16 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef CASHFLOWSTREAM_HPP
#define CASHFLOWSTREAM_HPP
#include "edginc/SimpleCashFlowStream.hpp"
#include "edginc/CashFlowStreamCreditSupport.hpp"
#include "edginc/ITaxableInst.hpp"

DRLIB_BEGIN_NAMESPACE

/** CashFlow instrument - fixed amounts at known future dates */

class PRODUCTS_DLL CashFlowStream: public SimpleCashFlowStream,
                      public CreditSupport::Interface,
                      public ITaxableInst::Basic,
                      public ITaxableInst::WithCoupons {
public:
    static CClassConstSP const TYPE;
    friend class CashFlowStreamHelper;
    friend class CashFlowStreamClosedForm;
    friend class CashFlowStreamCreditSupport;

    /** Constructor only needed for demo purposes */
    CashFlowStream(CashFlowArray* cfl,
                   const string&  discountName);
              
    virtual CreditSupportSP createCreditSupport(CMarketDataSP market){
                return CreditSupportSP(new CashFlowStreamCreditSupport(this, market));}

    /** Support for ITaxableInstBasic */
    virtual const DateTime getFinalPaymentDate() const;

    /** Support for ITaxableInstWithCoupons */
    virtual CashFlowArrayConstSP getCoupons() const;

private:
    CashFlowStream();
    CashFlowStream(const CashFlowStream& rhs);
    CashFlowStream& operator=(const CashFlowStream& rhs);
    void price(Control* control, CResults* results) const;
    void requests(Control* control, CResults* results) const;

};

typedef smartPtr<CashFlowStream> CashFlowStreamSP;

DRLIB_END_NAMESPACE
#endif
