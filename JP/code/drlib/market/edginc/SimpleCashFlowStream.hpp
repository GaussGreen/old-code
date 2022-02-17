//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SimpleCashFlowStream.hpp
//
//   Description : CashFlow instrument with no credit or tax support
//                 Taken from CashFlowStream.hpp
//
//   Author      : Gordon Stephens
//
//   Date        : 13 April 2005
//
//
//----------------------------------------------------------------------------

#ifndef SIMPLECASHFLOWSTREAM_HPP
#define SIMPLECASHFLOWSTREAM_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** SimpleCashFlow instrument - fixed amounts at known future dates */

class MARKET_DLL SimpleCashFlowStream: public CInstrument, 
                            public virtual ClosedForm::IIntoProduct,
                            public virtual LastSensDate,
                            public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;
    friend class SimpleCashFlowStreamHelper;
    friend class SimpleCashFlowStreamClosedForm;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    SimpleCashFlowStream(CashFlowArray* cfl, const string&  discountName);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

protected:
    SimpleCashFlowStream();
    SimpleCashFlowStream(CClassConstSP clazz);
    SimpleCashFlowStream(CClassConstSP clazz, CashFlowArray* cfl, const string&  discountName);

    SimpleCashFlowStream(const SimpleCashFlowStream& rhs);
    SimpleCashFlowStream& operator=(const SimpleCashFlowStream& rhs);
    void price(Control* control, CResults* results) const;
    void requests(Control* control, CResults* results) const;

    CashFlowArraySP   cfl;
    YieldCurveWrapper discount;
    DateTime          valueDate;
};

typedef smartPtr<SimpleCashFlowStream> SimpleCashFlowStreamSP;

DRLIB_END_NAMESPACE
#endif
