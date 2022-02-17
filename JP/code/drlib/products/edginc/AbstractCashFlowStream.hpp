//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : AbstractCashFlowStream.hpp
//
//   Description : Price a stream of arbitrary cashflows
//
//   Author      : Gordon Stephens
//
//   Date        : 15 June 2005
//

//
//----------------------------------------------------------------------------

#ifndef ABSTRACTCASHFLOWSTREAM_HPP
#define ABSTRACTCASHFLOWSTREAM_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/AbstractCashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

/** SimpleCashFlow instrument - fixed amounts at known future dates */

class PRODUCTS_DLL AbstractCashFlowStream: public CInstrument, 
                            public virtual ClosedForm::IIntoProduct,
                            public virtual LastSensDate,
                            public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity*) const;

    virtual bool sensShift(Theta* shift);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    AbstractCashFlowStream();
    AbstractCashFlowStream(AbstractCashFlowArray cfl, const string&  discountName);

    AbstractCashFlowStream(const AbstractCashFlowStream& rhs);
    AbstractCashFlowStream& operator=(const AbstractCashFlowStream& rhs);
    void price(Control* control, CResults* results, IForwardRatePricerSP model) const;
    void requests(Control* control, CResults* results, IForwardRatePricerSP model) const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultAbstractCashFlowStream();

    AbstractCashFlowArray   cfl;
    YieldCurveWrapper discount;
    DateTime          valueDate;
};

typedef smartPtr<AbstractCashFlowStream> AbstractCashFlowStreamSP;

DRLIB_END_NAMESPACE
#endif
