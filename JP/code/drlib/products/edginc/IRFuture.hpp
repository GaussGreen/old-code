//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRFuture.hpp
//
//   Description : Interest rate future
//
//   Author      : Andrew J Swain
//
//   Date        : 11 January 2002
//
//
//----------------------------------------------------------------------------

#ifndef _IRFUTURE_HPP
#define _IRFUTURE_HPP
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interest rate future */
class PRODUCTS_DLL IRFuture: public CInstrument, 
                public CClosedFormLN::IIntoProduct,
                public LastSensDate,
                public Theta::Shift {
public:
    static CClassConstSP const TYPE;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual bool sensShift(Theta* shift);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

private:
    friend class IRFutureHelper;
    friend class IRFutureClosedForm;

    IRFuture();
    IRFuture(const IRFuture& rhs);
    IRFuture& operator=(const IRFuture& rhs);
    void price(CResults* results) const;

    DateTime          intStartDate;
    DateTime          intEndDate;
    DateTime          lastTradeDate;
    string            dcc;

    IRVolBaseWrapper  vol;
    YieldCurveWrapper discount;
    DateTime          valueDate;
};

DRLIB_END_NAMESPACE
#endif
