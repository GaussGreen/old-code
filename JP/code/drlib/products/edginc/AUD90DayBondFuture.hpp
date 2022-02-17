//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AUD90DayBondFuture.hpp
//
//   Description : Interest rate future for AUD 90 days bond future.
//
//   Author      : Keiji Kitazaw
//
//   Date        : 20 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef _AUD90DayBondFuture_HPP
#define _AUD90DayBondFuture_HPP
#include "edginc/ClosedFormLN.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/ClosedFormIRLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interest rate future */
class PRODUCTS_DLL AUD90DayBondFuture: public GenericSimpleIR, 
                          public virtual CClosedFormLN::IIntoProduct,
                          public virtual ClosedFormIRLN::IIntoProduct,
                          public virtual ISensitiveIRVolPoints,                     
                          public virtual LastSensDate{
public:
    static CClassConstSP const TYPE;

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
    virtual ClosedFormIRLN::IProduct* createProduct(ClosedFormIRLN* model) const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    // for IR vega
    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP outputName,
        const IModel*      model) const;

private:
    friend class AUD90DayBondFutureHelper;
    friend class AUD90DayBondFutureClosedForm;

    AUD90DayBondFuture();
    AUD90DayBondFuture(const AUD90DayBondFuture& rhs);
    AUD90DayBondFuture& operator=(const AUD90DayBondFuture& rhs);
    void price(Control* control, CResults* results) const;

    DateTime          intStartDate;
    DateTime          intEndDate;
    DateTime          lastTradeDate;
    string            dcc;

};

DRLIB_END_NAMESPACE
#endif
