//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadBoy.hpp
//
//   Description : Class that always leaks, has access errors
//                 and differences to prove that testcmp and
//                 purify are working
//
//   Author      : Andrew J Swain
//
//   Date        : 17 August 2001
//
//----------------------------------------------------------------------------

#ifndef QLIB_BADBOY_HPP
#define QLIB_BADBOY_HPP

#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"

DRLIB_BEGIN_NAMESPACE

class ADDINS_DLL BadBoy: public CInstrument, 
              public ClosedForm::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual void GetMarket(const IModel*, 
                           const CMarketDataSP) {};

    virtual void Validate() {};

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    /** Returns the value date (aka today) the instrument is currently
        pricing for */
    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency */
    virtual string discountYieldCurveName() const;

private:
    BadBoy();
    void price(CResults* results) const;
    friend class BadBoyHelper;
    friend class BadBoyClosedForm;

};

DRLIB_END_NAMESPACE

#endif
