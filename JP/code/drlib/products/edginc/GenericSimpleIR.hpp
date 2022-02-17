//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericSimpleIR.hpp
//
//   Description : Base class for simple interest rate generic instruments
//
//   Author      : Andrew J Swain
//
//   Date        : 18 January 2002
//
//
//----------------------------------------------------------------------------

#ifndef _GENERICSIMPLEIR_HPP
#define _GENERICSIMPLEIR_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL GenericSimpleIR: public CInstrument,
                       virtual public Theta::Shift {
public:
    static CClassConstSP const TYPE;

    virtual ~GenericSimpleIR();

    virtual DateTime getValueDate()const;

    /** Get the asset and discount market data */
    virtual void GetMarket(const IModel*          model, 
                           const CMarketDataSP    market);
    
    bool sensShift(Theta* theta);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

protected:
    GenericSimpleIR(CClassConstSP clazz);

    DateTime          valueDate;
    YieldCurveWrapper discount;
    YieldCurveSP      coupon; // transient - set to 'projection' curve
    IRVolBaseWrapper  vol;

private:
    friend class GenericSimpleIRHelper;

    GenericSimpleIR(); // not implemented
    GenericSimpleIR(const GenericSimpleIR& rhs);
    GenericSimpleIR& operator=(const GenericSimpleIR& rhs);
};

typedef smartPtr<GenericSimpleIR> GenericSimpleIRSP;

DRLIB_END_NAMESPACE
#endif

