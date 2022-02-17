//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubSimple.hpp
//
//   Description : Simple stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef STUBSIMPLE_HPP
#define STUBSIMPLE_HPP

#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

/** Simple stub */

class MARKET_DLL StubSimple: public Stub {
public:
    static CClassConstSP const TYPE;
    friend class StubSimpleHelper;

    StubSimple();
    virtual ~StubSimple();
    
    /** how big is the stub payment ? */
    virtual double payment(const DateTime&           prevCouponDate,
                           const DateTime&           nextCouponDate,
                           const DateTime&           stubStart,
                           const DateTime&           stubEnd,
                           double                    rate,
                           const DayCountConvention* dcc) const;
    
    /** returns a string description e.g. SIMPLE */
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
