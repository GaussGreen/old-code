//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubNone.hpp
//
//   Description : No stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef STUBNONE_HPP
#define STUBNONE_HPP

#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

/** Simple stub */

class MARKET_DLL StubNone: public Stub {
public:
    static CClassConstSP const TYPE;
    friend class StubNoneHelper;

    StubNone();
    virtual ~StubNone();
    
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
