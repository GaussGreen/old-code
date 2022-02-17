//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StubBond.hpp
//
//   Description : Bond stub
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef STUBBOND_HPP
#define STUBBOND_HPP

#include "edginc/Stub.hpp"

DRLIB_BEGIN_NAMESPACE

/** Bond stub */

class MARKET_DLL StubBond: public Stub {
public:
    static CClassConstSP const TYPE;
    friend class StubBondHelper;

    StubBond();
    virtual ~StubBond();
    
    /** how big is the stub payment ? */
    virtual double payment(const DateTime&           prevCouponDate,
                           const DateTime&           nextCouponDate,
                           const DateTime&           stubStart,
                           const DateTime&           stubEnd,
                           double                    rate,
                           const DayCountConvention* dcc) const;
    
    /** returns a string description e.g. Bond */
    virtual string toString() const;
};

DRLIB_END_NAMESPACE


#endif
