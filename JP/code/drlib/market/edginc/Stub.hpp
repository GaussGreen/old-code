//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Stub.hpp
//
//   Description : Stub payment interface
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDGSTUB_HPP
#define EDGSTUB_HPP

#include "edginc/DateTime.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/config.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE
// defines an interface to be implemented by concrete classes
class MARKET_DLL Stub:public CObject {
public:
    static CClassConstSP const TYPE;

    virtual ~Stub();

    /** how big is the stub payment ? */
    virtual double payment(const DateTime&           prevCouponDate,
                           const DateTime&           nextCouponDate,
                           const DateTime&           stubStart,
                           const DateTime&           stubEnd,
                           double                    rate,
                           const DayCountConvention* dcc) const = 0;

    /** returns a string description e.g. SIMPLE */
    virtual string toString() const = 0;

protected:
    friend class StubHelper;
    Stub(CClassConstSP clazz);
};

typedef smartConstPtr<Stub> StubConstSP;
typedef smartPtr<Stub> StubSP;

DRLIB_END_NAMESPACE
#endif
