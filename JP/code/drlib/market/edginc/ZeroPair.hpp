#ifndef EDR_ZEROPAIR_HPP
#define EDR_ZEROPAIR_HPP

#include "edginc/ZeroCurve.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

// keep (potentially) a pair of zero curves in cache
class MARKET_DLL ZeroPair
{
public:
    // default constructor
    ZeroPair();

    // useful constructor
    ZeroPair(ZeroCurveSP discZC, ZeroCurveSP growZC);

    ZeroCurveSP discZC;
    ZeroCurveSP growZC;

    // is there anybody in there?
    bool empty() const;

    const ZeroCurve* get(bool useProjectionCurve) const;
};
DECLARE_REF_COUNT(ZeroPair);

DRLIB_END_NAMESPACE
#endif
