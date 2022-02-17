//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : StubPlacement.hpp
//
//   Description : Provides single point for defining stub placement.
//
//   Author      : Richard Appleton
//
//   Date        : 31st August 2005
//
//----------------------------------------------------------------------------

#ifndef STUB_PLACEMENT_HPP
#define STUB_PLACEMENT_HPP

#include "edginc/config.hpp"
#include "edginc/MaturityPeriod.hpp"
#include <string>


using namespace std;    // string
DRLIB_BEGIN_NAMESPACE


class MARKET_DLL StubPlacement
{
public:
    StubPlacement(const string& spec);  // spec = FRONT, BACK or AUTO
    ~StubPlacement();

    bool isLongStub() const;
    bool isShortStub() const;

    /** 
     * Returns TRUE if a stub should be at the end. The decision is
     * based on the following.
     * 1. First the default stub position is checked. If this is a back
     *    stub then return TRUE - if a front stub return FALSE.
     * 2. If the default stub position is auto then the function will
     *    return TRUE, unless there is a stub in which case FALSE is returned.
     */
    bool isEndStub(
        const MaturityPeriod& interval, 
        const DateTime&       startDate, 
        const DateTime&       endDate) const;

private:
    bool   frontStub;
    bool   backStub;
};


DRLIB_END_NAMESPACE
#endif // STUB_PLACEMENT_HPP
