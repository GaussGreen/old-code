//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : StubPlacement.cpp
//
//   Description : Provides single point for defining stub placement.
//
//   Author      : Richard Appleton
//
//   Date        : 31st August 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/SwapTool.hpp"

DRLIB_BEGIN_NAMESPACE


StubPlacement::StubPlacement(const string& spec)
{
    frontStub = CString::equalsIgnoreCase(spec,"Front");
    backStub = CString::equalsIgnoreCase(spec,"Back");
}


StubPlacement::~StubPlacement()
{
}


// for now only short stubs supported
bool StubPlacement::isShortStub() const
{
    return true;
}


bool StubPlacement::isLongStub() const
{
    return !isShortStub();
}


// ALIB: stub.c#941
bool StubPlacement::isEndStub(
    const MaturityPeriod& tenor,
    const DateTime&       startDate,
    const DateTime&       endDate) const
{
    if (frontStub)
    {
        return false;
    }
    else if (backStub)
    {
        return true;
    }
    else
    {
        string ivl;
        int numIvls;
        tenor.decompose(numIvls, ivl);

        int numIntervals;
        int extraDays;
        SwapTool::countDates(startDate, endDate, numIvls, ivl, &numIntervals, &extraDays);
        return !(extraDays > 0);
    }
}


DRLIB_END_NAMESPACE
