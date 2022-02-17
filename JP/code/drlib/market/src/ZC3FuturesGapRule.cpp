//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3FuturesGapRule.cpp
//
//   Description : Implementation for futures gap rules.
//
//   Author      : Richard Appleton
//
//   Date        : 17th May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3FuturesGapRule.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/ZC3CurveSegment.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CompoundBasis.hpp"


DRLIB_BEGIN_NAMESPACE


ZC3FuturesGapRule::ZC3FuturesGapRule(const string& type, int pGapDays)
  : gapDays(pGapDays)
{
    if (CString::equalsIgnoreCase(type, "N")
     || CString::equalsIgnoreCase(type, "None", 4))
    {
        none = true;
        extend = false;
    }
    else if (CString::equalsIgnoreCase(type, "I")
     || CString::equalsIgnoreCase(type, "Ignore", 6))
    {
        none = false;
        extend = false;
    }
    else if (CString::equalsIgnoreCase(type, "E")
     || CString::equalsIgnoreCase(type, "Extend", 6))
    {
        none = false;
        extend = true;
    }
    else
    {
        string msg = Format::toString("Unknown futures gap rule %s", type.c_str());
        throw ModelException("ZC3FuturesGapRule::ZC3FuturesGapRule", msg);
    }
}


// ALIB: szcswaps.c#718
bool ZC3FuturesGapRule::ignore(
        const ZC3SwapData&      swapData, 
        ZC3ZeroCurve&           szcCurrent,
        const ZC3ZeroCurve&     target,
        const BadDayConvention& badDayConv,
        const Holiday&          holidays) const
{
    if (!none && swapData.getEndDate() > target.getBaseDate())
    {
        DateTime lastIgnored = target.endDate().rollDate(gapDays);
        DateTime adjusted = badDayConv.adjust(swapData.getEndDate(), &holidays);

        if (adjusted <= lastIgnored && adjusted > target.endDate())
        {
            if (extend)
            {
                szcCurrent.extendLastRate(target.endDate(), adjusted);
            }

            return true;
        }
    }

    return false;
}


// ALIB: zcrecurse.c#623
DateTime ZC3FuturesGapRule::adjust(const DateTime& lastCashDate) const
{
    if (!none && !extend)   // ie. 'Ignore'
    {
        return lastCashDate.rollDate(gapDays);
    }

    return lastCashDate;
}


DRLIB_END_NAMESPACE
