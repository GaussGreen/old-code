//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3FuturesGapRule.hpp
//
//   Description : Defines interface for gap rules used when adding swaps 
//                 during curve bootstrapping
//
//   Author      : Richard Appleton
//
//   Date        : 17th May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_FUTURES_GAP_RULE_HPP
#define ZC3_FUTURES_GAP_RULE_HPP

#include "edginc/config.hpp"


DRLIB_BEGIN_NAMESPACE
class DateTime;
class BadDayConvention;
class Holiday;

class ZC3SwapData;
class ZC3ZeroCurve;


class MARKET_DLL ZC3FuturesGapRule
{
public:
    ZC3FuturesGapRule(const string& type, int gapDays);

    bool ignore(
        const ZC3SwapData&      swapData, 
        ZC3ZeroCurve&           szcCurrent,
        const ZC3ZeroCurve&     target,
        const BadDayConvention& badDayConv,
        const Holiday&          holidays) const;

    DateTime adjust(const DateTime& lastCashDate) const;

private:
    int  gapDays;
    bool extend;
    bool none;
};


DRLIB_END_NAMESPACE

#endif // ZC3_FUTURES_GAP_RULE_HPP
