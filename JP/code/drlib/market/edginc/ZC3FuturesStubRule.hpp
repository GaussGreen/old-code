//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3FuturesStubRule.hpp
//
//   Description : Defines interface for stub rules used when adding futures 
//                 during curve bootstrapping
//
//   Author      : Richard Appleton
//
//   Date        : 18th May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_FUTURES_STUB_RULE_HPP
#define ZC3_FUTURES_STUB_RULE_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ZC3CurveInstrument.hpp"


DRLIB_BEGIN_NAMESPACE

class BadDayConvention;
class Holiday;
class IRVolBase;
class ZC3ZeroCurve;


class ZC3FuturesStubRule;
typedef refCountPtr<ZC3FuturesStubRule> ZC3FuturesStubRuleSP;


class MARKET_DLL ZC3FuturesStubRule
{
public:
    static ZC3FuturesStubRuleSP make(const string& type);

    virtual ~ZC3FuturesStubRule() {}

    ZC3RateDataArraySP ratesUsage(
        const ZC3RateDataArray& ratesData,
        const ZC3TurnArray&     turnsData,
        const DateTime*         futuresMMDate,
        const ZC3ZeroCurve&     stubCurve, /* base date, adjustments, adjDcc, hasUTurns */
        bool                    mtmAdjustment,
        const IRVolBase*        volModelIR);

protected:
    bool rangesOverlap(const DateTime& date) const;
    bool rangesTouch() const;
    bool gapExists() const;
    void retainRate();
    void removeDeleted();
    void include(const ZC3RateData* rate);

private:
    /*
     * First scan all the rates and see whether or not we have any futures and
     * also whether we have any money market rates that we may wish to exclude.
     */
    void ratesInitUsage();

    virtual void process() = 0;

protected:
    const ZC3ZeroCurve* stubCurve;
    ZC3RateDataArraySP  ratesData;
    const ZC3TurnArray* turnsData;
    bool                mtmAdjustment;
    const IRVolBase*    volModelIR;
    const DateTime*     futuresMMDate;
    DateTime            futuresStartDate;
    DateTime            futuresEndDate;
    DateTime            lastMMDate;
    DateTime            badMMDate;
    ZC3RateDataArray    deleted;
};


// do nothing special
class MARKET_DLL ZC3FuturesStubRuleNone : public ZC3FuturesStubRule
{
private:
    void process();
};


// interpolate linearly on 2 closest MM rates
class MARKET_DLL ZC3FuturesStubRuleSimple : public ZC3FuturesStubRule
{
private:
    void process();
    void futuresStubMarketInterp();
};


// interpolate using curve interpolation on 2 closest MM rates
class MARKET_DLL ZC3FuturesStubRuleInterp : public ZC3FuturesStubRule
{
private:
    void process();
    void futuresStubInterp();
};


// interpolate using flat fwds on the futures MM date
class MARKET_DLL ZC3FuturesStubRuleFlatFwds : public ZC3FuturesStubRule
{
private:
    void   process();
    void   futuresStubFlatFwds();
    double futuresStubRateFlatForwards(
        const DateTime&           date,
        double                    mmRate,
        const DayCountConvention& mmDcc,
        const ZC3RateData&        future,
        const ZC3TurnArray&       turns
        );
    double turnEffectInPeriod(
        const DateTime&           periodStartDate,
        const DateTime&           periodEndDate,
        const ZC3TurnArray&       turns,
        const DayCountConvention& dcc,
        double*                   turnYF,
        int*                      turnDays
        ) const;
};



DRLIB_END_NAMESPACE

#endif // ZC3_FUTURES_STUB_RULE_HPP
