//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : SpreadEquityFunc.hpp
//
//   Description : spread as function of equity. base class
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#ifndef EDG_SPREAD_EQUITY_FUNC_H
#define EDG_SPREAD_EQUITY_FUNC_H

#include "edginc/Class.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class DateTime;
FORWARD_DECLARE(TimeMetric);

/** A base class spread as function of equity. */
class MARKET_DLL SpreadEquityFunc: public CObject
{
public:
    // return continuous compounded spread
    // allow user specific daycount calculated dt as input 
    inline double getSpreadCC(double spot,
                              const DateTime &start,
                              const DateTime &end,
                              double dt=0) const
    {
        double spread;
        getSpreadCC(start, end, 1, &spot, &spread, dt);
        return spread;
    }
                                                
    virtual void getSpreadCC(const DateTime &start,
                             const DateTime &end,
                             int     nbSpot, 
                             const double *spots,
                             double *spreads,
                             double dt=0) const =0;

    virtual void setTimeMetric(const TimeMetric *metric);

    static CClassConstSP const TYPE;
        
    static void load(CClassSP& clazz);

protected:
    SpreadEquityFunc(const CClassConstSP& clazz);

    TimeMetricSP metric;
};

typedef smartConstPtr<SpreadEquityFunc> SpreadEquityFuncConstSP;
typedef smartPtr<SpreadEquityFunc>      SpreadEquityFuncSP;

DRLIB_END_NAMESPACE
#endif
