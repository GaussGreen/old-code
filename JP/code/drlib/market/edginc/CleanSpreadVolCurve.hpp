//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : CleanSpreadVolCurve.hpp
//
//   Description : vol curve for clean spread curve
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CLEAN_SPREAD_VOL_CURVE_H
#define EDG_CLEAN_SPREAD_VOL_CURVE_H

#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TimeMetric);

/** A clean spread vol curve is associated with a clean spread curve */
class MARKET_DLL CleanSpreadVolCurve: public CObject
{
public:
// not done yet
    static CClassConstSP const TYPE;
    friend class CleanSpreadVolCurveHelper;
    friend class CleanSpreadVolCurveAddin;
    friend class DDEModule;

    CleanSpreadVolCurve();

    CleanSpreadVolCurve(const DateTime&     valueDate,
                        const ExpiryArray*  expiries,
                        const DoubleArray*  vols,
                        double meanReversion,
                        const TimeMetric *timeMetric);

    DateTimeArray getDates() const;
    DoubleArray * getVols() const;
    double  getMeanReversion() const;
    const TimeMetric *getTimeMetric() const;

    // calc the variance of the spot clean spread. dates must be after value date and is increasing array
    void getSpotSprdVars(const DateTimeArray &dates, DoubleArray &vars) const;

protected:
    DateTime        valueDate;
    ExpiryArraySP   expiries;
    CDoubleArraySP  vols;
    double          meanReversion;
    TimeMetricSP    timeMetric;
};

typedef smartConstPtr<CleanSpreadVolCurve> CleanSpreadVolCurveConstSP;
typedef smartPtr<CleanSpreadVolCurve>      CleanSpreadVolCurveSP;

DRLIB_END_NAMESPACE
#endif
