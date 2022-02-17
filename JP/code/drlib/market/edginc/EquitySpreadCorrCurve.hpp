//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : EquitySpreadCorrCurve.hpp
//
//   Description : correlation between equity and clean spread curve
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EQUITY_SPREAD_CORR_CURVE_H
#define EDG_EQUITY_SPREAD_CORR_CURVE_H

#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/CleanSpreadCurve.hpp"


DRLIB_BEGIN_NAMESPACE

/** A equity spread corr curve is associated with a clean spread vol curve */
class MARKET_DLL EquitySpreadCorrCurve: public CObject
{
public:
// not done yet
    static CClassConstSP const TYPE;
    friend class EquitySpreadCorrCurveHelper;
    friend class EquitySpreadCorrCurveAddin;
	friend class DDEModule;

    EquitySpreadCorrCurve();

    EquitySpreadCorrCurve(const DateTime&     valueDate,
                     const ExpiryArray*  expiries,
                     const DoubleArray*  corrs);

    DateTimeArray getDates() const;
    DoubleArray *getCorrs() const;

protected:
    DateTime        valueDate;
    ExpiryArraySP   expiries;
    CDoubleArraySP  corrs;
};

typedef smartConstPtr<EquitySpreadCorrCurve> EquitySpreadCorrCurveConstSP;
typedef smartPtr<EquitySpreadCorrCurve>      EquitySpreadCorrCurveSP;

DRLIB_END_NAMESPACE
#endif
