//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : CleanSpreadCurve.hpp
//
//   Description : a clean spread spread curve (aka. hazard rate curve)
//
//   Author      : André Segger
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CLEAN_SPREAD_CURVE_H
#define EDG_CLEAN_SPREAD_CURVE_H

#include "edginc/Class.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/CreditSpreadRhoPointwise.hpp"
#include "edginc/AdjCreditSpreadRhoParallel.hpp"
#include "edginc/AdjCreditSpreadRhoPointwise.hpp"
#include "edginc/CreditSpreadLevel.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/DefaultRates.hpp"



DRLIB_BEGIN_NAMESPACE

/** A clean spread curve should be derived from a CDS par spread curve - it should not be an external class
    to be used as part of an interface, at least in the near future */

FORWARD_DECLARE(CleanSpreadCurve);

class MARKET_DLL CleanSpreadCurve: public CObject,
                        public ICreditCurve
{
public:
// not done yet
    static CClassConstSP const TYPE;
    friend class CleanSpreadCurveHelper;
    friend class CleanSpreadCurveAddin;
	friend class DDEModule;

    CleanSpreadCurve();

    CleanSpreadCurve(const DateTime&     valueDate,
                     const ExpiryArray*  expiries,
                     const DoubleArray*  rates);

    /** returns the hazard rate between two dates conditioned on no default up to the start date */
    double getCleanSpread(const DateTime& startDate,
                          const DateTime& endDate) const;

    /** returns the hazard rate between two dates conditioned on no default up to the start date */
    double getCleanSpread(const DateTime& endDate) const;

    /** returns the hazard rate between two dates conditioned on no default up to the start date */
    double getCleanSpread(const DateTime& endDate, const int basis) const;

    /** returns the default probability between two dates conditioned on no default up to the start date */
    double getDefaultProb(const DateTime& startDate,
                          const DateTime& endDate) const;

    /** returns the default probability between two dates conditioned on no default up to the start date */
    double getDefaultProb(const DateTime& endDate) const;

    /** returns 1 - the default probability */
    double getDefaultPV(const DateTime& startDate,
                        const DateTime& endDate) const;

    /** returns 1 - the default probability */
    double getDefaultPV(const DateTime& endDate) const;

    double getCurrentSpread(const DateTime& valueDate,
                            const DateTime& maturityDate) const;

    double getCurrentSpread(const DateTime& valueDate,
                            const DateTime& maturityDate,
                            const BadDayConvention* bdc,
                            const Holiday* hols) const;

	static CleanSpreadCurveSP cleanSpreadCurveFromDF(
		DateTime valueDate,
		DateTimeArray dates,
		DoubleArray   dfs);

	static CleanSpreadCurveSP convertFromDefaultRates(
				DefaultRatesConstSP curveIn);

    DateTimeArray getDates() const;

	DateTime getValueDate() const
	{ return valueDate; }

	/**
	 * Create risky yield curve.
	 */
	virtual IYieldCurveSP makeRiskyCurve(
		const IYieldCurve& yc, 
		const DateTime*    maturityDate=NULL) const;

protected:
    DateTime        valueDate;       // today
    ExpiryArraySP   expiries;
    CDoubleArraySP  cleanSpreads;
};

//typedef smartConstPtr<CleanSpreadCurve> CleanSpreadCurveConstSP;
//typedef smartPtr<CleanSpreadCurve>      CleanSpreadCurveSP;


DECLARE(CleanSpreadCurveArray);

DRLIB_END_NAMESPACE
#endif
