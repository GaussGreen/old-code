//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadCurve.hpp
//
//   Description : a (credit) spread curve
//
//   Author      : Ning Shen
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_SPREAD_CURVE_HPP
#define CREDIT_SPREAD_CURVE_HPP

#include "edginc/Class.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/CreditSpreadRhoPointwise.hpp"
#include "edginc/AdjCreditSpreadRhoParallel.hpp"
#include "edginc/AdjCreditSpreadRhoPointwise.hpp"
#include "edginc/CreditSpreadLevel.hpp"
#include "edginc/CreditSpreadParallelShift.hpp"
#include "edginc/CreditSpreadPropShift.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/CreditCurve.hpp"


DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CreditSpreadCurve : public MarketObject,
                          virtual public ICreditCurve,
                          virtual public CreditSpreadRhoParallel::IRestorableShift,
                          virtual public CreditSpreadRhoPointwise::IRestorableShift,
                          virtual public AdjCreditSpreadRhoParallel::IRestorableShift,
                          virtual public AdjCreditSpreadRhoPointwise::IRestorableShift,
                          virtual public CreditSpreadLevel::IShift,
                          virtual public CreditSpreadParallelShift::IShift,
                          virtual public CreditSpreadPropShift::IShift

{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    friend class CashSwapCurve;

    CreditSpreadCurve();

    CreditSpreadCurve(double spread);   // flat curve

    CreditSpreadCurve(const string&      credName,
                      const ExpiryArray* expiries,
                      const DoubleArray& rates,
					  double recoveryPct=0);

    virtual void validatePop2Object();

    virtual string getName() const;

    double getCurrentSpread(const DateTime& valueDate,
                            const DateTime& maturityDate)const;

    double getCurrentSpread(const DateTime& valueDate,
                            const DateTime& maturityDate,
                            const BadDayConvention* bdc,
                            const Holiday* hols) const;

    CreditSpreadCurve* getScaledCurve(double scalingFactor) const;

    void addSpread(double spreadOffset);

	/**
	 * Create risky yield curve.
	 */
	virtual IYieldCurveSP makeRiskyCurve(
		const IYieldCurve& yc, 
		const DateTime*    maturityDate=NULL) const;

    /** CreditSpreadRhoParallel sensitivity support */
    string sensName(CreditSpreadRhoParallel* shift) const;
    bool sensShift(CreditSpreadRhoParallel* shift);
    void sensRestore(CreditSpreadRhoParallel* shift);

    /** CreditSpreadRhoPointwise sensitivity support */
    string sensName(CreditSpreadRhoPointwise* shift) const;
    ExpiryArrayConstSP sensExpiries(CreditSpreadRhoPointwise* shift) const;
    bool sensShift(CreditSpreadRhoPointwise* shift);
    void sensRestore(CreditSpreadRhoPointwise* shift);

    /** AdjCreditSpreadRhoParallel support */
    string sensName(AdjCreditSpreadRhoParallel* shift) const;
    bool sensShift(AdjCreditSpreadRhoParallel* shift);
    void sensRestore(AdjCreditSpreadRhoParallel* shift);

    /** AdjCreditSpreadRhoPointwise support */
    string sensName(AdjCreditSpreadRhoPointwise* shift) const;
    ExpiryArrayConstSP sensExpiries(AdjCreditSpreadRhoPointwise* shift) const;
    bool sensShift(AdjCreditSpreadRhoPointwise* shift);
    void sensRestore(AdjCreditSpreadRhoPointwise* shift);

    // CreditSpreadLevel support
    virtual string sensName(CreditSpreadLevel* shift) const;
    virtual bool sensShift(CreditSpreadLevel* shift);

    // CreditSpreadParallelShift scenario support
    virtual string sensName(CreditSpreadParallelShift* shift) const;
    virtual bool sensShift(CreditSpreadParallelShift* shift);

    // CreditSpreadPropShift support
    virtual string sensName(CreditSpreadPropShift* shift) const;
    virtual bool sensShift(CreditSpreadPropShift* shift);

    ExpiryArrayConstSP getExpiries()const;

    DateTimeArray getExpiryDates(const DateTime& today)const;

    double recovery() const;

protected:
    CreditSpreadCurve(const CClassConstSP& clazz);
    // may add bench mark lables if market rates have it
    string          name;
    ExpiryArraySP   expiries;
    CDoubleArray    spreads;

	double recoveryPct;
};

typedef smartConstPtr<CreditSpreadCurve> CreditSpreadCurveConstSP;
typedef smartPtr<CreditSpreadCurve> CreditSpreadCurveSP;
#ifndef QLIB_CREDITSPREADCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CreditSpreadCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CreditSpreadCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CreditSpreadCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CreditSpreadCurve>);
#endif
// support for wrapper class
typedef MarketWrapper<CreditSpreadCurve> CreditSpreadCurveWrapper;
#ifndef QLIB_CREDITSPREADCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CreditSpreadCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CreditSpreadCurve>);
#endif

DRLIB_END_NAMESPACE
#endif
