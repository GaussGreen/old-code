//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CREDIT_CURVE_HPP
#define EDG_CREDIT_CURVE_HPP

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/YieldCurve.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ICreditCurve: virtual public IObject {
public:

     virtual double getCurrentSpread(const DateTime& valueDate,
                                     const DateTime& maturityDate) const = 0;

     virtual double getCurrentSpread(const DateTime& valueDate,
                                     const DateTime& maturityDate,
                                     const BadDayConvention* bdc,
                                     const Holiday* hols) const = 0;

    /**
     * An object for representing getCurrentSpread() in the Results: either
     * a CDouble, or an Untweakable if the calculation failed
     *
     * Now that QuantoCDSParSpreads::getCurrentSpread() is fully implemented,
     * there is a risk of it failing, so various places where
     * getCurrentSpread() is invoked need guarding.  (Maybe rather than use
     * this specific "convenience method", it would be better to make a more
     * general macro?)
     */

    IObjectSP getCurrentSpreadOrUntweakable(const DateTime& valueDate,
                                            const DateTime& maturity) const;

	/**
	 * Create risky yield curve.
	 */
	virtual IYieldCurveSP makeRiskyCurve(
		const IYieldCurve& yc, 
		const DateTime*    maturityDate=NULL) const = 0;

    static CClassConstSP const TYPE; 

private:
    static void load(CClassSP& clazz);
};

// typedef for smart pointers to ISensitiveStrikes
typedef smartConstPtr<ICreditCurve> ICreditCurveConstSP;
typedef smartPtr<ICreditCurve>      ICreditCurveSP;
#ifndef QLIB_CREDITCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ICreditCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ICreditCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ICreditCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ICreditCurve>);
#endif

// support for wrapper class
typedef MarketWrapper<ICreditCurve> CreditCurveWrapper;
#ifndef QLIB_CREDITCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<ICreditCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<ICreditCurve>);
#endif

DRLIB_END_NAMESPACE

#endif // EDG_CREDIT_CURVE_HPP
