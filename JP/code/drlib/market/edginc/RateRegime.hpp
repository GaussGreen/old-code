#ifndef EDR_RATEREGIME_HPP
#define EDR_RATEREGIME_HPP

#include "edginc/RegimeFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class RateRegime: public RegimeFactor {
public:
    static CClassConstSP const TYPE;
	friend class RateRegimeAddin;
	friend class MultiRegimeFactor;

	void validatePop2Object(); 

	double calcRisklessBond(double tradYear) const;

	double calcYieldCurve(double tradYear) const;

	DoubleArray calcYieldCurve(DoubleArray tradYears) const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    RateRegime();

    static IObject* defaultCtor();
};
typedef smartPtr<RateRegime> RateRegimeSP;
typedef smartConstPtr<RateRegime> RateRegimeConstSP;

DRLIB_END_NAMESPACE

#endif
