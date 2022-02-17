#ifndef EDR_VolRegime_HPP
#define EDR_VolRegime_HPP

#include "edginc/RegimeFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class VolRegime: public RegimeFactor {
public:
    static CClassConstSP const TYPE;
	friend class VolRegimeAddin;
	friend class MultiRegimeFactor;

	void validatePop2Object(); 

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VolRegime();

    static IObject* defaultCtor();
};
typedef smartPtr<VolRegime> VolRegimeSP;
typedef smartConstPtr<VolRegime> VolRegimeConstSP;

DRLIB_END_NAMESPACE

#endif
