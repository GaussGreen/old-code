/**
 * @file IRSmile2QTweak.hpp
 */

#ifndef DRLIB_IRSmile2QTweak_H
#define DRLIB_IRSmile2QTweak_H

#include "edginc/Void.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/Additive.hpp"


DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "interest rate mean reversion" as a property of market names.
 *
 */

class IRSmile2QTweak: public ScalarRiskPropertySensitivity,
                               public virtual Additive {

    static void load(CClassSP& clazz);
    virtual ScalarRiskPropertySensitivity::Deriv deriv() const;

public:

    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    static const string NAME;

    struct RISKMGR_DLL Property: public CObject {
        static CClassConstSP const TYPE;
        Property() : CObject(TYPE) {}
    
       typedef Void Qualifier;
       enum { discrete = 0 };
    private:
        static void load(CClassSP& clazz);
    };
    
    IRSmile2QTweak(double shiftSize = DEFAULT_SHIFT);
};



DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
