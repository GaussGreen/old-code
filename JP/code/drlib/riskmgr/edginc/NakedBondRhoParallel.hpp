/**
 * @file NakedBondRhoParallel.hpp
 */

#ifndef QLIB_NakedBondRhoParallel_H
#define QLIB_NakedBondRhoParallel_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(NakedBondRhoParallel)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries simultaneously.
 *
 * NakedBondRhoParallel has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL NakedBondRhoParallel: public ScalarRiskPropertySensitivity,
                    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    NakedBondRhoParallel(double shiftSize = DEFAULT_SHIFT);
    ~NakedBondRhoParallel();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
