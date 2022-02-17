/**
 * @file NakedBondRhoPointwise.hpp
 */

#ifndef QLIB_NakedBondRhoPointwise_H
#define QLIB_NakedBondRhoPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(NakedBondRhoPointwise)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries simultaneously.
 *
 * NakedBondRhoPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview. */


class RISKMGR_DLL NakedBondRhoPointwise: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    NakedBondRhoPointwise(const string& name, double shiftSize);
    NakedBondRhoPointwise(double shiftSize = DEFAULT_SHIFT);
    ~NakedBondRhoPointwise();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
