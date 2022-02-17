/**
 * @file OptionStrikeRhoPointwise.hpp
 */

#ifndef QLIB_OptionStrikeRhoPointwise_H
#define QLIB_OptionStrikeRhoPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OptionStrikeRhoPointwise)

/**
 * A greek calculated by tweaking each market name's volatility at all its
 * defined expiries simultaneously.
 *
 * OptionStrikeRhoPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview. */

class RISKMGR_DLL OptionStrikeRhoPointwise: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    OptionStrikeRhoPointwise(const string& name, double shiftSize);
    OptionStrikeRhoPointwise(double shiftSize = DEFAULT_SHIFT);
    ~OptionStrikeRhoPointwise();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
