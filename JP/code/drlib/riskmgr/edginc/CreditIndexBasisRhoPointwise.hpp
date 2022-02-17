/**
 * @file ParSpreadRhoParallel.hpp
 */

#ifndef QLIB_CreditIndexBasisRhoPointwise_H
#define QLIB_CreditIndexBasisRhoPointwise_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CreditIndexBasisRhoPointwise)

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * This implementation uses the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL CreditIndexBasisRhoPointwise: 
    public PerNameRiskPropertySensitivity<ExpiryWindow>,
    public virtual Additive,
    public virtual CreditTweak {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

    DoubleArrayConstSP adaptiveTweakSizes(MultiTweakGroupConstSP tweakGroup, 
                                          OutputNameConstSP name,
                                          ExpiryArrayConstSP expiries) const;

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const string NAME;

    CreditIndexBasisRhoPointwise(double shiftSize = DEFAULT_SHIFT);
    ~CreditIndexBasisRhoPointwise();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
