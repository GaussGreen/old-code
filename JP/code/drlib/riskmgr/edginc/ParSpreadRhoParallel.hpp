/**
 * @file ParSpreadRhoParallel.hpp
 */

#ifndef QLIB_ParSpreadRhoParallel_H
#define QLIB_ParSpreadRhoParallel_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ParSpreadRhoParallel)

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * ParSpreadRhoParallel has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class RISKMGR_DLL ParSpreadRhoParallel: public ScalarRiskPropertySensitivity,
                            public virtual Additive,
                            public virtual CreditTweak {

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    ParSpreadRhoParallel(double shiftSize = DEFAULT_SHIFT,
                         const string& packetName = NAME);
    ~ParSpreadRhoParallel();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
