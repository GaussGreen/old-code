/**
 * file LossMappingSens.hpp
   Author   : Juan Carlos Porras      Dec 06
   Definitions for std dev tweak to loss mapping
 */

#ifndef QLIB_LossMappingSens_H
#define QLIB_LossMappingSens_H

#include "edginc/Additive.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(LossMappingSens)

class RISKMGR_DLL LossMappingSens: public ScalarRiskPropertySensitivity,
                            public virtual Additive,
							public virtual CreditTweak	{

    static void load(CClassSP& clazz);
    Deriv deriv() const;

public:

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    LossMappingSens(double shiftSize = DEFAULT_SHIFT,
                    const string& packetName = NAME);
    ~LossMappingSens();
};

DRLIB_END_NAMESPACE


#endif
