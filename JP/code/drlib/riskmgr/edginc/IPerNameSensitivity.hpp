/**
 * @file IPerNameSensitivity.hpp
 */

#ifndef QLIB_IPerNameSensitivity_H
#define QLIB_IPerNameSensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/TweakTypeID.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IPerNameSensitivity)
FORWARD_DECLARE(OutputName)

/**
 * Interface for sensitivities which only apply to certain names.
 *
 * I introduced this in order to provide interim support some model/instrument
 * mechanisms based on interrogation of Control::getCurrentSensitivity(),
 * without having to inherit PerNameRiskPropertySensitivity from ScalarShift.
 * Next step is to clean up those mechanisms, at which point this can be
 * dropped.
 */

class RISKMGR_DLL IPerNameSensitivity: public virtual IObject,
                           public virtual ITweakTypeID  {

public:

    static CClassConstSP const TYPE;

    IPerNameSensitivity();
    ~IPerNameSensitivity();

    virtual OutputNameConstSP getMarketDataName() const = 0;

    virtual OutputNameArrayConstSP allNames(const IObject* object) const = 0;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
