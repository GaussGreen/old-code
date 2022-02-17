/**
 * @file SimpleTweakNameResolver.hpp
 */

#ifndef QLIB_SimpleTweakNameResolver_H
#define QLIB_SimpleTweakNameResolver_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/TweakNameResolver.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(SimpleTweakNameResolver)

/**
 * A ITweakNameResolver which implements nameMatches() just as equality.
 */

class RISKMGR_DLL SimpleTweakNameResolver: public CObject,
                               public virtual ITweakNameResolver {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    OutputNameConstSP name;

public:

    SimpleTweakNameResolver(OutputNameConstSP name);
    ~SimpleTweakNameResolver();

    OutputNameConstSP getMarketDataName() const;

    bool nameMatches(const OutputName& name, IObjectConstSP obj);
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
