//
// C++ Interface: SuperCubeHelper
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SUPERCUBEHELPER
#define SUPERCUBEHELPER

#include "edginc/coreConfig.hpp"
#include "edginc/ISuperCubeHooks.h"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

// Mix-in Helper class for transformations that redirect all SuperCube hooks to generators.
class RNG_DLL SuperCubeHelper : public virtual ISuperCubeHooks,
            public virtual IRNG
{
    ISuperCubeHooks * getHook() const
    {
        return dynamic_cast<ISuperCubeHooks*>(getGenerator().get());
    }
public:
    virtual void superCubeAdjustment()
    { // FIXME: delete after SupCub is finished
        getHook()->superCubeAdjustment();
    }

    virtual void seekToPath(int iPath)
    { // using current state, seek to path -- SuperCube
        getHook()->seekToPath(iPath);
    }
    virtual long getSeed() const
    {
        return getHook()->getSeed();
    }
};

CORE_END_NAMESPACE

#endif
