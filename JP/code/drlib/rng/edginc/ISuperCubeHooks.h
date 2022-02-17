//
// C++ Interface: ISuperCubeHooks
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ISUPERCUBEHOOKS_H
#define ISUPERCUBEHOOKS_H

#include "edginc/coreConfig.hpp"

CORE_BEGIN_NAMESPACE

/** SuperCube  requires some extra functionality from RNG.
 * We capture these requirements in a separate interface.
 */

class  RNG_DLL ISuperCubeHooks
{
public:
    virtual void superCubeAdjustment() = 0; ///< gives chance to do any superCube specific adjustments after clone(). The original code flushed internal pipeline of gasdev2 transformation.
    virtual void seekToPath(int iPath) = 0; ///< resets generator to a repeatable state unique for the i-th path. Original code did that via setting seed to a function of iPath and resetting generator.
    virtual ~ISuperCubeHooks()
    {}
    virtual long getSeed() const = 0; ///< used in some superCube printings. no meaning outside a few generators.
};

CORE_END_NAMESPACE
#endif
