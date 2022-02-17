//
// C++ Interface: IRNG
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef IRNG_H
#define IRNG_H

#include "edginc/coreConfig.hpp"
#include "edginc/DECLARESP.h"
#include "edginc/IRNGGenerator.h"

#include "edginc/ISuperCubeHooks.h"

CORE_BEGIN_NAMESPACE

class IRNG;
DECLARESP(IRNG);

/** Random Number Generator Object:
    IRNG is simply a transformation atop of a lower level IRNGGenerator object.
    It will also require clone() functionality.
 */
class  RNG_DLL IRNG
{
public:
    virtual double fetch() = 0; ///< get next tranformed RNG
    virtual IRNGGeneratorSP getGenerator() const = 0; ///< this allows different transformations to be applied to the same generator
    virtual IRNGSP  clone() const = 0; ///< clone current state of the transformation
    virtual ~IRNG()
    {}  ///< we have smart pointers of this type
    double operator () (void)
    {
        return fetch();
    }  ///< IRNG is a Generator in STL sense
};

/** ISuperCubeRNG is the interface for all transformations that support SuperCube's extended functionality
*/

class  RNG_DLL ISuperCubeRNG : public virtual IRNG,
            public virtual ISuperCubeHooks
{
};

DECLARESP(ISuperCubeRNG);

/** All transformations that return Normal (Gaussian) RNGs */
class RNG_DLL INormalRNG : public virtual IRNG
{
};

DECLARESP(INormalRNG);

/** All transformations that return Normal distribution and support SuperCube additional functionality.
 */
class  RNG_DLL INormalSuperCubeRNG : public virtual INormalRNG,
            public virtual ISuperCubeRNG
{
};

DECLARESP(INormalSuperCubeRNG);

/** IUniformRNG interface is for transformations that return uniform distribution
 */

class  RNG_DLL IUniformRNG : public virtual IRNG
    {}
;

DECLARESP(IUniformRNG);
/** IUniformSuperCubeRNG  is for transformations that return uniform distribution and support SuperCube additional functionality
 */
class  RNG_DLL IUniformSuperCubeRNG : public virtual IUniformRNG,
            public virtual ISuperCubeRNG
    {}
;

DECLARESP(IUniformSuperCubeRNG);

CORE_END_NAMESPACE
#endif
