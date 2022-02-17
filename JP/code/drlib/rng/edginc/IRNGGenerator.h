//
// C++ Interface: IRNGGenerator
//
// Description: Base classes for primitive generators (uniform, grid, etc.)
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef IRNGGENERATOR_H
#define IRNGGENERATOR_H

#include "edginc/coreConfig.hpp"
#include "edginc/DECLARESP.h"
#include "edginc/ISuperCubeHooks.h"

#include <cassert>

CORE_BEGIN_NAMESPACE

/** Base class for all RNG Generators.
 
    The users shouldn't use generators directly: generators are meant to be piped through a transformation class to obtain RNG. (Transformation can be identity (see TrivialRNG class) as well. The generators can produce some numbers (via fetch()) and be cloned (that supposedly copies internal state).
 
 Generators can have a finer type. So far, most generators produce uniform numbers, but other types are possible (ex: LowDisc. unlikely to pass uniformity test; Generators based on files of pre-recorded numbers can give another example of a non-uniform generator.
 
 Where do we get generators?
 - Our own code (see SC_ran2Gen), mainly a few generators from the Numerical Rec.
 - Libraries:
  a. GSL (GNU Scientific Library) -- a huge selection of generators. Wrappings available in <GSLGen.h>
  b. Intel's MKL  -- Intel specific library; supposedly very fast; wrappings not available yet
  c. SPRNG -- works very fast on 64 bit platforms; wrappings not available yet
  d. QuantLib -- C++ libary; bindings not available yet
 
 Currently, generators that are to be used with SuperCube have to support a bit more of functionality (like seekToPath). This functionality is split into ISuperCubeHooks interface and (uniform) generators that implement it are marked with ISuperCubeRNGGen interface.
 
    The distinction of generators and transformation follows design of other libraries and allows one to add additional type information to both generators/transformations. For example, it may be a mistake to use rejection based transformation on a low-disc. sequence. Types will allow this type of constraints.
 
*/
class IRNGGenerator;
DECLARESP(IRNGGenerator);

class  RNG_DLL IRNGGenerator
{
public:

    virtual double fetch() = 0; ///< Fetch next random number
    virtual IRNGGeneratorSP  clone() const = 0; ///< Clone the existing state
    virtual ~IRNGGenerator()
    {} ///< Allow smart pointers of this type
    virtual void debug()
    {} ///< name says it all
    double operator() (void)
    {
        return fetch();
    } ///< allow STL algoriths that use Generator objects
};

/** Interface for uniform generators */
class  RNG_DLL IUniformRNGGen : public virtual IRNGGenerator
    {}
;
DECLARESP(IUniformRNGGen);

/** SuperCube expects some more functionality from a generator, fact captured by the ISuperCubeRNGGen interface
*/
class  RNG_DLL ISuperCubeRNGGen : public virtual IUniformRNGGen,
            public virtual ISuperCubeHooks
    {}
;

DECLARESP(ISuperCubeRNGGen);

CORE_END_NAMESPACE

#endif

