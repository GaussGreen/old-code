
//
// C++ Interface: Trivial<IRNGGen>
//
// Description: provides trivial (identity) transformation for RNG generators.
// The template is simply to avoid the same code for IUniform, IUniformSuperCube, etc.
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef TRIVIAL_TRANSFORMATION_H
#define TRIVIAL_TRANSFORMATION_H

#include "edginc/coreConfig.hpp"
#include "edginc/DECLARESP.h"
#include "edginc/IRNGGenerator.h"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

/** Empty is an aux class from which we derive in case we don't want to have a default impelmentation of SuperCube additional information.
 */
class RNG_DLL Empty
    {}
;


/** Template to create trivial transformations
 * \param RNGGenInterface minimal Generator interface
 * \param RNGInterface desired IRNG interface tom implement
 * \param Helper a helper class (usually Empty, but later we'll encounder \class SuperCubeHelper class
    Trivial transformation extends "Helper" and implements RNGInterface
    It has create(smartPtr) method that accepts RNGGenInterfaceSP as an argument.
 */
template <class RNGGenInterface = IRNGGenerator, class RNGInterface = IRNG, class Helper = Empty >
class RNG_DLL Trivial : public virtual RNGInterface,
            public Helper // mix-in class to ease SuperCubeHooks for example
{

public:
    DECLARESP(RNGGenInterface);
    DECLARESP(RNGInterface);
    DECLARESP(Trivial);
private:
    RNGGenInterfaceSP gen;
protected:
    Trivial(RNGGenInterfaceSP _gen) : gen(_gen)
    {}

public:
    virtual IRNGGeneratorSP getGenerator() const
    {
        return gen;
    }
    virtual double fetch()
    {
        return gen->fetch();
    }
    virtual IRNGSP  clone() const
    {
        RNGGenInterfaceSP  genClone = DYNAMIC_POINTER_CAST<RNGGenInterface> (gen->clone());
        return TrivialSP(new Trivial(genClone));
    } // clone current state of the transformation
    static  TrivialSP create(RNGGenInterfaceSP _gen)
    {
        return TrivialSP(new Trivial(_gen));
    }
};


// Take any generator, and convert it into same generic RNG
typedef Trivial<IRNGGenerator, IRNG> TrivialRNG;
DECLARESP(TrivialRNG); // or we could typedef TrivialRNG::TrivialSP TrivialRNGSP, same as typedef just creates an alias.



/** Take any Uniform Generator and convert it to a Uniform RNG
*  Examples:
 
        SC_ran2SP gen = SC_ran2::create(-1L); // create generator
        TrivialUniformRNGSP rng = TrivialUniformRNG::create(gen); // transform generator into a uniform stream
        INormalRNGSP normal = SC_gasdev2::create(gen); // transfor generator into a normal stream
 
        rng->fetch(); // get one random number
        normal->fetch(); // get another advancing the same generator
        ....
 */


typedef Trivial<IUniformRNGGen, IUniformRNG> TrivialUniformRNG;
DECLARESP(TrivialUniformRNG);

class SuperCubeHelper;
typedef Trivial<ISuperCubeRNGGen, IUniformSuperCubeRNG, SuperCubeHelper>  TrivialUniformSuperCubeRNG;
DECLARESP(TrivialUniformSuperCubeRNG);

CORE_END_NAMESPACE

#endif
