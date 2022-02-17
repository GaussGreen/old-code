//
// C++ Interface: GSLRNG
//
// Description: Some of the transformations implemented by GSL for
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GSLRNG_H
#define GSLRNG_H

#include "edginc/coreConfig.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <functional>
#include <iostream>
#include "edginc/IRNG.h"
#include "edginc/SuperCubeHelper.h"
#include "edginc/DECLARESP.h"
#include "edginc/GSLRNGGen.h"
#include "edginc/TrivialTransform.h"
CORE_BEGIN_NAMESPACE


/** This template is to automate creation of transformations (RNG) based on GSL functionality.
 See http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distributions.html
 for the list of available transformations. We instantiate only 3: uniform (i.e. transformation is the identity), normal and fast normal (based on the Ziggurat method).
 One technical twist is that we need to wrap functions into classes.
 
 Template takes 2 or 3  arguments.
 The third argument is to facilitate instantiation compatible with SuperCube.
 The first argument specifies the desired transformation (wrapped in a class).
 The second argument specifies the desired type of the RNG generator (so we can be more specific than IRNG, ex: INormalRNG)
 TODO: do it via traits;
*/
template <class func, class RNGInterface = IRNG, class Helper1 = Empty>
class  RNG_DLL GSLTransformation : virtual public RNGInterface,
            public Helper1 // mainly to reuse SuperCubeHelper
{
    DECLARESP(RNGInterface);
    DECLARESP( GSLTransformation );
    GSLSuperCubeGenSP gen; // FIXME: templatize this

protected:
    GSLTransformation(GSLSuperCubeGenSP _gen) : gen(_gen)
    {}

public:
    virtual double fetch()
    { // get next tranformed RNG
        return func::eval(gen->getGen());
    }
    virtual IRNGSP  clone() const
    { // clone current state of the transformation
        return GSLTransformationSP(new GSLTransformation(DYNAMIC_POINTER_CAST<GSLSuperCubeGen>(gen->clone())));
    }
    virtual IRNGGeneratorSP getGenerator() const
    {
        return gen;
    }
    static GSLTransformationSP create(GSLSuperCubeGenSP _gen)
    {
        GSLTransformationSP rng =  GSLTransformationSP(new GSLTransformation(_gen));
        return rng;
    }
};


class  RNG_DLL GSL_trivial
{
public:
    static double eval (const gsl_rng * gen)
    {
        return gsl_rng_uniform(gen);
    }
};

class  RNG_DLL  GSL_ran_ugaussian_ratio_method
{
public:
    static double eval (const gsl_rng * gen)
    {
        return gsl_ran_ugaussian_ratio_method(gen);
    }
};
class  RNG_DLL  GSL_ran_ugaussian
{
public:
    static double eval (const gsl_rng * gen)
    {
        return gsl_ran_ugaussian(gen);
    }
};
// Some typedefs
// so one can use it like GSLNormalRNG::create(GSLRNGGen::create(gsl_rng_ranlux389))
typedef GSLTransformation<GSL_trivial, IUniformRNG> GSLUniformRNG;
DECLARESP( GSLUniformRNG);

typedef GSLTransformation<GSL_ran_ugaussian_ratio_method, INormalRNG> GSLZigguratNormalRNG;
DECLARESP( GSLZigguratNormalRNG);

typedef GSLTransformation<GSL_ran_ugaussian, INormalRNG> GSLNormalRNG;
DECLARESP( GSLNormalRNG);


/////////////////////////// SuperCube analogs



typedef GSLTransformation<GSL_trivial,
IUniformSuperCubeRNG,
SuperCubeHelper> GSLUniformSuperCubeRNG;
DECLARESP( GSLUniformSuperCubeRNG);

typedef GSLTransformation<GSL_ran_ugaussian_ratio_method,
INormalSuperCubeRNG,
SuperCubeHelper> GSLZigguratNormalSuperCubeRNG;
DECLARESP(    GSLZigguratNormalSuperCubeRNG);

typedef GSLTransformation<GSL_ran_ugaussian,
INormalSuperCubeRNG,
SuperCubeHelper> GSLNormalSuperCubeRNG;

DECLARESP( GSLNormalSuperCubeRNG);

#if 0
static IRNGSP rngs [] = {
                            GSLUniformSuperCubeRNG::create(GSLRNGGen::create()),
                            GSLUniformSuperCubeRNG::create(GSLRNGGen::create()),
                            GSLZigguratNormalSuperCubeRNG::create(GSLRNGGen::create()),
                            GSLNormalSuperCubeRNG::create(GSLRNGGen::create()) };
#endif

CORE_END_NAMESPACE

#endif
