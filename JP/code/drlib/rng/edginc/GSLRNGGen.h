//
// C++ Interface: GSLQRNG
//
// Description: wrap up uniform RNG generators exported by GSL
// Clients will have to link with libgsl and libgslcblas
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef GSLQRNGGEN_H
#define GSLQRNGGEN_H

#include "edginc/coreConfig.hpp"
#include <gsl/gsl_rng.h>

#include "edginc/IRNGGenerator.h"
#include "edginc/DECLARESP.h"

CORE_BEGIN_NAMESPACE

class IGSLGen
{
public:
    virtual gsl_rng * getGen() const = 0; // return GSL generator

    // All GSL generators support notion of a seed; setting it resets generator into a state that is fully determined by the seed. The "0" is generator specific ("default") seed; see <gsl_rng.h> doc
    virtual void   setSeed(unsigned long seed)
    {
        gsl_rng_set(getGen(), seed); //
    }
};

class GSLRNGGen;
DECLARESP(GSLRNGGen);
/** For GSL, see documentation at
    http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html
 http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
*/
class  RNG_DLL GSLRNGGen : public virtual IUniformRNGGen,
            public virtual IGSLGen
{
public:
    // Fetch one uniformly distributed number
    virtual double fetch()
    {
        return gsl_rng_uniform(gen);
    }


    virtual IRNGGeneratorSP clone() const
    {
        gsl_rng * cl = gsl_rng_clone(gen);
        return GSLRNGGenSP(new GSLRNGGen(cl));
    }


    /** Recommended ones:
    gsl_rng_mt19937 Mersenne Twister
    gsl_rng_ranlux  “luxury random numbers”
    ....
    see http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
    */
    static GSLRNGGenSP create(const gsl_rng_type *_type = NULL, unsigned long _seed = 0)
    {
        return GSLRNGGenSP(new GSLRNGGen(_type, _seed));
    }
    virtual ~GSLRNGGen()
    {
        gsl_rng_free(gen);
        gen = NULL;
    }
    gsl_rng * getGen() const
    {
        return gen;
    }

protected:
    // Create object using the current state of the generator
    GSLRNGGen(gsl_rng * _gen) : gen(_gen)
    {}
    // Note that you can tweak RNG type and the initial SEED via env. variables (see web doc above)
    GSLRNGGen(const gsl_rng_type *_type = NULL, unsigned long _seed = 0) : gen(NULL)
    {
        static bool once = (gsl_rng_env_setup(), true); // make sure gsl_rng_env_setup is called >0 times
        gen = gsl_rng_alloc(_type ? _type : gsl_rng_default);
        if (_seed != 0L)  // according to GSL, 0L is a special value, meaning default one.
            setSeed(_seed);
    }
private:
    
    GSLRNGGen(const GSLRNGGen&); // disable ability to copy these objects, as that would alias "gen" more than once.
    GSLRNGGen& operator=(const GSLRNGGen&);

    gsl_rng * gen;
}
;

class GSLSuperCubeGen;
DECLARESP (GSLSuperCubeGen);

class  RNG_DLL GSLSuperCubeGen : public GSLRNGGen,
            public virtual ISuperCubeRNGGen
{
protected:
    GSLSuperCubeGen(gsl_rng * _gen) : GSLRNGGen(_gen)
    {}
    GSLSuperCubeGen(const gsl_rng_type *_type = NULL, unsigned long _seed = 0) : GSLRNGGen(_type, _seed)
    {}
public:
    IRNGGeneratorSP  clone() const
    {
        gsl_rng * cl = gsl_rng_clone(getGen());

        return GSLSuperCubeGenSP(new GSLSuperCubeGen(cl));
    }
    static GSLSuperCubeGenSP create(const gsl_rng_type *_type = NULL, unsigned long _seed = 0)
    {
        return GSLSuperCubeGenSP(new GSLSuperCubeGen(_type, _seed));
    }

    // SuperCubeHooks interface
    virtual void superCubeAdjustment()
    {}

    // This is a naive attempt to duplicate functionality from SC_ran2Gen class
    // where seekToPath(iPath) should restart generator at the new point.
    // of course, the "seed" is meaningless outside of initialization, so we set it to 1 and reinitialize.
    // FIXME: switch to streams in the future;

    virtual void seekToPath(int iPath)
    {
        const long seedStride = 1000;
        const long seedMod = 16061968;
        long seed =1;
        seed = (seed +  seedStride* (long)(1+iPath) )% seedMod + iPath;
        if(seed > 0)
            seed *= -1;
        setSeed( seed);
    }
    virtual long getSeed() const
    {
        throw "can't do that";
    }
};

CORE_END_NAMESPACE

#endif
