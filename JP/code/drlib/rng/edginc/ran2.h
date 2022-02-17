// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2002 J.P. Morgan-Chase & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*
  ran2.h
  Author: M.Huq, Derivatives Research

  Header for ran2 related functions. Include as extern "C" from C++.

  NOTE: There are two new classes, UniformRandomSequence and
  GaussianRandomSequence. These are intended to replace ran2 and ran2.
  It is recommended that the two new classes are used in lieu of ran2 and
  ran2.


 */

#ifndef _RAN2__H_
#define _RAN2__H_

#include "edginc/coreConfig.hpp"
#include "edginc/DECLARESP.h"
#include "edginc/IRNGGenerator.h"
#include "edginc/IRNG.h"
#include "edginc/TrivialTransform.h"

CORE_BEGIN_NAMESPACE
//////////////////// Generators //////////////////////////
class Ran2Gen;
DECLARESP(Ran2Gen);

/** Ran2Gen class is a Uniform RNG Generator that uses ran2 Num.Receip. implementation.
    Instances of this class do not share state with other Ran2Gen instances. If you need this functionality, for example to match numbers from the original SuperCube library, please use \class Ran2StaticGen.
*/
class RNG_DLL Ran2Gen : public virtual IUniformRNGGen
{
protected:
    static const int NTAB = 32;
    // FIXME: old implementation assumed theres only 1 generator, so all fields were static
    // double check if anyone creates more that one generator
    long seed;
    long idum2;
    long iy;
    long iv[NTAB];

    Ran2Gen(long _seed) :  seed(_seed), idum2(123456789L), iy(0L)
    {}
public:
    static Ran2GenSP create(long _seed = -1L)
    {
        return Ran2GenSP(new Ran2Gen(_seed));
    }
    virtual IRNGGeneratorSP  clone() const
    {
        return Ran2GenSP(new Ran2Gen(*this));
    }
    virtual double fetch();
    virtual void debug()
    {
        std::fprintf(stdout, "debug: seed= %ld\n", seed);
    }
    long    getSeed() const
    {
        return seed;
    }

};

class SC_ran2Gen;
DECLARESP(SC_ran2Gen);

class RNG_DLL SC_ran2Gen : public Ran2Gen,
            public virtual ISuperCubeRNGGen
{
    protected:
        SC_ran2Gen(long _seed) : Ran2Gen(_seed)
        {}

    public:
        static SC_ran2GenSP    create(long _seed = -1L)
        {
            return SC_ran2GenSP(new SC_ran2Gen(_seed));
        }
        virtual IRNGGeneratorSP    clone() const
        {
            return SC_ran2GenSP(new SC_ran2Gen(*this));
        }

        virtual void superCubeAdjustment(); // FIXME: delete after SupCub is finished
        virtual void seekToPath(int iPath); // using current state, seek to path -- SuperCube
        virtual long getSeed() const
        {
            return Ran2Gen::getSeed();
        }
};



/*******************
* Generator to replicate the original SuperCube numbers. Note that
*/
class Ran2StaticGen;
DECLARESP(Ran2StaticGen);

/** Ran2StaticGen class shares its state between instances. It has the same functionality as the old C-base code.
 * Fetching numbers with a negative seed will reset RNG. Different clients of this class should reinitialize RNG before using it. Otherwise, inner state is changed in some random way.
 */
class RNG_DLL Ran2StaticGen : public virtual IUniformRNGGen
{
protected:
    static const int NTAB = 32;
    // FIXME: old implementation assumed theres only 1 generator, so all fields were static
    // double check if anyone creates more that one generator
    long seed;
    static long idum2;
    static long iy;
    static long iv[NTAB];

    // NB: a negative seed will reinitialize static variables above
    Ran2StaticGen(long _seed) :  seed(_seed) /*, idum2(123456789L), iy(0L)*/ {}
public:
    static  Ran2StaticGenSP create(long _seed)
    {
        return Ran2StaticGenSP(new Ran2StaticGen(_seed));
    }
    virtual IRNGGeneratorSP  clone() const
    {
        return Ran2StaticGenSP(new Ran2StaticGen(*this));
    }
    virtual double fetch();
    virtual void debug()
    {
        std::fprintf(stdout, "debug: seed= %ld\n", seed);
    }
    long    getSeed() const
    {
        return seed;
    } ///< get the current state of "seed"

};



class SC_ran2StaticGen;
DECLARESP(SC_ran2StaticGen);

class RNG_DLL SC_ran2StaticGen : public /*Ran2Gen*/ Ran2StaticGen,
            public virtual ISuperCubeRNGGen
{
protected:
    SC_ran2StaticGen(long _seed) : /*Ran2Gen*/ Ran2StaticGen(_seed)
    {}

public:
    static SC_ran2StaticGenSP    create(long _seed = -1L)
    {
        return SC_ran2StaticGenSP(new SC_ran2StaticGen(_seed));
    }
    virtual IRNGGeneratorSP    clone() const
    {
        return SC_ran2StaticGenSP(new SC_ran2StaticGen(*this));
    }

    virtual void superCubeAdjustment(); // FIXME: delete after SupCub is finished
    virtual void seekToPath(int iPath); // using current state, seek to path -- SuperCube
    virtual long getSeed() const
    {
        return Ran2StaticGen::getSeed();
    }
};


////////////////////// Usable RNG //////////////////////////////

// Legacy class:
class SC_ran2;
DECLARESP(SC_ran2);

class RNG_DLL SC_ran2 : public TrivialUniformRNG,
            public virtual IUniformSuperCubeRNG
{
    ISuperCubeHooks * hooks;
protected:
    SC_ran2(ISuperCubeRNGGenSP _rng) : TrivialUniformRNG(_rng), hooks(_rng.get())
    {}
public:
    static  SC_ran2SP create(ISuperCubeRNGGenSP _rng)
    {
        return SC_ran2SP(new SC_ran2(_rng));
    }
    virtual void superCubeAdjustment()
    {
        hooks->superCubeAdjustment();
    } // FIXME: delete after SupCub is finished
    virtual void seekToPath(int iPath)
    {
        hooks->seekToPath(iPath);
    } // using current state, seek to path -- SuperCube
    virtual long getSeed() const
    {
        return hooks->getSeed();
    }
};

CORE_END_NAMESPACE

#endif // _RAN2__H_
