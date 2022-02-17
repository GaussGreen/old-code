// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2002 J.P. Morgan-Chase & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*
  gasdev2.h
  Author: M.Huq, Derivatives Research
 
  Header for gasdev2 related functions. Include as extern "C" from C++.
 
  NOTE: There are two new classes, UniformRandomSequence and
  GaussianRandomSequence. These are intended to replace ran2 and gasdev2.
  It is recommended that the two new classes are used in lieu of ran2 and
  gasdev2.
 
 
 */

#ifndef _GASDEV2__H_
#define _GASDEV2__H_

#include "edginc/DECLARESP.h"

//#include "edginc/ran2.h"
#include "edginc/IRNGGenerator.h"
#include "edginc/IRNG.h"
#include "edginc/ISuperCubeHooks.h"
#include "edginc/SuperCubeHelper.h"

CORE_BEGIN_NAMESPACE

// Num recepies GasDev2 algorithm
// To create an instance, call create() method
class GasDev2;
DECLARESP(GasDev2);

/** GasDev2 is a gaussian random number generator based on the gasdev2 algorithm.
 * It uses acceptance/rejection and consumes pairs of uniform numbers. Each time a pair of normal numbers is generated and excess is buffered.
*/

class RNG_DLL GasDev2 : public virtual INormalRNG
{
    IUniformRNGGenSP uniform;
protected:
    bool iset;
    double gset;
    GasDev2(IUniformRNGGenSP rng, bool _iset, double _gset);

public:
    virtual IRNGGeneratorSP getGenerator() const
    {
        return uniform;
    }
    virtual IRNGSP  clone() const
    {
        return GasDev2SP(new GasDev2(DYNAMIC_POINTER_CAST<IUniformRNGGen>(uniform->clone()), iset, gset));
    }
    static GasDev2SP create (IUniformRNGGenSP rng, bool _iset = false, double _gset = 0.);
    virtual double fetch(); // get next tranformed RNG
};

class SC_gasdev2;
DECLARESP(SC_gasdev2);

/** GasDev2 class with support for SuperCube needes.
 */
class RNG_DLL SC_gasdev2  : public GasDev2,
            public SuperCubeHelper,
            public virtual INormalSuperCubeRNG
{
    ISuperCubeHooks *hooks;
    void init();
protected:
    SC_gasdev2(ISuperCubeRNGGenSP rng, bool _iset, double _gset);
public:
    static SC_gasdev2SP     create(ISuperCubeRNGGenSP rng, bool _iset = false, double _gset = 0.);
    virtual IRNGSP  clone() const;

};

CORE_END_NAMESPACE

#endif // _GASDEV2__H_
