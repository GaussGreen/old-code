// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* PseudoRandomSequence.h                                                 *
 * Author: M.Huq                                                          *
 *                                                                        *
 * PseudoRandomSequence.h: interface for the PseudoRandomSequence class.  */
// -----------------------------------------------------------------------

#ifndef _SC_PSEUDORANDOMSEQUENCE__H
#define _SC_PSEUDORANDOMSEQUENCE__H

#include "edginc/Sequence.h"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

// Transformer of RNG into vector based ISequence objects
// Transformation is essentially leap-frogging with step = dimensions
class  RNG_DLL PseudoRandomSequence : public Sequence // FIXME deprecate later
{
    IRNGSP      rng;
public:

    PseudoRandomSequence(IRNGSP _rng, size_t _dimension) :
            Sequence(_dimension),
            rng(_rng)
    {}
    virtual void populateVector()
    {
        size_t n = getDimension();
        double * v = getVector();
        for(size_t i= 0; i!= n; ++i)
            v[i]=rng->fetch();
    }
    virtual void populateVector(int iPath)
    {
        populateVector();
    } // FIXME: code smell
    IRNGSP    getRNG() const
    {
        return rng;
    }
};

DECLARESP(PseudoRandomSequence);

/**
    PseudoRandom template allows one to create PseudoRandomSequence with a specific transformatiuon (RNG) and a minimum level of Generator (RNGGen).
    The resulting object should be initialized with dimension and a  generator of type atleast specified in the template.
    If generator has default create() method, one can pass only number of dimensions.
*/

template <class RNG , class RNGGen >
class  RNG_DLL PseudoRandom : public PseudoRandomSequence
{
    DECLARESP(RNG);
    DECLARESP(RNGGen);
    DECLARESP(PseudoRandom);

    RNGSP   rng;

public:

    PseudoRandom(RNGGenSP gen, size_t _dimension)
            : PseudoRandomSequence(_dimension, RNG::create(gen))
    {}
    // If Generator has crate() method that accepts no params, one is free to use it.
    // Note that we cannot join the two constructors, as some types, like IUniformRNGGen have no create() meth.
    PseudoRandom(size_t _dimension)
            : PseudoRandomSequence(_dimension, RNG::create(RNGGen::create()))
    {}
}
;

#if 0
// Example how we could create new Sequences by matching generators
#include "edginc/ran2.h"
typedef PseudoRandom<SC_ran2, SC_ran2Gen> UniformRandomSequenceNew;
typedef PseudoRandom<SC_gasdev2, SC_ran2Gen> GaussianRandomSequenceNew;
#endif

#if 0
class  RNG_DLL PseudoRandomSequence : public Sequence
{
public:
    //  PseudoRandomSequence();    // Default constructor
    PseudoRandomSequence(IRNGSP rng, const int dimensions); // Constructor that should be used.
    virtual ~PseudoRandomSequence();
    virtual void populateVector();
    virtual void populateVector(int iPath);

private:
    IRNGSP rng;
};
#endif

CORE_END_NAMESPACE

#endif // !defined(_SC_PSEUDORANDOMSEQUENCE_H)
