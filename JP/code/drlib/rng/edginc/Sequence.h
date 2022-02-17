// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* Sequence.h: interface for the Sequence class.                       *
 * Author: M.Huq                                                       */
// -----------------------------------------------------------------------

#ifndef _SC_SEQUENCE_H
#define _SC_SEQUENCE_H

#include "edginc/coreConfig.hpp"
#include "edginc/DECLARESP.h"

#include <cstdio>
#include <vector>

CORE_BEGIN_NAMESPACE

/** Sequence is the base class of the SuperCube library that produces multi-dimensional (Q)RNG.
 * The class hardcodes storage for the n-dimensional vector and supports re-initialization with a different seed.
 */
class  RNG_DLL Sequence
{
public:

    virtual ~Sequence();
    virtual void populateVector() = 0;
    virtual void populateVector(int iPath);
    virtual double * getVector();
    virtual size_t getDimension() const;
    virtual void advanceSequence(const int numberElements); // FIXME: has unexpected semantics!!!
    virtual void reinitialize(const long _seed); // FIXME : move out of base class

protected:
    int nDimensions;    // Number of dimensions of random numbers to return.
    std::vector<double> vector;     // Array returning the sequence.

    Sequence();         // FIXME: get rid of constructor()/initialize(....) approach
    Sequence(size_t nDim);

private:
    // These two are mis-implemented (copy pointer to FILE), so disable them for now
    Sequence(const Sequence & t);
    Sequence & operator=(const Sequence &t);

};

DECLARESP(Sequence);

// Sequence that is backed by a file
// File format:
// line1: nDimensions nItems #nItems is ignored
// line2: float1 float2 ... float_nDimensions
// line3: float float ...

class RNG_DLL FileSequence : public Sequence
{
    std::FILE *filePointer;
public:
    FileSequence(const char *fileName); // Constructor that takes a file. FIXME: move to a new class
    virtual void populateVector();
    virtual ~FileSequence();
};

DECLARESP(FileSequence);

CORE_END_NAMESPACE

#endif // !defined(_SC_SEQUENCE_H)
