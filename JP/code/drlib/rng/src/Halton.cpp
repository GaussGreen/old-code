// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// Halton.cpp
// Author: M. Huq
// Implementation of the basic halton sequence class.
// -----------------------------------------------------------------------
#include "edginc/Halton.h"
#include "edginc/SCException.h"

#include <stdlib.h>
#include <iostream>
#include <math.h>

CORE_BEGIN_NAMESPACE

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Halton::Halton(const long setSeed, const int setDimensions) : LowDiscrepancySequence(setSeed, setDimensions)
{
    // Build prime number list...
    setInitialPrime(2);
    allocatePrimesList(nDimensions+1);
    populatePrimesList(nDimensions+1);
#ifdef SC_DEBUG

    dumpPrimesList("Halton");
#endif
}

Halton::~Halton()
{
}

void Halton::populateVector()
{
    // Loop over each of the dimensions and compute the sequence
    for(int j = 0; j < nDimensions; j++) {
        double p_inv = 1./PrimeNumberList[j];
        double pow_p_inv = p_inv;
        /* compute digit expansion of n in the current base */
        getDigitExpansion(seed, PrimeNumberList[j]);

        /* Construct each radical inverse function */
        vector[j] = 0.0;
        for(int k = 0 ; k < digitCount; k++) {
            vector[j] += (double)digitExpansion[k] * pow_p_inv;
            pow_p_inv *= p_inv;
        }
    }
    ++seed;
}


void Halton::advanceSequence(const int numberElements)
{
    //NOTE: each time we increment seed, we are jumping over a
    //vector of size "nDimensions" worth of Halton numbers.
    //This is what happens in the original advanceSequence()
    //of the base class (check Sequence.cpp)
    seed+=numberElements;
}

CORE_END_NAMESPACE

