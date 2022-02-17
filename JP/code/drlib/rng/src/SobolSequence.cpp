// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// SobolSequence.cpp: implementation of the SobolSequence class.
// Author: M.Huq
// -----------------------------------------------------------------------

#include "edginc/SobolSequence.h"
#include "edginc/SCException.h"

#include <stdlib.h>
#include <iostream>

#include "sobol.h" // DIMENSION, MAXBIT and MAXORDPOLY defined here


CORE_BEGIN_NAMESPACE


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
SobolSequence::SobolSequence()
{
    maximumSobolDim = 0;
    Sobol_gen = NULL;
    Sobol_xn = NULL;
    Sobol_n = NULL;
    Sobol_nstream = 0;
}

SobolSequence::SobolSequence(const long setSeed, const int maxSobolDimension) : LowDiscrepancySequence(setSeed, maxSobolDimension)
{
    nDimensions = maxSobolDimension;
    // seed = setSeed;

    // parameter checking
    if(nDimensions > DIMENSION) {
        throw SCException( __FILE__, __LINE__,
                           "Too high dimension for initializing Sobol generator");
    }

    // Allocate memory for Sobol sequence generator
    Sobol_gen = (int **)malloc(sizeof(int *) * MAXBIT);
    for(int i = 0; i < MAXBIT; i++) {
        Sobol_gen[i] = (int *)malloc(sizeof(int) * DIMENSION);
        if(!Sobol_gen[i]) {
            throw SCException( __FILE__, __LINE__,
                               "Could not allocate memory for Sobol_gen");
        }
    } //i
    Sobol_xn = new int [DIMENSION];
    Sobol_n  = new int [DIMENSION];

    // Initialize the Sobol Sequence code.
    if(sobinit(nDimensions,
               Sobol_gen,
               Sobol_xn,
               Sobol_n,
               &Sobol_nstream)!=0) {
        throw SCException( __FILE__, __LINE__,
                           "Could not initialize Sobol sequence generator.\n Call to sobinit returned failure.");
    } // sobinit

    // Generate first 50 sobol as a rule as initialization
    //  for(int it = 0; it < 50; it++){
    //    if(sobvect(nDimensions,
    //        vector,
    //        Sobol_gen,
    //        Sobol_xn,
    //        Sobol_n,
    //        &Sobol_nstream)!= 0){
    //    throw SCException( __FILE__, __LINE__,
    //         "Could not generate Sobol sequence: see stderr");
    //    }
    //  }
}

SobolSequence::~SobolSequence()
{
    if(Sobol_gen) {
        for(int i = 0; i < MAXBIT; i++) {
            free(Sobol_gen[i]);
            Sobol_gen[i] = 0;
        }
        free(Sobol_gen);
    }
    if(Sobol_xn) {
        delete [] Sobol_xn;
        Sobol_xn = 0;
    }
    if(Sobol_n) {
        delete [] Sobol_n;
        Sobol_n = 0;
    }
}


void SobolSequence::populateVector()
{
    // Check basic parameters
    if(nDimensions > DIMENSION || nDimensions > Sobol_nstream) {
        throw SCException( __FILE__, __LINE__,
                           "Dimension too large or Sobol generator not initialized with sufficiently high dimension.");

    }
    if(sobvect(nDimensions,
               getVector(),
               Sobol_gen,
               Sobol_xn,
               Sobol_n,
               &Sobol_nstream)!= 0) {
        throw SCException( __FILE__, __LINE__,
                           "Call to sobvect returned failure. Could not populate vector.");
    }
}
void SobolSequence::populateVector(int iPath)
{
    populateVector();
}



void SobolSequence::advanceSequence(const int numberElements)
{
    //NOTE: each time we increment seed, we are jumping over a
    //vector of size "nDimensions" worth of Sobol numbers.
    //This is what happens in the original advanceSequence()
    //of the base class (check Sequence.cpp)

    int b, j, i;

    // basic idea:
    // update Sobol_xn[] so that it reflects the prior element and
    // then invoke populateVector to get the desired elements.
    // The trick is to NOT use the Antonov-Saleev to skip ahead since
    // it is sequential. Instead use the original formulation by Sobol:
    // S_n = b1.VC1 XOR b2.VC2 XOR ... bk.VCk  where
    // [bk...b2b1] = Gray code for n.
    // in our case, n = Sobol_n + numberElements

    for (i=0;i<nDimensions; i++) {
        Sobol_n[i] += numberElements;
        // we need to pick the direction numbers corresponding to the
        // Gray code of Sobol_n[i]; Easy trick to get Gray code of Sobol_n[i]:
        j = Sobol_n[i] ^ (Sobol_n[i]/2);

        // forget the past; we are jumping ahead from the beginning
        Sobol_xn[i] = 0;

        for (b=0; j; j>>=1, b++) {
            if (j & 0x1) {
                Sobol_xn[i] ^= Sobol_gen[b][i];
            }
        }
        vector[i] = Sobol_xn[i]/((double)(1<<MAXBIT));
    }

}

CORE_END_NAMESPACE
