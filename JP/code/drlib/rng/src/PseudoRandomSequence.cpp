// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
// PseudoRandomSequence.cpp: implementation of the PseudoRandomSequence class.
// Author: M.Huq, Credit DR
// -----------------------------------------------------------------------
#include "edginc/PseudoRandomSequence.h"
#include "edginc/SCException.h"

CORE_BEGIN_NAMESPACE

#if 0
// Numerical recipies includes
// extern "C" {
// double SC_gasdev2(long *idum);
// }

//////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////
// PseudoRandomSequence::PseudoRandomSequence() : seed(-1)
// {
// }

PseudoRandomSequence::PseudoRandomSequence(IRNGSP _rng, const int setDimensions) :
        Sequence(setDimensions),
        rng(_rng)
{
    if(rng.get() == 0) {
        throw SCException(__FILE__, __LINE__,
                          "Bad value for seed (=0)");
    }
    //setSeed( _setSeed);
    /*  nDimensions = setDimensions;
      if(!vector)
     vector = new double [nDimensions];
      if(!vector){
        throw SCException(__FILE__, __LINE__,
            "Failed to allocate memory for vector");
      }*/
}


void PseudoRandomSequence::populateVector()
{
    if(vector.empty()) {
        throw SCException(__FILE__, __LINE__,
                          "vector not allocated for population");
    }
    for(int i=0; i<nDimensions; i++) {
        vector[i] = rng->fetch(); // SC_gasdev2(&seed);
    }
}

void PseudoRandomSequence::populateVector(int iPath)
{
    populateVector();
}



#endif

CORE_END_NAMESPACE
