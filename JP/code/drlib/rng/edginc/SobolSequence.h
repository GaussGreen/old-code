// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// SobolSequence.h: interface for the SobolSequence class.
// Author: M.Huq
//
// Wed Jun 20 15:47:44 EDT 2001 - M.Huq
//    Created new version of SobolSequence class that utilizes SPRNG version
//    of Sobol sequence generator.
// -----------------------------------------------------------------------

#ifndef _SC_SOBOLSEQUENCE_H
#define _SC_SOBOLSEQUENCE_H

#include "edginc/LowDiscrepancySequence.h"

CORE_BEGIN_NAMESPACE

class RNG_DLL SobolSequence : public LowDiscrepancySequence
{
public:
    SobolSequence();  // Default constructor. Does nothing
    SobolSequence(const long setSeed, const int maxSobolDimension);
    virtual ~SobolSequence();
    virtual void populateVector();
    virtual void populateVector(int iPath);
    void advanceSequence(const int numberElements);


protected:
    int maximumSobolDim;
    // Arrays required for modified SPRNG routines
    int **Sobol_gen; // pre-calculated value of vc in Fox' paper. sobinit inits
    int *Sobol_xn;   // stores the previous xn values in each dimension
    int *Sobol_n;    // stores the current n-value in each dimension
    int Sobol_nstream; // stores the number of sobol streams being used.
};

DECLARESP(SobolSequence);

CORE_END_NAMESPACE
#endif // !defined(_SC_SOBOLSEQUENCE_H)
