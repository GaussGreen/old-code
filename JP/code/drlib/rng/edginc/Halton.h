// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* Halton.h                                                                *
 * Author: M.Huq                                                           *
 * Class definition for basic Halton sequence.                             *
 *                                                                         *
 * Reference:                                                              *
 * "Computational Investigations of Low-Discrepancy Sequences",            *
 * L.Kocis and W.J.Whiten, ACM Transactions on Mathematical Software,      *
 * Vol. 23, No.2 June 1997, Pages 266-294.                                 *
 ----------------------------------------------------------------------- */
#ifndef _SC_HALTON_H
#define _SC_HALTON_H

#include "edginc/LowDiscrepancySequence.h"
#include "edginc/DECLARESP.h"

CORE_BEGIN_NAMESPACE

class  RNG_DLL Halton : public LowDiscrepancySequence
{
public:
    Halton(const long setSeed, const int setDimensions);
    virtual ~Halton();
    virtual void populateVector();
    virtual void advanceSequence(const int numberElements);
private:
    long seed;

};

DECLARESP(Halton);

CORE_END_NAMESPACE

#endif // !defined(_SC_HALTON_H)
