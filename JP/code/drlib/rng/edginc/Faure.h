// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* Faure.h                                                                 *
 * Author: M.Huq                                                           *
 * Class definition for basic Faure sequence.                              *
 *                                                                         *
 * Reference:                                                              *
 * L. Finschi, "Quasi-Monte-Carlo: An empirical study on Low-Discrepancy   *
 *              Sequences"                                                 *
   ----------------------------------------------------------------------- */
#ifndef _SC_FAURE_H
#define _SC_FAURE_H

#include "edginc/LowDiscrepancySequence.h"

CORE_BEGIN_NAMESPACE

class  RNG_DLL Faure : public LowDiscrepancySequence
{
public:
    Faure();
    Faure(const long setSeed, const int setDimensions);
    virtual ~Faure();
    virtual void populateVector();
    virtual void advanceSequence(const int numberElements);

protected:
    int closestPrime;
    long seed;
};

CORE_END_NAMESPACE

#endif // !defined(_SC_FAURE_H)
