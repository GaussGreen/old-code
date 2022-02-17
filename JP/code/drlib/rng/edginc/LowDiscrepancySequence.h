// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* LowDiscrepancySequence.h                                                 *
 * Author: M.Huq                                                            *
 * Class definition for LowDiscrepancySequence.                             *
 * This class contains methods that are common to Halton and Faure sequences*
 * and hence is an intermediate class from which Halton and Faure classes   *
 * will be derived.                                                         *
 * This class itself will not return any valid sets of numbers.             */
// -----------------------------------------------------------------------
#ifndef _SC_LOWDISCREPANCYSEQUENCE_H
#define _SC_LOWDISCREPANCYSEQUENCE_H

#include "edginc/Sequence.h"

CORE_BEGIN_NAMESPACE

class  RNG_DLL LowDiscrepancySequence : public Sequence
{
protected:
    double discrepancy;  // Some measure of discrepancy
    int digitCount;
    int digitBufferSize; // Size in constructor if not passed
    int *digitExpansion;  //

    // The following will be useful for Halton and Faure-like sequences
    // where a list of prime numbers may be necessary.
    // We allocate memory for these primes from within the derived classes
    // and provide methods here for that allocation.
    int nPrimes;         // Number of primes
    int initialPrime;    // Initial prime set in constructor and changeable.
    int *PrimeNumberList; // List of prime numbers used by some sequences.
    int primeArraySize;


    // Methods internal to this class and derived class
    virtual void getDigitExpansion(const int number, const int base);
    virtual void allocatePrimesList(const int primeSize);
    virtual void populatePrimesList(const int thisMany);
    virtual void setInitialPrime(const int thisInitialPrime);
    virtual void dumpPrimesList(const char *tag);

public:
    LowDiscrepancySequence();
    LowDiscrepancySequence(const long seed, const int nDimensions);
    LowDiscrepancySequence(const long setSeed, const int setDimensions, const int setDigitBufferSize);
    virtual ~LowDiscrepancySequence();
private:

};

DECLARESP(LowDiscrepancySequence);

CORE_END_NAMESPACE

#endif // !defined(_SC_LOWDISCREPANCYSEQUENCE_H)
