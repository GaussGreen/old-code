// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* LowDiscrepancySequence.cpp                                                   *
 * Author: M. Huq                                                               *
 * Methods for the LowDiscrepancySequence class.                                *
----------------------------------------------------------------------- */
#include "edginc/LowDiscrepancySequence.h"
#include "edginc/SCException.h"

#include <iostream>
#include <stdlib.h>

CORE_BEGIN_NAMESPACE


using namespace std;

/////////////////////////////////////
// Constructors
/////////////////////////////////////
LowDiscrepancySequence::LowDiscrepancySequence()
{
    //  vector = 0;
    discrepancy = 0.0;
    digitCount = 0;
    digitBufferSize = 0;
    digitExpansion = NULL;
    nPrimes = 0;
    initialPrime = 0;
    PrimeNumberList = NULL;
    primeArraySize = 0;
}

/* Constructor that internally sets the size of the digit expansion array */
LowDiscrepancySequence::LowDiscrepancySequence( const long _setSeed, const int setDimensions )
{
    //    setSeed(_setSeed);
    nDimensions = setDimensions;

    vector.assign( nDimensions, 0.0 );

    // Handle digit count arrays here
    digitCount = 0;
    digitBufferSize = 5000; // Hardcode to an internal maximum.
    digitExpansion = new int [ digitBufferSize ];
    if ( !digitExpansion ) {
        throw SCException( __FILE__, __LINE__,
                           "Could not allocate memory for digitExpansion" );
    }

    // Initialize other variables
    discrepancy = 99999999.99;
    nPrimes = 0;
    PrimeNumberList = NULL;
    primeArraySize = 0;
    initialPrime = 2;
}

/* Constructor that takes in the size of the digit expansion array */
LowDiscrepancySequence::LowDiscrepancySequence( const long _setSeed,
        const int setDimensions,
        const int setDigitBufferSize ) :
        Sequence( setDimensions )
{
    //    setSeed(_setSeed);

    // Handle digit count arrays here
    digitCount = 0;
    digitBufferSize = setDigitBufferSize;
    digitExpansion = new int [ digitBufferSize ];
    if ( !digitExpansion ) {
        throw SCException( __FILE__, __LINE__,
                           "Could not allocate memory for digitExpansion" );
    }

    // Initialize other variables
    discrepancy = 99999999.99;
    nPrimes = 0;
    initialPrime = 2;
    PrimeNumberList = NULL;
    primeArraySize = 0;
    initialPrime = 2;
}

/////////////////////////////////////
// Destructor
/////////////////////////////////////
LowDiscrepancySequence::~LowDiscrepancySequence()
{
    // Free up all the arrays that are allocated
    if ( digitExpansion ) {
        delete [] digitExpansion;
        digitExpansion = 0;
    }
    if ( PrimeNumberList ) {
        delete [] PrimeNumberList;
        PrimeNumberList = 0;
    }
}

/////////////////////////////////////
// Method to get digit expansion.
/////////////////////////////////////
void LowDiscrepancySequence::getDigitExpansion( const int number, const int base )
{
    static char cRoutine[] = "LowDiscrepancySequence::getDigitExpansion(int,int)";
    int modVal;
    modVal = number;
    digitCount = 0;
#ifdef MYDEBUG

    cout << cRoutine << ": Digit expansion for %d in base %d: " << number << " " << base << endl;
#endif

    while ( modVal != 0 ) {
        digitExpansion[ digitCount ] = modVal % base;
#ifdef MYDEBUG

        cout << digitExp[ ( *digitCount ) ] << " ";
#endif

        digitCount++;
        modVal = ( int ) ( modVal / base );
    }
#ifdef MYDEBUG
    cout << "\n digit count = " << digitCount << endl;
#endif

} /* End getDigitExpansion */


/////////////////////////////////////
// Method to allocate primes list
/////////////////////////////////////
void LowDiscrepancySequence::allocatePrimesList( const int primeSize )
{
    if ( !PrimeNumberList ) {
        PrimeNumberList = new int [ primeSize ];
        if ( !PrimeNumberList ) {
            throw SCException( __FILE__, __LINE__,
                               "Could not allocate memory for PrimeNumberList" );
        }
    }
    primeArraySize = primeSize;
}

/////////////////////////////////////
// Method to populate primes list
/////////////////////////////////////
void LowDiscrepancySequence::populatePrimesList( const int thisMany )
{
    if ( !PrimeNumberList ) {
        PrimeNumberList = new int [ thisMany ];
        if ( !PrimeNumberList ) {
            throw SCException( __FILE__, __LINE__,
                               "Could not allocate memory for PrimeNumberList" );
        }
        primeArraySize = thisMany;
    }
    if ( primeArraySize < thisMany ) {
        delete [] PrimeNumberList;
        PrimeNumberList = new int [ thisMany ];
        if ( !PrimeNumberList ) {
            throw SCException( __FILE__, __LINE__,
                               "Could not allocate memory for PrimeNumberList" );
        }
        primeArraySize = thisMany;
    }

    // Now populate list
    int number = initialPrime;
    int i, j;
    int remainder;
    PrimeNumberList[ 0 ] = initialPrime;
    nPrimes = 1;

    number++;
    /* loop over number of primes desired */
    while ( nPrimes < thisMany ) {
        i = 0;
        while ( number < 99999999 ) {
            remainder = 1;
            j = 0;
            while ( remainder != 0 && j < nPrimes ) {
                remainder = number % PrimeNumberList[ j ];
                j++;
            }
            if ( remainder != 0 ) {
                /* found a prime */
                break;
            } else {
                number++;
            }
        }
        PrimeNumberList[ nPrimes ] = number;
        number++;
        nPrimes++;
    }
}

/////////////////////////////////////
// Method to set initial prime
/////////////////////////////////////
void LowDiscrepancySequence::setInitialPrime( const int thisInitialPrime )
{
    initialPrime = thisInitialPrime;
    if ( PrimeNumberList ) {
        PrimeNumberList[ 0 ] = initialPrime;
    }
}

/////////////////////////////////////
// Method to dump prime number list
/////////////////////////////////////
void LowDiscrepancySequence::dumpPrimesList( const char *tag )
{
    if ( !PrimeNumberList ) {
        cout << "dumpPrimesList : Zero prime numbers - populate first!" << endl;
    } else {
        cout << "=======Prime Number List : " << tag << " =======" << endl;
        for ( int i = 0 ; i < nPrimes; i++ ) {
            cout << " prime[" << i << "] = " << PrimeNumberList[ i ] << endl;
        }
        cout << "================================================" << endl;
    }
}
CORE_END_NAMESPACE
