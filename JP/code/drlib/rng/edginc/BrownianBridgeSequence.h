// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* BrownianBridgeSequence.h                                                 *
 * Author: M. Huq                                                           *
 *                                                                          *
 * Class for basic brownian bridge based on Arnon Levys code.               *
 * Essentially this is a front-end to that code. The class constructor      *
 * initializes the memory and data for the brownian bridge. The             *
 * populateVector method samples the brownian bridge on a per path basis.   *
 * populateVector() returns a vector that contains the entire list of       *
 * deviates. By "space" dimension and time. Brownian bridge applied in the  *
 * "time" direction.                                                        *
--------------------------------------------------------------------------  */
#ifndef _SC_BROWNIANBRIDGESEQUENCE_H
#define _SC_BROWNIANBRIDGESEQUENCE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "edginc/Sequence.h"
#include "edginc/IRNG.h"
#include "VectorMatrix.h"

CORE_BEGIN_NAMESPACE

/// Specifiers for dimension type per date.
const long BB_PSEUDO_RANDOM = -1L;
const long BB_GRID = 0L ;
/** A one-dimensional grid. Can be used on exactly
    one dimension among all
    dimensions (across all factors). Strongly
    recommended
    that this be used
    on the most important dimension
    - usually last time point in single factor
    implementation.
    - usually last time point in most important
    dimension in a multi-factor
    implementation.
*/
const long BB_SOBOL = 1L;   ///< Constructor will count up sobol sequences for Arnon's code.


class RNG_DLL BrownianBridgeSequence : public Sequence
{
public:
    BrownianBridgeSequence();
    BrownianBridgeSequence( int dimension,   // dimensionality of deviate that we are after.
                            int numDates,    // number of dates or stratifications.
                            long **dimTypeList, // type specifier for dimension by date.
                            int numberPaths,   // number of paths that will be sampled.
                            INormalSuperCubeRNGSP rng );   // Seed for pseudo-random part.

    // Following constructor takes in as a last argument a pointer to an array of doubles.
    // The array size should be the same as the number of paths to be used. The grid in
    // brownian bridge can be specified through this interface. If this pointer is NULL
    // an internal grid is set up uniform in probability.
    // If the array is specified the values specified should be gaussian deviates. For example
    // they can have the values from -5 to +5. Use this feature at your own risk.
    BrownianBridgeSequence( int dimension,   // dimensionality of deviate that we are after.
                            int numDates,    // number of dates or stratifications.
                            long **dimTypeList, // type specifier for dimension by date.
                            int numberPaths,   // number of paths that will be sampled.
                            INormalSuperCubeRNGSP rng ,     // Seed for pseudo-random part.
                            double *gaussianGridPtr ); // Pointer to gaussian grid


    // Method only to be used when using the default constructor.
    void initializeBrownianBridge( int dimensions,   // dimensionality of deviate that we are after.
                                   int numDates,    // number of dates or stratifications.
                                   long **dimTypeList, // type specifier for dimension by date.
                                   int numberPaths,   // number of paths that will be sampled.
                                   INormalSuperCubeRNGSP rngIn ,     // Seed for pseudo-random part.
                                   double *gaussianGridPtr ); // pointer to grid, gaussianized
    //The following constructor takes as an argument a pointer to a Sequence
    // to be used for the supports of the bridge.
    BrownianBridgeSequence( int dimensions, // dimensionality of deviate that we are after.
                            int numDates,   // number of dates or stratifications.
                            long **dimTypeList, // type specifier for dimension by date.
                            int numberPaths,   // number of paths that will be sampled.
                            INormalSuperCubeRNGSP rngIn ,     // Seed for pseudo-random part.
                            double *gaussianGridPtr,  // pointer to grid, gaussianized
                            SequenceSP ldsPtr ); // Pointer to low-discrepancy sequence to be used
    void initializeBrownianBridge( int dimensions,   // dimensionality of deviate that we are after.
                                   int numDates,    // number of dates or stratifications.
                                   long **dimTypeList, // type specifier for dimension by date.
                                   int numberPaths,   // number of paths that will be sampled.
                                   INormalSuperCubeRNGSP rngIn ,     // Seed for pseudo-random part.
                                   double *gaussianGridPtr,  // pointer to grid, gaussianized
                                   SequenceSP ldsPtr ); // Pointer to low-discrepancy sequence to be used


    virtual ~BrownianBridgeSequence();
    virtual void populateVector();         // Default populateVector. Not to be used.
    virtual void populateVector( int thisPath );  // populateVector per path. Path parameter determines grid.
    virtual double *getDeviatesbyPathDate( int iPath, int iDate ); // Use this method to get deviates.
    virtual int getNumberBBFactors() const
    {
        return factorDimension;
    }
    virtual int getNumberBBDates() const
    {
        return numberDates;
    }
    virtual SequenceSP getLowDiscrepancySequence() const
    {
        return lowDiscSeqPtr;
    }
    virtual void setLowDiscrepancySequence( SequenceSP theLDSSequence )
    {
        lowDiscSeqPtr = theLDSSequence;
    }
    virtual void getDeviatesbyFactorDate( int iFactor,
                                          int iDate,
                                          double *&returnVector );

    virtual void getBrownianBridgeStructure( tnti_mat &theMat );


private:
    INormalSuperCubeRNGSP rng;

    int factorDimension;
    int numberDates;
    long **factorTypeList;
    int nPaths;
    double *gaussianVect;
    double *probabilities;
    long *listDates;
    int gridLocation;
    int gridFactor;
    double deltaGrid;
    double *ptr2GaussianGrid;
    SequenceSP lowDiscSeqPtr;   ///< Pointer for low-discprepancy sequence to be used within BB.
    SequenceSP internalLowDiscSeqPtr; ///< internal pointer for memory allocation
};

DECLARESP( BrownianBridgeSequence );

CORE_END_NAMESPACE

#endif // !defined(_SC_BROWNIANBRIDGESEQUENCE_H)
