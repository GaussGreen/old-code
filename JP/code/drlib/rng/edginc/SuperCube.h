// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* SuperCube.h                                                              *
 * Author: M.Huq, Credit DR                                                 *
 * Class definition for SuperCube.                                          *
 -----------------------------------------------------------------------*/
#ifndef _SC_SUPERCUBE_H
#define _SC_SUPERCUBE_H

#include "edginc/Sequence.h"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

const size_t MAX_SOBOL = 395;      // Maximum allowable sobol for brownian bridge.

enum RngTypes { GET_UNIFORM=0, GET_GAUSSIAN= 1, GET_GAUSSIAN_PADDING= 2, GET_BB_STYLE= 3}; // types for the UniformityFlag -- should not be needed eventually
//==============================================
//SuperCube definition begins here.
//==============================================
/**
 * SuperCube class allows one to put ideas behind LatinHyperCybe into multi dimensional case.
 * I.e. it allows one to generate more random samples by permuting the existing ones.
 * \author Mijan Huq (see also http://stowe/supercube/ )
*/
class RNG_DLL SuperCube : public Sequence
{
public:
    // Parameters to specify whether uniform random deviates or
    // gaussian deviates are to be returned.


    SuperCube(); // Default constructor.
    //========================================================
    // The following are constructor - initializeObject pairs
    //========================================================

#if 0
    // Must specify primary partition sequence to use through
    // the Sequence pointer given by useSequence

    SuperCube(
        INormalSuperCubeRNGSP rng,              // seed for random shuffle
        int numberPaths,    // Number of paths
        int maxDimension,        // max dimensions to set up for
        int sizePartition,       // number of partitions to break up in.
        SequenceSP useSequence); // Sequence to be used. Base class pointer.

    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              // seed for random shuffle
        int numberPaths,    // Number of paths
        int maxDimension,        // max dimensions to set up for
        int sizePartition,       // size of partitions to break up in.
        SequenceSP useSequence); // Sequence to be used. Base class pointer.

#endif
    //========================================================
    /// Must specify primary partition sequence to use through
    /// the Sequence pointer given by useSequence
    /// Method internally figures out the number of partitions
    /// to use. Adds additional padding if maxDimension not large
    /// enough.

    SuperCube(
        INormalSuperCubeRNGSP rng,              ///< Gaussian RNG capable to report its underlying uniform Generator
        int numberPaths,    ///< Number of paths
        int maxDimension,        ///< max dimensions to set up for
        int sizePartition,       ///< number of partitions to break up in.
        SequenceSP useSequence,  ///< Sequence to be used (to fill primary partition?)
        int UniformityFlag = GET_UNIFORM);    ///< Flag to set uniform sequence(GET_UNIFORM = 0) or gaussian (GET_GAUSSIAN = 1). Used for padding (?)
    /** Object initialization is done not via constructor, but through a family of initializeObject() functions
     * Mainly to allow "reinitialize" possibility (?)
     */
    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              ///< seed for random shuffle
        int numberPaths,   ///< Number of paths
        int maxDimension,        ///< max dimensions to set up for
        int sizePartition,       ///< number of partitions to break up in.
        SequenceSP useSequence,  ///< Sequence to be used. Base class pointer.
        int UniformityFlag);     ///<  RNG dist to use for padding (?) GET_UNIFORM | GET_GAUSSIAN

    //========================================================
    // Must specify primary partition sequence to use through
    // the Sequence pointer given by useSequence pointer.
    // number of partitions passed in as opposed to above
    // constructor.

    SuperCube(INormalSuperCubeRNGSP rng,              // seed for random shuffle
              int numberPaths,    // Number of paths
              int maxDimension,        // max dimensions to set up for
              int sizePartition,       // size of partitions to break up in.
              int numPartitions /*= -1*/,       // number of partitions
              int seqDimension,        // dimensionality of low-disc sequence : < sizePartition.
              SequenceSP useSequence,  // Sequence to be used. Base class pointer.
              int UniformityFlag);     // Flag to set uniformity

    virtual void initializeObject(INormalSuperCubeRNGSP rng,            // seed for random shuffle
                                  int numberPaths,   // Number of paths
                                  int maxDimension,        // max dimensions to set up for
                                  int sizePartition,       // size of partitions
                                  int numPartitions,       // number of partitions to break up in.
                                  int seqDimension,        // Dimension of low-descrepency sequence
                                  SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                                  int UniformityFlag);

    //========================================================
    // Recommended version to use that avoids the complexity
    // of using the above constructors.

    SuperCube(INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
              int numberPaths,         // number of paths of data desired.
              int maxDimension,        // Overall dimensionality of the SuperCube.
              int sizePartition,       // Dimensionality of partition within SuperCube.
              int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
              int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
              int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
              // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
              int numBBFactors);       // number of dimensions for a single BB deviate at one date.
    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
        int numberPaths,         // number of paths of data desired.
        int maxDimension,        // Overall dimensionality of the SuperCube.
        int sizePartition,       // Dimensionality of partition within SuperCube.
        int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
        int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
        int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
        // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
        int numBBFactors);       // number of dimensions for a single BB deviate at one date.

    //========================================================
    // UAYOR version.
    // Allows for the passage of an user-specified final grid
    // for the Brownian bridge.

    SuperCube(
        INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
        int numberPaths,         // number of paths of data desired.
        int maxDimension,        // Overall dimensionality of the SuperCube.
        int sizePartition,       // Dimensionality of partition within SuperCube.
        int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
        int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
        int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
        // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
        int nBBfactors,          // number of dimensions for a single BB deviate at one date.
        double *BBptr2GaussianGrid);    // Pointer to gaussian grid for Brownian bridge. If NULL internal grid is constructed.

    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
        int numberPaths,         // number of paths of data desired.
        int maxDimension,        // Overall dimensionality of the SuperCube.
        int sizePartition,       // Dimensionality of partition within SuperCube.
        int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
        int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
        int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
        // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
        int nBBfactors,        // number of dimensions for a single BB deviate at one date.
        double *BBptr2GaussianGrid);    // Pointer to gaussian grid for Brownian bridge. If NULL internal grid is constructed.




    //========================================================
    // Destructor and other methods
    //========================================================

    virtual ~SuperCube();

    virtual void populatePrimaryPartition();

    virtual void populateVector();

    virtual void populateVector(int nPaths);

    //====================================================================
    // The following method advances the low-discrepancy sequence used in
    // the SuperCube brownian bridge if used. Basically steps up by the
    // given number of iterations.
    virtual void advanceBBLDSequence(const int numberIterations);
    //====================================================================
    // The following method advances the pseudo-random sequence used in
    // the SuperCube brownian bridge. Basically steps up by the given number of
    // iterations.
    virtual void advancePSRSequence(const int numberIterations);
    //====================================================================
    // The following method reinitializes the Supercube with a new seed.
    // Internally it retains memory allocated and reinitializes to use
    // the new seed being passed in.
    virtual void reinitialize(INormalSuperCubeRNGSP rng);
    //====================================================================
    // The following method reinitializes the Supercube with a new seed.
    // Internally it retains memory allocated and reinitializes to use
    // the new seed being passed in. This method also advances the Sobol
    // sequence used internally.
    virtual void reinitialize(INormalSuperCubeRNGSP rng,
                              const int numberLDSIterations);
    //====================================================================
    // The following method reinitializes the Supercube with a new seed,
    // number of LDS iterations to advance by and a with a prespecified
    // gaussian grid.
    virtual void reinitialize(INormalSuperCubeRNGSP rng,
                              const int numberLDSIterations,
                              const double *thisGaussianGrid);

protected:
    int *bridgeLocations; // Number of bridge locations
    int numPaths;        // Number of paths
    int MaxDimensions;   // Size of entire supercube
    int nPartitions;  // Number of partitions.
    int partitionSize;  // Size of partition set within the constructor.
    int sequenceDimension; // Size of vector to be sampled from Sequence passed in.

    int uniformityFlag;  // Flag to set uniform sequence(GET_UNIFORM = 0) or gaussian (GET_GAUSSIAN = 1)
    int numBBdates;      // Number of Brownian bridge dates if used.
    int numBBFactors;      // Number of Brownian bridge factors if used.
    double *BBgaussianGrid; // Pointer to grid for brownian bridge
    long **dimTable;     // Table needed for brownian bridge factors if used.
    int ldsBBcount;      // Count of number of LDS dimensions used in BB

    SequenceSP internalBasePtr; // Pointer to base sequence if set up internally.
    SequenceSP baseSequence; // Pointer to base class of sequence to be used.
    SequenceSP BBLDSequence; // Pointer to sequence for ld sequence for Br.Bridge
    double **primary;
    int **indexTable;

    virtual void initializeIndexTable(IUniformSuperCubeRNGSP uniform);
    double *deviatesByPath;

    //  long m_endSeed;            // Storage for seed at the end of primary partition computation
    /*  void setSeed( long thisSeed) {seed =thisSeed;}
      long seed;
    public:
        long getSeed() const {return seed;} // FIXME*/
protected:
    // FIXME: this is not right, but the old code switches between uniform/normal in a way that makes it dangerous to assume that we can guess type from the uniformityFlag
    INormalSuperCubeRNGSP getNormalRNG(); // returns rng or creates a normal one based on the RNG
    IUniformSuperCubeRNGSP getUniformRNG(); // returns a uniform RNG extracted from the rgn.
    void resetRNGGen(INormalSuperCubeRNGSP u);

private:
    INormalSuperCubeRNGSP nrng;

    INormalSuperCubeRNGSP normal;
    IUniformSuperCubeRNGSP uniform;


    void storeRNGGen(void);
    INormalSuperCubeRNGSP getStoredRNGGen() const;
    INormalSuperCubeRNGSP storedRNGGen; // clone of "rng" at the   end of primary partition computation (ex-m_endSeed)

};

DECLARESP(SuperCube);

CORE_END_NAMESPACE

#endif // !defined(_SC_SUPERCUBE_H)
