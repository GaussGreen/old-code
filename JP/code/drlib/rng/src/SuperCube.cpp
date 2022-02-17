// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                     *
 * SuperCube.h: implementation for the SuperCube class.                *
 * Author: M.Huq, Credit DR                                            *
 *                                                                     *
 *  Basic idea here is create a block of data of the size of samples   *
 *  that we need. In our case this is a vector of doubles. Take this   *
 *  vector and partition it into numPartitions. In the first partition *
 *  use your favorite sequence. In the subsequent partitions permute   *
 *  the order of the numbers and shuffle.                              *
 *  $Name$
 -----------------------------------------------------------------------*/


#include "edginc/SuperCube.h"
#include "edginc/SCException.h"

#include "edginc/PseudoRandomSequence.h"
#include "edginc/SobolSequence.h"
#include "edginc/BrownianBridgeSequence.h"
#include "VectorMatrix.h" // tnti_* types
#include "edginc/SuperCubeHelper.h"
#include "edginc/TrivialTransform.h"
#include "InvertCumNormal.h"

#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>

double SC_InvertCumNormal (double x);

CORE_BEGIN_NAMESPACE

using namespace std;


//////////////////////////////////////////////////////////////////////
// Constructors
// The following set of constructors are listed in pairs.
// For each constructor there is an equivalent initializeObject
// method. This done so as to allow for a default constructor
// declaration.
//////////////////////////////////////////////////////////////////////

//====================
// Default constructor
//====================
SuperCube::SuperCube()
{
    //        uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = SequenceSP();
    baseSequence = SequenceSP();
    BBLDSequence = SequenceSP();
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    //        m_endSeed = 0;
}

#if 0
//==============================================
//Constructor to support original functionality.
//Grows the super-cube to fit partitions that
//dont fit.
//==============================================
SuperCube::SuperCube(INormalSuperCubeRNGSP rng,              // seed for random shuffle
                     int numberPaths,         // Number of paths
                     int maxDimension,        // max dimensions to set up for
                     int sizePartition,       // size of partitions to break up in.
                     SequenceSP useSequence) // Sequence to be used. Base class pointer.

{
    uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = 0;
    baseSequence = 0;
    BBLDSequence = 0;
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    internalBasePtr = 0;  // set internal base pointer to null since passed in.
    uniformityFlag = GET_UNIFORM;
    //        m_endSeed = 0;

    initializeObject(uniform,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     useSequence);
} //SuperCube::SuperCube
void SuperCube::initializeObject(
    INormalSuperCubeRNGSP rng,              // seed for random shuffle
    int numberPaths,          // Number of paths
    int maxDimension,        // max dimensions to set up for
    int sizePartition,       // size of partitions
    SequenceSP useSequence) // Sequence to be used. Base class pointer.

{
    uniformityFlag = GET_UNIFORM; // Assume original functionality
    partitionSize = sizePartition;
    MaxDimensions = maxDimension;
    //        m_endSeed = 0;

    if(maxDimension % sizePartition != 0)
    {
        cout << "SuperCube: maxDimension is not a multiple of size of partitions. Resetting.\n";
        // make the maximum Dimension divisible by partitionSize with no remainder.
        MaxDimensions = MaxDimensions - (MaxDimensions % sizePartition) + sizePartition;
        cout << "SuperCube::SuperCube : Changed MaxDimensions to " << MaxDimensions << endl;
    }
    nPartitions = MaxDimensions / partitionSize;
    sequenceDimension = useSequence->getDimension();
    // call the core method.
    initializeObject(uniform,
                     numberPaths,
                     MaxDimensions,
                     partitionSize,
                     nPartitions,
                     sequenceDimension,
                     useSequence,
                     uniformityFlag);
}

#endif
//==============================================
SuperCube::SuperCube(INormalSuperCubeRNGSP rng,              // seed for random shuffle
                     int numberPaths,         // Number of paths
                     int maxDimension,        // max dimensions to set up for
                     int sizePartition,       // size of partitions to break up in.
                     SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                     int UniformityFlag)     // Flag to set uniformity
{
    //        uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = SequenceSP();
    baseSequence = SequenceSP();
    BBLDSequence = SequenceSP();
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    internalBasePtr = SequenceSP();  // set internal base pointer to null since passed in.
    //        m_endSeed = 0;
    initializeObject(rng,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     useSequence,
                     UniformityFlag);
} // SuperCube::SuperCube
void SuperCube::initializeObject(INormalSuperCubeRNGSP rng,              // seed for random shuffle
                                 int numberPaths,         // Number of paths
                                 int maxDimension,        // max dimensions to set up for
                                 int sizePartition,       // size of partitions
                                 SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                                 int UniformityFlag)
{
    uniformityFlag = UniformityFlag;
    //        m_endSeed = 0;
    partitionSize = sizePartition;
    MaxDimensions = maxDimension;
    if(maxDimension % sizePartition != 0) {
        cout << "SuperCube: maxDimension is not a multiple of size of partitions. Resetting.\n";
        // make the maximum Dimension divisible by partitionSize with no remainder.
        MaxDimensions = MaxDimensions - (MaxDimensions % sizePartition) + sizePartition;
        cout << "SuperCube::SuperCube : Changed MaxDimensions to " << MaxDimensions << endl;
    }
    if(partitionSize == 0) {
        nPartitions = 0;
    } else {
        nPartitions = MaxDimensions / partitionSize;
    }
    sequenceDimension = useSequence->getDimension();
#ifdef SCBB_DEBUG

    cout << "sequenceDimension = " << sequenceDimension << endl;
#endif // SCBB_DEBUG
    // call the core method.
    initializeObject(rng,
                     numberPaths,
                     MaxDimensions,
                     partitionSize,
                     nPartitions,
                     sequenceDimension,
                     useSequence,
                     UniformityFlag);
} //SuperCube::initializeObject

//==========================================================
// Following constructor is a new interface to the SuperCube
// that internally specifies whether or not to use
//  0) All pseudo-random by setting low-discrepancy dimension
//     to zero
//  1) Have a combination of low-discrepancy and
//     pseudo-random numbers in a primitive.
//==========================================================
SuperCube::SuperCube(INormalSuperCubeRNGSP rng, // seed for pseudo-random deviates used in both shuffling and PSR part.
                     int numberPaths,         // number of paths of data desired.
                     int maxDimension,        // Overall dimensionality of the SuperCube.
                     int sizePartition,       // Dimensionality of partition within SuperCube.
                     int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
                     int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
                     int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
                     // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
                     int nBBFactors)         // number of dimensions for a single BB deviate at one date.
{
    //        uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = SequenceSP();
    baseSequence = SequenceSP();
    BBLDSequence = SequenceSP();
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    //        m_endSeed = 0;

    if (numPartitions == -1)
        nPartitions = (partitionSize == 0) ? 0 : MaxDimensions / partitionSize;

    initializeObject(rng,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     numPartitions,
                     lowDiscDimension,
                     brownianBridgeFlag,
                     nBBFactors);
}



void SuperCube::initializeObject(INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
                                 int numberPaths,         // number of paths of data desired.
                                 int maxDimension,        // Overall dimensionality of the SuperCube.
                                 int sizePartition,       // Dimensionality of partition within SuperCube.
                                 int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
                                 int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
                                 int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
                                 // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
                                 int nBBfactors)        // number of dimensions for a single BB deviate at one date.
{
    //        BBgaussianGrid = 0;
    //        m_endSeed = 0;
    initializeObject(rng,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     numPartitions,
                     lowDiscDimension,
                     brownianBridgeFlag,
                     nBBfactors,
                     BBgaussianGrid = NULL);
}
//==========================================================
// This version of the constructor takes in a pointer to
// an array containing the grid.
//==========================================================
SuperCube::SuperCube(INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
                     int numberPaths,         // number of paths of data desired.
                     int maxDimension,        // Overall dimensionality of the SuperCube.
                     int sizePartition,       // Dimensionality of partition within SuperCube.
                     int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
                     int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
                     int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
                     // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
                     int nBBfactors,          // number of dimensions for a single BB deviate at one date.
                     double *BBptr2GaussianGrid)  // Pointer to gaussian grid for Brownian bridge. If NULL internal grid is constructed.
{
    //        uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = SequenceSP();
    baseSequence = SequenceSP();
    BBLDSequence = SequenceSP();
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    //        m_endSeed = 0;
    initializeObject(rng,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     numPartitions,
                     lowDiscDimension,
                     brownianBridgeFlag,
                     nBBfactors,
                     BBptr2GaussianGrid);
}
void SuperCube::initializeObject(INormalSuperCubeRNGSP _nrng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
                                 int numberPaths,         // number of paths of data desired.
                                 int maxDimension,        // Overall dimensionality of the SuperCube.
                                 int sizePartition,       // Dimensionality of partition within SuperCube.
                                 int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
                                 int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
                                 int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
                                 // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
                                 int nBBfactors,        // number of dimensions for a single BB deviate at one date.
                                 double *BBptr2GaussianGrid)    // Pointer to gaussian grid for Brownian bridge. If NULL internal grid is constructed.
{
    int iFacts, iDates;
    // Check to make sure dimensions are correct
    if(maxDimension < sizePartition * numPartitions)
    {
        throw SCException(__FILE__, __LINE__,
                          "maxDimension must be greater than sizePartition * numPartitions");
    }
    if(lowDiscDimension > sizePartition)
    {
        throw SCException(__FILE__, __LINE__,
                          "lowDiscDimension cannot exceed sizePartition");
    }
    nrng = _nrng;
    numPaths = numberPaths;
    MaxDimensions = maxDimension;
    nPartitions = numPartitions;
    partitionSize = sizePartition;
    nDimensions = maxDimension;
    numBBFactors = nBBfactors;
    BBgaussianGrid = BBptr2GaussianGrid;
    //        m_endSeed = 0;

    // Check if we are to set up a brownian bridge.
    if(brownianBridgeFlag != 0)
    {
        // Check input parameter needed for brownian bridge
        if(numBBFactors < 1 && lowDiscDimension % numBBFactors != 0) {
            throw SCException(__FILE__, __LINE__,
                              "numBBFactors must be an integer > 0 which is a factor of lowDiscDimension");
        }
        // Algorithm for specifying bridge location is based upon a binary-tree construct.
        // Given T slices split into 2. Store T/2
        // Then split remaining two parts into sets of two more. Start at higher end in storing.
        // Thus ratios laid out as follows:
        // T, T/2, 3T/4, T/4, 7T/8, 5T/8 and so forth until count reaches the brownian bridge count.
        if(brownianBridgeFlag > lowDiscDimension && brownianBridgeFlag < 2 ) {
            throw SCException(__FILE__, __LINE__,
                              "Number of brownian bridges must be between 2 and lowDiscDimension");
        }

        bridgeLocations = new int [brownianBridgeFlag];
        if(!bridgeLocations) {
            throw SCException(__FILE__, __LINE__,
                              "Could not allocate memory for bridgeLocations");
        }

        int nBBdates = (int)(lowDiscDimension/numBBFactors);
        numBBdates = nBBdates;

        bridgeLocations[0] = numBBdates-1; // T point
        int bcount = 1; // T point is already secured.
        double ratios[1000];
        if(brownianBridgeFlag > 1) {
            int level =0;
            double inv_pow_val = 1.;
            int pow_val = 1;
            while(bcount < brownianBridgeFlag) {
                level++;
                inv_pow_val = inv_pow_val * 0.5;
                for(int m = pow_val; m>0; m--) {
                    ratios[bcount-1]=((2.0*m -1.0) * inv_pow_val);
                    bcount++;
                }
                pow_val *= 2;
            }//bcount
            for(int iBridge=1; iBridge < brownianBridgeFlag; iBridge++) {

                bridgeLocations[iBridge] = (int)(ratios[iBridge-1] * (numBBdates-1));
                //cout << ratios[iBridge] << " " << bridgeLocations[iBridge] << endl;
            }//iBridge
        }

        // Insert sobol sequences for each of the bridges
        // ensuring that we don't exceed the maximum allowable sobol
        dimTable = (long **)malloc(sizeof(long *) * numBBFactors);
        for(int i=0;i<numBBFactors; i++) {
            dimTable[i] = (long *)malloc(sizeof(long)*nBBdates);
            if(!dimTable[i]) {
                throw SCException(__FILE__, __LINE__,
                                  "Could not allocate memory for dimTable");
            }
        } // i
        // Load up pseudo-random for all by default
        for(iFacts = 0 ; iFacts < numBBFactors; iFacts++) {
            for(iDates = 0; iDates < numBBdates; iDates++)
                dimTable[iFacts][iDates] = BB_PSEUDO_RANDOM;
        } // iFact


        iFacts = 0;
        int sobolCount = 0;
        // Loop over factors on the outside
        while(sobolCount < MAX_SOBOL && iFacts < numBBFactors) {
            iDates = 0;
            while(iDates < numBBdates && sobolCount < MAX_SOBOL) {
                int iBridge = 0;
                while(iBridge < brownianBridgeFlag && sobolCount < MAX_SOBOL) {
                    if(iDates == bridgeLocations[iBridge]) {
                        dimTable[iFacts][iDates] = BB_SOBOL;
                        sobolCount++;
                    }
                    iBridge++;
                } //iBridge
                iDates++;
            } // iDates
            iFacts++;
        } // iFacts

        // Apply grid at end for factor 0
#ifndef SC_NO_GRID
        dimTable[0][nBBdates-1] = BB_GRID;
#endif // SC_NO_GRID

        // Count up number of Sobol dimensions to be used.
        // store this for when we reinitialize.
#ifdef SC_NO_GRID

        ldsBBcount = sobolCount;
#else

        ldsBBcount = sobolCount-1;
#endif  // SC_NO_GRID
        //cout << "LDS count = " << ldsBBcount << endl;


        // Create a brownian bridge sequence
        if(ldsBBcount == 0)
            ldsBBcount = 1;
        SobolSequenceSP bridgeSobSeq(new SobolSequence(-1, ldsBBcount)); // Sobol doesn't use seed
        BBLDSequence = bridgeSobSeq;



        BrownianBridgeSequenceSP bSeq (new BrownianBridgeSequence(numBBFactors,
                                       numBBdates,
                                       dimTable,
                                       numPaths,
                                       getNormalRNG(),
                                       BBgaussianGrid,
                                       BBLDSequence));
        //seed = bSeq->getSeed(); // FIXME not needed now, as we pass rng by ref.
        //************************
        // Since base Sequence is created internally we need to actively
        // delete it as well in the destructor. This pointer points to the
        // appropriate area of memory to delete.
        internalBasePtr = bSeq;
        //************************
        // Initialize the supercube
        initializeObject(nrng,
                         numberPaths,
                         maxDimension,
                         sizePartition,
                         numPartitions,
                         lowDiscDimension,
                         bSeq,
                         GET_BB_STYLE);
        //Free bridgeLocations
        delete [] bridgeLocations;
        bridgeLocations = 0;
    } else
    {
        numBBdates = 0;
        numBBFactors = 0;
        // Create a generic Sobol sequence.
        long lowDiscSeed = -1; // Right now not used in Sobol.
        if(lowDiscDimension != 0) {
            SobolSequenceSP sSeq(new SobolSequence(lowDiscSeed, lowDiscDimension));
            // Initialize the super-cube. This generates the sobol sequence and pseudo-random sequences
            // as needed.
            initializeObject(nrng,
                             numberPaths,
                             maxDimension,
                             sizePartition,
                             numPartitions,
                             lowDiscDimension,
                             sSeq,
                             GET_GAUSSIAN);
        } else {
            PseudoRandomSequenceSP sSeq(new PseudoRandomSequence(getNormalRNG() /*lowDiscSeed*/, 1));  // Do nothing sequence
            // Initialize the super-cube. This generates the sobol sequence and pseudo-random sequences
            // as needed.
            initializeObject(nrng,
                             numberPaths,
                             maxDimension,
                             sizePartition,
                             numPartitions,
                             lowDiscDimension,
                             sSeq,
                             GET_GAUSSIAN);
        }
    }// brownianBridgeFlag

} // SuperCube::initializeObject


//==============================================
//Lower level constructor in new framework
//initializeObject here is the core method called
//by all above.
//==============================================
SuperCube::SuperCube(INormalSuperCubeRNGSP rng,              // seed for random shuffle
                     int numberPaths,         // Number of paths
                     int maxDimension,        // max dimensions to set up for
                     int sizePartition,       // size of partitions to break up in.
                     int numPartitions,       // number of partitions
                     int seqDimension,        // dimensionality of low-disc sequence : < sizePartition.
                     SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                     int UniformityFlag)      // Flag to set uniformity
{
    //         uniform = NULL;
    bridgeLocations = 0;
    numPaths = 0;
    MaxDimensions = 0;
    nPartitions = 0;
    partitionSize = 0;
    sequenceDimension = 0;
    internalBasePtr = SequenceSP();
    baseSequence = SequenceSP();
    BBLDSequence = SequenceSP();
    uniformityFlag = 0;
    dimTable = 0;
    ldsBBcount = 0;
    BBgaussianGrid = 0;
    numBBdates = 0;
    numBBFactors = 0;
    deviatesByPath = 0;
    indexTable = 0;
    //        m_endSeed = 0;
    initializeObject(rng,
                     numberPaths,
                     maxDimension,
                     sizePartition,
                     numPartitions,
                     seqDimension,
                     useSequence,
                     UniformityFlag);
} //SuperCube::SuperCube
// core method.
void SuperCube::initializeObject(INormalSuperCubeRNGSP _nrng,              // seed for random shuffle
                                 int numberPaths,         // Number of paths
                                 int maxDimension,        // max dimensions to set up for
                                 int sizePartition,       // Size of partitions.
                                 int numPartitions,      // number of partitions to break up in.
                                 int seqDimension,       // dimension of low discrepancy sequence.
                                 SequenceSP useSequence,
                                 int UniformityFlag) // Sequence to be used. Base class pointer.

{
    nrng = _nrng;
    numPaths = numberPaths;
    MaxDimensions = maxDimension;
    partitionSize = sizePartition;
    nPartitions = numPartitions;
    sequenceDimension = seqDimension;
    baseSequence = useSequence;
    uniformityFlag = UniformityFlag;
    //        m_endSeed = 0;

    // Set vector size used by populateVector/getVector.
    nDimensions = MaxDimensions;

    if (numPaths < 0)
        throw SCException(__FILE__, __LINE__,
                          "Negative number of paths not allowed");

    deviatesByPath = new double [numPaths];
    if(deviatesByPath == 0)
    {
        throw SCException (__FILE__, __LINE__,
                           "Could not allocate memory for devaitesByPath");
    }

    // Allocate memory for the primary partition.
    primary = new double * [numPaths];
    for(int i = 0; i < numPaths; i++)
    {
        primary[i] = new double [partitionSize];
    }

    //     if(!vector)
    //  vector = new double [nDimensions];
    //     if(!vector){
    //  throw SCException(__FILE__, __LINE__, "Could not allocate memory for vector");
    //     }
    vector.assign(nDimensions, 0);
    // Initialize index table
    if(nPartitions !=0)
    {
        initializeIndexTable(getUniformRNG());

        // Set up partition 0 ; shuffling done during population of vectors.
        populatePrimaryPartition();
    }
}


//================================
//Destructor
//================================
SuperCube::~SuperCube()
{
    int iPath;
    //if(vector){
    //  delete vector;
    //  vector = NULL;
    //}
    // Empty arrays
    if(nPartitions != 0) {
        if(indexTable) {
            for(iPath = 0; iPath < numPaths; iPath++) {
                delete [] indexTable[iPath];
            }//iPath
            delete [] indexTable;
            indexTable = 0;
        }//indexTable
        if(primary) {
            for(int iPart=0; iPart < numPaths; iPart++) {
                delete [] primary[iPart];
            }//iPart
            delete [] primary;
            primary = 0;
        }//primary
    }
    if(BBLDSequence.get() != NULL) {
        //delete static_cast<SobolSequence *>(BBLDSequence);
        //                delete BBLDSequence;
        BBLDSequence = SequenceSP();
    }


    if(baseSequence.get() != NULL && internalBasePtr.get() != NULL ) {
        //static_cast<BrownianBridgeSequence *>(baseSequence)->~BrownianBridgeSequence();
        // delete baseSequence;
        internalBasePtr = SequenceSP();// Points to same area of memory.
        baseSequence = SequenceSP();
    }

    if(dimTable) {
        for(int i=0;i<numBBFactors; i++) {
            if(dimTable[i]) {
                free(dimTable[i]);
                dimTable[i] = 0;
            }
        }
        free(dimTable);
        dimTable = 0;
    }
    if(deviatesByPath) {
        delete [] deviatesByPath;
        deviatesByPath = 0;
    }
    //if(bridgeLocations){
    //      delete [] bridgeLocations;
    //      bridgeLocations = 0;
    //  }


}

//===================================
//Populates the primary partition
//from which all scrambled partitions
//are derived.
//===================================
void SuperCube::populatePrimaryPartition()
{
    int iPath;

    if(uniformityFlag == GET_GAUSSIAN) {
        // For each path populate the primary partition.
        for(iPath = 0 ; iPath < numPaths; iPath++) {
            // Populate the first sub-grid with the sequence passed in.
            if(sequenceDimension!=0) {
                baseSequence->populateVector(iPath);
                double *RNvector = baseSequence->getVector();
                for(int i = 0; i < sequenceDimension; i++) {
                    primary[iPath][i] = SC_InvertCumNormal(RNvector[i]);
                } //i
            }
            // For the rest of the sequence put in pseudo-random.
            for(int j = sequenceDimension; j < partitionSize; j++) {
                primary[iPath][j] = getNormalRNG()->fetch(); //SC_gasdev2(&seed); // NORMAL
            }
        } // iPath

    } else if(uniformityFlag == GET_GAUSSIAN_PADDING) { // Leaves low-discrepancy part alone
        // For each path populate the primary partition.
        for(iPath = 0 ; iPath < numPaths; iPath++) {
            // Populate the first sub-grid with the sequence passed in.
            if(sequenceDimension!=0) {
                baseSequence->populateVector(iPath);
                double *RNvector = baseSequence->getVector();
                for(int i = 0; i < sequenceDimension; i++) {
                    primary[iPath][i] = (RNvector[i]);
                } //i
            }
            // For the rest of the sequence put in pseudo-random.
            for(int j = sequenceDimension; j < partitionSize; j++) {
                primary[iPath][j] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);          // NORMAL
            }
        } // iPath


    } else if(uniformityFlag == GET_BB_STYLE) {

        // Brownian bridge special case
        BrownianBridgeSequenceSP BBptr = DYNAMIC_POINTER_CAST<BrownianBridgeSequence > (baseSequence);
        //BBptr = (BrownianBridgeSequence *)baseSequence;
        numBBFactors = BBptr->getNumberBBFactors();
        numBBdates = BBptr->getNumberBBDates();
        if(numBBdates == 0) {
            throw SCException(__FILE__, __LINE__,
                              "Number of brownian bridge dates is zero.");
        }
        if(numBBdates * numBBFactors != sequenceDimension) {
            cout << __FILE__ << " "
            <<__LINE__
            << " WARNING: numBBdates * numBBFactors not equal to sequenceDimension" << endl;
            cout << "----numBBdates = " << numBBdates << endl;
            cout << "----numBBFactors = " << numBBFactors << endl;
            cout << "----sequenceDimension = " << sequenceDimension << endl;
        }
        // Assume first numBBdates by numBBFactors are the low-discrepancy region
        int primIndex;
        for(int iFact = 0; iFact < numBBFactors; iFact++) {
            for(int iDate = 0; iDate < numBBdates; iDate++) {
                BBptr->getDeviatesbyFactorDate(iFact, iDate, deviatesByPath);
                // Now have an array full of path-wise deviates for given date-factor
                // Translate
                primIndex = iFact + numBBFactors * iDate;
                for(int k=0; k<numPaths; k++) {
                    primary[k][primIndex] = deviatesByPath[k];
                }
            }//iFact
        }//iDate


        // Now assume that brownian bridge is filled in. Continue to fill in rest if there is
        // For the rest of the sequence put in pseudo-random.
        for(iPath=0; iPath < numPaths; iPath++) {
            for(int j = sequenceDimension; j < partitionSize; j++) {
                primary[iPath][j] = getNormalRNG()->fetch(); //SC_gasdev2(&seed); // NORMAL
            }
        } //iPath

    } else if (uniformityFlag == GET_UNIFORM) {
        // For each path populate the primary partition.
        for(iPath = 0 ; iPath < numPaths; iPath++) {
            // Populate the first sub-grid with the sequence passed in.
            double *RNvector;
            if(sequenceDimension!=0) {
                baseSequence->populateVector(iPath);
                RNvector = baseSequence->getVector();
                std::copy(RNvector, RNvector+ sequenceDimension, primary[iPath]);
                //                                 for(int i = 0; i < sequenceDimension; i++) {
                //                                         primary[iPath][i] = (RNvector[i]);
                //                                 } //i
            }
            // For the rest of the sequence put in pseudo-random.
            for(int j = sequenceDimension; j < partitionSize; j++) {
                primary[iPath][j] = getUniformRNG()->fetch(); //SC_ran2(&seed);                 // UNIFORM
            }

        } // iPath

    } else {
        assert(! "Unhandled branch"); // uniformityFlag
    }
    // m_endSeed= seed;
    storeRNGGen();
} // End populatePrimaryPartition

void SuperCube::storeRNGGen(void)
{
    // endRNG = INormalSuperCubeRNGSP(rng->clone()); // store RNG at the end of this stage, to allow per path fetching later.
    /// FIXME: in the original code the seed is always -1 here ?!
    storedRNGGen = DYNAMIC_POINTER_CAST<INormalSuperCubeRNG>(nrng->clone()); // (new SC_ran2(-1L));
    //    nrng->debug();
    //storedRNGGen = SC_ran2Gen::create(-1L); // move this logic to seekToPath() impl.
    //storedRNGGen->debug();

}


INormalSuperCubeRNGSP SuperCube::getStoredRNGGen(void) const
{
    return storedRNGGen;
}
//==================================
//Overload base-class populateVector
//==================================
void SuperCube::populateVector()
{
    throw SCException(__FILE__, __LINE__,
                      "Method will not work with SuperCube. Use version with iPath passed in.");
} // End populateVector


/*********************************************************
 Function to populate the internal vector accessible by
 getVector().
 It uses the indexTable to figure out the shuffling based
 upon the primary partition.
*********************************************************/

void SuperCube::populateVector(int iPath)
{
    /*
            long seedStride = 1000;
            long seedMod = 16061968;*/

    // Fetch deviates that belong in the primary or secondary partitions.
    int vIndex = 0;
    for(int iPartition = 0 ; iPartition < nPartitions; iPartition++) {
        int shuffle_index = indexTable[iPath][iPartition];
        std::copy(primary[shuffle_index],
                  primary[shuffle_index] + partitionSize,
                  & vector[iPartition*partitionSize]);
        //                 for(int i = 0; i < partitionSize; i++) {
        //                         //vIndex = i + (partitionSize * iPartition);
        //                         vector[vIndex] =primary[shuffle_index][i];
        //                         vIndex++;
        //                 }
    }  // iPartition
    vIndex = nPartitions * partitionSize;
    // If there is a padded region fetch it here.
    if(vIndex < MaxDimensions) {
#if 0
        // Set up seed for each path using deterministic formula
        seed = (m_endSeed +  seedStride* (long)(1+iPath) )% seedMod + iPath;
        if(seed > 0)
            seed *= -1;
        // reinitialize gasdev2
        SC_init_gasdev2();
        // For the rest of the super-cube pad with pseudo-random
        //vIndex = nPartitions * partitionSize;
        if(uniformityFlag == GET_UNIFORM) {
            while(vIndex < MaxDimensions) {
                vector[vIndex] = SC_ran2(&seed);
                vIndex++;
            }

        } else {
            while(vIndex < MaxDimensions) {
                vector[vIndex] = SC_gasdev2(&seed);
                vIndex++;
            }
        } //uniformityFlag
#endif
        /*INormalSuperCubeRNGSP*/
        // nrng = INormalSuperCubeRNGSP(getStoredRNGGen()->clone());
        // nrng->debug();
        // ISuperCubeRNGSP rng;
        // FIXME we should store RNG, and pre
        /** original code again;
            to seek to a path we reset the current generator and derive RNG according to a magic flag.
        */
        resetRNGGen( DYNAMIC_POINTER_CAST<INormalSuperCubeRNG>(getStoredRNGGen()->clone()) );
        ISuperCubeRNGSP rng;
        if (uniformityFlag == GET_UNIFORM)
            rng = getUniformRNG(); // SC_ran2::create(nrng);
        else // NORMAL
            rng = getNormalRNG(); // SC_gasdev2::create(nrng);

        rng->seekToPath(iPath);
        for(vIndex = nPartitions * partitionSize;
                vIndex < MaxDimensions;
                ++vIndex)
            vector[vIndex] = rng->fetch();
    }//if < MaxDimensions : that is if there is a padded region.
} // End populateVector(int iPath)

/*****************************************************************
 Function to initialize the index table.
 Idea here being that we preinitialize an index table to determine
 shuffling order. These indices are picked such that there is no
 repetition of any index along the same path for all partitions.
*****************************************************************/
void SuperCube::initializeIndexTable(IUniformSuperCubeRNGSP uniform)
{
    int iPath;
    int iPartition;
    // memory allocation
    if(!indexTable) {
        indexTable = new int * [numPaths];
        for(iPath = 0; iPath < numPaths; iPath++) {
            indexTable[iPath] = new int [nPartitions];
        }//iPath
    }

    // Base vector : which we will shuffle
    tnti_vec baseIndices(numPaths);
    // Assign first partition indices
    iPartition = 0;
    for(iPath = 0; iPath < numPaths; iPath++) {
        indexTable[iPath][iPartition] = iPath;
        baseIndices[iPath] = iPath;
    }

    // Shuffle the partitions that follow.
    for(iPartition = 1; iPartition < nPartitions; iPartition++) {
        // randomly reorder current partition. Repeat numPaths times.
        //FIXME wouldn't std::random_shuffle(start, end) be simpler ?!
        for(int pRep = 1; pRep < numPaths; pRep++) {
            // OLDCODE: int rPath = (int)(pRep*SC_ran2(&seed));
            int rPath = (int)(pRep*uniform->fetch());
            std::swap(baseIndices[rPath], baseIndices[pRep]);
#if 0

            int tmp = baseIndices[rPath];
            baseIndices[rPath] = baseIndices[pRep];
            baseIndices[pRep] = tmp;
#endif

        } // pRep
    } // iPartition

    for(iPath = 0 ; iPath < numPaths; iPath++) {
        int tmp = baseIndices[iPath];
        for(iPartition = 1; iPartition < nPartitions; iPartition++) {
            // Set the baseindex to be current
            indexTable[iPath][iPartition] = tmp;
        } // iPath
    } // iPartition

    // clean up base storage
#ifndef USE_TNT
    baseIndices.empty();
#endif

} // End InitializeIndexTable

//====================================================================
// The following method advances the low-discrepancy sequence used in
// the SuperCubeBBridge. Basically steps up by the given number of
// iterations.
void SuperCube::advanceBBLDSequence(const int numberIterations)
{
    if(!BBLDSequence) {
        throw SCException(__FILE__, __LINE__,
                          "BBLDSequence = NULL : Check to see if you are using a Brownian bridge construct or if the Brownian bridge is initialized");
    }
    BBLDSequence->advanceSequence(numberIterations);
}// void SuperCube::advanceBBLDSequence(const int numberIterations)

//====================================================================
// The following method advances the pseudo-random sequence used in
// the SuperCubeBBridge. Basically steps up by the given number of
// iterations.
// Right now done by sampling N times with gasdev
void SuperCube::advancePSRSequence(const int numberIterations)
{
    //double tmpx ;
    for(int it = 0; it < numberIterations; it++) {
        nrng->fetch(); // tmpx = SC_gasdev2(&seed);
    }
}// void SuperCube::advancePSRSequence(const int numberIterations)

//=======================================================================
// Method that reinitializes the SuperCubeBBridge object to take a new
// seed.
void SuperCube::reinitialize(INormalSuperCubeRNGSP rng)
{
    // First check that this SuperCubeBBridge object is instantiated.
    // At least check that the primary partition is allocated.
    if(numPaths == 0 )
        throw SCException(__FILE__,__LINE__,
                          "reinitialize can only be used after an instance of SuperCube has been created using initializeObject");
    reinitialize(rng, 0);

}// void SuperCube::reinitSuperCubeBBridge(const long seed)
//====================================================================
//=======================================================================
// Method that reinitializes the SuperCubeBBridge object to take a new
// seed and advances the low-discrepancy sequence by numberIteartions
void SuperCube::reinitialize(INormalSuperCubeRNGSP gen,
                             const int numberLDSIterations,
                             const double *thisGaussianGrid)
{

    BBgaussianGrid = const_cast<double *>(thisGaussianGrid);
    // First check that this SuperCubeBBridge object is instantiated.
    // At least check that the primary partition is allocated.
    if(numPaths == 0 )
        throw SCException(__FILE__,__LINE__,
                          "reinitialize can only be used after an instance of SuperCube has been created using initializeObject");
    // Set the seed


    gen->superCubeAdjustment();
    resetRNGGen(gen);
    //        gen->debug();
    /*  _uniform->init();
       if( seed > 0 ){
        seed *= -1;
        }
        if(seed == 0)
        throw SCException(__FILE__,__LINE__,
                  "zero seed passed in. Will lead to problems in pseudo-random number generator");*/
    // check to see if BBLDSequence is defined - if it is then a BrownianBridge is being used.
    if(nPartitions != 0 && BBLDSequence.get() != NULL) {

        //delete static_cast<SobolSequence *>(BBLDSequence);
        //                static_cast<SobolSequence*>(BBLDSequence)->~SobolSequence();
        //                delete BBLDSequence;
        BBLDSequence = SequenceSP();
        //delete static_cast<BrownianBridgeSequence *>(baseSequence);
        //                static_cast<BrownianBridgeSequence*>(baseSequence)->~BrownianBridgeSequence();
        //                delete baseSequence;
        baseSequence = SequenceSP();

        // Have to call SC_gasdev2 with -ve seed to reinitialize.
        // Global variables.................
        // double tmp = SC_gasdev2(&seed);

        if(ldsBBcount == 0)
            ldsBBcount = 1;
        BBLDSequence = SobolSequenceSP(new SobolSequence(-1, ldsBBcount)); // SobolSequence doesn't use seed
        if(!BBLDSequence)
            throw SCException(__FILE__, __LINE__,
                              "Could not allocate SobolSequence");
        if(numberLDSIterations!=0)
            BBLDSequence->advanceSequence(numberLDSIterations);

        BrownianBridgeSequenceSP bSeq (new BrownianBridgeSequence(numBBFactors,
                                       numBBdates,
                                       dimTable,
                                       numPaths,
                                       getNormalRNG(),
                                       BBgaussianGrid,
                                       BBLDSequence));
        //               seed = bSeq->getSeed();
        // Set the seed in the brownian bridge sequence
        //                setSeed(seed);
        baseSequence = bSeq;
    }

    // reinitialize the SuperCube datastructures
    //======================
    // Initialize index table
    if(nPartitions !=0) {
        initializeIndexTable(getUniformRNG());

        // Set up partition 0 ; shuffling done during population of vectors.
        populatePrimaryPartition();
    }
}// void SuperCube::reinitialize(const long aSeed, const int numberLDSIterations, const double *grid);


//====================================================================
// The following method reinitializes the Supercube with a new seed,
// number of LDS iterations to advance by and a with a prespecified
// gaussian grid. Uses previous BBgrid

void SuperCube::reinitialize(INormalSuperCubeRNGSP rng,
                             const int numberLDSIterations)
{
    reinitialize(rng, numberLDSIterations,BBgaussianGrid);

}/* void SuperCube::reinitialize(const long aSeed,
    const int numberLDSIterations )
 */


IUniformSuperCubeRNGSP SuperCube::getUniformRNG()
{
    if (! uniform) {
        IRNGGeneratorSP gen = getNormalRNG()->getGenerator();
        ISuperCubeRNGGenSP g =  DYNAMIC_POINTER_CAST<ISuperCubeRNGGen>(gen);
        uniform = TrivialUniformSuperCubeRNG::create(g);
    }
    //        uniform = SC_ran2::create(nrng); // FIXME: convert to abstract fact.
    return uniform;
}
/** Convert input RNG into a normal oneAQR1
    If rng is already Normal -- do nothing;
    Otherwise, find RNG's uniform generator and convert to Normal via GasDev2
*/

INormalSuperCubeRNGSP SuperCube::getNormalRNG()
{
    if (! normal)
        normal = nrng;

    return normal;
}

void SuperCube::resetRNGGen(INormalSuperCubeRNGSP u)
{
    nrng = u;
    normal = INormalSuperCubeRNGSP();
    uniform = IUniformSuperCubeRNGSP();
}
CORE_END_NAMESPACE
