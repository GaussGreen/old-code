// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/* BrowianBridgeSequence.cpp                                                     *
 * Author: M. Huq                                                                *
 * Methods for the BrownianBridgeSequence class.                                 *
   ----------------------------------------------------------------------- */

#include "edginc/BrownianBridgeSequence.h"
#include "edginc/SobolSequence.h"

// Exception handling
#include "edginc/SCException.h"

#include "sobmc_grid_multi.h"

#include "InvertCumNormal.h"
#include "VectorMatrix.h"

#include <cstdlib>
#include <iostream>

CORE_BEGIN_NAMESPACE

using namespace std;

//************
// Default Constructor
//************
BrownianBridgeSequence::BrownianBridgeSequence()
{
    factorDimension = 0;
    numberDates = 0;
    nPaths = 0;
    deltaGrid = 0.0;
    gridLocation = 0;
    gridFactor = 0;
    listDates=0;
    gaussianVect=0;
    probabilities=0;
    factorTypeList =0;
    ptr2GaussianGrid = 0;
    lowDiscSeqPtr = SequenceSP();
    internalLowDiscSeqPtr = SequenceSP();

}

BrownianBridgeSequence::BrownianBridgeSequence(int dimensions,  // dimensionality of deviate that we are after.
        int numDates,   // number of dates or stratifications.
        long **dimTypeList,// type specifier for dimension by date.
        int numberPaths,  // number of paths that will be sampled.
        INormalSuperCubeRNGSP rngIn )    // Seed for pseudo-random part.
{
    factorDimension = 0;
    numberDates = 0;
    nPaths = 0;
    deltaGrid = 0.0;
    gridLocation = 0;
    gridFactor = 0;
    listDates=0;
    gaussianVect=0;
    probabilities=0;
    factorTypeList =0;
    ptr2GaussianGrid = 0;
    lowDiscSeqPtr = SequenceSP();
    internalLowDiscSeqPtr = SequenceSP();

    ptr2GaussianGrid = 0;
    initializeBrownianBridge(dimensions,  // dimensionality of deviate that we are after.
                             numDates,   // number of dates or stratifications.
                             dimTypeList,// type specifier for dimension by date.
                             numberPaths,  // number of paths that will be sampled.
                             rngIn ,    // Seed for pseudo-random part.
                             ptr2GaussianGrid);
} // BrownianBridgeSequence::BrownianBridgeSequence

BrownianBridgeSequence::BrownianBridgeSequence(int dimensions,  // dimensionality of deviate that we are after.
        int numDates,   // number of dates or stratifications.
        long **dimTypeList,// type specifier for dimension by date.
        int numberPaths,  // number of paths that will be sampled.
        INormalSuperCubeRNGSP rngIn ,    // Seed for pseudo-random part.
        double *gaussianGridPtr) // pointer to grid, gaussianized
{
    factorDimension = 0;
    numberDates = 0;
    nPaths = 0;
    deltaGrid = 0.0;
    gridLocation = 0;
    gridFactor = 0;
    listDates=0;
    gaussianVect=0;
    probabilities=0;
    factorTypeList =0;
    ptr2GaussianGrid = 0;
    lowDiscSeqPtr = SequenceSP();
    internalLowDiscSeqPtr = SequenceSP();

    initializeBrownianBridge(dimensions,  // dimensionality of deviate that we are after.
                             numDates,   // number of dates or stratifications.
                             dimTypeList,// type specifier for dimension by date.
                             numberPaths,  // number of paths that will be sampled.
                             rngIn ,    // Seed for pseudo-random part.
                             gaussianGridPtr);

} // BrownianBridgeSequence::BrownianBridgeSequence (int,int,long,int,long,double*)

//--------------------------------------------------------
// Method to initialize the Brownian bridge.
//--------------------------------------------------------
void BrownianBridgeSequence::initializeBrownianBridge(int dimensions,  // dimensionality of deviate that we are after.
        int numDates,   // number of dates or stratifications.
        long **dimTypeList,// type specifier for dimension by date.
        int numberPaths,  // number of paths that will be sampled.
        INormalSuperCubeRNGSP rngIn ,    // Seed for pseudo-random part.
        double *gaussianGridPtr) // pointer to grid, gaussianized

{
    // This version of the constructor assumes that a Sobol sequence will
    // be used as the low-discprency sequence part of the brownian bridge.
    // Thus we initialize a SobolSequence and appropriately assign it.
    // Count up the number of required Sobol dimensions first.
    int sobolCount = 0;
    for(int i = 0; i < dimensions; i++)
    {
        for(int j = numDates-1; j>=0; j--) {
            if(dimTypeList[i][j] == BB_SOBOL) {
                sobolCount++;
            }
        }
    }
    if(sobolCount == 0)
        sobolCount = 1;
    SobolSequenceSP theSobolSequence(new SobolSequence(-1, sobolCount));
    internalLowDiscSeqPtr = theSobolSequence;

    initializeBrownianBridge(dimensions,
                             numDates,
                             dimTypeList,
                             numberPaths,
                             rngIn ,
                             gaussianGridPtr,
                             theSobolSequence);
}

BrownianBridgeSequence::BrownianBridgeSequence(int dimensions,  // dimensionality of deviate that we are after.
        int numDates,   // number of dates or stratifications.
        long **dimTypeList,// type specifier for dimension by date.
        int numberPaths,  // number of paths that will be sampled.
        INormalSuperCubeRNGSP rngIn ,    // Seed for pseudo-random part.
        double *gaussianGridPtr, // pointer to grid, gaussianized
        SequenceSP ldsPtr) // Pointer to low-discrepancy sequence to be used

{
    factorDimension = 0;
    numberDates = 0;
    nPaths = 0;
    deltaGrid = 0.0;
    gridLocation = 0;
    gridFactor = 0;
    listDates=0;
    gaussianVect=0;
    probabilities=0;
    factorTypeList =0;
    ptr2GaussianGrid = 0;
    lowDiscSeqPtr = SequenceSP();
    internalLowDiscSeqPtr = SequenceSP();

    initializeBrownianBridge(dimensions,  // dimensionality of deviate that we are after.
                             numDates,   // number of dates or stratifications.
                             dimTypeList,// type specifier for dimension by date.
                             numberPaths,  // number of paths that will be sampled.
                             rngIn ,    // Seed for pseudo-random part.
                             gaussianGridPtr,
                             ldsPtr);

} // BrownianBridgeSequence::BrownianBridgeSequence (int,int,long,int,long,double*)

//--------------------------------------------------------
// Method to initialize the Brownian bridge.
//--------------------------------------------------------
void BrownianBridgeSequence::initializeBrownianBridge(int dimensions,  // dimensionality of deviate that we are after.
        int numDates,   // number of dates or stratifications.
        long **dimTypeList,// type specifier for dimension by date.
        int numberPaths,  // number of paths that will be sampled.
        INormalSuperCubeRNGSP rngIn ,    // Seed for pseudo-random part.
        double *gaussianGridPtr, // pointer to grid, gaussianized
        SequenceSP ldsPtr) // Pointer to low-discrepancy sequence to be used


{
    // Number of dimensions per date per path
    factorDimension = dimensions;
    numberDates     = numDates;
    nPaths = numberPaths;
    // set the date delta
    deltaGrid = 1.0 / (double)numberPaths;
    // This version of the constructor sets up the pointer to the gaussian grid as NULL.
    // This means that the next layer of brownian bridge code sets up the grid internally.
    ptr2GaussianGrid = gaussianGridPtr;
    // Set the low-discrepancy sequence pointer.
    lowDiscSeqPtr = ldsPtr;

    // Allocate memory for factorTypeList
    factorTypeList = new long * [factorDimension];
    for(int d = 0 ; d < factorDimension; d++)
    {
        factorTypeList[d] = new long  [numberDates];
        if(!factorTypeList[d]) {
            throw SCException( __FILE__, __LINE__,
                               "Failed to allocate memory for factorTypeList");
        }
    } // d
    long sobolCount = 0L;
    gridLocation = -99;
    gridFactor = -99;
    for(int i = 0; i < factorDimension; i++)
    {
        for(int j = numDates-1; j>=0; j--) {
            if(dimTypeList[i][j] == BB_SOBOL) {
                sobolCount = sobolCount + 1L;
                factorTypeList[i][j] = sobolCount;
            } else if(dimTypeList[i][j] == BB_PSEUDO_RANDOM) {
                factorTypeList[i][j] = BB_PSEUDO_RANDOM;
            } else if(dimTypeList[i][j] == BB_GRID) {
                factorTypeList[i][j] = BB_GRID;
                if(gridFactor != -99 && gridFactor != i) {
                    throw SCException( __FILE__, __LINE__,
                                       "Can only apply grid factor to one dimension only");
                }
                gridFactor   = i;
                gridLocation = j;
            } else {
                throw SCException( __FILE__, __LINE__,
                                   "Invalid dimType.\n Can only be one of BB_GRID, BB_SOBOL, BB_PSEUDO_RANDOM");
            }
        } // j
    } // i
    rng = rngIn;
    //     setSeed(rngIn);
    nDimensions = factorDimension * numberDates;
    // Initialize return vector
    //  if(!vector)
    //      vector = new double [factorDimension];
    vector.assign(factorDimension, 0);
    if(!gaussianVect)
        gaussianVect = new double [factorDimension];
    if(!gaussianVect)
    {
        throw SCException( __FILE__, __LINE__,
                           "Could not allocate memory for gaussianVect");
    }
    if(!listDates)
        listDates = new long [factorDimension];
    if(!listDates)
    {
        throw SCException( __FILE__, __LINE__,
                           "Could not allocate memory for listDates");
    }

    for(int iFactor = 0; iFactor < factorDimension; iFactor++)
    {
        listDates[iFactor]=(long)numberDates;
        //cout << "listDates[" <<iFactor <<"]=" << listDates[iFactor] << endl;
    } // iFactor

    if(!probabilities)
        probabilities = new double [factorDimension];
    if(!probabilities)
    {
        throw SCException( __FILE__, __LINE__,
                           "Could not allocate memory for probabilities");
    }

    // Initialize the Brownian bridge
#ifdef SCBB_DEBUG
    cout <<"factor dimension " << factorDimension << endl;
    for(int ab = 0; ab < factorDimension; ab++)
    {
        cout << "listDates[" << ab << "] = " << listDates[ab] << endl;
    }
    cout << "nPaths " << nPaths << endl;
    for(int b = 0; b < factorDimension; b++)
    {
        for(int ld =0 ; ld< listDates[b]; ld++) {
            cout << "factorTypeList[" << b << "][" << ld <<"] = " << factorTypeList[b][ld] << endl;
        }
    }
#endif //SCBB_DEBUG
    if(SGM_Initialize_Sobol_Bridge(
                (long)factorDimension,   // number of factors in a sim
                listDates,       //number of dimensions per factor (time steps)
                (long)nPaths,     // number of paths
                1,            // number of Runs
                (long)nPaths,     // number of paths to fetch
                factorTypeList,      // dimensionType array
                rng,        // seed
                FALSE,        // Sample tails flag
                0,            // number of points for sampling dist
                0,            // sampling interval bound
                1,            // max Prob ratio
                ptr2GaussianGrid,
                lowDiscSeqPtr)!=0)
    {
        throw SCException(__FILE__, __LINE__,
                          "Could not initialize Brownian bridge. SGM_Initialize_Sobol_Bridge failed");
    }

    if( SGM_Prepare_Sobol_Bridge() !=0 )
    {
        throw SCException(__FILE__, __LINE__,
                          "Could not prepare Brownian bridge. SGM_Prepare_Sobol_Bridge failed");
    }


} //BrownianBridgeSequence::initializedBrownianBridge
//***********
// Destructor
//***********
BrownianBridgeSequence::~BrownianBridgeSequence()
{
    // Free up memory
    if(listDates) {
        delete [] listDates;
        listDates = 0;
    }
    if(gaussianVect) {
        delete [] gaussianVect;
        gaussianVect = 0;
    }
    if(probabilities) {
        delete [] probabilities;
        probabilities = 0;
    }
    if(internalLowDiscSeqPtr.get()) {
        //                delete internalLowDiscSeqPtr;
        internalLowDiscSeqPtr = SequenceSP();
    }

    // Free Brownian bridge related data.
    SGM_Free_Sobol_Bridge();

    if(factorTypeList) {
        for(int i=0;i<factorDimension;i++) {
            delete [] factorTypeList[i];
            factorTypeList[i] = 0;
        }
        delete [] factorTypeList;
        factorTypeList = 0;
    }

} // BrownianBridgeSequence::~BrownianBridgeSequence
//***********************************************
// Default method to populate vector. Turned off.
//***********************************************
void BrownianBridgeSequence::populateVector()
{
    throw SCException( __FILE__, __LINE__,
                       "Must use populateVector(int thisPath) to access data on a per path basis");
}//populateVector()
//************************************************
// Method to populate vector. Pass in current path
// Needed for determining grid on final time.
//************************************************
void BrownianBridgeSequence::populateVector(int thisPath)
{
    // Layout of return vector factors by dates.
    // (factor,date)
    // (0,0),(1,0),(2,0) ... (0,1),(1,1),(2,1),...
    int count =0;

    for(long iDate = 0; iDate < numberDates; iDate++) {
        SGM_Fetch_Sobol_Bridge (SGM_FETCH_ONE_DIM_ALL_FACTORS,
                                0,
                                iDate,
                                (long)thisPath,
                                //                                         rng,
                                gaussianVect,
                                probabilities
                               );
        for(int iFactor=0; iFactor < factorDimension; iFactor++) {
            vector[count] = gaussianVect[iFactor];
            count++;
        } // iFactor
    } //iDate
    //cout << "GRID " << thisPath << " " << vector[numberDates-1] << endl;
} // populateVector(int thisPath)
//************************************************
// Work around to problem with BBridge sampling.
// Have to sample in iDate-iPath order. Inner loop
// on iPath and outer on iDate or else problems.
//************************************************
double * BrownianBridgeSequence::getDeviatesbyPathDate(int iPath, int iDate)
{
    //  double *retVect = new double [factorDimension];
    //  if(!retVect){
    //    throw SCException(__FILE__, __LINE__,
    //            "Could not allocate memory for retVect");
    //  }

    if(SGM_Fetch_Sobol_Bridge (SGM_FETCH_ONE_DIM_ALL_FACTORS,
                               0,
                               (long)iDate,
                               (long)iPath,
                               //                                    rng,
                               gaussianVect,
                               probabilities
                              )!=0) {
        throw SCException(__FILE__, __LINE__,
                          "Could not fetch Brownian bridge sequence. SGM_Fetch_Sobol_Bridge failed");
    }
    //  for(int iFactor=0; iFactor < factorDimension; iFactor++){
    //    retVect[iFactor] = gaussianVect[iFactor];
    //  } // iFactor

    //  return retVect;
    return gaussianVect;
}
//=============================================================
// Method to get deviates by date and factor. This is a possibly
// more efficient approach to fetching deviates. The innermost
// loop is over paths.
void BrownianBridgeSequence::getDeviatesbyFactorDate(int iFactor,
        int iDate,
        double *&returnVector)
{
    if(SGM_Fetch_Sobol_Bridge (SGM_FETCH_ONE_DIM_ONE_FACTOR,
                               (long)iFactor,
                               (long)iDate+1,
                               nPaths-1,
                               //                                    rng,
                               returnVector,
                               probabilities
                              )!=0) {
        throw SCException(__FILE__, __LINE__,
                          "Could not fetch Brownian bridge sequence. SGM_Fetch_Sobol_Bridge failed");
    }
}/*
void BrownianBridgeSequence::getDeviatesbyFactorDate(int iFactor,
                             int iDate,
                             double *&returnVector){
 */

//=============================================================
void BrownianBridgeSequence::getBrownianBridgeStructure(tnti_mat &theMat)
{
    int iFactor, iDate;
    if(theMat.size() ==0)
        theMat.assign(factorDimension, std::vector<int>(numberDates, 0));
    //    theMat.newsize(factorDimension,numberDates);
    if(!factorTypeList)
        throw SCException(__FILE__, __LINE__,
                          "factorTypeList not allocated. initialize Brownian bridge first.");
    for(iFactor = 0 ; iFactor < factorDimension; iFactor++) {
        for(iDate = 0 ; iDate < numberDates; iDate++) {
            theMat[iFactor][iDate] = factorTypeList[iFactor][iDate]; // simply copies matrix between two formats
        }//iDaet
    }//iFactor
}// long **BrownianBridgeSequence::getBrownianBridgeStructure()

CORE_END_NAMESPACE
