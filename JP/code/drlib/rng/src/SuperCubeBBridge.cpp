// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                        *
 * SuperCubeBBridge.cpp:                                                  *
 *      Implementation for the SuperCube-Brownian bridge class.           *
 * Author: M.Huq, Credit DR                                               *
 *                                                                        *
 *  Implementation of Supercube for SRM. The default here is that the     *
 *  Brownian bridge is contained within the SuperCube. There are two      *
 *  arrays that are passed in. The first is an importance array that      *
 *  contains the order of importance of each factor. The second is an     *
 *  array that contains a flag telling SupercubeBBridge whether each      *
 *  factor has path-depedence or not. These two arrays determine how to   *
 *  place the Brownian bridge supports.                                   *
 *                                                                        *
 -----------------------------------------------------------------------*/

//#include "edginc/VectorMatrix.h"
#include "edginc/SuperCubeBBridge.h"
#include "edginc/SCException.h"
#include "edginc/BrownianBridgeSequence.h"
#include "edginc/PseudoRandomSequence.h"
#include "edginc/SobolSequence.h"
#include "Version.h"
//#include "edginc/SC_stl.h"
#include "VectorMatrix.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

CORE_BEGIN_NAMESPACE

using namespace std;

double  SC_InvertCumNormal(double x);
int     SGM_Prepare_Sobol_Bridge();


//=======================================================
// Default constructor.
SuperCubeBBridge::SuperCubeBBridge() : SuperCube()
{
    numberAllFactors = 0;
    numberBridgeSupports = 0;
    pathDimTableSize = 0;
    pathDimMin = 0;
    pathDimDelta = 0;
    pathDimTable = 0;
    internalGaussianGrid = 0;
    pathDateIndexMap = 0;
    LDS_multiplier = 1;

    scbbPathDepList = 0;    // Vector containing path dependence information (0=false, 1=true)
    origPathDepList = 0;    // Vector to contain index mapping for factors
    importance = 0;        // Vector to contain original importance array
    invImportance = 0;        // Vector to contain inverse index mapping
    pathDepMapping =0; // Array to store path-dependence for SC
    theRandomGrid = RandomGridSP();
    whichRun = 0;
    numberRuns = 1;
}

//========================================================
// SRM-version
// Interface assumes brownian bridge type structure laid out
// in factors and dates.
// Look at SuperCube.h and Sequence.h for information on
// inherited variables.
SuperCubeBBridge::SuperCubeBBridge (
    INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permutations of partitions
    int numberPaths,         // Number of paths to extract
    int numberFactors,       // Number of factors
    int numberDates,         // Number of dates to consider.
    int *theImportance,         // Array containing importance ordering of each factor
    int *pathDependence) : SuperCube()     // Array indicating path-dependence (1) or not (0).
{
    numberAllFactors = 0;
    numberBridgeSupports = 0;
    pathDimTableSize = 0;
    pathDimMin = 0;
    pathDimDelta = 0;
    pathDimTable = 0;
    internalGaussianGrid = 0;
    pathDateIndexMap = 0;
    LDS_multiplier = 1;

    scbbPathDepList = 0;    // Vector containing path dependence information (0=false, 1=true)
    origPathDepList = 0;    // Vector to contain index mapping for factors
    importance = 0;        // Vector to contain original importance array
    invImportance = 0;        // Vector to contain inverse index mapping
    pathDepMapping =0; // Array to store path-dependence for SC
    theRandomGrid = RandomGridSP();
    whichRun = 0;
    numberRuns = 1;
    initializeObject ( rng,
                       numberPaths,
                       numberFactors,
                       numberDates,
                       theImportance,
                       pathDependence);

} // SuperCubeBBridge::SuperCubeBBridge

//===========================================================
// Method for initializing SuperCubeBBridge object
// Important inputs:
//   int *importance : array of importance rankings per factor
//   int *pathDependence : array of pathDependence per factor
//   These two arrays are fundamental in setting up the structure
//   of the underlying SuperCube.
//
// * importance array determines which factors are more important
//   in terms of laying down supports for the Brownian bridge.
// * path depdence array allows values of -1, 0 and 1
//   -1 => no bridge supports. Then this factor is shoved outside of the partitions
//    0 => support only at the final time.
//    1 => supports laid out in time dependent on LDS availability
//   SuperCube set up as follows:
//   1. Count up the number of factors that have path dependence not -1
//      This forms the number of factors on which we wish to create a
//      SuperCube out of. Factors with path dep -1 are set aside to be
//      sampled from a pseudo-random number generator.
//   2. Construct now the SuperCube structure on the factors that remain.
//   3. Initialize the Supercube
//   4. When sampling sample factors with path dep -1 from pseudo-random
//      and others from the SuperCube according to importance and path
//      dependence.
void SuperCubeBBridge::initializeObject (
    INormalSuperCubeRNGSP _rng,              // Seed for pseudo-random generator and random-permuations of partition
    int numberPaths,         // Number of paths to extract
    int numberFactors,       // Number of factors
    int numberDates,         // Number of dates to consider.
    int *theImportance,      // Array containing importance ordering of each factor
    int *pathDependence)     // Array indicating path-dependence (1) or not (0).
{
    internalGaussianGrid = 0;
    initializeObject (_rng,
                      numberPaths,
                      numberFactors,
                      numberDates,
                      theImportance,
                      pathDependence,
                      internalGaussianGrid);    // grid determined internally. Set to null.
}

//========================================================
// SRM-version II
// The following version takes as an argument a parameter
// that allows setting up an uniform grid in probability space
// to be used by the Brownian bridge (if there is to be one).
// The parameter is an offset from 0 or 1 probability. If epsilon is
// chosen then we set up a grid of N points that is centered in values from
// epsilon to 1-epsilon.
SuperCubeBBridge::SuperCubeBBridge (
    INormalSuperCubeRNGSP rng,         // Seed for pseudo-random generator and random-permuations of partitions
    int numberPaths,    // Number of paths to extract
    int numberFactors,  // Number of factors
    int numberDates,    // Number of dates to consider.
    int *theImportance,    // Array containing importance ordering of each factor
    int *pathDependence,// Array indicating path-dependence (1) or not (0).
    double epsilon) : SuperCube()    // Parameter for setting up grid.
{
    numberAllFactors = 0;
    numberBridgeSupports = 0;
    pathDimTableSize = 0;
    pathDimMin = 0;
    pathDimDelta = 0;
    pathDimTable = 0;
    internalGaussianGrid = 0;
    pathDateIndexMap = 0;
    LDS_multiplier = 1;

    scbbPathDepList = 0;    // Vector containing path dependence information (0=false, 1=true)
    origPathDepList = 0;    // Vector to contain index mapping for factors
    importance = 0;        // Vector to contain original importance array
    invImportance = 0;        // Vector to contain inverse index mapping
    pathDepMapping =0; // Array to store path-dependence for SC
    theRandomGrid = RandomGridSP();
    whichRun = 0;
    numberRuns = 1;
    initializeObject (rng,
                      numberPaths,
                      numberFactors,
                      numberDates,
                      theImportance,
                      pathDependence,
                      epsilon);
} //SuperCubeBBridge::SuperCubeBBridge


// This function is called by the above constructor. It is provided to get around issues
// of scope. The recommended approach to initializing a SuperCubeBBridge is first declare
// the object with the default constructor and then use initializeObject.
void SuperCubeBBridge::initializeObject (
    INormalSuperCubeRNGSP rng,       // Seed for pseudo-random generator and random-permuations of partitions
    int numberPaths,    // Number of paths to extract
    int numberFactors,  // Number of factors
    int numberDates,    // Number of dates to consider.
    int *theImportance,    // Array containing importance ordering of each factor
    int *pathDependence,// Array indicating path-dependence (1) or not (0).
    double epsilon)     // Parameter for setting up grid.
{
    int iPath;
    if(numberPaths == 0)
        throw SCException (__FILE__, __LINE__,
                           "number of paths found to be zero! Check constructor.");
    // Set up the grid.
    double min = epsilon;
    double max = 1.0 - epsilon;
    // Center the grid about the min and max so that the limit as epsilon-> 0 leads the min, max->0,1
    // uniformly with the grid properly centered.
    double step_size = (max-min)/(double)numberPaths;
    // Allocation memory for the grid
    internalGaussianGrid = new double [numberPaths];
    if(!internalGaussianGrid)
        throw SCException (__FILE__, __LINE__,
                           "Could not allocate memory for internalGaussianGrid!");
    // Compute values now.
    for(iPath = 0 ; iPath < numberPaths; iPath++)
    {
        internalGaussianGrid[iPath] = min + step_size *0.5 + (double)iPath * step_size;
    }//iPath
    for(iPath = 0 ; iPath < numberPaths; iPath++)
    {
        internalGaussianGrid[iPath] = SC_InvertCumNormal(internalGaussianGrid[iPath]);
    }//iPath

    // initialize the supercube brownian bridge with this grid now.
    initializeObject (rng,
                      numberPaths,
                      numberFactors,
                      numberDates,
                      theImportance,
                      pathDependence,
                      internalGaussianGrid);
}// void SuperCubeBBridge::initializeObject
//========================================================
// SRM-version III
// Multi-run version that uses a randomly sampled grid
// at the final time for each run. For M runs of a N-path
// integration the grid size will be MxN long. All M runs
// will fully sample the entire grid. This is a useful way
// to get higher resolution with smaller memory impact.
SuperCubeBBridge::SuperCubeBBridge (
    INormalSuperCubeRNGSP rng, // Seed for pseudo-random generator and random-permuations of partitions
    int numberPaths,         // Number of paths to extract
    int numberFactors,       // Number of factors
    int numberDates,         // Number of dates to consider.
    int *theImportance,         // Array containing importance order ing of each factor
    int *pathDependence,     // Array indicating path-dependence (1) or not (0).
    int randomGridFlag) : SuperCube()
{     // Flag to toggle random grid. If > 1 t he number of runs
    // assumed for random sampling the grid is
    //= randomGridFlag. Upon each reinitial ization then
    // a new grid of size numberPaths is sa mpled.
    numberAllFactors = 0;
    numberBridgeSupports = 0;
    pathDimTableSize = 0;
    pathDimMin = 0;
    pathDimDelta = 0;
    pathDimTable = 0;
    internalGaussianGrid = 0;
    pathDateIndexMap = 0;
    LDS_multiplier = 1;

    scbbPathDepList = 0;    // Vector containing path dependence information (0=false, 1=true)
    origPathDepList = 0;    // Vector to contain index mapping for factors
    importance = 0;        // Vector to contain original importance array
    invImportance = 0;        // Vector to contain inverse index mapping
    pathDepMapping =0; // Array to store path-dependence for SC
    theRandomGrid = RandomGridSP();
    whichRun = 0;
    numberRuns = 1;
    initializeObject (rng,
                      numberPaths,
                      numberFactors,
                      numberDates,
                      theImportance,
                      pathDependence,
                      randomGridFlag);
} // SuperCube::SuperCube(...)

void SuperCubeBBridge::initializeObject (
    INormalSuperCubeRNGSP rng,
    int numberPaths,
    int numberFactors,
    int numberDates,
    int *theImportance,
    int *pathDependence,
    int randomGridFlag)
{
    if(randomGridFlag < 1) {
        // Set up a static grid
        theRandomGrid = RandomGridSP(); // Use this later to detect this condition
        internalGaussianGrid = 0;
        initializeObject (rng,
                          numberPaths,
                          numberFactors,
                          numberDates,
                          theImportance,
                          pathDependence,
                          internalGaussianGrid);
    } else {
        // initialize the random grid
        whichRun  =0;
        numberRuns = randomGridFlag;
        theRandomGrid = RandomGridSP(new RandomGrid);
        theRandomGrid->initializeRandomGrid(this->getUniformRNG(),
                                            numberPaths,
                                            randomGridFlag);
        theRandomGrid->populateVector();
        internalGaussianGrid = theRandomGrid->getVector();
        initializeObject (rng, //theRandomGrid->getSeed(),
                          numberPaths,
                          numberFactors,
                          numberDates,
                          theImportance,
                          pathDependence,
                          internalGaussianGrid);
    }
}/*
void SuperCubeBBridge::initializeObject (
                INormalSuperCubeRNGSP rng,
                int numberPaths,
                int numberFactors,
                int numberDates,
                int *theImportance,
                int *pathDependence,
                int randomGridFlag);
*/

//========================================================
// Internal version that takes a gaussian grid as an argument.
// This is used by the first two constructors above.
SuperCubeBBridge::SuperCubeBBridge (
    INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
    int numberPaths,         // Number of paths to extract
    int numberFactors,       // Number of factors
    int numberDates,         // Number of dates to consider.
    int *theImportance,         // Array containing importance ordering of each factor
    int *pathDependence,     // Array indicating path-dependence (1) or not (0).
    double *GaussianGridIn)  : SuperCube()// Gaussian grid.
{
    numberAllFactors = 0;
    numberBridgeSupports = 0;
    pathDimTableSize = 0;
    pathDimMin = 0;
    pathDimDelta = 0;
    pathDimTable = 0;
    internalGaussianGrid = 0;
    pathDateIndexMap = 0;
    LDS_multiplier = 1;

    scbbPathDepList = 0;    // Vector containing path dependence information (0=false, 1=true)
    origPathDepList = 0;    // Vector to contain index mapping for factors
    importance = 0;        // Vector to contain original importance array
    invImportance = 0;        // Vector to contain inverse index mapping
    pathDepMapping =0; // Array to store path-dependence for SC
    theRandomGrid = RandomGridSP();
    whichRun = 0;
    numberRuns = 1;
    initializeObject (
        rng,
        numberPaths,
        numberFactors,
        numberDates,
        theImportance,
        pathDependence,
        GaussianGridIn);
}
void SuperCubeBBridge::initializeObject (
    INormalSuperCubeRNGSP _gen,              // Seed for pseudo-random generator and random-permuations of partitions
    int numberPaths,         // Number of paths to extract
    int numberFactors,       // Number of factors
    int numberDates,         // Number of dates to consider.
    int *theImportance,         // Array containing importance ordering of each factor
    int *pathDependence,     // Array indicating path-dependence (1) or not (0).
    double *GaussianGridIn)  // Gaussian grid
{
    int iFactor;
    int iDates;
    int iBridge;

    // This version gets grid passed in
    BBgaussianGrid = GaussianGridIn;

    //==================INPUT CHECKING===============================
    /*        if(aSeed > 0)
                    aSeed *= -1;*/
    //         assert(_gen->getSeed() == -1);
    _gen->superCubeAdjustment();
    resetRNGGen( _gen);

    if(numberFactors < 0 || numberDates < 0)
        throw SCException(__FILE__, __LINE__,
                          "Invalid value for number of dates or factors : < 0");
    if(numberPaths < 0)
        throw SCException(__FILE__, __LINE__,
                          "Invalid value for number of paths : < 0");
    if(numberFactors < 0)
        throw SCException(__FILE__, __LINE__,
                          "Invalid value for total number of factors : < 0");
    if(numberDates < 0)
        throw SCException(__FILE__, __LINE__,
                          "Invalid value for total number of Dates : < 0");
    //Check to ensure that importance and pathDependence arrays are defined
    if(!theImportance || !pathDependence)
    {
        throw SCException(__FILE__, __LINE__,
                          "theImportance or pathDependence array not defined.");
    }
    // Check to ensure that importance array contains all factor numbers
    // and none appear more than once
    int curIndex = numberFactors-1;
    while(curIndex >= 0)
    {
        if(theImportance[curIndex] < 0 ||
                theImportance[curIndex] >= numberFactors)
            throw SCException(__FILE__, __LINE__,
                              "Error in importance array : Values must be within the range 0 to number_factors-1");
        int count = 0;
        for(int i = 0; i < numberFactors; i++) {
            if(theImportance[i] == curIndex) {
                count++;
            }
        }
        if(count != 1) {
            throw SCException (__FILE__, __LINE__,
                               "Error in importance array : Check to ensure that each value appears once only.");
        } //count
        curIndex--;
    }// curIndex
    for(iFactor = 0; iFactor < numberFactors; iFactor++)
    {
        if(theImportance[iFactor] < 0 || theImportance[iFactor] > numberFactors-1) {
            throw SCException (__FILE__, __LINE__,
                               "Error in importance array : Values must be in range 0 to number Factors -1");
        }
    }//iFactor

    // Check path-dependence array has values -1, 0 and 1 only
    for(curIndex =0; curIndex < numberFactors; curIndex++)
    {
        if(pathDependence[curIndex] <-1 || pathDependence[curIndex] > 1) {
            throw SCException (__FILE__, __LINE__,
                               "Error in pathDependence array : Can have value >= -1 and <=1 only.");
        }
    }//curIndex
    //==================END INPUT CHECKING===========================




    //==================INITIALIZE SUPERCUBEBBRIDGE AND PARENT CLASS PARAMETERS=======================
    // Basic parameters unaffected by importance and path-dependence
    numPaths = numberPaths;
    numberAllFactors = numberFactors;
    numBBdates = numberDates;
    //        seed = aSeed;
    nDimensions = numberAllFactors * numBBdates; // dimensionality of entire SC
    /*  if(!vector)
          vector = new double [nDimensions];*/
    vector.assign(nDimensions, 0.0);
    /*        initType = DEFAULTTYPE;
            filePointer = 0;*/
    //  fileIndex = 0;
    //    nSets = 0;
    // Initialize Sobol dimension - path number table
    initializePathDimTable();

    // Set SuperCube uniformity flag
    uniformityFlag = GET_BB_STYLE;

    //=========================================================================
    // Set up importance and path-dependence arrays
    //=================NOTE============================================
    // The importance array provides the mapping from
    // ****SuperCube arrays which are importance-rank-ordered***
    // This provides the book-keeping necesary for reallocation of
    // factors when Sampling the SuperCubeBBridge.
    // The following goes through the path-dependence parts and assigns
    // path-depedence flags to the primary partition. If there are no
    // path dependent deviates then add one at least.
    //
    // Primary partition thus contains deviates in the following order:
    //   0 1 2 ... n-3 n-2 n-1    (importance-ranking)
    //   0 0 0 ...  0   1   1     (path-dependence : 2 out of n factors)
    // When sampling a partition (secondary or primary)
    //
    // Note the importance array takes as its argument an external index
    // and returns the SC index.
    // invImportance takes as its argument a SC index and returns an external
    // index.
    //================================================================='
    // Data for path-dependence in two different orderings
    origPathDepList = new int[numberAllFactors];

    if(deviatesByPath == 0)
    {
        deviatesByPath = new double [numPaths];
    }

    //=========================================================================
    // copy over importance array.
    importance = new int [numberAllFactors];
    for(iFactor=0; iFactor<numberAllFactors; iFactor++)
    {
        importance[iFactor] = theImportance[iFactor];
    }//iFactor
    //==================SET UP INVERSE MAPPING=======================
    invImportance = new int [numberAllFactors];
    for(iFactor=0;iFactor<numberAllFactors;iFactor++)
    {
        invImportance[importance[iFactor]] = iFactor;
    }//iFactor

    for(iFactor=0;iFactor<numberAllFactors;iFactor++)
    {
        if(invImportance[iFactor] < 0 || invImportance[iFactor] > numberAllFactors-1)
            throw SCException(__FILE__,__LINE__,
                              "Inverse index map out of range - check importance values");
    }
    //=========================================================================
    // Compute the number of path-dependent factors.
    // We will construct the underlying SuperCube only for factors which are
    // not pseudo-random.
    // Thus if number of factors with path-dep of 0 or 1 is 5 and the number
    // of pseudo-random factors is 2 - then the Supercube size is 5 and overall
    // size is 7.
    // The number of path-dependent factors together with the number of LDS
    // dimensions available will determine the number of partions, supports etc.
    //======================
    // Count up path dependent factors
    //  int numberAllBridgeFactors=0;
    for(iFactor = 0; iFactor<numberAllFactors; iFactor++)
    {
        origPathDepList[iFactor] = pathDependence[iFactor];
        //    if(pathDependence[iFactor]>=0 )
        //      numberAllBridgeFactors++;
    }//iFactor
    // SuperCube size now in numberAllBridgeFactors

    //========================================================================
    // Now determine the partition layout
    // Variables to store the number of available Low-discrepancy variables
    // at the final time and for the bridge supports.
    int numberLDSBridge;
    int numberLDSFinal;
    // Given the number of paths figure out the number of allowable
    // low-discrepancy dimensions.
    // The following sets:
    // nPartition
    // partitionSize
    // MaxDimension
    // numBBfactors  (number of bridge factors in primary partition)
    //getPartitionData(numberLDSBridge, numberLDSFinal);
    initializeIndexMaps(numberLDSBridge, numberLDSFinal);
    // Set the SuperCube base sequence dimension - size of primary
    sequenceDimension = partitionSize;

    //========================================================
    // Memory allocation for Supercube primary partition
    /*  primary = new double * [partitionSize];
        for(int iPart = 0; iPart < partitionSize; iPart++){
          primary[iPart] = new double [numPaths];
        }//iPart*/
    // Allocate memory for the primary partition.
    primary = new double * [numPaths];
    for(int i = 0; i < numPaths; i++)
    {
        primary[i] = new double [partitionSize];
    }


    //======================
    // Set up Brownian bridge in primary partition.
    //======================
    // First allocate memory for arrays needed in Brownian bridge
    dimTable = new long * [numBBFactors];
    for(iFactor=0;iFactor<numBBFactors; iFactor++)
    {
        dimTable[iFactor] = new long [numBBdates];
        if(!dimTable[iFactor]) {
            throw SCException(__FILE__, __LINE__,
                              "Could not allocate memory for dimTable");
        }
    } // iFactor
    // Load up pseudo-random for all by default
    for(iFactor = 0 ; iFactor < numBBFactors; iFactor++)
    {
        for(iDates = 0; iDates < numBBdates; iDates++)
            dimTable[iFactor][iDates] = BB_PSEUDO_RANDOM;
    } // iFactor
    //==========================================================
    // Count up path dependent factors within primary partition.
    // Map to Supercube indexing scheme.
    // That is, reorder according to importance.
    //  for(iFactor = 0; iFactor < numberAllFactors; iFactor++){
    //    scbbPathDepList[iFactor] = pathDependence[invImportance[iFactor]];
    //}//iFactor
    //==========================================================
    // Now count up primary partition path dependence
    int numPrimaryPathDep = 0;
    for(iFactor=0; iFactor<numBBFactors; iFactor++)
    {
        if(scbbPathDepList[iFactor] >= 0)
            numPrimaryPathDep++;
    }//iFactor
    //  int minPathDeps;
    // #define TURNOFF20PERCENT
    // #ifndef TURNOFF20PERCENT
    //  if(numPrimaryPathDep == 0 && numberAllBridgeFactors != 0 && nPartitions >1){
    //    //Make sure that at least 20% of factors are path-dependent.
    //    // This way if no path-dep factors in primary but do exist in
    //    // secondary then they are covered.
    //    minPathDeps = MAX(1, (int)ceil(0.2 * numBBFactors));
    //    for(int iprim =0; iprim < minPathDeps; iprim++){
    //      scbbPathDepList[iprim] = 1; // Set from most important
    //    }//iprim
    //  }else{
    // #endif //TURNOFF20PERCENT
    //    minPathDeps = numPrimaryPathDep;
    //  // At this point minPathDeps contains the number of path depedent vars
    // #ifndef TURNOFF20PERCENT
    //  }
    // #endif
    //=========================================================================
    // Compute the support location dates in bisection order
    // Use the numBBdates as the maximum numbers of supports.
    // This merely contains the list of possible support dates in increasing
    // order.
    //=========================================================================
    bridgeLocations = new int [numBBdates];  // Allocate max memory
    if(!bridgeLocations)
    {
        throw SCException(__FILE__, __LINE__,
                          "Could not allocate memory for bridgeLocations");
    }

    bridgeLocations[0] = numBBdates-1; // T point
    numberBridgeSupports = 1; // T point is already secured.
    double *ratios = new double [numBBdates+1];
    for(int iratio = 0 ; iratio < numBBdates+1; iratio++)
        ratios[iratio] = 0;
    int level =0;
    while(numberBridgeSupports < numBBdates )
    {
        level++;
        for(int m = (int)pow(2.0,(double)(level-1)); m>0; m--) {
            ratios[numberBridgeSupports-1]=((2.0*m -1.0)/pow(2.0,(double)level));
            numberBridgeSupports++;
            if(numberBridgeSupports > numBBdates-1)
                break;
        }//m
    }//numberBridgeSupports


    for(iBridge=1; iBridge < numberBridgeSupports; iBridge++)
    {
        bridgeLocations[iBridge] =(int)(ratios[iBridge-1] * (numBBdates-1));
    }//iBridge
    delete [] ratios;


    //================================================================
    // At this point basic locations of the bridges are secured.
    // Have to figure out which bridges to actually place based upon
    // which factors are path dependent and how many LDS are available.
    // For each factor  - starting from most important lay down supports
    //================================================================
    // At this point we are in SuperCubeBBridge internal representation
    int LDScount = 0;   // LDS counter initialized.
    int gridCount =0;
    int iSupport = 0;
    //
    // On the first instance of a support where we have not placed a grid place it.
    for(iDates=numBBdates-1; iDates>=0; iDates--)
    {
        int supportFound=0;
        for(iFactor = 0 ; iFactor < numBBFactors; iFactor++) {
            if(scbbPathDepList[iFactor] >= 0 && iDates == numBBdates-1 && LDScount < numberLDSFinal+numberLDSBridge) {
                if(gridCount == 0) {
#ifndef SCBB_NO_GRID
                    dimTable[iFactor][iDates] = BB_GRID; // Set up the grid
                    supportFound++;
                    gridCount++;
#else

                    dimTable[iFactor][iDates] = BB_SOBOL;
                    supportFound++;
                    LDScount++;
#endif // SCBB_NOGRID

                } else {
                    dimTable[iFactor][iDates] = BB_SOBOL;
                    supportFound++;
                    LDScount++;
                } //over gridCount
            } else if(scbbPathDepList[iFactor] ==1 &&
                      LDScount<numberLDSFinal+numberLDSBridge &&
                      iSupport < numberBridgeSupports) {
                if(bridgeLocations[iSupport] < 0 || bridgeLocations[iSupport] > numBBdates-1)
                    cout << "Warning : " << iSupport << " " << bridgeLocations[iSupport]<< endl;
                dimTable[iFactor][bridgeLocations[iSupport]] = BB_SOBOL;
                LDScount++;
                supportFound++;
                if(LDScount >= numberLDSFinal+numberLDSBridge || iSupport >= numberBridgeSupports)
                    break;
            }// over if
        }//iFactor
        if(supportFound != 0)
            iSupport++;
    }//iDates
    ldsBBcount = std::max(LDScount,1);
    //cout << "LDS count = " << ldsBBcount << endl;


#ifdef SCBB_DEBUG

    cout << "============SuperCubeBBridge dump============" << endl;
    cout << "Number of partitions " << nPartitions << endl;
    cout << "Partition size       " << partitionSize << endl;
    cout << "Number of BB factors " << numBBFactors<< endl;
    cout << "Number of BB supports "<< numberBridgeSupports << endl;
    cout << "=============================================" << endl;
#endif // SCBB_DEBUG

    //======================
    // Initialize the BrownianBridgeSequence object
    // Create a brownian bridge sequence
    // BrownianBridgeSequence *bSeq ;
    // PseudoRandomSequence *pSeq =

    if(numBBFactors == 0 || nPartitions == 0)
    {
        // Set up dummy pseudo-random getter even though not used
        baseSequence = SequenceSP(new PseudoRandomSequence(this->getNormalRNG(), numberAllFactors));
        //    seed = baseSequence->getSeed();
    } else
    {
        BBLDSequence = SobolSequenceSP(new SobolSequence(-1, ldsBBcount));
        baseSequence = BrownianBridgeSequenceSP(new BrownianBridgeSequence(numBBFactors,
                                                numBBdates,
                                                dimTable,
                                                numPaths,
                                                this->getNormalRNG(),
                                                BBgaussianGrid,
                                                BBLDSequence));
        //                seed = bSeq->getSeed();
        //======================
        // Assign the primary partition sequence pointer.
        // baseSequence = bSeq;
        // delete pSeq;
    }

    // Set up internal base pointer for memory management. See SuperCube
    // destructor.
    internalBasePtr = baseSequence;

#ifdef SCBB_DEBUG

    cout << "==============" << endl;
    //        cout << seed << endl;
    cout << "Number paths " << numberPaths << endl;
    cout << "MaxDimensions" << MaxDimensions << endl;
    cout << "Partition Size "<< partitionSize << endl;
    cout << "nPartitions " << nPartitions << endl;
    cout << "Sequencedimensioon = " << sequenceDimension << endl;
    cout << "nDimensions = " << nDimensions << endl;
    cout << "numBBdates = " << numBBdates << endl;
    cout << "numBBFactors = " << numBBFactors << endl;
    cout << "==============" << endl;
#endif // SCBB_DEBUG
    //======================
    // Initialize index table
    if(nPartitions !=0)
    {
        initializeIndexTable(this->getUniformRNG()); /* (1) */

        // Set up partition 0 ; shuffling done during population of vectors.
        populatePrimaryPartition(); /* (2) */
        /* The problem in the above code is that both functions will use the same sequence of RNG!. Indeed, the BrownianBridgeSequence is created with seed -1; which is then stored in context.seed, so the sequence will be reset in populatePrimaryPartition(), and the same will happen in initializeIndexTable().
        */

    }


    // Clean up
    delete [] bridgeLocations;

}// SuperCubeBBridge::initializeObject


SuperCubeBBridge::~SuperCubeBBridge()
{
    if(dimTable) {
        for(int iFactor=0;iFactor<numBBFactors; iFactor++) {
            delete [] dimTable[iFactor];
            dimTable[iFactor] = 0;
        }
        delete [] dimTable;
        dimTable = 0;
    }
#if 0
    if(theRandomGrid) {
        delete theRandomGrid; // uses SP now
        theRandomGrid=0;
        // Also dereference internalGaussianGrid from the randomgrid
        internalGaussianGrid = 0;
    }
#endif
    if(internalGaussianGrid) {
        delete [] internalGaussianGrid;
        internalGaussianGrid = 0;
    }
    if(pathDimTable) {
        delete [] pathDimTable;
        pathDimTable = 0;
    }
    if(pathDateIndexMap) {
        delete [] pathDateIndexMap;
        pathDateIndexMap = 0;
    }

    if(scbbPathDepList) {
        delete [] scbbPathDepList;
        scbbPathDepList = 0;
    }
    if(origPathDepList) {
        delete [] origPathDepList;
        origPathDepList = 0;
    }
    if(importance) {
        delete [] importance;
        importance = 0;
    }
    if(invImportance) {
        delete [] invImportance;
        invImportance = 0;
    }
    if(pathDepMapping) {
        delete [] pathDepMapping;
        pathDepMapping =0;
    }
    if(deviatesByPath) {
        delete [] deviatesByPath;
        deviatesByPath = 0;
    }
    if(impPathDepList) {
        delete [] impPathDepList;
        impPathDepList = 0;
    }
}


//=======================================================================
// Override SuperCube method for populating vector. This takes into
// account indexing for path-dependent factors.
// At the end of this vector must contain factors by dates
// Wrapper to preserve the SuperCube structure.
void SuperCubeBBridge::populateVector(int iPath)
{
    getDeviateByPath(iPath, getVector());
}//SuperCubeBBridge::populateVector


//==========================================================================
// SuperCubeBBridge::getPartitionData()
// Method to determine distribution of low-discrepancy dimensions amongst
// factors. Proceed as follows:
// 1. Look up number of allowable sobol and grid dimensions. This from a table.
// 2. Given a maximum number of low-discrepancies to work with divide the number
//    into two parts. Those that will be applied at the final time and those
//    applied to the rest of the bridge for path-dependent factors. If there are
//    no path dependent factors then apply all to the final time - as many as needed.
// This function sets
//   - size of primary partition
//   - number of partitions
//   - number of bridge locations
void SuperCubeBBridge::getPartitionData(int &numberLDSBridge,
                                        int &numberLDSFinal)
{

    // Look up number of low-discrepancy sequences available
    int numberAvailLDS = getMaxLDSDim(numPaths);


    // Divide into two parts
    numberLDSBridge = (int)(numberAvailLDS / 2);
    numberLDSFinal = numberAvailLDS - numberLDSBridge;

    // If there are no path-dependent factors then do not use SC type structure
    if(numberAllBridgeFactors == 0) {
        numberLDSFinal = 0;
        numberLDSBridge = 0;
        nPartitions = 0;
        MaxDimensions = numberAllFactors * numBBdates;
        partitionSize = 0;
        numBBFactors = 0;
    } else {
        // Assign final time slice dimensions.
        if(numberAllBridgeFactors < numberLDSFinal) {
            // If the number of factors exceeds the number of available LDS dimensions
            // for the final timeslice then repartition the LDSs.
            numberLDSFinal = numberAllBridgeFactors;
            numberLDSBridge = numberAvailLDS - numberLDSFinal;
            // Have only one partition : set the partition data now.
            nPartitions = 1;
            MaxDimensions = numberAllBridgeFactors * numBBdates;
            partitionSize = MaxDimensions;
            numBBFactors = numberAllBridgeFactors;

        } else {
            // Otherwise, construct partitions based upon availability
            nPartitions = (int)ceil((double)numberAllBridgeFactors / (double)numberLDSFinal); // round up
            partitionSize = numberLDSFinal * numBBdates;
            MaxDimensions = partitionSize * nPartitions;
            numBBFactors = numberLDSFinal;   // Number of BB factors. This is the number
            // of factors for the Brownian bridge
            // constructed in the primary partition.
        }
    } // if over zero numberAllBridgeFactors
}//SuperCubeBBridge::getPartitionData()
//===================================================================================
// SuperCubeBBridge::initializePathDimTable()
// Method to initialize pathDimTable. The data in this table are calibrated externally
// and hard-coded here. This way they are private to the library.
//===================================================================================
void SuperCubeBBridge::initializePathDimTable()
{
    // 8192 16384 32768, 65536, 131072, 262144, 524288 1048576
    //   2    4     10     20     20      40     40      80
    //
    pathDimTableSize = 8;
    pathDimMin = 8192;
    pathDimDelta = 2;
    pathDimTable = new int [pathDimTableSize];

    // Values
    pathDimTable[0] = 4;   // 8192
    pathDimTable[1] = 4;   // 16384
    pathDimTable[2] = 4;   // 32768
    pathDimTable[3] = 4;   // 65536
    pathDimTable[4] = 8;   // 131072
    pathDimTable[5] = 8;   // 262144
    pathDimTable[6] = 16;   // 524288
    pathDimTable[7] = 32;   // 1048576

    for(int iRow = 0; iRow < pathDimTableSize; iRow++) {
        pathDimTable[iRow] *= LDS_multiplier;
        if(pathDimTable[iRow] > MAX_SOBOL) {
            pathDimTable[iRow] = MAX_SOBOL;
        }
    }//iRow
}//SuperCubeBBridge::initializePathDimTable

//===================================================================================
// SuperCubeBBridge::getMaxLDSDim(const int iPath)
// Method to return the maximum number of LDS dimensions for a given number
// of paths.
//===================================================================================
int SuperCubeBBridge::getMaxLDSDim(const int thisNumPath)
{
    // Sanity check
    if (thisNumPath < 0 ) {
        throw SCException(__FILE__, __LINE__,
                          "thisNumPath out of range in getMaxLDSDim(const int thisNumPath)");
    }
    int pathIndex = 0;
    int pIndex = 0;
    // FIXME: what this loop is doing? finds dimMin*2^i < thisNum <= dimMin*2^{i+1} ?
    while(pathIndex < pathDimTableSize) {
        int pIndexNext = pathDimMin*pow(2.0,pathIndex);
        if(thisNumPath > pIndex && thisNumPath <= pIndexNext) {
            break;
        } else {
            pIndex = pIndexNext;
            pathIndex++;
        }
    }//pIndex

    if (pathIndex < 0 ) {
        throw SCException(__FILE__, __LINE__,
                          "pathIndex out of range in getMaxLDSDim(const int thisNumPath)");
    }
    //  cout << "pathDimTableSize " << pathDimTableSize << endl;
    //  cout << "pathIndex " << pathIndex << endl;
    //  cout << "pathDimTable[pathIndex " << pathDimTable[pathIndex] << endl;
    return pathDimTable[pathIndex];
}//SuperCubeBBridge::getMaxLDSDim

//===================================================================================
// Method to return a block of deviates for a given factor and date index pair.
// This method is more efficient to use than say looping over paths and fetching
// factor-path or date-path pairs.
//===================================================================================
void SuperCubeBBridge::getDeviateByFactorDate(const long iFactor,
        const long iDate,
        double *returnArray)
{

    // Check that returnArray is allocated.
    if(!returnArray) {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByFactorDate : returnArray not allocated");
    }


    // If a pseudorandom factor then return gasdev
    if(origPathDepList[iFactor] == -1) {
        for(int iPath = 0; iPath < numPaths; iPath++) {
            returnArray[iPath] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
        }//iPath
    } else {
        // Find the SuperCube index this iFactor belongs in
        int scFactor = global2scFactorMap[iFactor];
        int scPartition = global2scPartitionMap[iFactor];

        int shuffle_index;

        // Figure out the parition index
        int partIndex = scFactor + iDate*numBBFactors;

        // Now populate the return array
        if(scPartition != 0) {
            for(int iPath = 0; iPath < numPaths; iPath++) {
                shuffle_index = indexTable[iPath][scPartition];
                returnArray[iPath] = primary[shuffle_index][partIndex];
            }//iPath
        } else {
            for(int i = 0; i < numPaths ; ++i)
                returnArray[i] = primary[i][partIndex];
            //memcpy(returnArray, primary[partIndex],sizeof(double)*numPaths);
        }
    }//if over path-dependence of factor

}/*
int SuperCubeBBridge::getDeviateByFactorDate(const long iFactor,
                         const long iDate,
                         double *&returnArray);
*/

//========================================================================
// SuperCubeBBridge::getDeviateByPathDate(long iPath, long iDate)
// Method to get deviates by path and date index. This is called
// by all methods that require deviates. Returns a vector of factors.
//========================================================================
void SuperCubeBBridge::getDeviateByPathDate(long iPath, long iDate,
        double *returnArray)
{
    //    int tmpIndex;
    int iFactor;

#ifndef SC_TURN_ON_CHECK

    if(!returnArray) {
        throw SCException(__FILE__, __LINE__,
                          "Memory not allocated for returnArray");
    }

    // Check that we have a valid date index
    if(iDate < 0 || iDate > numBBdates-1) {
        throw SCException(__FILE__, __LINE__,
                          "Date index passed in is out of range in getDeviateByPathDate");
    }
    // Check that we have a valid path index
    if(iPath < 0 || iPath > numPaths-1) {
        throw SCException(__FILE__, __LINE__,
                          "Path index passed in is out of range in getDeviateByPathDate");
    }
#endif // SC_TURN_ON_CHECK


    // Check to see if not all gasdev
    if(numberAllBridgeFactors != 0) {
        for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++) {
            // If a pseudorandom factor then return gasdev
            if(origPathDepList[iFactor] == -1) {
                returnArray[iFactor] = getNormalRNG()->fetch(); // SC_gasdev2(&seed);
            } else {
                // Find the SuperCube index this iFactor belongs in
                int scFactor = global2scFactorMap[iFactor];
                int scPartition = global2scPartitionMap[iFactor];

                int shuffle_index;

                // Figure out the parition index
                int partIndex = scFactor + iDate*numBBFactors;

                // Now populate the return array
                shuffle_index = indexTable[iPath][scPartition];
                returnArray[iFactor] = primary[shuffle_index][partIndex];
            } //if over origPathDepList
        } //iFactor
    } else {
        for(iFactor = 0; iFactor < numberAllFactors; iFactor++) {
            returnArray[iFactor] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
        }//iFactor
    }//if over numberAllBridgeFactors
}//SuperCubeBBridge::getDeviateByPathDate
//=======================================================================
// Method to obtain deviates given a path index. Returns an array of
// deviates by factors and dates
void SuperCubeBBridge::getDeviateByPath(long iPath, double *returnArray)
{
    int iDate;
    int iFactor;

    if(!returnArray) {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByPath - returnArray not allocated");
    }

    // Check that we have a valid date index
    if(iPath < 0 || iPath > numPaths-1) {
        throw SCException(__FILE__, __LINE__,
                          "Path index is out of range in getDeviateByPath");
    }

    //
    // Check to see if not all gasdev
    int indexByPath;
    if(numberAllBridgeFactors != 0) {
        for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++) {
            // If a pseudorandom factor then return gasdev
            if(origPathDepList[iFactor] == -1) {
                for(iDate = 0; iDate < numBBdates; iDate++) {
                    indexByPath = iFactor + iDate * numberAllFactors;
                    returnArray[indexByPath] = getNormalRNG()->fetch(); // SC_gasdev2(&seed);
                }//iDate
            } else {
                // Find the SuperCube index this iFactor belongs in
                int scFactor = global2scFactorMap[iFactor];
                int scPartition = global2scPartitionMap[iFactor];

                int shuffle_index;

                // Figure out the parition index
                for(iDate = 0; iDate < numBBdates; iDate++) {
                    indexByPath = iFactor + iDate * numberAllFactors;
                    int partIndex = scFactor + iDate*numBBFactors;

                    // Now populate the return array
                    shuffle_index = indexTable[iPath][scPartition];
                    returnArray[indexByPath] = primary[shuffle_index][partIndex];
                }//iDate
            }
        } //iFactor
    } else {
        for(iFactor = 0; iFactor < numberAllFactors; iFactor++) {
            for(iDate = 0; iDate < numBBdates; iDate++) {
                indexByPath = iFactor + iDate * numberAllFactors;
                returnArray[indexByPath] = getNormalRNG()->fetch(); // SC_gasdev2(&seed);
            }//iDate
        }//iFactor
    }//if over numberAllBridgeFactors
}//SuperCubeBBridge::getDeviateByPath
//=======================================================================
// Method to obtain deviates given a date index. Returns an array of
// deviates by factors and paths.
void SuperCubeBBridge::getDeviateByDate(long iDate,
                                        double *returnArray)
{
    int iPath;
    int iFactor;

    if(!returnArray) {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByDate - returnArray not allocated");
    }

    // Check that we have a valid date index
    if(iDate < 0 || iDate > numBBdates-1) {
        throw SCException(__FILE__, __LINE__,
                          "Date index passed in is out of range in getDeviateByDate");
    }

    //
    int indexByDate;
    if(numberAllBridgeFactors != 0) {
        for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++) {
            // If a pseudorandom factor then return gasdev
            if(origPathDepList[iFactor] == -1) {
                for(iPath = 0; iPath < numPaths; iPath++) {
                    indexByDate = iFactor + iPath * numberAllFactors;
                    returnArray[indexByDate] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
                }//iDate
            } else {
                // Find the SuperCube index this iFactor belongs in
                int scFactor = global2scFactorMap[iFactor];
                int scPartition = global2scPartitionMap[iFactor];

                int shuffle_index;

                // Figure out the parition index
                int partIndex = scFactor + iDate*numBBFactors;

                // Now populate the return array
                for(iPath = 0; iPath < numPaths; iPath++) {
                    indexByDate = iFactor + iPath * numberAllFactors;
                    shuffle_index = indexTable[iPath][scPartition];
                    returnArray[indexByDate] = primary[shuffle_index][partIndex];
                }//iPath
            }//if on origPathDepList
        } //iFactor
    } else {
        for(iFactor = 0; iFactor < numberAllFactors; iFactor++) {
            for(iPath=0; iPath < numPaths; iPath++) {
                indexByDate = iFactor + iPath * numberAllFactors;
                returnArray[indexByDate] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
            }//iPath
        }//iFactor
    }//if over numberAllBridgeFactors
}//SuperCubeBBridge::getDeviateByDate
//=======================================================================
// Method to obtain deviates given a date index. Returns an array of
// deviates by factors and paths.
void SuperCubeBBridge::getDeviateByDate(long iDate,
                                        double **returnArray)
{
    int iPath;
    int iFactor;

    if(!returnArray) {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByDate - returnArray not allocated");
    }

    // Check that we have a valid date index
    if(iDate < 0 || iDate > numBBdates-1) {
        throw SCException(__FILE__, __LINE__,
                          "Date index passed in is out of range in getDeviateByDate");
    }

    //
    if(numberAllBridgeFactors != 0) {
        for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++) {
            // If a pseudorandom factor then return gasdev
            if(origPathDepList[iFactor] == -1) {
                for(iPath = 0; iPath < numPaths; iPath++) {
                    returnArray[iFactor][iPath] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
                }//iDate
            } else {
                // Find the SuperCube index this iFactor belongs in
                int scFactor = global2scFactorMap[iFactor];
                int scPartition = global2scPartitionMap[iFactor];

                int shuffle_index;

                // Figure out the parition index
                int partIndex = scFactor + iDate*numBBFactors;

                // Now populate the return array
                for(iPath = 0; iPath < numPaths; iPath++) {
                    shuffle_index = indexTable[iPath][scPartition];
                    returnArray[iFactor][iPath] = primary[shuffle_index][partIndex];
                }//iPath
            }//if on origPathDepList
        } //iFactor
    } else {
        for(iFactor = 0; iFactor < numberAllFactors; iFactor++) {
            for(iPath=0; iPath < numPaths; iPath++) {
                returnArray[iFactor][iPath] = getNormalRNG()->fetch(); //SC_gasdev2(&seed);
            }//iPath
        }//iFactor
    }//if over numberAllBridgeFactors
}//SuperCubeBBridge::getDeviateByDate(2darray version)
//=============================================
// Method to output SuperCubeBBridge structure
void SuperCubeBBridge::DumpSCStructure(const char *fileName)
{
    int iFactor,iDate;

    fstream fileStream;
    ostream *theStream;
    if(!fileName) {
        theStream = &cerr;
    } else {
        fileStream.open(fileName,ios::out);
        theStream = &fileStream;
        if(!theStream)
            theStream = &cerr;
    }
    // First the Basic properties
    *theStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    *theStream << "SuperCube version " << getSCVersion() << endl;
    *theStream << "# of paths : " << numPaths << endl;
    *theStream << "# of overall factors " << numberAllFactors << endl;
    *theStream << "# of factors in SuperCube structure:  "
    << numberAllBridgeFactors << endl;
    *theStream << "# of factors sampled with pure pseudo-random : ie with no path dep "
    << numberAllFactors - numberAllBridgeFactors << endl;
    *theStream << "# of Partitions in SuperCube " << nPartitions << endl;
    *theStream << "# of factors per partition " << numBBFactors << endl;
    *theStream << "# of dates overall (no supercubing in dates)" << numBBdates << endl;
    *theStream << "# of Available sobol dimensions "
    << getMaxLDSDim(numPaths) << endl;
    //   *theStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    //   *theStream << "Path-Dimension table" << endl;
    //   (*theStream).setf(ios::right);
    //   *theStream << setw(10)<< "Path#: ";
    //   int thedim = pathDimMin/pathDimDelta;
    //   int iDim;
    //   for(iDim = 0; iDim < pathDimTableSize; iDim++){
    //     thedim *= pathDimDelta;
    //     *theStream << setw(10) << thedim;
    //   }//iDim
    //   *theStream << endl;
    //   (*theStream).setf(ios::right);
    //   *theStream << setw(10)<< "Dim: ";
    //   for(iDim = 0; iDim < pathDimTableSize; iDim++)
    //     *theStream << setw(10) << pathDimTable[iDim];
    //   *theStream << endl;



    *theStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    *theStream << "Importance (I), Path dpendence (P) dump by factor" << endl;
    (*theStream).setf(ios::right);
    *theStream << "F: ";
    for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++)
        *theStream << setw(6) << iFactor;
    *theStream << endl;
    *theStream << "I: ";
    for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++)
        *theStream << setw(6) << importance[iFactor];
    *theStream << endl;
    *theStream << "P: ";
    for(iFactor = 0 ; iFactor < numberAllFactors; iFactor++)
        *theStream << setw(6) << origPathDepList[iFactor];
    *theStream << endl;
    // Brownian bridge structure dump if there
    if(numBBFactors !=0) {
        tnti_mat BBstructure;
        BrownianBridgeSequenceSP BBptr = DYNAMIC_POINTER_CAST<BrownianBridgeSequence>(baseSequence);
        BBptr->getBrownianBridgeStructure(BBstructure);
        *theStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        *theStream << "Labels: -1 = PSR ; 0 = grid  ; > 0 sobol number used " << endl;
        *theStream << "Primary partition factors and layout" << endl;
        *theStream <<setw(10) << "PrimFact" << ":";
        for(iFactor = 0 ; iFactor < numBBFactors; iFactor++)
            *theStream << setw(6) << invImportance[iFactor];
        *theStream << endl;
        for(iDate = numBBdates-1; iDate >=0; iDate--) {
            *theStream <<setw(10) << iDate << ":" ;
            for(iFactor = 0 ; iFactor < numBBFactors; iFactor++) {
                *theStream << setw(6) << BBstructure[iFactor][iDate];
            }//iFactor
            *theStream << endl;
        }//iDate
    } // if over numBBFactors
}// void SuperCubeBBridge::DumpSCStructure(ostream *theStream);

//**************************************************************
// Method to computing index maps for getting deviates for
// getDeviateByPathDate(...).
//
// New version of initPathIndexMaps()..
// Populates pathDateIndexMap
void SuperCubeBBridge::initPathDateIndexMaps()
{
    int iFactor, iPartition, iStart;
    //Original SC representation
    tnti_vec baseIndex;
    baseIndex.resize(numBBFactors*nPartitions);

    // Given an initial index array in the full supercube
    for(iFactor = 0; iFactor < numBBFactors*nPartitions; iFactor++) {
        baseIndex[iFactor] = iFactor;
        //cout << "baseIndex[" << iFactor << "]="<<baseIndex[iFactor]<<endl;
    }//iFactor
    // Intermediate representation where the secondary partitions are reshuffled
    tnti_vec interIndex;
    interIndex.resize(numBBFactors*nPartitions);
    tnti_vec shuffleTable;
    shuffleTable.resize(numBBFactors);
    // Primary partition as is
    iPartition = 0;
    for(iFactor = 0; iFactor < numBBFactors; iFactor++) {
        interIndex[iFactor] = baseIndex[iFactor];
        //cout << "interIndex[" << iFactor << "]="<<interIndex[iFactor]<<endl;
    }//iFactor
    // Create the shuffle table from the primary partition
    int path1Index = 0;
    int path0Index = numBBFactors-1;
    for(iFactor = 0; iFactor < numBBFactors; iFactor++) {
        if(scbbPathDepList[pathDepMapping[iFactor]] == 1) {
            shuffleTable[path1Index] = iFactor;//relative to primary
            path1Index++;
        } else {
            shuffleTable[path0Index] = iFactor;//relative to primary
            path0Index--;
        }
    }//iFactor
    // Loop over secondary partitions
    //
    for(iPartition = 1; iPartition < nPartitions; iPartition++) {
        iStart = iPartition * numBBFactors;
        // Given the path-dependence of this partition reassign factors
        // based on primary partition shuffle table
        path1Index =0;
        path0Index =numBBFactors-1;
        for(iFactor = 0; iFactor < numBBFactors; iFactor++) {
            if(iFactor + iStart < numberAllBridgeFactors) {
                if(scbbPathDepList[pathDepMapping[iFactor+iStart]] == 1) {
                    interIndex[iFactor+iStart] = shuffleTable[path1Index] + iStart; // translate to current partition
                    path1Index++;
                } else {
                    interIndex[iFactor+iStart] = shuffleTable[path0Index] + iStart; // translate to current partition
                    path0Index--;
                }
                //cout << "interIndex[" << iFactor +iStart<< "]="<<interIndex[iFactor+iStart]<<endl;
            }// if over iFactor+iStart
        }//iFactor
    }//iPartition
    //
    // At this point we have a well defined interIndex array with nPartitions by numBBFactors elements.
    // This is in the SuperCube representation SC2 (remapped by pathDepMapping)
    // Map to intermediate SuperCube representation
    tnti_vec fullIndex;
    fullIndex.resize(numberAllFactors);
    for(int iAll=0; iAll < numberAllFactors;iAll++) {
        fullIndex[iAll] = numberAllFactors-1;
    }//iAll
    for(iFactor = 0; iFactor < numberAllBridgeFactors; iFactor++) {
        fullIndex[pathDepMapping[iFactor]] = interIndex[iFactor];
        //cout << "fullIndex[" << pathDepMapping[iFactor] << "]="<<interIndex[iFactor]<<endl;
    }
    // Memory allocation check
    if(!pathDateIndexMap) {
        pathDateIndexMap = new int [numberAllFactors];
        if(!pathDateIndexMap) {
            throw SCException(__FILE__, __LINE__,
                              "Could not allocate memory in initPathDateIndexMaps");
        }
    }
    // Initialize array
    for(int i = 0; i < numberAllFactors; i++) {
        pathDateIndexMap[i] = i;
    }
    for(iFactor = 0; iFactor < numberAllFactors; iFactor++) {
        pathDateIndexMap[iFactor] = fullIndex[iFactor];
    }
}// void SuperCubeBBridge::initPathDateIndexMaps()


//====================================================================
// Method to set multiplier for Low Discrepancy sequences
// Should be set prior to initializeObject(...).
void SuperCubeBBridge::setLDSMultiplier(const int thisValue)
{
    LDS_multiplier = thisValue;
}// void setLDSMultiplier(const int thisValue){

//====================================================================
// Method to get multiplier for Low Discrepancy sequences
int SuperCubeBBridge::getLDSMultiplier()
{
    return LDS_multiplier;
}// int getLDSMultiplier(){
//====================================================================
// Method to dump the path-dimension table to a given ostream
void SuperCubeBBridge::dumpLDSDimensionTable(ostream *thisStream)
{
    int factor = 1;
    *thisStream << "============== Table of #paths to #LDS dimensions " << endl;
    for(int iRow = 0 ; iRow < pathDimTableSize; iRow++) {
        *thisStream << "#paths = " << pathDimMin * factor
        << " : #dimensions = " << pathDimTable[iRow] << endl;
        factor *= 2;
    }//iRow
}// void dumpLDSDimensionTable(ostream *thisStream);
//=====================================================================
//
// Method for index map creation.
// This is a new version geared towards fetching deviates by factor and
// date for all paths.
// Takes as input:
//    importance array
//    path-dependence array
// Creates the following maps
//    1. Map of global factor to SC factor if path-dep not -1
//       First map global factor according to importance.
//       Second, strip out -1 factors - importance is in natural ordering
//       at this point.
//       Third, create a map containing global factor on one hand and SC
//       factor on the other. Add -1 factors to pseudo-random map perhaps.
//
//    2. Mapping of SC-factor to partition level factor.
//       First, figure out which partition the factor belongs in.
//       Second, if a secondary partition rearrange factor assignment
//       based upon path-dependence. Do the creation of a list of 1-factors
//       and 0-factors within the primary partition and then dole those out
//       for the secondary partition in question.
//       Third, create a map that maps SC-factor to actual scfactor in the
//       Supercube.
//
//    3. Create the mapping which is the convolution of maps 1 and 2. This
//       Gives the mapping from a global factor to an actual supercube factor
//
//=====================================================================
void SuperCubeBBridge::initializeIndexMaps(int &numberLDSBridge, int &numberLDSFinal)
{
    int iFactor,iPartition;
    SC_int_map ext2scIndexMap; // Index map from external to supercube index mapping keeping out PSR factors
    SC_int_map sc2partitionIndexMap; // Index map from supercube index to partition
    SC_int_map sc2partFactorIndexMap; // Index map from supercube index to partition level factor
    // importance array contains the global index map
    // natural ordering is 0,1,2,3,...
    // Create new path-dependence mapping for all factors
    impPathDepList = new int [numberAllFactors];
    if(!impPathDepList) {
        throw SCException(__FILE__,__LINE__,
                          "Could not allocate memory for impPathDepList in SuperCubeBBridge::initializeIndexMaps");
    }//!impPathDepList
    // Reassign path-dependence in natural representation - ie increasing importance 0,1,...,n
    for(iFactor=0; iFactor<numberAllFactors; iFactor++) {
        impPathDepList[importance[iFactor]] = origPathDepList[iFactor];
    }//iFactor
    // Correct till now
    // Now construct an index array for path-dependence without the -1 factors.
    scbbPathDepList = new int [numberAllFactors];
    if(!scbbPathDepList) {
        throw SCException(__FILE__,__LINE__,
                          "Could not allocate memory for scbbPathDepList in SuperCubeBBridge::initializeIndexMaps");
    }//!scbbPathDepList
    numberAllBridgeFactors = 0;
    for(iFactor=0; iFactor<numberAllFactors; iFactor++) {
        if(impPathDepList[iFactor] !=-1) {
            scbbPathDepList[numberAllBridgeFactors] = impPathDepList[iFactor];
            // add entry for external index to SC-index - keeping out -1 factors
            ext2scIndexMap[invImportance[iFactor]] = numberAllBridgeFactors;
            //       cout << "ext2scIndexMap : " << iFactor << " "
            //     << invImportance[iFactor] << " "
            //     << numberAllBridgeFactors << endl;
            numberAllBridgeFactors++;
        }// not -1 factor
    }//iFactor
    //  cout << "---------------------------" << endl;

    //
    // Set up the partition structure
    getPartitionData(numberLDSBridge, numberLDSFinal);
    //
    // With the given partition structure find out how supercube factors map
    // in given multiple partition.
    // First construct list of path-dep (1) and non-path-dep factors (0) in
    // primary partition.
    int *pathDepList = new int[numBBFactors];
    // Go through primary partition and tally up path-dep (1) factors and
    // non-path-dep (0) factors.
    int numPathDepFactors = 0;
    int numNonPathDepFactors = numBBFactors-1;
    for(iFactor=0;iFactor <numBBFactors; iFactor++) {
        if(scbbPathDepList[iFactor] == 1) {
            pathDepList[numPathDepFactors] = iFactor;
            numPathDepFactors++;
        } else {
            pathDepList[numNonPathDepFactors] = iFactor;
            numNonPathDepFactors--;
        }
    }//iFactor
    //   for(iFactor=0;iFactor <numBBFactors; iFactor++){
    //     cout << "pathDepList: "
    //   << iFactor << " "
    //   << pathDepList[iFactor] << " "
    //   << scbbPathDepList[pathDepList[iFactor]] << endl;
    //   }//iFactor

    int oldPartition = -9;
    for(iFactor=0; iFactor<numberAllBridgeFactors; iFactor++) {
        // For the given factor obtain the partiton id
        iPartition = iFactor / numBBFactors; // floor
        //    cout << "iPartition = " << iPartition << endl;
        if(iPartition > nPartitions-1) {
            throw SCException(__FILE__, __LINE__,
                              "SuperCubeBBridge::initializeIndexMaps: factor-partition error: iPartition>nPartitions-1") ;
        }//iPartition>nPartitions-1
        // Add to map for partition
        sc2partitionIndexMap[iFactor] = iPartition;
        // determine the factor within the partition in question
        int partIndex = iFactor - numBBFactors * iPartition;
        if(iPartition == 0) {
            sc2partFactorIndexMap[iFactor] = partIndex;
        } else {
            // If not primary partition then ration out path-dependent and non-path dependent
            // factors accordingly.
            if(iPartition!=oldPartition) {//reset path-dep lists
                numPathDepFactors = 0;
                numNonPathDepFactors = numBBFactors-1;
                oldPartition = iPartition;
            }
            if(scbbPathDepList[iFactor] == 1) {
                sc2partFactorIndexMap[iFactor] = pathDepList[numPathDepFactors];
                numPathDepFactors++;
            } else {
                sc2partFactorIndexMap[iFactor] = pathDepList[numNonPathDepFactors];
                numNonPathDepFactors--;
            }
        }
        //     cout << "Internal factor : " << iFactor
        //   << " outside factor : " << invImportance[iFactor]
        //   << " partition      : " << sc2partitionIndexMap[iFactor]
        //   << " sc2partFactor  : " << sc2partFactorIndexMap[iFactor] << endl;
    }//iFactor
    //   cout << "-----------------------" << endl;

    // Now convolve the two sets of maps
    // Note that global2sc* are in external representation.
    // sc2* is in SCBB representation without PSR
    for(iFactor=0;iFactor<numberAllFactors; iFactor++) {
        // remember to map to global index space.
        if(origPathDepList[iFactor] == -1) {
            global2scFactorMap[iFactor]    = -1;
            global2scPartitionMap[iFactor] = -1;
        } else {
            global2scFactorMap[iFactor] = sc2partFactorIndexMap[ext2scIndexMap[iFactor]];
            global2scPartitionMap[iFactor] = sc2partitionIndexMap[ext2scIndexMap[iFactor]];
        }
    }//iFactor

    //   for(iFactor=0; iFactor<numberAllFactors; iFactor++){
    //     cout << iFactor << " " << global2scFactorMap[iFactor]
    //       << " " << global2scPartitionMap[iFactor] << " "
    //       << scbbPathDepList[iFactor] << endl;
    //   }//iFactor

    // Clean up local variables
    sc2partFactorIndexMap.erase(sc2partFactorIndexMap.begin(),
                                sc2partFactorIndexMap.end()   );
    sc2partitionIndexMap.erase(sc2partitionIndexMap.begin(),
                               sc2partitionIndexMap.end()     );
    ext2scIndexMap.erase(ext2scIndexMap.begin(),
                         ext2scIndexMap.end()     );
    delete [] pathDepList;
}// void SuperCubeBBBridge::initializeIndexMaps()
//====================================================================
// The following method reinitializes the Supercube with a new seed.
// Internally it retains memory allocated and reinitializes to use
// the new seed being passed in. This method also advances the Sobol
// sequence used internally.
void SuperCubeBBridge::reinitialize(const INormalSuperCubeRNGSP rng,
                                    const int numberLDSIterations)
{
    if(theRandomGrid.get() == NULL) {
        SuperCube::reinitialize(rng,numberLDSIterations,internalGaussianGrid);
    } else {
        if(whichRun == numberRuns-1) {
            cout << whichRun << " " << numberRuns << endl;
            throw SCException(__FILE__, __LINE__,
                              "SuperCubeBBridge::reinitialize(...) - grid already fully sampled - number of runs exceeded");
        }
        // Fetch another random sample of the grid
        theRandomGrid->setRNG(this->getUniformRNG());
        theRandomGrid->populateVector();
        internalGaussianGrid = theRandomGrid->getVector();
        whichRun++;
        SuperCube::reinitialize(rng /*-1*theRandomGrid->getSeed()*/,
                                numberLDSIterations,
                                internalGaussianGrid);
    }
}
CORE_END_NAMESPACE
