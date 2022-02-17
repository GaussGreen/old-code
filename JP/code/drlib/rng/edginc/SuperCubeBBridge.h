// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                        *
 * SuperCubeBBridge.h:                                                    * 
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

#ifndef _SC_SUPERCUBEBBRIDGE_H
#define _SC_SUPERCUBEBBRIDGE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "edginc/coreConfig.hpp"
#include "edginc/SCException.h"
#include "edginc/SuperCube.h"
#include "edginc/RandomGrid.h"

#include "edginc/IRNG.h"
#include <map>

CORE_BEGIN_NAMESPACE

using namespace std;

class RNG_DLL SuperCubeBBridge : public SuperCube
{
public:
    //========================
    // Default constructor
    SuperCubeBBridge();

    //========================================================
    // SRM-version
    // Interface assumes brownian bridge type structure laid out
    // in factors and dates.
    //  Arrays importance and pathDependence must have size numFactors.
    SuperCubeBBridge (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence);     // Array indicating path-dependence (1) or not (0).


    // This function is called by the above constructor. It is provided to get around issues
    // of scope. The recommended approach to initializing a SuperCubeBBridge is first declare
    // the object with the default constructor and then use initializeObject.
    virtual void initializeObject (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence);    // Array indicating path-dependence (1) or not (0).
    //========================================================
    // SRM-version II
    // The following version takes as an argument a parameter
    // that allows setting up an uniform grid in probability space
    // to be used by the Brownian bridge (if there is to be one).
    // The parameter is an offset from 0 or 1 probability. If epsilon is
    // chosen then we set up a grid of N points that is centered in values from
    // epsilon to 1-epsilon.
    SuperCubeBBridge (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        double epsilon);         // Parameter for setting up grid.


    // This function is called by the above constructor. It is provided to get around issues
    // of scope. The recommended approach to initializing a SuperCubeBBridge is first declare
    // the object with the default constructor and then use initializeObject.
    virtual void initializeObject (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        double epsilon);         // Parameter for setting up grid.

    //========================================================
    // SRM-version III
    // Multi-run version that uses a randomly sampled grid
    // at the final time for each run. For M runs of a N-path
    // integration the grid size will be MxN long. All M runs
    // will fully sample the entire grid. This is a useful way
    // to get higher resolution with smaller memory impact.
    SuperCubeBBridge (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        int randomGridFlag);     // Flag to toggle random grid. If > 1 the number of runs
    // assumed for random sampling the grid is
    //= randomGridFlag. Upon each reinitialization then
    // a new grid of size numberPaths is sampled.


    // This function is called by the above constructor. It is provided to get around issues
    // of scope. The recommended approach to initializing a SuperCubeBBridge is first declare
    // the object with the default constructor and then use initializeObject.
    virtual void initializeObject (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        int randomGridFlag);     // Flag to toggle random grid. If > 1 the number of runs
    // assumed for random sampling the grid is
    //= randomGridFlag. Upon each reinitialization then
    // a new grid of size numberPaths is sampled.


    //========================================================
    // Internal version that takes a gaussian grid as an argument.
    // This is used by the first three constructors above.
    SuperCubeBBridge (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        double *BBgaussianGrid); // Gaussian grid.


    // This function is called by the above constructor. It is provided to get around issues
    // of scope. The recommended approach to initializing a SuperCubeBBridge is first declare
    // the object with the default constructor and then use initializeObject.
    virtual void initializeObject (
        INormalSuperCubeRNGSP rng,              // Seed for pseudo-random generator and random-permuations of partitions
        int numberPaths,         // Number of paths to extract
        int numberFactors,       // Number of factors
        int numberDates,         // Number of dates to consider.
        int *theImportance,         // Array containing importance ordering of each factor
        int *pathDependence,     // Array indicating path-dependence (1) or not (0).
        double *BBgaussianGrid); // Gaussian grid


    // On SOLARIS compiler warnings about hiding following versions of
    // SuperCube::initializeObject. Thus declare them and call parent
    // class method. Should not have to do this?
    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              // seed for random shuffle
        int numberPaths,    // Number of paths
        int maxDimension,        // max dimensions to set up for
        int sizePartition,       // size of partitions to break up in.
        SequenceSP useSequence)
    { // Sequence to be used. Base class pointer.

        SuperCube::initializeObject(
            rng,
            numberPaths,
            maxDimension,
            sizePartition,
            useSequence,
            GET_UNIFORM);
    };
    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,// seed for pseudo-random deviates used in both shuffling and PSR part.
        int numberPaths, // number of paths of data desired.
        int maxDimension,  // Overall dimensionality of the SuperCube.
        int sizePartition,       // Dimensionality of partition within SuperCube.
        int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
        int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
        int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
        // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
        int nBBfactors,        // number of dimensions for a single BB deviate at one date.
        double *BBptr2GaussianGrid)
    {
        SuperCube::initializeObject(rng,
                                    numberPaths,
                                    maxDimension,
                                    sizePartition,
                                    numPartitions,
                                    lowDiscDimension,
                                    brownianBridgeFlag,
                                    nBBfactors,
                                    BBptr2GaussianGrid);

    };
    virtual void initializeObject(
        INormalSuperCubeRNGSP rng,              // seed for pseudo-random deviates used in both shuffling and PSR part.
        int numberPaths,         // number of paths of data desired.
        int maxDimension,        // Overall dimensionality of the SuperCube.
        int sizePartition,       // Dimensionality of partition within SuperCube.
        int numPartitions,       // Number of partitions within SuperCube. Anything outside is treated as PSR.
        int lowDiscDimension,    // Dimensionality of low discrepancy sequence to be used. If zero uses PSR.
        int brownianBridgeFlag,  // 0=off, >0=on. If on, put in a brownian bridge with this many bridges evenly spaced.
        // Users job to specify the number of bridges relative to the lowDescrepancy dimension.
        int numSBBFactors)
    {       // number of dimensions for a single BB deviate at one date.
        SuperCube:: initializeObject(rng,
                                     numberPaths,
                                     maxDimension,
                                     sizePartition,
                                     numPartitions,
                                     lowDiscDimension,
                                     brownianBridgeFlag,
                                     numSBBFactors);
    };
    virtual void initializeObject(INormalSuperCubeRNGSP rng,         // seed for random shuffle
                                  int numberPaths, // Number of paths
                                  int maxDimension,    // max dimensions to set up for
                                  int sizePartition,   // size of partitions
                                  int numPartitions,   // number of partitions to break up in.
                                  int seqDimension,   // Dimension of low-descrepency sequence
                                  SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                                  int UniformityFlag)
    {
        SuperCube::initializeObject( rng,
                                     numberPaths,
                                     maxDimension,
                                     sizePartition,
                                     numPartitions,
                                     seqDimension,
                                     useSequence,
                                     UniformityFlag);
    };
    virtual void initializeObject(INormalSuperCubeRNGSP rng,          // seed for random shuffle
                                  int numberPaths, // Number of paths
                                  int maxDimension,    // max dimensions to set up for
                                  int sizePartition,   // number of partitions to break up in.
                                  SequenceSP useSequence,  // Sequence to be used. Base class pointer.
                                  int UniformityFlag)
    {
        SuperCube::initializeObject( rng,
                                     numberPaths,
                                     maxDimension,
                                     sizePartition,
                                     useSequence,
                                     UniformityFlag);
    };



    //=================================
    // Desctructor and methods
    virtual ~SuperCubeBBridge();

    //=================================
    // The following method populates an
    // internal vector for a given path.
    // Calls getDeviateByPath(...)
    virtual void populateVector(int nPaths);
    virtual void populateVector()
    {
        throw SCException(__FILE__, __LINE__,
                          "Method will not work with SuperCubeBBridge. Use version with iPath passed in.");
    };

    //====================================================================
    // The following method returns deviates for the date index passed in.
    // iDate runs from 0 to numberDates-1.
    // The array returned contains for each factor all paths. That is,
    // F(0,0),F(0,1),F(0,2),...,F(0,numberPaths-1), F(1,0),...,F(numberFactors-1,numberPaths-1)
    virtual void getDeviateByDate(long iDate, double *returnArray);
    virtual double * getDeviateByDate(long iDate)
    {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByDate - use version getDeviateByDate(long iDate, double *&returnArray) instead");
    };
    //====================================================================
    // The following method returns deviates for the date index passed in.
    // This version returns a two-dimensional array [factor] by [path]
    // The return array must be allocated externally.
    virtual void getDeviateByDate(long iDate, double **returnArray);

    //====================================================================
    // The following method returns deviates for a date and path index
    // passed in. The array returned contains only the factors for the
    // given path-date combination. That is,
    // F(0),F(1),F(2),...,F(numberFactors-1)
    virtual void getDeviateByPathDate(long iPath, long iDate,
                                      double *returnArray);
    virtual double * getDeviateByPathDate(long iPath, long iDate)
    {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByPathDate - use version getDeviateByPathDate(long iPath, long iDate, double *&returnArray) instead");
    };
    //====================================================================
    // The following method returns deviates for a path index
    // passed in. The array returned contains the factors by number of dates.
    // That is,
    // F(0,0),F(1,0),F(2,0),...,F(numberFactors-1,0), F(0,1),...,F(numberFactors-1,numberDates-1)
    virtual void getDeviateByPath(long iPath, double *returnArray);
    virtual double * getDeviateByPath(long iPath)
    {
        throw SCException(__FILE__, __LINE__,
                          "SuperCubeBBridge::getDeviateByPath - use version getDeviateByPath(long iPath,  double *&returnArray) instead");
    };
    //====================================================================
    // The following method fetches deviates for all paths given a factor
    // and date index. The results are returned in a preallocated array
    // shown below as returnArray.
    void getDeviateByFactorDate(const long iFactor,
                                const long iDate,
                                double *returnArray);
    //====================================================================
    // Method to dump bridge properties
    virtual void DumpSCStructure(const char *fileName);

    //====================================================================
    // Methods to set and get multiplier for Low Discrepancy sequences
    // Should be set prior to initializeObject(...).
    virtual void setLDSMultiplier(const int thisValue);
    virtual int getLDSMultiplier();
    // Method to get path to LDS dimension table
    virtual void dumpLDSDimensionTable(ostream *thisStream);

    //====================================================================
    // The following method reinitializes the Supercube with a new seed.
    // Internally it retains memory allocated and reinitializes to use
    // the new seed being passed in. This method also advances the Sobol
    // sequence used internally.
    virtual void reinitialize(INormalSuperCubeRNGSP rng,
                              const int numberLDSIterations);

    //====================================================================
    // The following method reinitializes the Supercube with a new seed.
    // Internally it retains memory allocated and reinitializes to use
    // the new seed being passed in.
    virtual void reinitialize(INormalSuperCubeRNGSP rng)
    {
        SuperCube::reinitialize(rng);
    };

    //====================================================================
    // The following method reinitializes the Supercube with a new seed,
    // number of LDS iterations to advance by and a with a prespecified
    // gaussian grid.
    virtual void reinitialize(INormalSuperCubeRNGSP rng,
                              const int numberLDSIterations,
                              const double *thisGaussianGrid)
    {
        SuperCube::reinitialize(rng, numberLDSIterations, thisGaussianGrid);
    };




protected:
    int *scbbPathDepList; //Vector containing path dependence information (-1=PSR, 0=1support,1=multisupport)
    int *impPathDepList; // Vector containing path dependence information for SuperCube only (0,1 values)
    int *origPathDepList; // Vector to contain index mapping for factors
    int *importance;        // Vector to contain original importance array
    int *invImportance;        // Vector to contain inverse index mapping
    int *pathDepMapping; // Array to store path-dependence for SC
    int  *pathDateIndexMap;  // Array to store index map for getDeviateByPathDate
    int numberAllFactors;       // Number of factors being passed through importance/pathDependence
    int numberAllBridgeFactors; // SuperCube BBridge factors
    int numberBridgeSupports;   // Number of bridge supports. Determined internally.
    // Method to get partition data based on number of available LDS dimensions.
    // Sets internally the number of partitions, partition size and
    // the number of factors within the partitions.
    void getPartitionData(int &numberLDSBridge, int &numberLDSFinal );

    void initializePathDimTable(); // Method to set up path dim table.
    void initializeIndexMaps(int &numberLDSBridge,
                             int &numberLDSFinal); // Method to initialize internal index maps
    int getMaxLDSDim(const int thisNumPath); // Method to get max LDS dimensions allowed.
    int pathDimTableSize;       // Size of path-dimensions table
    int pathDimMin;    // Minimum value inside table
    int pathDimDelta;           // delta-paths for each entry in table
    int *pathDimTable;          // The table of dimensions for a given number of paths
    double *internalGaussianGrid; // Internal grid for Brownian bridge so as to make memory management easier.
    //====================================================================
    // Method to initialize index map for getDeviateByPathDate.
    // Use of such index maps internally improves performance drastically.
    void initPathDateIndexMaps();
    // multiplier for increasing Low-discrepancy dimensions
    int LDS_multiplier;
    // Maps for indexing
    typedef std::map<int, int> SC_int_map;
    SC_int_map global2scFactorMap;
    SC_int_map global2scPartitionMap;

    // Following grid pointers added for random grid sampling
    int    whichRun;         // run index
    int    numberRuns;       // Number of runs expected - full grid size = numberRuns*numPaths
    RandomGridSP  theRandomGrid; // Random grid object


};

DECLARESP(SuperCubeBBridge);
CORE_END_NAMESPACE
#endif // !defined(_SC_SUPERCUBEBBRIDGE_H)
