// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                     *
 * GridUtility.cpp:                                                    *
 * Author: M.Huq, Credit DR                                            *
 *                                                                     *
 *  $Name$
 -----------------------------------------------------------------------*/

#include "edginc/SCException.h"

#include "InvertCumNormal.h"
#include "edginc/GridUtility.h"
#include <cmath>

CORE_BEGIN_NAMESPACE

using namespace std;



//----------------------------------------------------------------------
// Function to initialize a grid
// Inputs:
// min  =  min of grid - actual min used is min+0.5/Npaths
// max  =  max of grid - actual max used is max-0.5/Npaths
// Npaths = number of paths
//
// Returns a pointer to the grid. Memory is allocated internally.
// User is expected to free this memory.
//
double *SC_initializeGrid(double min, double max, int Npaths)
{
    double *theGrid;
    if(Npaths < 1) {
        throw SCException(__FILE__,__LINE__,
                          "SC_initializeGrid : Npaths < 1");
    }
    // allocate memory for grid
    theGrid = new double [Npaths];
    if(!theGrid) {
        throw SCException(__FILE__,__LINE__,
                          "SC_initializeGrid : Could not allocate memory for Grid");
    }
    // compute the grid parameters.
    // We always offset by dx/2 from either end of the grid.
    // This way with increasing numbers of paths we are domain filling.
    double deltaGrid = (max-min)/(double)Npaths;
    double Gmin = min + deltaGrid*0.5;
    int iPath;
    for(iPath = 0; iPath < Npaths; iPath++) {
        theGrid[iPath] = Gmin + iPath * deltaGrid;
    }//iGrid
    for(iPath = 0; iPath < Npaths; iPath++) {
        theGrid[iPath] = SC_InvertCumNormal(theGrid[iPath]);
    }//iGrid

    return theGrid;
}// double *SC_initializeGrid(double min, double max, int Npaths);

//----------------------------------------------------------------------
// Functions for randomly sampling a given grid
// Inputs:
// Npaths = #npaths on grid
// min    = minimum value on grid. Actual min used is min + 0.5/Npaths
// max    = maximum value on grid. Actual max used is max - 0.5/Npaths
// Outputs:
// theFullGrid = double array containing full grid.
// theGridMap  = characteristic array storing
// Memory allocated within routine
void SC_initializeRandomGrid(int Npaths,
                             double min,
                             double max,
                             double *&theFullGrid,
                             char   *&theGridMap)
{
    int iPath;

    // Set up the full grid.
    theFullGrid = SC_initializeGrid(min, max, Npaths);
    // allocate memory for the characteristic map
    theGridMap = new char [Npaths];
    if(theGridMap==0) {
        throw SCException(__FILE__,
                          __LINE__,
                          "initializeRandomGrid : Could not allocate memory for theGridMap");
    }

    for(iPath = 0; iPath < Npaths; iPath++) {
        theGridMap[iPath] = 0;
    }//iPath

}/*
void initializeRandomGrid(int Npaths,
     double min,
     double max,
     double *&theFullGrid,
     char   *&theGridMap)
 */

//----------------------------------------------------------------------------
// Function to randomly sample full grid.

// NpathsFull = Size of complete grid
// NpathsSampleGrid = size of sampling grid- must be < NpathsFull
// fullGrid   = full grid      - already initialized
// sampleGrid = sampling grid  - already preallocated
// theGridMap = characteristic map - already preinitialized.
// fullGrid and theGridMap should be initialized with initializeRandomGrid(..)
//
void SC_getRandomSampleGrid(IUniformSuperCubeRNGSP uniform, /*long &aSeed*/
                            int NpathsFull,
                            int NpathsSampleGrid,
                            double *fullGrid,
                            double *&sampleGrid,
                            char   *&theGridMap)
{
    int iSample;
    int random_index;
    // Check that there are enough points left to sample.
    int iFull = 0;
    int freeCount = 0;
    while(iFull < NpathsFull) {
        if(theGridMap[iFull]==0) {
            freeCount++;
        }
        iFull++;
    }//iFull
    if(freeCount < NpathsSampleGrid) {
        throw SCException(__FILE__,
                          __LINE__,
                          "SC_getRandomSampleGrid : Insufficient gridpoints left to sample. Nruns by Nsamples should be <= Nfull");
    }
    for(iSample = 0; iSample < NpathsSampleGrid; iSample++) {
        do {
            random_index = (int)(((NpathsFull)* uniform->fetch() /*SC_ran2(&aSeed)*/));
            if(random_index < 0 || random_index > NpathsFull-1) {
                throw SCException(__FILE__, __LINE__,
                                  "SC_getRandomSampleGrid : random index out of range");
            }
        } while(theGridMap[random_index]==1);
        sampleGrid[iSample] = fullGrid[random_index];
        theGridMap[random_index] = 1;
    }//iSample
}
CORE_END_NAMESPACE

