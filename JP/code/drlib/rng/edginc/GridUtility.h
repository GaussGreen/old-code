// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                     *
 * GridUtility.h:                                                    *
 * Author: M.Huq, Credit DR                                            *
 *                                                                     *
 *  $Name$
 -----------------------------------------------------------------------*/
#ifndef _SC_GRIDUTILITY__H
#define _SC_GRIDUTILITY__H

#include "edginc/coreConfig.hpp"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

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
double *SC_initializeGrid(double min, double max, int Npaths);

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
                             char   *&theGridMap);

//----------------------------------------------------------------------------
// Function to randomly sample full grid.

// NpathsFull = Size of complete grid
// NpathsSampleGrid = size of sampling grid- must be < NpathsFull
// fullGrid   = full grid      - already initialized
// sampleGrid = sampling grid  - already preallocated
// theGridMap = characteristic map - already preinitialized.
// fullGrid and theGridMap should be initialized with initializeRandomGrid(..)

void SC_getRandomSampleGrid(IUniformSuperCubeRNGSP uniform, /*long &aSeed*/
                            int NpathsFull, ///< Size of complete grid
                            int NpathsSampleGrid, ///< size of sampling grid- must be < NpathsFull
                            double *fullGrid, ///< already initialized, should be initialized with initializeRandomGrid(..)
                            double *&sampleGrid, ///< sampling grid  - already preallocated
                            char   *&theGridMap); ///< should be initialized with initializeRandomGrid(..)

CORE_END_NAMESPACE

#endif // _SC_GRIDUTILITY__H
