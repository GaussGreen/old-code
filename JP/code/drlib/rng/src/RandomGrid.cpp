// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                        *
 * RandomGrid.cpp                                                         *
 *      Implementation for RandomGrid class.                              *
 * Author: M.Huq, Credit DR                                               *
 *                                                                        *
 -----------------------------------------------------------------------*/


#include "edginc/RandomGrid.h"
#include "edginc/SCException.h"

#include "InvertCumNormal.h"
#include <iostream>

CORE_BEGIN_NAMESPACE

using namespace std;

// Default constructor
RandomGrid::RandomGrid() : Sequence()
{
    m_numberSamples = 0;
    m_randomGridCharMap = 0;
    m_theFullGrid = 0;
    m_theRandomSampleGrid = 0;
}

/*
RandomGrid::initializeRandomGrid(const int numberPoints,
                 const int numberSamples){
  */
void RandomGrid::initializeRandomGrid(IUniformRNGSP _uniform,
                                      const int numberPoints,
                                      const int numberSamples,
                                      const double min,
                                      const double max)
{
    if(numberPoints < 1) {
        throw SCException(__FILE__, __LINE__,
                          "numberPoints <1 in RandomGrid::initializeRandomGrid");
    }
    if(numberSamples < 1) {
        throw SCException(__FILE__, __LINE__,
                          "numberSamples <1 in RandomGrid::initializeRandomGrid");
    }
    if(max > 1.0 || max < 0.0) {
        throw SCException(__FILE__, __LINE__,
                          "max must have value in [0,1] in RandomGrid::initializeRandomGrid.");
    }
    if(min < 0.0 || min > 1.0) {
        throw SCException(__FILE__, __LINE__,
                          "min must have value in [0,1] in RandomGrid::initializeRandomGrid.");
    }
    m_gridSpacing = (max - min)/(double)(numberPoints*numberSamples);
    if(m_gridSpacing == 0.0) {
        throw SCException(__FILE__, __LINE__,
                          "max = min in RandomGrid::initializeRandomGrid.");
    }
    nDimensions = numberPoints;
    uniform = _uniform;
    //  seed = theSeed;
    m_numberSamples = numberSamples;
    m_sampleCount = 0;
    m_randomGridCharMap = new char [numberPoints * numberSamples];
    if(!m_randomGridCharMap) {
        throw SCException(__FILE__, __LINE__,
                          "RandomGrid::initializeRandomGrid - could not allocate memory for m_randomGridCharMap");
    }
    m_theFullGrid = new double [numberPoints * numberSamples];
    if(!m_theFullGrid) {
        throw SCException(__FILE__, __LINE__,
                          "RandomGrid::initializeRandomGrid - could not allocate memory for m_theFullGrid");
    }
    //   vector = new double [numberPoints];
    //   if(!vector){
    //     throw SCException(__FILE__, __LINE__,
    //            "RandomGrid::initializeRandomGrid - could not allocate memory for vector");
    //   }
    vector.assign(numberPoints, 0.0);

    m_fullGridMin = min + m_gridSpacing*0.5;
    m_fullGridMax = max - m_gridSpacing*0.5;

    int iPath;
    //Initialize the full grid
    for(iPath = 0; iPath < nDimensions*m_numberSamples; iPath++) {
        m_theFullGrid[iPath] = m_fullGridMin + iPath * m_gridSpacing;
    }//iPath
    for(iPath = 0; iPath < nDimensions*m_numberSamples; iPath++) {
        m_theFullGrid[iPath] = SC_InvertCumNormal(m_theFullGrid[iPath]);
        //cout <<iPath << " "<< m_theFullGrid[iPath] << endl;
    }//iPath
    // initialize the characteristic map
    for(iPath = 0; iPath < nDimensions*m_numberSamples; iPath++) {
        m_randomGridCharMap[iPath] = 0;
    }//iPath

}

RandomGrid::~RandomGrid()
{
    if(m_randomGridCharMap) {
        delete [] m_randomGridCharMap;
        m_randomGridCharMap = 0;
    }
    if(m_theFullGrid) {
        delete [] m_theFullGrid;
        m_theFullGrid = 0;
    }
    //   if(vector){
    //     delete [] vector;
    //     vector = 0;
    //   }
} // RandomGrid::~RandomGrid()

void RandomGrid::populateVector()
{
    int iSample;
    int random_index;

    if(m_sampleCount > m_numberSamples-1) {
        throw SCException(__FILE__, __LINE__,
                          "RandomGrid::populateVector - already sampled entire grid.  m_sampleCount > m_numberSamples");
    }
    // Check that there are enough points left to sample.
    int iFull = 0;
    int freeCount = 0;
    while(iFull < nDimensions*m_numberSamples) {
        if(m_randomGridCharMap[iFull]==0) {
            freeCount++;
        }
        iFull++;
    }//iFull
    if(freeCount < nDimensions) {
        throw SCException(__FILE__,
                          __LINE__,
                          "RandomGrid::populateVector : Insufficient gridpoints left to sample. nPoints by nSamples should be <= nFullGridSize");
    }
    for(iSample = 0; iSample < nDimensions; iSample++) {
        do {
            random_index = (int)(((nDimensions*m_numberSamples)*uniform->fetch()));
            if(random_index < 0 || random_index > nDimensions*m_numberSamples-1) {
                throw SCException(__FILE__, __LINE__,
                                  "RandomGrid::populateVector : random index out of range");
            }
        } while(m_randomGridCharMap[random_index]==1);
        vector[iSample] = m_theFullGrid[random_index];
        m_randomGridCharMap[random_index] = 1;
    }//iSample
    m_sampleCount ++;
}// void RandomGrid::populateVector();
CORE_END_NAMESPACE
