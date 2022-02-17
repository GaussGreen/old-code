// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2001 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
/*                                                                        *
 * RandomGrid.h                                                           *
 *      Implementation for RandomGrid class.                              *
 * Author: M.Huq, Credit DR                                               *
 *                                                                        *
 -----------------------------------------------------------------------*/

#ifndef _SC_RANDOMGRID_H
#define _SC_RANDOMGRID_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "edginc/Sequence.h"
#include "edginc/SCException.h"
#include "edginc/DECLARESP.h"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

//#include "edginc/ran2.h"
class  RNG_DLL RandomGrid : public Sequence
{
protected:
    // Number of samples to be obtained
    int m_numberSamples;

    // Sample counter
    int m_sampleCount;

    // Characteristic map for storing which points on the grid have been
    // used.
    char *m_randomGridCharMap;
    // Double pointer to full grid array
    double *m_theFullGrid;
    // Double pointer to random sample grid array
    double *m_theRandomSampleGrid;
    // Grid min value
    double m_fullGridMin;
    // Grid max value
    double m_fullGridMax;
    // Grid spacing
    double m_gridSpacing;

public:
    // Default constructor
    RandomGrid();

    void initializeRandomGrid(IUniformRNGSP uniform,
                              const int numberPoints,
                              const int numberSamples,
                              const double min = 0.0,
                              const double max = 1.0);
    virtual ~RandomGrid();
    virtual double getGridSpacing()
    {
        return m_gridSpacing;
    };
    virtual int getSampleCount()
    {
        return m_sampleCount;
    };
    virtual void populateVector();
    virtual void populateVector(int iPath)
    {
        throw SCException(__FILE__, __LINE__,
                          "RandomGrid::populateVector(int) not implemented - use RandomGrid::populateVector(void) instead");
    };

    void setRNG(IUniformRNGSP _uniform)
    {
        uniform = _uniform;
    }
private:
    IUniformRNGSP uniform;
};

DECLARESP(RandomGrid);

CORE_END_NAMESPACE

#endif // _SC_RANDOMGRID_H
