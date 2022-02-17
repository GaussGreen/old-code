/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/gridwithcv.h
// Purpose:     Grid with change of variable class
// Author:      Nabil
// Created:     2003/11/19
// RCS-ID:      $Id: gridwithcv.h,v 1.13 2004/11/10 16:37:14 afrolov Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/mesh/gridwithcv.h
    @brief Grid with change of variable class

    Implementation of the grid with change of variable class.
    
 */

#ifndef _ITO33_NUMERIC_MESH_GRIDWITHCV_H_
#define _ITO33_NUMERIC_MESH_GRIDWITHCV_H_

#include <stddef.h>
#include "ito33/common.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{

class GridWithCV
{
public:

  GridWithCV():
    m_bChangeOfVar(false),
    m_bSpaceGridByRef(false),
    m_dFixedBoundary(-1),
    m_pdS(0), m_pdLogS(0),
    m_pdTimes(0), m_pdCorrectionS(0)
  {
  }

  // default copy ctor is ok

  ~GridWithCV() 
  {
    if(!m_bSpaceGridByRef)
    {
      delete [] m_pdS;
      delete [] m_pdLogS;
      m_pdS = 0;
    }
    
    delete [] m_pdCorrectionS;
    m_pdCorrectionS = 0;
  }
  
  size_t GetIndex(double dS);
  size_t GetIndexLog(double dLogS);
  
  bool HasFixedboundary() 
  { 
    return (m_nNbCurrentS != m_nNbS);
  }

  bool
    m_bChangeOfVar,
    m_bSpaceGridByRef;

  size_t
    m_nNbS,
    m_nNbTimes,
    m_nNbCurrentS,

    /// global time index of the first time point
    m_nIdxTimeBegin;

  double
    m_dFixedBoundary;

  double
    *m_pdS,      // size = m_nNbS.
    *m_pdLogS;    // size = m_nNbS.
 
  // Grid[i].m_pdTimes points on the beginning of the i^th grid in the 
  // global time mesh table. 
  double
    *m_pdTimes;      
 
  double  
    *m_pdCorrectionS;  

};

} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_GRIDWITHCV_H_
