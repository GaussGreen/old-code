/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/arraystructasgrid.h
// Purpose:     class for array structured as grid with change of variable.
// Author:      Nabil
// Created:     2003/12/04
// RCS-ID:      $Id: arraystructasgrid.h,v 1.5 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/mesh/arraystructasgrid.h
    @brief class for array structured as grid with change of variable.

    Implementation of the ArrayStructAsGrid template class.
    
 */

#ifndef _ITO33_NUMERIC_MESH_ARRAYSTRUCTASGRID_H_
#define _ITO33_NUMERIC_MESH_ARRAYSTRUCTASGRID_H_

#include "ito33/array.h"
#include "ito33/numeric/mesh/gridwithcv.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{

template <class T>
class ArrayStructAsGrid
{
public:

  ArrayStructAsGrid(Array<GridWithCV>& pGrids, size_t nNbGrids)
  {
    size_t
      nI,
      nNbTotalTime;

    nNbTotalTime = 0;
    for(nI = 0; nI < nNbGrids; nI++)
      nNbTotalTime += pGrids[nI].m_nNbTimes;

    m_ppArrayStruct = new T *[nNbGrids];
    m_ppArrayStruct[0] = new T [nNbTotalTime];
    for(nI = 1; nI < nNbGrids; nI++)
      m_ppArrayStruct[nI] = m_ppArrayStruct[nI - 1] + 
        pGrids[nI - 1].m_nNbTimes;
  }
  
  ~ArrayStructAsGrid()
  {
    delete [] m_ppArrayStruct[0];
    delete [] m_ppArrayStruct;
  }
  
  T& operator()(size_t nIdxGrid, size_t nIdxTime)
  {
    return m_ppArrayStruct[nIdxGrid][nIdxTime];
  }

protected:
  
  T **m_ppArrayStruct;

};

} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_ARRAYSTRUCTASGRID_H_
