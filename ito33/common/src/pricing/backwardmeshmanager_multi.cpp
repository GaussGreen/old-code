/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/backwardmeshmanager_multi.cpp
// Purpose:     base multi-grid mesh manager for backward PDE problems
// Author:      Nabil
// Created:     2004/03/17
// RCS-ID:      $Id: backwardmeshmanager_multi.cpp,v 1.10 2004/10/05 09:13:46 pedro Exp $
// Copyright:   (c) 1999-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/pricing/backwardmeshmanager_multi.h"


namespace ito33
{

namespace pricing
{

void BackwardMeshManagerMulti::SetupMe()
{
  BackwardMeshManager::SetupMe();

  ConstructSpaceMesh();
}

void BackwardMeshManagerMulti::GoAhead()
{
  if(m_nIdxTimeGrid == 0) // We are at a time of change of grids.
  {
    // Save data of the mesh before the change of grid.
    m_nOldNbS = m_nNbS;
    for(size_t nIdxS = 0; nIdxS < m_nOldNbS; nIdxS++)
    {
      m_pdOldS[nIdxS] = m_pdS[nIdxS];
      m_pdOldLogS[nIdxS] = m_pdLogS[nIdxS];
    }
    
    m_nIdxGrid--;
    m_nIdxTimeGrid = m_pGrids[m_nIdxGrid].m_nNbTimes - 1;
    m_nNbS = m_pGrids[m_nIdxGrid].m_nNbS;
    UpdateSpaceMeshes();
    m_bIsEndOfGrid = true;
  }    
  else if(m_bIsInitialTime) 
  {
    m_bIsInitialTime = false;
    // don't need to change grid
    m_nOldNbS = m_nNbS;
    for(size_t nIdxS = 0; nIdxS < m_nOldNbS; nIdxS++)
      m_pdOldS[nIdxS] = m_pdS[nIdxS];

    m_bIsEndOfGrid = true;
  }
  else
  {
    BackwardMeshManager::GoAhead();
    
    m_nIdxTimeGrid--;
    UpdateSpaceMeshes();
    m_bIsEndOfGrid = false;
  }
}

void BackwardMeshManagerMulti::UpdateSpaceMeshes()
{
  size_t
    nIdxS;

  if( IsChangeOfVar() )
    for(nIdxS = 0; nIdxS < m_pGrids[m_nIdxGrid].m_nNbS; nIdxS++)
      m_pdS[nIdxS] = m_pGrids[m_nIdxGrid].m_pdS[nIdxS] 
        * m_pGrids[m_nIdxGrid].m_pdCorrectionS[m_nIdxTimeGrid];
  else
    for(nIdxS = 0; nIdxS < m_pGrids[m_nIdxGrid].m_nNbS; nIdxS++)
      m_pdS[nIdxS] = m_pGrids[m_nIdxGrid].m_pdS[nIdxS]; 
  
  //log S
  for(nIdxS = 0; nIdxS < m_pGrids[m_nIdxGrid].m_nNbS; nIdxS++)
    m_pdLogS[nIdxS] = log(m_pdS[nIdxS]); 
}

} // namespace pricing

} // namespace ito33

