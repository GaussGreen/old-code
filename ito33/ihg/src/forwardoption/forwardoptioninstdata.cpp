/////////////////////////////////////////////////////////////////////////////
// Name:        forwardoption/forwardoptioninstdata.cpp
// Purpose:     implementation of option instdata class using forward PDE 
// Author:      Wang
// RCS-ID:      $Id: forwardoptioninstdata.cpp,v 1.7 2004/10/04 18:04:08 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/payoff.h"

#include "ihg/forwardoptioninstdata.h"

using ito33::ihg::ForwardOptionInstData;

void ForwardOptionInstData::Alloc(size_t nNbS)
{
  InstData::Alloc(nNbS);
}

void ForwardOptionInstData::Init()
{
  // Get the space mesh
  m_pdS = m_forwardOptionMeshes.GetS();

  m_nNbS = m_forwardOptionMeshes.GetNbS();
										                		
  m_pdLogS = m_forwardOptionMeshes.GetLogS();

  /*
  m_BoundaryCondition.SetLeft(ito33::numeric::BCType_Dirichlet, 
                              m_forwardOptionParams.GetSpotSharePrice() - (m_forwardOptionMeshes.GetS())[0] );    
  m_BoundaryCondition.SetRight(ito33::numeric::BCType_Dirichlet, 0.0 );    
  */

  Alloc(m_nNbS);    
}


void ForwardOptionInstData::SetInitialValue()
{
  // Setup initial conditions.    
  m_forwardOptionParams.GetPayoff()->Get(m_pdS, m_pdPrices.Get(), m_nNbS);
}


void ForwardOptionInstData::UpdateBeforeStep()
{
  InstData::UpdateBeforeStep();

  m_bIsHRTimeOnly = m_model.IsHazardRateTimeOnly();

}
