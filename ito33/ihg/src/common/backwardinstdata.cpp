/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/backwardinstdata.h
// Purpose:     implementation for ihg backward instdata class
// Author:      Wang
// Created:     2004/02/13
// RCS-ID:      $Id: backwardinstdata.cpp,v 1.10 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/event.h"
#include "ito33/pricing/params.h"

#include "ihg/model.h"
#include "ihg/backwardinstdata.h"

using ito33::ihg::BackwardInstData;

void BackwardInstData::Alloc(size_t nNbS)
{
  InstData::Alloc(nNbS);

  // allocate the vol array for vega computation
  if (m_bComputeVega)
    m_pdVols = Array<double>(nNbS);
  
  // Allocate the vega arrays
  if (m_bComputeVega)
  {
    m_pdVegas = Array<double>(nNbS);
    m_pdOldVegas = Array<double>(nNbS);
    m_pdOldOldVegas = Array<double>(nNbS);
  }
  
  // Allocate the fugit arrays
  if (m_bComputeFugit)
  {
    m_pdFugits = Array<double>(nNbS);
    m_pdOldFugits = Array<double>(nNbS);
    m_pdOldOldFugits = Array<double>(nNbS);
  }
}

void BackwardInstData::Swap()
{
  InstData::Swap();

  if (m_bComputeVega)
  {
    swap(m_pdOldOldVegas, m_pdVegas);
    swap(m_pdOldOldVegas, m_pdOldVegas);
  }
  
  if (m_bComputeFugit)
  {
    swap(m_pdOldOldFugits, m_pdFugits);
    swap(m_pdOldOldFugits, m_pdOldFugits);
  }
}

void BackwardInstData::UpdateBeforeStep()
{
  InstData::UpdateBeforeStep();

  // Get the volatilities, needed for vega compute
  if (m_bComputeVega)
    m_model.GetVols(m_pdS, m_pdVols.Get(), m_nNbS);
}

void BackwardInstData::ApplyEvent(const ito33::pricing::Event *pEvent)
{
  InstData::ApplyEvent(pEvent); 

  if (m_bComputeVega)
    pEvent->ApplyToGreek(m_pdS, m_pdVegas.Get(), m_nNbS);
  
  if (m_bComputeFugit)
    pEvent->ApplyToGreek(m_pdS, m_pdFugits.Get(), m_nNbS);
}

double BackwardInstData::GetInitialSpot() const
{
  return m_params.GetSpotSharePrice(); 
}
