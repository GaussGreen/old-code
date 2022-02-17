/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/instdata.h
// Purpose:     implementation for ihg instdata class
// Author:      David Pooley
// Created:     2003/12/10
// RCS-ID:      $Id: instdata.cpp,v 1.10 2004/10/04 18:04:07 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/event.h"

#include "ihg/model.h"
#include "ihg/instdata.h"

using ito33::ihg::InstData;

void InstData::Alloc(size_t nNbS)
{
  InstData::BaseClass::Alloc(nNbS);
  
  // Allocate space for ihg specific arrays 
  m_pdVolsSquared = Array<double>(nNbS);
  m_pdHazardRates = Array<double>(nNbS);
}

void InstData::UpdateBeforeStep()
{
  InstData::BaseClass::UpdateBeforeStep();

  // Get the squared vols, needed by the price PDE
  m_model.GetVolsSquared(m_pdS, m_pdVolsSquared.Get(), m_nNbS);

  // Get the hazard rates
  m_model.GetHazardRates(m_pdS, m_pdHazardRates.Get(), m_nNbS); 
}

void InstData::ApplyEvent(const ito33::pricing::Event *pEvent)
{
  pEvent->ApplyToPrice(m_pdS, m_pdPrices.Get(), m_nNbS);
}

