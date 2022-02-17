/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/instdata.cpp
// Author:      Wang
// Created:     2004/02/10
// RCS-ID:      $Id: instdata.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdata.h"

using namespace ito33::pricing;

void InstData::Alloc(size_t nNbS)
{
  // Allocate the price arrays
  m_pdPrices = Array<double>(nNbS);
  m_pdOldPrices = Array<double>(nNbS);
  m_pdOldOldPrices = Array<double>(nNbS);
 
  // Allocate storage for delta and gamma (used by Greek PDEs)
  m_pdDeltas = Array<double>(nNbS);
  m_pdGammas = Array<double>(nNbS);
}

void InstData::Swap()
{
  swap(m_pdOldOldPrices, m_pdPrices);
  swap(m_pdOldOldPrices, m_pdOldPrices);
}

void InstData::UpdateBeforeStep()
{
  InstDataTimeOnly::UpdateBeforeStep();

  Swap();
}
