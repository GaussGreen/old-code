/////////////////////////////////////////////////////////////////////////////
// Name:        forwardcds/forwardcdsnumoutput.cpp
// Purpose:     implementation of CDS NumOutput class using forward PDE 
// Author:      David
// RCS-ID:      $Id: forwardcdsnumoutput.cpp,v 1.13 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"


#include "ihg/forwardcdsinstdata.h"
#include "ihg/forwardcdsnumoutput.h"

#include "ito33/ihg/modeloutput.h"

#include "ito33/numeric/predicatetime.h"

using namespace ito33::ihg;

// implement the AutoPtrDeleter for ForwardCDSNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::ForwardCDSNumOutput);
}

void ForwardCDSNumOutput::Init(ForwardCDSInstData& /* instdata */)
{
  // Nothing to do
}

void ForwardCDSNumOutput::UpdateMe
     (ForwardCDSInstData& instdata, double dTime)
{
  // Save the price and the terms that determine the price at this timestep.
  // Check if this function has been called more than once at the same time.
  if ( m_Times.size() > 0 && numeric::AreTimesEqual(dTime, m_Times.back() ) )
  {
    // Same as previous time (e.g. an event occured)
    // Over-write the previous data
    m_CDSPrices.back() = instdata.m_dCDSPrice;

    m_RecoveryTerms.back() = instdata.m_dRecoveryTerm;
    m_SpreadTerms.back() = instdata.m_dSpreadTerm;
    m_AccruedTerms.back() = instdata.m_dAccruedTerm;
  }
  else
  {
    // New timestep
    m_Times.push_back(dTime);

    m_CDSPrices.push_back(instdata.m_dCDSPrice);

    m_RecoveryTerms.push_back(instdata.m_dRecoveryTerm);
    m_SpreadTerms.push_back(instdata.m_dSpreadTerm);
    m_AccruedTerms.push_back(instdata.m_dAccruedTerm);
  }

  ASSERT_MSG(m_Times.size() == instdata.m_nCurrentIndex + 1, 
    "Data in numoutput is out of synch with the time mesh");
}


void ForwardCDSNumOutput::Finalize(ForwardCDSInstData& /* instdata */)
{
  // Convert the lists into vectors
  size_t nNbValues = m_Times.size();

  m_pdTimes.resize(nNbValues);
  m_pdCDSPrices.resize(nNbValues);
  m_pdAccruedTerms.resize(nNbValues);
  m_pdSpreadTerms.resize(nNbValues);
  m_pdRecoveryTerms.resize(nNbValues);

  size_t nIdx;
  std::list<double>::const_iterator iterTimes = m_Times.begin();
  std::list<double>::const_iterator iterPrices = m_CDSPrices.begin();
  std::list<double>::const_iterator iterAccrued = m_AccruedTerms.begin();
  std::list<double>::const_iterator iterSpread = m_SpreadTerms.begin();
  std::list<double>::const_iterator iterRecovery = m_RecoveryTerms.begin();

  for (nIdx = 0; nIdx < nNbValues; nIdx++)
  {
    m_pdTimes[nIdx] = *iterTimes;
    m_pdCDSPrices[nIdx] = *iterPrices;
    m_pdRecoveryTerms[nIdx] = *iterRecovery;
    m_pdSpreadTerms[nIdx] = *iterSpread;
    m_pdAccruedTerms[nIdx] = *iterAccrued;

    iterTimes++;
    iterPrices++;
    iterRecovery++;
    iterSpread++;
    iterAccrued++;
  }
}

ito33::shared_ptr<ModelOutput> ForwardCDSNumOutput::GetModelOutput()
{ 
  // Construct the output class
  ito33::shared_ptr<ModelOutput> pOutput (new ModelOutput());

  pOutput->SetPrice(m_pdCDSPrices.back());

  return pOutput;
}
