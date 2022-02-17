/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/parbondmeshmanager.cpp
// Purpose:     parbond mesh manager for backward PDE problems
// Created:     2005/05/20
// RCS-ID:      $Id: parbondmeshmanager.cpp,v 1.7 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

namespace ito33
{

  using finance::CashFlowStream;

namespace pricing
{


ParBondMeshManager::ParBondMeshManager(ParBondParams &params, Model &model)
                             : BackwardMeshManagerFix(params, model),
                               m_parbondParams(params)    
{ 
}

void ParBondMeshManager::ConstructSpaceMesh()
{  
  std::vector<double> vecGrid( m_params.GenerateSpaceMesh(m_model) );

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);
  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx];
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
}

void ParBondMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  ParBond& parbond = m_parbondParams.GetParBond();
  
  shared_ptr<CashFlowStream> spreadStream = parbond.GetSpreadStream();

  CashFlowStream::const_iterator pPayment = spreadStream->begin();
 
  double dValuationTime = m_parbondParams.GetValuationTime();
  double dPaymentTime, dOldPaymentTime;

  // find out the first payment that's significant for us. 
  // What's the convention when there is a payment at valuation time?
  dOldPaymentTime = GetDoubleFrom( spreadStream->GetContractingDate() );
  while ( numeric::IsBefore(GetDoubleFrom(pPayment->first), dValuationTime) ) 
  { 
    dOldPaymentTime = GetDoubleFrom( pPayment->first );
    ++pPayment;
  }

  size_t nIdxT = 0;

  while ( pPayment != spreadStream->end() )
  {
    dPaymentTime = GetDoubleFrom( pPayment->first );

    double dTang = pPayment->second / (dPaymentTime - dOldPaymentTime);

    while ( numeric::IsBefore(m_pdTimes[nIdxT], dPaymentTime) )
    {
      // (1. + accrued) * R
      m_pdRecoveryValues[nIdxT]
          = parbond.GetRecoveryRate()
          * ( 1 + dTang * (m_pdTimes[nIdxT] - dOldPaymentTime));

      nIdxT++;
    }

    // Take the right limit at a spread payment date for backward problem
    m_pdRecoveryValues[nIdxT++] = parbond.GetRecoveryRate()
                                * ( 1 + pPayment->second);

    dOldPaymentTime = dPaymentTime;

    ++pPayment;
  }
}

void ParBondMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();
}

void ParBondMeshManager::SetupMeTimeOnly()
{
  BackwardMeshManager::SetupMe();

  ComputeRecoveryValues();
}


} // namespace pricing

} // namespace ito33

