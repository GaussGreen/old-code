/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cdsmeshmanager.cpp
// Purpose:     cds mesh manager for backward PDE problems
// Created:     2004/03/02
// RCS-ID:      $Id: cdsmeshmanager.cpp,v 1.35 2006/07/31 16:00:28 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

namespace ito33
{

  using finance::CashFlowStreamUniform;

namespace pricing
{


CDSMeshManager::CDSMeshManager(CDSParams &params, Model &model)
                             : BackwardMeshManagerFix(params, model),
                               m_cdsParams(params)    
{ 
}

void CDSMeshManager::ConstructSpaceMesh()
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

void CDSMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  CDS& cds = m_cdsParams.GetCDS();
  
  double dRecoveryRate = cds.GetRecoveryRate();

  double dRecoveryValue = 1. - dRecoveryRate;

  const CashFlowStreamUniform& spreadStream = cds.GetSpreadStream();

  CashFlowStreamUniform::const_iterator pPayment = spreadStream.begin();
 
  double dValuationTime = m_cdsParams.GetValuationTime();
  double dPaymentTime, dOldPaymentTime;

  // find out the first payment that's significant for us. 
  // What's the convention when there is a payment at valuation time?
  dOldPaymentTime = GetDoubleFrom( spreadStream.GetContractingDate() );
  while ( numeric::IsBefore(GetDoubleFrom(pPayment->first), dValuationTime) ) 
  { 
    dOldPaymentTime = GetDoubleFrom( pPayment->first );
    ++pPayment;
  }

  size_t nIdxT = 0;

  while ( pPayment != spreadStream.end() )
  {
    dPaymentTime = GetDoubleFrom( pPayment->first );

    double dTang = pPayment->second / (dPaymentTime - dOldPaymentTime);

    while ( numeric::IsBefore(m_pdTimes[nIdxT], dPaymentTime) )
    {
      // (1. + accrued) * (1. - R) - accured = 1. - R - R * accrued  
      m_pdRecoveryValues[nIdxT] = dRecoveryValue 
                                //- dRecoveryRate * dTang * (m_pdTimes[nIdxT] - dOldPaymentTime);
                                - dTang * (m_pdTimes[nIdxT] - dOldPaymentTime);

      nIdxT++;
    }

    // Take the right limit at a spread payment date for backward problem
    m_pdRecoveryValues[nIdxT++] = dRecoveryValue;

    dOldPaymentTime = dPaymentTime;

    ++pPayment;
  }
}

void CDSMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();
}

void CDSMeshManager::SetupMeTimeOnly()
{
  BackwardMeshManager::SetupMe();

  ComputeRecoveryValues();
}


} // namespace pricing

} // namespace ito33

