/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/edsmeshmanager.cpp
// Purpose:     EDS mesh manager for backward PDE problems
// Created:     2005/01/26
// RCS-ID:      $Id: edsmeshmanager.cpp,v 1.8 2006/04/21 09:25:45 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/array.h"
#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/genraf.h"

#include "ito33/numeric/mesh/barrierspacemesh.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/edsparams.h"
#include "ito33/pricing/edsmeshmanager.h"

namespace ito33
{

  using finance::CashFlowStreamUniform;
  
  using namespace numeric;
  using namespace numeric::mesh;

namespace pricing
{


EDSMeshManager::EDSMeshManager(EDSParams& params, Model& model)
                             : BackwardMeshManagerFix(params, model),
                               m_edsParams(params)    
{ 
}

// todo: all these magic numbers should be documented and shared with others
void EDSMeshManager::ConstructSpaceMesh()
{  
/*
  BarrierSpaceMesh barrierMesh;

  barrierMesh.SetLowerBarrier(m_edsParams.GetEDS().GetBarrier());

  double dSqrTotalVol = m_model.GetSquaredTotalVolForMesh
                        ( 
                          m_params.GetStoppingTime(),
                          m_params.GetSpotSharePrice() 
                        );
  
  double dDiff = m_params.GetDiffusionSize(dSqrTotalVol);
  if ( dDiff < 0.4 )
    dDiff = 0.4;

  double dConv = m_params.GetBackwardConvectionSize(dSqrTotalVol);
  if ( dConv < 0. )
    dConv = 0.;

  barrierMesh.SetDiffusionSize(dDiff);
  barrierMesh.SetConvectionSize(dConv);

  barrierMesh.SetSpecialPoint( m_edsParams.GetSpotSharePrice(), false);

  std::vector<double> pdMesh;
  barrierMesh.Build( m_params.GetNumParams()->GetNbSpaceSteps(),
                     m_edsParams.GetSpotSharePrice(),
                     pdMesh);


  m_nNbS = pdMesh.size();
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);

  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdS[nIdx] = pdMesh[nIdx];
    m_pdLogS[nIdx] = log( pdMesh[nIdx] );
  }

  return;
*/

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                            ( 
                              m_params.GetStoppingTime(),
                              m_params.GetSpotSharePrice() 
                            );

  double dDiffusion = m_params.GetMeshParams()->GetGridSpan()
                    * m_params.GetDiffusionSize(dSquaredTotalVol);
  if ( dDiffusion < 0.6)
    dDiffusion = 0.6;

  double dConvectionSize = m_params.GetBackwardConvectionSize(dSquaredTotalVol);
  if ( dConvectionSize < 0.)
    dConvectionSize = 0.;

  const numeric::MeshParams* pMeshParams = m_params.GetMeshParams();

  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();
    
  /*
  // Original parameters. Do not work well if the barrier is very small
  // compared to the spot price (eg. < 10%)
  double dA = log( m_edsParams.GetEDS().GetBarrier() );
  double dB = dA + dDiffusion + dConvectionSize;  
  double dMax = 10.0;
  double dDelta = dDiffusion / nNbPoints;
  */
  
  
  double dA = log( m_edsParams.GetEDS().GetBarrier() );
  double dB = log( m_edsParams.GetSpotSharePrice() * 2.0 );
  double dMax = dB + dDiffusion + dConvectionSize;
  if ( dMax < log( m_edsParams.GetSpotSharePrice() * 100.0 ) )
    dMax = log( m_edsParams.GetSpotSharePrice() * 100.0 );
  double dDelta = (dB - dA) / nNbPoints;
  

  double dDeltaMin = dDelta / pMeshParams->GetSpaceAccumulation();
  double dStretch = pMeshParams->GetSpaceStretch();
  double dDeltaMax = 0.5;

  // Use an upper bound for the normal delta, otherwise the result can be very
  // bad
  if ( dDelta > 0.15 )
    dDelta = 0.15;

  m_nNbS = GenRafSize( dB - dA, dDeltaMin, dDelta, dStretch);

  int iNbS2 = 1;
  
  if ( dB < dMax )
    iNbS2 = GenRafSize( dMax - dB, dDelta, dDeltaMax, dStretch);
  
  m_pdLogS = Array<double>(m_nNbS);
  
  int iNbS;

  // Generate the mesh from dA to dB and refine around the center which for now
  // is just the barrier
  numeric::mesh::GenRaf(dA, dB, dDeltaMin, dDelta, 
                        dStretch, 1, 
                        m_pdLogS.Get(), iNbS);

  

  // Generate the mesh from dB to a large number, beginning with previous 
  // mesh spacing, and growing to a max delta.
  //dDelta = m_pdLogS[iNbS-1] - m_pdLogS[iNbS-2];
  if ( dB < dMax )
    numeric::mesh::GenRaf(dB, dMax, dDelta, dDeltaMax, 
                          dStretch, 1, 
                          m_pdLogS.Get() + iNbS - 1, iNbS2);

  m_nNbS = iNbS + iNbS2 - 1;
  m_pdS = Array<double>(m_nNbS);

  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
  
}

void EDSMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  EDS& eds = m_edsParams.GetEDS();
  
  double dRecoveryRate = eds.GetRecoveryRate();

  double dRecoveryValue = 1. - dRecoveryRate;

  const CashFlowStreamUniform& spreadStream = eds.GetSpreadStream();

  finance::CashFlowStream::const_iterator pPayment = spreadStream.begin();
 
  double dValuationTime = m_params.GetValuationTime();
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
      // 1. - R - accrued  
      m_pdRecoveryValues[nIdxT] = dRecoveryValue 
                                - dTang * (m_pdTimes[nIdxT] - dOldPaymentTime);

      nIdxT++;
    }

    // Take the right limit at a spread payment date for backward problem
    m_pdRecoveryValues[nIdxT++] = dRecoveryValue;

    dOldPaymentTime = dPaymentTime;

    ++pPayment;
  }
}


} // namespace pricing

} // namespace ito33

