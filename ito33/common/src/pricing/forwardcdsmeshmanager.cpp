/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/forwardcdsmeshmanager.cpp
// Purpose:     cds mesh manager for forward PDE problems
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: forwardcdsmeshmanager.cpp,v 1.23 2006/08/21 16:10:14 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/optionliketype.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/specialtimes.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/forwardcdsparams.h"
#include "ito33/pricing/forwardcdsmeshmanager.h"

using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

void ForwardCDSMeshManager::ConstructSpaceMesh()
{
 
  if (m_params.GetMeshParams()->GetUniformSpaceGrid() == true)
  {
    ConstructUniformSpaceMesh();
    return;
  }

  OptionSpaceMesh osgGrid;

  osgGrid.SetOptionType( Option_Call );

  // Get the diffusion size from the used numerical model
  double 
    dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                       ( 
                         m_forwardCDSParams.GetStoppingTime(),
                         m_forwardCDSParams.GetSpotSharePrice() 
                       );

  osgGrid.SetDiffusionSize
          ( m_forwardCDSParams.GetDiffusionSize(dSquaredTotalVol) );

  // Compute the convection size for this problem
  osgGrid.SetConvectionSize
          ( m_forwardCDSParams.GetForwardConvectionSize(dSquaredTotalVol) );

  // Generate a space mesh centered at zero
  std::vector<double> vecGrid;
  
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_forwardCDSParams.GetNumParams()->GetNbSpaceSteps();

  osgGrid.Build(nNbPoints,
                1., // log( m_params.GetSpotSharePrice() / m_params.GetSpotSharePrice() ),
                vecGrid);

  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  // Re-center the space mesh on the spot
  double dLogSpot = log( m_forwardCDSParams.GetSpotSharePrice() );
  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx] + dLogSpot;
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }
}

void ForwardCDSMeshManager::SetupMe()
{
  ForwardMeshManagerFix::SetupMe();


  // Pre-compute the recovery values and spreads
  m_pdAccruedFraction = Array<double>(m_nNbTimes);
  m_pdSpreads = Array<double>(m_nNbTimes);

  SpecialTimes paymentDateList;

  const finance::CashFlowStream& 
    spreadStream = m_forwardCDSParams.GetCDSes().GetSpreadStream();

  finance::CashFlowStream::const_iterator 
    paymentDates = spreadStream.begin();

  double dValuationTime = m_forwardCDSParams.GetValuationTime();

  double dPaymentTime, dOldPaymentTime;

  size_t nIdxT = 0;

  // find out the first payment date that's significant for us. 
  // What's the convention to be used when pricing date is a payment date?
  dOldPaymentTime = m_forwardCDSParams.GetValuationTime();
  while (numeric::IsBefore(GetDoubleFrom(paymentDates->first), dValuationTime))  
  { 
    dOldPaymentTime = GetDoubleFrom(paymentDates->first);
    ++paymentDates;
  }

  while ( paymentDates != spreadStream.end() )
  {
    dPaymentTime = GetDoubleFrom(paymentDates->first);

    double dInversePeriod = 1.0 / (dPaymentTime - dOldPaymentTime);

    while ( numeric::IsBefore(m_pdTimes[nIdxT], dPaymentTime) )
    {
      m_pdAccruedFraction[nIdxT] = dInversePeriod
                                 * (m_pdTimes[nIdxT] - dOldPaymentTime);
      m_pdSpreads[nIdxT] = paymentDates->second;
      nIdxT++;
    }

    // Set the accrued fraction to 1.0 at the spread dates. This is assumed
    // by the stepper to properly handle the discontinuity
    m_pdAccruedFraction[nIdxT] = 1.0;
    m_pdSpreads[nIdxT] = paymentDates->second;
    nIdxT++;
    
    dOldPaymentTime = dPaymentTime;

    ++paymentDates;
  }

}

void ForwardCDSMeshManager::ConstructUniformSpaceMesh()
{

  // Force the points 0.0 and the spot.  Continue until approximately
  // ten times the strike.  This means about 1/5 of the points will be
  // between 0 and the strike.

  m_nNbS = m_params.GetNumParams()->GetNbSpaceSteps();
  size_t nNbTmp = m_nNbS / 10;
  double dDeltaT = m_params.GetSpotSharePrice() / nNbTmp;
  
  // Construct the actual mesh
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS[0] = dDeltaT/10000.0;
  //m_pdS[0] = 1.e-16;
  m_pdLogS[0] = log( m_pdS[0] );
  size_t nIdx;
  for (nIdx = 1; nIdx < m_nNbS - 3; nIdx++)
  {    
    m_pdS[nIdx] = nIdx * dDeltaT;
    m_pdLogS[nIdx] = log(m_pdS[nIdx]);
  }

  m_pdS[m_nNbS - 3] = m_pdS[m_nNbS - 4] + dDeltaT*10;
  m_pdLogS[m_nNbS - 3] = log( m_pdS[m_nNbS-3] );

  m_pdS[m_nNbS - 2] = m_pdS[m_nNbS - 3] + dDeltaT*100;
  m_pdLogS[m_nNbS - 2] = log( m_pdS[m_nNbS-2] );

  m_pdS[m_nNbS - 1] = m_pdS[m_nNbS - 2] + dDeltaT*1000;
  m_pdLogS[m_nNbS - 1] = log( m_pdS[m_nNbS-1] );

}
