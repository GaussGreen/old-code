/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/optionmeshmanager.cpp
// Purpose:     option mesh manager for backward PDE problems
// Created:     2004/02/11
// RCS-ID:      $Id: optionmeshmanager.cpp,v 1.35 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/yieldcurve.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/mesh/optionspacemesh.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/optionparams.h"
#include "ito33/pricing/optionmeshmanager.h"

using namespace ito33::finance;
using namespace ito33::numeric;
using namespace ito33::numeric::mesh;
using namespace ito33::pricing;

// Added for testing uniform mesh refinement
// extern size_t nGlobalRefine;

OptionMeshManager::OptionMeshManager(OptionParams& params, Model& model)
                                   : BackwardMeshManagerFix(params, model),
                                     m_optionParams(params)  { }

void OptionMeshManager::ConstructSpaceMesh()
{
  if (m_params.GetMeshParams()->GetUniformSpaceGrid() == true)
  {
    ConstructUniformSpaceMesh();
    return;
  }

  const Option &option = m_optionParams.GetOption();
  
  // Note that for path dependent each path has its own
  // strike who is different from the actual option
  // contract. Allow to have mesh refined around the strike
  // of each path
  double dStrike = m_optionParams.GetStrike();

  ASSERT_MSG( dStrike >= 0, "Negative strike.");    

  OptionSpaceMesh osgGrid;

  osgGrid.SetOptionType( option.GetOptionType() );

  // Get the diffusion size from the used numerical model
  double dSquaredTotalVol = m_model.GetSquaredTotalVolForMesh
                            ( m_params.GetStoppingTime(), dStrike );

  osgGrid.SetDiffusionSize
          ( m_params.GetDiffusionSize(dSquaredTotalVol) );

  double dConvection = m_model.GetConvection
                       ( m_params.GetStoppingTime(), dStrike );

  // Compute the convection size for this problem
  osgGrid.SetConvectionSize
          ( m_params.GetBackwardConvectionSize(dConvection) );

  double dExtra = m_params.GetSpotSharePrice() / dStrike;
  osgGrid.SetExtraHorizon( dExtra );

  // Generate a space mesh centered at zero
  std::vector<double> vecGrid;
  
  // Number of requested points, get it from m_pNumParams 
  size_t nNbPoints = m_params.GetNumParams()->GetNbSpaceSteps();

  // Filter out case of strike being near zero (which happens, for example,
  // in put-call parity calculation)
  const double dStrikeCutoff = 1.e-10;
  if ( dStrike < dStrikeCutoff)
  {
    osgGrid.Build(nNbPoints,
                  log( m_params.GetSpotSharePrice() ),
                  vecGrid);
  }
  else
  {
    osgGrid.Build(nNbPoints,
                  log( m_params.GetSpotSharePrice() / dStrike ),
                  vecGrid);
  }


  m_nNbS = vecGrid.size();
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS = Array<double>(m_nNbS);

  // Re-center the space mesh on the strike
  double dLogStrike;
  if ( dStrike < dStrikeCutoff )
    dLogStrike = 0.0;
  else
    dLogStrike = log( dStrike );
    
  
  for (size_t nIdx = 0; nIdx < m_nNbS; nIdx++)
  {
    m_pdLogS[nIdx] = vecGrid[nIdx] + dLogStrike;
    m_pdS[nIdx] = exp(m_pdLogS[nIdx]);
  }

  /*
  // This section added for testing.
  // Refine the mesh
  for (size_t nTmp = 0; nTmp < nGlobalRefine; nTmp++)
  {
    Array<double> tmpMesh(m_nNbS*2-1);
    for (size_t nIdx = 0; nIdx < m_nNbS-1; nIdx++)
    {
      tmpMesh[nIdx*2] = m_pdS[nIdx];
      tmpMesh[nIdx*2 + 1] = (m_pdS[nIdx] + m_pdS[nIdx+1])/2.0;;
    }
    tmpMesh[m_nNbS*2-2] = m_pdS[m_nNbS-1];

    m_nNbS = m_nNbS*2-1;
    m_pdS = tmpMesh;
  }
  */
}


void OptionMeshManager::ConstructUniformSpaceMesh()
{

  // Force the points 0.0 and the strike.  Continue until approximately
  // three times the strike.  This means about 1/3 of the points will be
  // between 0 and the strike.

  m_nNbS = m_params.GetNumParams()->GetNbSpaceSteps();
  size_t nNbTmp = m_nNbS / 3;
  double dDeltaT = m_optionParams.GetStrike() / nNbTmp;
  
  // Construct the actual mesh
  m_pdS = Array<double>(m_nNbS);
  m_pdLogS = Array<double>(m_nNbS);
  m_pdS[0] = 0.0;
  m_pdLogS[0] = log(dDeltaT / 100.0);
  for (size_t nIdx = 1; nIdx < m_nNbS; nIdx++)
  {    
    m_pdS[nIdx] = nIdx * dDeltaT;
    m_pdLogS[nIdx] = log(m_pdS[nIdx]);
  }

}

void OptionMeshManager::ComputeRecoveryValues()
{
  // Pre-compute the recovery values
  m_pdRecoveryValues = Array<double>(m_nNbTimes);

  Option &option = m_optionParams.GetOption();
  
  size_t nIdxT;

  switch ( option.GetOptionType() )
  {
  case Option_Other:

  case Option_Digital:  
  
  case Option_Call: 
    
    for (nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
      m_pdRecoveryValues[nIdxT] = 0.;

    break;

  case Option_Put: // K e^{ - \int_t^T r(s) ds} 
  {  
    double dStrike = m_optionParams.GetStrike();
    
    if (option.GetExerciseType() == ExerciseType_American)
      for (nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
        m_pdRecoveryValues[nIdxT] = dStrike;
    else
    {
      Array<double> pdCF(m_nNbTimes);
      
      m_optionParams.GetYieldCurve()->GetCompoundFactor
        (m_pdTimes.Get(), pdCF.Get(), m_nNbTimes);

      double dTmp = 1. / pdCF[m_nNbTimes - 1];
      for (nIdxT = 0; nIdxT < m_nNbTimes; nIdxT++)
        m_pdRecoveryValues[nIdxT] = dStrike * pdCF[nIdxT] * dTmp;
    }

    break;
  }
  
  default:  

    FAIL("Option pricer doesn't handle current option type.");    
  }
}

void OptionMeshManager::SetupMe()
{
  BackwardMeshManagerFix::SetupMe();

  ComputeRecoveryValues();
}
