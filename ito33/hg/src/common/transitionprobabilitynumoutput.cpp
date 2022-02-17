/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/transitionprobabilitynumoutput.cpp
// Purpose:     implementation of HG numoutput class for transition probability
// Created:     2006/03/31
// RCS-ID:      $Id: transitionprobabilitynumoutput.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/hg/transitionprobabilityoutput.h"

#include "hg/numoutputanalytical.h"
#include "hg/transitionprobabilitynumoutput.h"
#include "hg/transitionprobabilityinstdata.h"

// implement the AutoPtrDeleter for TransitionProbabilityNumOutput
namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(hg::TransitionProbabilityNumOutput);
}

namespace ito33
{

namespace hg
{
  using namespace numeric;


shared_ptr<TransitionProbabilityOutput>
TransitionProbabilityNumOutput::Solve
(TransitionProbabilityInstData& instdata)
{
  m_nNbRegimes = instdata.m_nNbRegimes;
  
  // Get the space mesh. Also sets m_nNbS.
  m_pdS = instdata.GetSpaceMesh(m_nNbS);
  m_nNbX = m_nNbS * m_nNbRegimes;

  shared_ptr<TransitionProbabilityOutput>
    pOutput(new TransitionProbabilityOutput);

  std::vector<double> returns(m_nNbS);

  for (size_t n = 0; n < m_nNbS; n++)
    returns[n] = m_pdS[n] / instdata.GetInitialSpot();

  pOutput->SetSpaceMesh(returns);

  std::vector< shared_ptr<finance::ModelOutput> > pOutputs;

  Array<double> pdF(m_nNbX);
  std::vector<double> pdProbas(m_nNbX);
  std::vector<double> pdProbaToDefaults(m_nNbRegimes);

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  { 
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      pdF[nIdx] = 0.;

    pdF[instdata.GetIndexSpot() + nIdxR * m_nNbS] = 1;

    SolveTransposeSystem(instdata, pdF.Get(), &pdProbas[0]);

    // Modeloutput for now is just a place holder for numoutput
    // since we need usually values at all regimes
    shared_ptr<finance::ModelOutput> pMO(new finance::ModelOutput);

    // We are only interested at the values at the valuation date
    // so just use numoutput analytical to simplify the code
    shared_ptr<NumOutputAnalytical>
      pNumOutput( new NumOutputAnalytical(m_params) );
      
    pNumOutput->SetFinalValues(returns, pdProbas); 

    pMO->SetNumOutput(pNumOutput);
 
    pOutputs.push_back(pMO);   
    
    double dProbaToDefault = 1.;
    for (size_t nIdx = 0; nIdx < m_nNbX; nIdx++)
      dProbaToDefault -= pdProbas[nIdx];

    pdProbaToDefaults[nIdxR] = dProbaToDefault;
  }

  pOutput->SetProbabilities(pOutputs);
  
  pOutput->SetProbabilitiesToDefault(pdProbaToDefaults);

  return pOutput;
} 


} // namespace hg

} // namespace ito33
