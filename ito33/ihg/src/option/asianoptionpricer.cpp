/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/option/asianoptionpricer.cpp
// Purpose:     asian option pricer class
// Author:      ITO 33 Canada
// Created:     April 7, 2005
// RCS-ID:      $Id: asianoptionpricer.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/asianoptionparams.h"
#include "ito33/pricing/asianoption.h"
#include "ito33/pricing/pathdeppricer.h"
#include "ito33/pricing/asianoptionmeshmanager.h"

#include "ito33/numeric/predicatedouble.h"

#include "ihg/model.h"
#include "ihg/optionnumoutput.h"
#include "ihg/asianoptionpricer.h"
#include "ihg/optionpathdepstructure.h"

namespace ito33
{

namespace ihg
{
  
typedef pricing::Engine< pricing::AsianOptionParams, pricing::OptionMeshManager,
                         OptionInstData, OptionStepper, OptionNumOutput
                       > AsianOptionEngine; 


AsianOptionPricer::AsianOptionPricer(pricing::AsianOptionParams& params, 
                           Model& model,
                           const finance::ComputationalFlags& flags)
                         : m_params(params), 
                           m_model(model),
                           m_flags(flags),
                           m_pNumOutput(new OptionNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}

bool AsianOptionPricer::IsPathDependent()
{
  // It is path dependent if at least one of the following is true
  // - has a cash or pseudocash dividend
  // - stock dependent vol and hazard rate

  if (  m_params.IsPathDependent() )
    return true;

  // hazard rate must be independent of S
  if ( !m_model.IsHazardRateTimeOnly() )
      return true;

  // volatility must be independent of S
  if ( !m_model.IsVolatilityTimeOnly() )
    return true;

  return false;
}

ito33::AutoPtr<OptionNumOutput> AsianOptionPricer::Price()
{
  bool bHasSimilarityReduction = true;

  if ( IsPathDependent() )
    bHasSimilarityReduction = false;
    
  // Averaging grid
  std::vector<double> pdGridY;

  // Create the grid
  m_params.ConstructPathDepGrid(pdGridY, m_model, bHasSimilarityReduction);

  // Clone the params and create the paths
  std::vector< AutoPtr<pricing::OptionParams> > pAsianOptionParams;
  pAsianOptionParams = m_params.ConstructParam(pdGridY);

  // Construct event list
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
  pathDepEvents = m_params.ConstructEvents(bHasSimilarityReduction);

  // Get the path to save
  size_t nPathToSave = m_params.GetPathToSave(pdGridY);
 
  size_t nNbPaths = pdGridY.size();
  size_t nIdx;
  
  // Do the appropriate mesh construction
  Array< shared_ptr<pricing::OptionMeshManager> > ppMeshes(nNbPaths);

  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    pricing::AsianOptionParams* pAsianParams 
      = (pricing::AsianOptionParams*) pAsianOptionParams[nIdx].get();
    
    ppMeshes[nIdx] = make_ptr( new pricing::AsianOptionMeshManager
                                   (*pAsianParams, m_model) );
  }

  //Create the structure
  ihg::OptionPathDepStructure optionPath( ppMeshes, pdGridY, 
    pAsianOptionParams, m_model, m_flags, pathDepEvents, nPathToSave);
  
  // Solve
  pricing::PathDepPricer pathDepPricer;
  pathDepPricer.Price(optionPath, pathDepEvents); 

  AutoPtr<OptionNumOutput> pNumOutput = optionPath.GetOutput();

  return pNumOutput;
}


} // namespace ihg

} // namespace ito33
