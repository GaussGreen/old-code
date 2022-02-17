/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/asianoption/asianoptionpricer.cpp
// Purpose:     asian option pricer class
// Created:     2006/03/01
// RCS-ID:      $Id: asianoptionpricer.cpp,v 1.7 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/binarysearch.h"
#include "ito33/array.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/asianoptionparams.h"
#include "ito33/pricing/asianoptionmeshmanager.h"
#include "ito33/pricing/optionparams.h"

#include "hg/asianoptionpricer.h"
#include "hg/model.h"
#include "hg/pathdeppricer.h"
#include "hg/pathdepstructure.h"
#include "hg/optionnumoutput.h"
#include "hg/optionstepper.h"
#include "hg/optioninstdata.h"

namespace ito33
{

namespace hg
{
  
typedef OptionStepper<StepperFix> MyStepper;

typedef pricing::Engine< pricing::AsianOptionParams, pricing::OptionMeshManager,
                         OptionInstData, MyStepper, OptionNumOutput
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


ito33::AutoPtr<OptionNumOutput> AsianOptionPricer::Price()
{
  bool bHasSimilarityReduction = true;

  if ( m_params.IsPathDependent() )
    bHasSimilarityReduction = false;

  std::vector<double> pdGridY;

  // Construct Averaging grid
  m_params.ConstructPathDepGrid(pdGridY, m_model, bHasSimilarityReduction);
 
  // Clone the params and create the paths
  std::vector< AutoPtr<pricing::OptionParams> > pAsianOptionParams;
  pAsianOptionParams = m_params.ConstructParam(pdGridY);

  // Construct event list
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
  pathDepEvents = m_params.ConstructEvents(bHasSimilarityReduction);

  // Get Path to save
  size_t nPathToSave = m_params.GetPathToSave(pdGridY);

  size_t nNbPaths = pdGridY.size();
  size_t nIdx;
  
  // Do the appropriate mesh construction
  std::vector< shared_ptr<pricing::OptionMeshManager> > ppMeshes(nNbPaths);

  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    pricing::AsianOptionParams* pAsianParams 
      = (pricing::AsianOptionParams*) pAsianOptionParams[nIdx].get();
    
    ppMeshes[nIdx] = make_ptr( new pricing::AsianOptionMeshManager
                                   (*pAsianParams, m_model) );
  }

  std::vector< std::vector<double> > ppdGrids(1);
  ppdGrids[0] = pdGridY;

  PathDepStructure<OptionInstData, 
                   pricing::OptionMeshManager, pricing::OptionParams,
                   OptionNumOutput, MyStepper>
    pathDepStruct(ppMeshes, pAsianOptionParams, m_model, m_flags,
                  pathDepEvents, ppdGrids, nPathToSave);

  pathDepStruct.PrepareForTimestepping();

  PathDepPricer<OptionInstData, pricing::OptionMeshManager, pricing::OptionParams, 
                OptionNumOutput, MyStepper>
    pathDepPricer(pathDepStruct);

  pathDepPricer.Price(pathDepEvents); 

  AutoPtr<OptionNumOutput> pNumOutput = pathDepStruct.GetNumOutput();

  return pNumOutput;
}


} // namespace hg

} // namespace ito33
