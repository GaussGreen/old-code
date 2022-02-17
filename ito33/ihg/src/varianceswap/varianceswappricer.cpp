/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/varianceswap/varianceswappricer.h
// Purpose:     variance swap pricer class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswappricer.cpp,v 1.12 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/error.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/pathdeppricer.h"

#include "ihg/model.h"
#include "ihg/varianceswappathdepstructure.h"
#include "ihg/varianceswappricer.h"
#include "ihg/varianceswapnumoutput.h"

extern const ito33::finance::Error ITO33_VARIANCESWAP_IHG_CONDITIONAL;

namespace ito33
{

namespace ihg
{
  
VarianceSwapPricer::VarianceSwapPricer(pricing::VarianceSwapParams& params, 
                           Model& model,
                           const finance::ComputationalFlags& flags)
                         : m_params(params), 
                           m_model(model),
                           m_flags(flags)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}


bool VarianceSwapPricer::HasSimilarityReduction()
{
  // It is path dependent if at least one of the following is true
  // - has a cash or pseudocash dividend
  // - stock dependent vol and hazard rate
  // Check params for the first and the model for the latter
  if ( !m_params.HasSimilarityReduction() ) 
    return false;

  // hazard rate must be independent of S
  if ( !m_model.IsHazardRateTimeOnly() )
    return false;

  // volatility must be independent of S
  if ( !m_model.IsVolatilityTimeOnly() )
    return false;

  // check the similarity reduction computational flags setting
  if ( m_flags.GetUseSimilarityReductions() == false )
    return false;

  return true;
}

ito33::AutoPtr<VarianceSwapNumOutput> VarianceSwapPricer::Price()
{
  // Conditional payoff is not currently supported
  CHECK_COND( m_params.GetVarianceSwap().GetSwapPayoffType() 
              != finance::SwapPayoff_Conditional, 
              ITO33_VARIANCESWAP_IHG_CONDITIONAL);

  bool bHasSimilarityReduction = HasSimilarityReduction();

  // Create a grid for the average of the squared returns
  std::vector<double> 
    pdAvgSqrReturnGrid( m_params.ConstructAvgSqrReturnGrid() );

  // Create a grid for the previous spot prices
  std::vector<double>
    pdPreviousSpotGrid( m_params.ConstructPreviousSpotGrid
                                 (m_model, bHasSimilarityReduction) );

  // Create a master 'S' grid for all the paths
  std::vector<double> pdMasterGrid( m_params.ConstructMasterGrid(m_model) );

  // Clone the params and create the paths
  std::vector< AutoPtr<pricing::VarianceSwapParams> > pParams;
  pParams = m_params.ConstructParams(pdAvgSqrReturnGrid, pdPreviousSpotGrid);

  // Find the path to save
  size_t nIdxPathToSave = m_params.GetPathToSave(pdAvgSqrReturnGrid, 
                                                 pdPreviousSpotGrid);

  // Create the list of events   
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
  pathDepEvents = m_params.ConstructEvents(bHasSimilarityReduction);

  // If forward starting, save output surfaces, theta, etc.
  bool bSaveOutputAfterEvents = m_params.IsForwardStarting();

  // Create the structure  
  ihg::VarianceSwapPathDepStructure 
    pathDepStruct( pdAvgSqrReturnGrid, pdPreviousSpotGrid, pdMasterGrid,
                   pParams, m_model, m_flags, pathDepEvents, nIdxPathToSave);
  
  // Solve
  pricing::PathDepPricer pathDepPricer;
  pathDepPricer.Price(pathDepStruct, pathDepEvents, bSaveOutputAfterEvents); 

  AutoPtr<VarianceSwapNumOutput> pNumOutput = pathDepStruct.GetOutput();

  return pNumOutput;
}


} // namespace ihg

} // namespace ito33
