/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/varianceswap/varianceswappricer.cpp
// Purpose:     HG option pricer class
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswappricer.cpp,v 1.15 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"

#include "ito33/numeric/predicatetime.h"

#include "hg/model.h"
#include "hg/pathdepstructure.h"
#include "hg/pathdeppricer.h"
#include "hg/varianceswapstepper.h"
#include "hg/varianceswapnumoutput.h"
#include "hg/varianceswappricer.h"

namespace ito33
{

namespace hg
{

using namespace pricing;

VarianceSwapPricer::VarianceSwapPricer(VarianceSwapParams& params, 
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
  // The HG model is always homogeneous.  Just need to check the params 
  // and the similarity reduction computational flags setting
  if ( m_flags.GetUseSimilarityReductions() == false )
    return false;

  return m_params.HasSimilarityReduction();
}


ito33::AutoPtr<VarianceSwapNumOutput> VarianceSwapPricer::Price()
{
  
  bool bHasSimilarityReduction = HasSimilarityReduction();

  // Create a grid for the average of the squared returns  
  std::vector<double> 
    pdAvgSqrReturnGrid( m_params.ConstructAvgSqrReturnGrid() );

  // Create a grid for the previous spot prices
  std::vector<double>
    pdPreviousSpotGrid( m_params.ConstructPreviousSpotGrid
                                 (m_model, bHasSimilarityReduction) );

  // Create a master 'S' grid for all the paths
  std::vector<double> pdMasterSGrid( m_params.ConstructMasterGrid(m_model) );

  // Clone the params and create the paths
  std::vector< AutoPtr<VarianceSwapParams> > pParams;
  pParams = m_params.ConstructParams(pdAvgSqrReturnGrid, pdPreviousSpotGrid);

  // Find the path to save
  size_t nIdxPathToSave = m_params.GetPathToSave(pdAvgSqrReturnGrid, 
                                                 pdPreviousSpotGrid);

  // Create the list of events   
  std::list< shared_ptr<PathDepEvent> > pathDepEvents;
  pathDepEvents = m_params.ConstructEvents(bHasSimilarityReduction);

  // Do the appropriate mesh construction
  size_t nNbPaths = pdAvgSqrReturnGrid.size() * pdPreviousSpotGrid.size();
  std::vector< shared_ptr<VarianceSwapMeshManager> > 
    ppMeshes(nNbPaths);
  
  size_t nNbMasterS = pdMasterSGrid.size();
  std::vector<double> pdTmpSpaceMesh( nNbMasterS );
  size_t nNbPrevious = pdPreviousSpotGrid.size();
  double dPreviousMaster = m_params.GetVarianceSwap().GetPreviousSharePrice();

  for (size_t nIdxPath = 0; nIdxPath < nNbPaths; nIdxPath++)
  {
    // Get the previous spot of this path
    size_t nIdxPrevious = nIdxPath % nNbPrevious;
    double dPreviousSpot = pdPreviousSpotGrid[nIdxPrevious];

    // The master grid is centered at the previous share price of the
    // original problem.  Scale the grid using the previous spot of this path.
    // This makes the path accurate at the previous spot of this path, which
    // desireable when interpolating along the diagonal
    double dScale = dPreviousSpot / dPreviousMaster;

    for (size_t nIdx = 0; nIdx < nNbMasterS; nIdx++)
      pdTmpSpaceMesh[nIdx] = dScale * pdMasterSGrid[nIdx];

    ppMeshes[nIdxPath] = make_ptr( new VarianceSwapMeshManager
                                       ( *(pParams[nIdxPath]), m_model,
                                        pdTmpSpaceMesh ) );
  }

  // If forward starting, save output surfaces, theta, etc.
  bool bSaveOutputAfterEvents = m_params.IsForwardStarting();

  // Create the structure
  std::vector< std::vector<double> > ppdGrids(2);
  ppdGrids[0] = pdAvgSqrReturnGrid;
  ppdGrids[1] = pdPreviousSpotGrid;

  PathDepStructure<VarianceSwapInstData, 
                   VarianceSwapMeshManager, 
                   VarianceSwapParams, 
                   VarianceSwapNumOutput, 
                   VarianceSwapStepper>
    pathDepStruct(ppMeshes, pParams, m_model, m_flags, 
                  pathDepEvents, ppdGrids, nIdxPathToSave);
  
  pathDepStruct.PrepareForTimestepping();

  // Solve
  PathDepPricer<VarianceSwapInstData, 
                VarianceSwapMeshManager, 
                VarianceSwapParams, 
                VarianceSwapNumOutput, 
                VarianceSwapStepper>
    pathDepPricer(pathDepStruct);

  pathDepPricer.Price(pathDepEvents, bSaveOutputAfterEvents); 

  AutoPtr<VarianceSwapNumOutput> pNumOutput = pathDepStruct.GetNumOutput();
  
  // If a conditional swap, price the fixed leg separately.  Subtract the
  // fixed leg from the floating leg which was priced above
  if ( m_params.GetVarianceSwap().GetSwapPayoffType() 
         == finance::SwapPayoff_Conditional )
  {
    AutoPtr<VarianceSwapNumOutput> pFixedNumOutput = PriceFixedLeg();

    pNumOutput->SubractNumOutput( *pFixedNumOutput);
  }

  return pNumOutput;
}


AutoPtr<VarianceSwapNumOutput> VarianceSwapPricer::PriceFixedLeg()
{

  // return object
  AutoPtr<VarianceSwapNumOutput> pNumOutput;

  // handle conditional swaps
  finance::SwapPayoffType 
    swapPayoff = m_params.GetVarianceSwap().GetSwapPayoffType();

  ASSERT_MSG( swapPayoff == finance::SwapPayoff_Conditional, 
              "Can only price fixed leg for conditional swap");

  if ( swapPayoff == finance::SwapPayoff_Conditional )
  {
    // State variable Y represents a count of the number of times the
    // spot is within the corridor at sampling return dates

    // Create a master 'S' grid for all the paths
    std::vector<double> pdMasterSGrid( m_params.ConstructMasterGrid(m_model) );

    // Construct grid for Y state variable. 
    std::vector<double> pdYGrid( m_params.ConstructConditionalGrid() );

    // Clone the params and create the paths
    std::vector< AutoPtr<VarianceSwapParams> > pParams;
    pParams = m_params.ConstructConditionalParams(pdYGrid);

    // Find the path to save
    size_t nIdxPathToSave = m_params.GetConditionalPathToSave(pdYGrid);

    // Create the S meshes using the master grid
    size_t nNbPaths = pParams.size();
    std::vector< shared_ptr<VarianceSwapMeshManager> > ppMeshes(nNbPaths);

    for (size_t nIdxPath = 0; nIdxPath < nNbPaths; nIdxPath++)
    {
      ppMeshes[nIdxPath] = make_ptr( new VarianceSwapMeshManager 
                                         ( *(pParams[nIdxPath]), m_model, 
                                           pdMasterSGrid ) );
      ppMeshes[nIdxPath]->SetIsConditionalFixed(true);
    }

    // Create the list of events   
    std::list< shared_ptr<PathDepEvent> > pathDepEvents;
    pathDepEvents = m_params.ConstructConditionalEvents();

    // If forward starting, save output surfaces, theta, etc.
    bool bSaveOutputAfterEvents = m_params.IsForwardStarting();

    // Create the structure
    std::vector< std::vector<double> > ppdGrids(1);
    ppdGrids[0] = pdYGrid;    

    PathDepStructure<VarianceSwapInstData, 
                     VarianceSwapMeshManager, 
                     VarianceSwapParams, 
                     VarianceSwapNumOutput, 
                     VarianceSwapStepper>
    pathDepStruct(ppMeshes, pParams, m_model, m_flags, 
                  pathDepEvents, ppdGrids, nIdxPathToSave);
  
    pathDepStruct.PrepareForTimestepping();

    // Solve
    PathDepPricer<VarianceSwapInstData, 
                  VarianceSwapMeshManager, 
                  VarianceSwapParams, 
                  VarianceSwapNumOutput, 
                  VarianceSwapStepper>
      pathDepPricer(pathDepStruct);

    pathDepPricer.Price(pathDepEvents, bSaveOutputAfterEvents); 

    pNumOutput = pathDepStruct.GetNumOutput();
  }

  // Avoid warning
  return pNumOutput;
}


} // namespace hg

} // namespace ito33
