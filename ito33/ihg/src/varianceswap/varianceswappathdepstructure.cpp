/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/varianceswap/varianceswappathdepstructure.cpp
// Purpose:     path dependent structure implementation for variance swaps
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswappathdepstructure.cpp,v 1.10 2006/08/20 09:31:05 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/pricing/varianceswapparams.h"
#include "ito33/pricing/varianceswapmeshmanager.h"
#include "ito33/pricing/varianceswap.h"
#include "ito33/pricing/pathdepevent.h"

#include "ihg/varianceswapinstdata.h"
#include "ihg/varianceswapstepper.h"
#include "ihg/varianceswapnumoutput.h"
#include "ihg/model.h"
#include "ihg/varianceswappathdepstructure.h"

namespace ito33
{
  using namespace pricing;
  using namespace numeric::mesh;

namespace ihg
{

VarianceSwapPathDepStructure::VarianceSwapPathDepStructure
                    (  
                       const std::vector<double>& pdAvgSqrReturnGrid, 
                       const std::vector<double>& pdPreviousSpotGrid,
                       const std::vector<double>& pdMasterSGrid,
                       std::vector< AutoPtr<VarianceSwapParams> >& ppParams,
                       Model& model,
                       const finance::ComputationalFlags& flags,
                       std::list< shared_ptr<PathDepEvent> >& pathDepEvents,
                       size_t nPathToSave
                    )
                    : m_flags(flags)

{
  size_t nNbReturns = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousSpotGrid.size();
  size_t nNbPaths = nNbReturns * nNbPrevious;

  // Create the meshes. The master grid can be scaled to make the
  // space grids for each path
  Array< shared_ptr<VarianceSwapMeshManager> > ppMeshes(nNbPaths);

  // Keep the masterSGrid concept since it should
  // be quicker than building all the S meshes, and can also guarantee that
  // all meshes are the same (this could perhaps exploited during event
  // interpolation)
  size_t nNbMasterS = pdMasterSGrid.size();
  std::vector<double> pdTmpSpaceMesh( nNbMasterS );
  double dPreviousMaster = 
    ppParams[nPathToSave]->GetVarianceSwap().GetPreviousSharePrice();

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

    ppMeshes[nIdxPath] = make_ptr
        ( new VarianceSwapMeshManager
              ( *(ppParams[nIdxPath]), model, pdTmpSpaceMesh ) );
  }

  InitPathDepStructure(ppMeshes, pdAvgSqrReturnGrid, pdPreviousSpotGrid, ppParams, 
    model, pathDepEvents, nPathToSave);


}
  
void VarianceSwapPathDepStructure::InitPathDepStructure(
    Array< shared_ptr<VarianceSwapMeshManager> > &ppMeshes,
    const std::vector<double>& pdAvgSqrReturnGrid, 
    const std::vector<double>& pdPreviousSpotGrid,
    std::vector< AutoPtr<VarianceSwapParams> >& ppParams,
    Model& model,
    std::list< shared_ptr<PathDepEvent> >& pathDepEvents,
    size_t nPathToSave)
{

  // Path to save for output
  m_nPathToSave = nPathToSave;

  // Storage for the 1-D paths
  size_t nNbReturns = pdAvgSqrReturnGrid.size();
  size_t nNbPrevious = pdPreviousSpotGrid.size();
  size_t nNbPaths = nNbReturns * nNbPrevious;

  m_path.resize(nNbPaths);
  m_pbIsActive.resize(nNbPaths);
 
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    // initialize params here, but release below since they are still needed
    ppParams[nIdx]->Init();
  }

  // If any continuous events are present, we need each path to have the 
  // same mesh (technically, the paths need to have the same mesh during
  // the continuous window, but for convenience, make the entire time
  // mesh identical)
  bool bHasContinuousEvent = false;
  std::list< shared_ptr<PathDepEvent> >::iterator iterEvents;
  for (iterEvents = pathDepEvents.begin();
       iterEvents != pathDepEvents.end();
       ++iterEvents)
  {
    if ( (*iterEvents)->IsContinuous() )
    {
      bHasContinuousEvent = true;
      break;
    }
  }

  if ( bHasContinuousEvent )
  {
    ConstructIdenticalTimeMeshes( ppMeshes.Get(), ppParams, model, 
                                  pathDepEvents, nNbPaths );
  }
  
  // Setup and save the meshes in the path dep structure
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    // In case individual meshes are constructed, each param object
    // must have access to the path dep events. If identical meshes
    // were constructed above, the events were included, and setting
    // them in params does no harm
    ppParams[nIdx]->SetPathDepEvents(pathDepEvents);

    ppMeshes[nIdx]->SetupMe();
  
    m_path[nIdx].meshes = ppMeshes[nIdx];
  }

  // Create the rest of the 1-D problems
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {     
    // By default, leave all paths on
    m_pbIsActive[nIdx] = true;

    // Create the instadata and stepper objects.  Params and meshes were 
    // created before. In each case, create a temporary object with the 
    // exact type (for later use), but store it in the data structure as a 
    // base type.
    shared_ptr<VarianceSwapInstData> 
      instDataVS(new VarianceSwapInstData
                         ( *(ppParams[nIdx]), model, *(ppMeshes[nIdx]) ) );

    // Currently, path dep events are not applied to greek data.  Compute
    // vega by pertubation
    instDataVS->m_bComputeVega = false;

    instDataVS->m_bComputeFugit = m_flags.GetComputeFugit();

    m_path[nIdx].instdata = instDataVS;
    m_path[nIdx].instdata->Init();

    m_path[nIdx].stepper = make_ptr( new VarianceSwapStepper
                                         ( *instDataVS, m_flags ) );
    m_path[nIdx].stepper->Init();

    // Set initial conditions
    m_path[nIdx].meshes->SetInitialState();
    m_path[nIdx].instdata->SetInitialValue();

    // Release control of params memory to this object
    m_path[nIdx].params = make_ptr( ppParams[nIdx].release() );

  }

  // Save the grid. Events may need access to the second state variable values
  // Simple in this case, since only 2 dimensions
  std::vector< std::vector<double> > ppdGrids(2);
  ppdGrids[0] = pdAvgSqrReturnGrid;
  ppdGrids[1] = pdPreviousSpotGrid;

  m_pppdGrids.resize(1);
  m_pppdGrids[0] = ppdGrids;
}

void VarianceSwapPathDepStructure::PrepareForTimestepping(size_t nPathToSave)
{
  // Initialize all the paths
  size_t nNbPaths = m_path.size();
  m_nPathToSave = nPathToSave;

  for (size_t nIdx = 0; nIdx < nNbPaths; nIdx++)
  {     
    // By default, leave all paths on
    m_pbIsActive[nIdx] = true;

    // Init the objects. They have already been created
    m_path[nIdx].instdata->Init();
    m_path[nIdx].stepper->Init();

    // Set initial conditions
    m_path[nIdx].meshes->SetInitialState();
    m_path[nIdx].instdata->SetInitialValue();

  } // loop over paths
}


void VarianceSwapPathDepStructure::ConstructIdenticalTimeMeshes
                         (
                           shared_ptr<VarianceSwapMeshManager>* ppMeshes,                            
                           std::vector< AutoPtr<VarianceSwapParams> >& ppParams, 
                           const pricing::Model& model,
                           std::list< shared_ptr<PathDepEvent> >& pathDepEvents,
                           size_t nNbPaths 
                         )
{
  // @todo code copied from meshmanager.cpp. Need to put into a function
  // question: the master CB doesn't have all the special times at other 
  //  params?
  // answer: probably, yes. But the structure only has access to the cloned  
  //  params, which have potentially been modified. The master params could 
  //  be passed to the constructor, but copying from all cloned params is 
  //  arguably more flexible/safer, and the overhead is small (only 
  //  non-path-dependent events are duplicated, and the number of these 
  //  events should be small)

  // Merge the special times from the params into a master list
  SpecialTimes allSpecialTimes;
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    SpecialTimes specialTimes;
    SpecialTimes::const_iterator iter;

    ppParams[nIdx]->GetSpecialTimes(specialTimes);
    for (iter = specialTimes.begin(); iter != specialTimes.end(); ++iter)
      allSpecialTimes.push_back(*iter);
  }

  double dValuationTime = ppParams[0]->GetValuationTime();
  double dStoppingTime = ppParams[0]->GetStoppingTime();

  // Get the special times of model for the construction of the time mesh
  SpecialTimes specialTimesTmp;

  model.GetSpecialTimes(specialTimesTmp);

  // Add the useful special times in the model to the list
  SpecialTimes::const_iterator i;

  for (i = specialTimesTmp.begin(); i != specialTimesTmp.end(); ++i)
    if (i->GetTime() >= dValuationTime && i->GetTime() < dStoppingTime)
      allSpecialTimes.push_back(*i);

  // add path dependent event times
  std::list< shared_ptr<PathDepEvent> >::iterator iterEvents;
  for (iterEvents = pathDepEvents.begin();
       iterEvents != pathDepEvents.end();
       ++iterEvents)
  {
    allSpecialTimes.push_back( SpecialTime((*iterEvents)->GetTime()) );
  }

  allSpecialTimes.sort();

  SpecialTimes::iterator iter1 = allSpecialTimes.begin();
  
  // the final, well sorted, duplicated removed special times
  SpecialTimes specialTimes;
  
  specialTimes.push_back(*iter1++);

  for ( ; iter1 != allSpecialTimes.end(); ++iter1)
    if ( *iter1 > specialTimes.back() )
      specialTimes.push_back(*iter1);
    else
    {
      if (iter1->GetRefineLevel() > specialTimes.back().GetRefineLevel())
      {
        specialTimes.pop_back();
        specialTimes.push_back(*iter1);
      }
    }

  // Create a time mesh using the given special times for the master cb
  ppMeshes[0]->SetupTimeMesh(specialTimes);

  // Now create all the individual path meshes using the
  // master mesh time grid
  for (nIdx = 1; nIdx < nNbPaths; nIdx++)
    ppMeshes[nIdx]->SetupTimeMesh(ppMeshes[0]->GetTimes(),
                                  ppMeshes[0]->GetSchemeTypes(),
                                  ppMeshes[0]->GetNbTimes() );
}

void VarianceSwapPathDepStructure::UpdateOutput(size_t nIdx)
{
  // Only save output for the specified path to save. Use the
  // previously saved instdata (with full type)
  if (nIdx == m_nPathToSave)
    m_pNumOutput->UpdateMe( *m_pInstDataToSave, 
                            m_path[m_nPathToSave].meshes->GetTime() );
}

void VarianceSwapPathDepStructure::InitPathToSave()
{
  // Assume that the path to save variable has been correctly set
  ASSERT_MSG(m_nPathToSave < m_path.size(), "Invalid path to save");

  // Set the analysis date
  if ( m_flags.GetAnalysisDate().IsValid() )
  {
    m_path[m_nPathToSave].params->SetAnalysisTime( 
      GetDoubleFrom(m_flags.GetAnalysisDate()) );
  }

  // Create, initialize, etc the numoutput
  m_pInstDataToSave = (VarianceSwapInstData*) m_path[m_nPathToSave].instdata.get();

  m_pNumOutput = AutoPtr<VarianceSwapNumOutput>
    ( new VarianceSwapNumOutput( (VarianceSwapParams&) *m_path[m_nPathToSave].params) );

  m_pNumOutput->SetFinalSave(true);

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                      ( m_flags.GetComputeSurface() );

  m_pNumOutput->Init( *m_pInstDataToSave );

  m_pNumOutput->UpdateMe(*m_pInstDataToSave, 
                         m_path[m_nPathToSave].meshes->GetTime());

}

void VarianceSwapPathDepStructure::Finalize()
{
  m_pNumOutput->Finalize(*m_pInstDataToSave);
}

AutoPtr<VarianceSwapNumOutput> VarianceSwapPathDepStructure::GetOutput()
{
  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
