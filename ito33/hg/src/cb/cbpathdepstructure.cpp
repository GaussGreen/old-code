/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cbpathdepstructure.cpp
// Purpose:     path dependent strcuture class implementation for HG model
// Created:     2005/04/11
// RCS-ID:      $Id: cbpathdepstructure.cpp,v 1.8 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <algorithm>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cb.h"
#include "ito33/pricing/pathdepevent.h"

#include "hg/cbinstdata.h"
#include "hg/cbstepper.h"
#include "hg/cbnumoutput.h"
#include "hg/model.h"
#include "hg/cbpathdepstructure.h"

namespace ito33
{
  using namespace pricing;
  using namespace numeric::mesh;

namespace hg
{

CBPathDepStructure::CBPathDepStructure
                    (  
                       std::vector<double>& pdGridY,
                       std::vector< AutoPtr<CBLikeParams> >& ppParams,
                       Model& model,
                       const finance::ComputationalFlags& flags,
                       std::list< shared_ptr<PathDepEvent> >& pathDepEvents,
                       size_t nPathToSave
                    )
                    : m_flags(flags)
{

  // Path to save for output
  m_nPathToSave = nPathToSave;

  // Storage for the 1-D paths
  size_t nNbPaths = pdGridY.size();

  m_path.resize(nNbPaths);
  m_pbIsActive.resize(nNbPaths);
 
  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {
    // initilize params here, but release below since they are still needed
    ppParams[nIdx]->Init();
  }

  // If any continuous events are present, we need each path to have the 
  // same mesh (technically, the paths need to have the same mesh during
  // the continuous window, but for convenience, make the entire time
  // mesh identical)
  bool bHasContinuousEvent = false;
  std::list< shared_ptr<pricing::PathDepEvent> >::iterator iterEvents;
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

  // Do the appropriate mesh construction
  Array< shared_ptr<CBMeshManager> > ppMeshes(nNbPaths);

  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
    ppMeshes[nIdx] = make_ptr( new CBMeshManager(*(ppParams[nIdx]), model) );

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

  // Construct numerical output class for the path to save
  // (in case meshes are different)
  if (m_nPathToSave < nNbPaths)
  {
    m_pNumOutput = AutoPtr<CBNumOutput>
                   ( new CBNumOutput(*ppParams[m_nPathToSave]) );
 
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
    shared_ptr<CBInstData> 
      instDataCB( new CBInstData(*(ppParams[nIdx]), model, *(ppMeshes[nIdx])) );

    instDataCB->m_bComputeFugit = false;

    m_path[nIdx].instdata = instDataCB;
    m_path[nIdx].instdata->Init();

    m_path[nIdx].stepper = make_ptr( new CBStepper(*instDataCB, m_flags) );
    m_path[nIdx].stepper->Init();

    // Set initial conditions
    m_path[nIdx].meshes->SetInitialState();
    m_path[nIdx].instdata->SetInitialValue();

    // Release control of params memory to this object
    m_path[nIdx].params = make_ptr( ppParams[nIdx].release() );

    // Init the numerical output (if correct path). Also, save the instdata,
    // since we need access to the type to update numoutput
    if (nIdx == m_nPathToSave)
    {
      m_pInstDataToSave = instDataCB.get();
      
      m_pNumOutput->SetFinalSave(true);

      m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

      m_pNumOutput->Init(*instDataCB);
      m_pNumOutput->UpdateMe(*instDataCB, m_path[nIdx].meshes->GetTime());
    }
  }


  // Save the grid. Events may need access to the second state variable values
  // Simple in this case, since only 2 dimensions
  std::vector< std::vector<double> > ppdGrids(1);
  ppdGrids[0] = pdGridY;

  m_pppdGrids.resize(1);
  m_pppdGrids[0] = ppdGrids;
}



void CBPathDepStructure::PrepareForTimestepping(size_t nPathToSave)
{
  // Construct numerical output class for the path to save
  size_t nNbPaths = m_path.size();
  m_nPathToSave = nPathToSave;
  if (m_nPathToSave < nNbPaths)
  {
    m_pNumOutput = AutoPtr<CBNumOutput>
      ( new CBNumOutput( (CBLikeParams&) *m_path[m_nPathToSave].params) );
  }

  size_t nIdx;
  for (nIdx = 0; nIdx < nNbPaths; nIdx++)
  {     
    // By default, leave all paths on
    m_pbIsActive[nIdx] = true;

    // Init the objects. They have already been created
    m_path[nIdx].instdata->Init();
    m_path[nIdx].stepper->Init();

    // Set initial conditions
    m_path[nIdx].meshes->SetInitialState();
    m_path[nIdx].instdata->SetInitialValue();

    // Init the numerical output (if correct path). Also, save the instdata,
    // since we need access to the type to update numoutput
    if (nIdx == m_nPathToSave)
    {
      m_pInstDataToSave = (CBInstData*) m_path[nIdx].instdata.get();

      m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

      m_pNumOutput->Init(*m_pInstDataToSave);
      m_pNumOutput->UpdateMe(*m_pInstDataToSave, m_path[nIdx].meshes->GetTime());
    }
  } // loop over paths
}


void CBPathDepStructure::ConstructIdenticalTimeMeshes
                         (
                           shared_ptr<CBMeshManager>* ppMeshes,                            
                           std::vector< AutoPtr<CBLikeParams> >& ppParams, 
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

void CBPathDepStructure::UpdateOutput(size_t nIdx)
{
  // Only save output for the specified path to save. Use the
  // previously saved instdata (with full type)
  if (nIdx == m_nPathToSave)
    m_pNumOutput->UpdateMe( *m_pInstDataToSave, 
                            m_path[m_nPathToSave].meshes->GetTime() );
}

void CBPathDepStructure::SetPathToSave(size_t nIdx)
{
  // Make sure the path is within range
  ASSERT_MSG(nIdx < m_path.size(), "Invalid path to save");

  size_t nOldPathToSave = m_nPathToSave;

  // Update path to save and switch the pointer to instdata
  m_nPathToSave = nIdx;
  m_pInstDataToSave = (CBInstData*) m_path[m_nPathToSave].instdata.get();

  // If the path to save was not originally set, create the numoutput now
  if (nOldPathToSave >= m_path.size() )
  {
    m_pNumOutput = AutoPtr<CBNumOutput>
      ( new CBNumOutput( (CBLikeParams&) *m_path[m_nPathToSave].params) );

    m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

    m_pNumOutput->Init( *m_pInstDataToSave );
  } 
}

void CBPathDepStructure::Finalize()
{
  m_pNumOutput->Finalize(*m_pInstDataToSave);
}

AutoPtr<CBNumOutput> CBPathDepStructure::GetNumOutput()
{
  return m_pNumOutput;
}


} // namespace hg

} // namespace ito33
