/////////////////////////////////////////////////////////////////////////////
// Name:        hg/pathdepstructure.h
// Purpose:     HG path dependent structure class
// Created:     2006/03/01
// RCS-ID:      $Id: pathdepstructure.h,v 1.10 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/pathdepstructure.h
    @brief HG path dependent structure class (template).
*/

#ifndef _HG_PATHDEPSTRUCTURE_H_
#define _HG_PATHDEPSTRUCTURE_H_

#include "ito33/vector.h"
#include "ito33/list.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/pathdepstructure.h"
#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/steppertridiag.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "hg/model.h"
#include "hg/payoff.h"


namespace ito33
{

namespace hg
{

/// HG path dependent structure class.
template<class InstDataType, 
         class MeshManagerType, 
         class ParamsType,
         class NumOutputType, 
         class StepperType>
class PathDepStructure : public pricing::PathDepStructure
{

public:

  /**
      Constructor.

      All these classes are needed to construct the stepper, instdata, etc.
   */
  PathDepStructure(std::vector< shared_ptr<MeshManagerType> >& ppMeshes,
                   std::vector< AutoPtr<ParamsType> >& ppParams,
                   Model& model,
                   const finance::ComputationalFlags& flags,
                   std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
                   std::vector< std::vector<double> > ppdGrids,
                   size_t nPathToSave=0)
  : m_flags(flags)
  {
    // Save the number of regimes so the pricer has access
    m_nNbRegimes = model.GetNbRegimes();

    // The path to save
    m_nPathToSave = nPathToSave;

    // Storage for the 1-D paths
    size_t nNbPaths = ppParams.size();

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

    if ( bHasContinuousEvent )
    {
      ConstructIdenticalTimeMeshes( &ppMeshes[0], ppParams, model, 
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
      // Create the instdata and stepper objects.  Params and meshes were 
      // created before. In each case, create a temporary object with the 
      // exact type (for later use), but store it in the data structure as a 
      // base type.
      shared_ptr<InstDataType> instData
        ( new InstDataType( *(ppParams[nIdx]), model, *(ppMeshes[nIdx]) ) );

      m_path[nIdx].instdata = instData;

      m_path[nIdx].stepper = make_ptr( new StepperType(*instData, m_flags) );

      // Release control of params memory to this object
      m_path[nIdx].params = make_ptr( ppParams[nIdx].release() );

    }

    // Save the grid. 
    // TODO: Either pass in full grid structure, or don't allow, or make
    //       new m_pppdGrids.  This code forces each state variable to have
    //       a constant grid (ie. X state variable, X constant grids)
    m_pppdGrids.resize(1);
    m_pppdGrids[0] = ppdGrids;
  }

  /// virtual Destructor for base class
  virtual ~PathDepStructure() { }

  /// @name virtual functions from base class
  //@{

  void InitPathToSave()
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
    m_pInstDataToSave = (InstDataType*) m_path[m_nPathToSave].instdata.get();

    m_pNumOutput = AutoPtr<NumOutputType>
      ( new NumOutputType( (ParamsType&) *m_path[m_nPathToSave].params) );

    m_pNumOutput->SetFinalSave(true);

    m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

    m_pNumOutput->Init( *m_pInstDataToSave );

    m_pNumOutput->UpdateMe(*m_pInstDataToSave, 
                           m_path[m_nPathToSave].meshes->GetTime());
  }

  void UpdateOutput(size_t nIdx)
  {
    // Only save output for the specified path to save. Use the
    // previously saved instdata (with full type)
    if (nIdx == m_nPathToSave)
      m_pNumOutput->UpdateMe( *m_pInstDataToSave, 
                              m_path[m_nPathToSave].meshes->GetTime() );
  }

  void UpdateOutputEndOfGrid(size_t nIdx) 
  {
    // Only save output for the specified path to save. Use the
    // previously saved instdata (with full type)
    if ( nIdx == m_nPathToSave)
      m_pNumOutput->UpdateMeAtEndOfGrid(*m_pInstDataToSave, 
                                     m_path[m_nPathToSave].meshes->GetTime() );
  }

  void Finalize()
  {
    m_pNumOutput->Finalize(*m_pInstDataToSave); 
  }

  double* GetPriceData(size_t nIdxPath)
  {
    // Return the data for the current regime
    size_t nNbS = m_path[nIdxPath].meshes->GetNbS();
    size_t nOffset = m_nCurrentRegime * nNbS;
    return &(m_path[nIdxPath].instdata->m_pdPrices[nOffset]);
  }
 
  //@}

  /**
      Gets the numerical output.  

      @return AutoPtr to numoutput for the requested path to save
   */
  AutoPtr<NumOutputType> GetNumOutput()
  {
    return m_pNumOutput;
  }

  /** 
      Prepares for timestepping by initializing the objects for each path.

      Must be called before pricing begins.  The initialization is 
      done here (and not in the consructor) in case an external payoff
      has been specified.
   */
  void PrepareForTimestepping()
  {
    // Initialize all the paths
    size_t nNbPaths = m_path.size();
    for (size_t nIdx = 0; nIdx < nNbPaths; nIdx++)
    {     
      // By default, leave all paths on
      m_pbIsActive[nIdx] = true;

      // Init the objects. They have already been created
      m_path[nIdx].instdata->Init();
      m_path[nIdx].stepper->Init();

      // Set initial conditions
      if ( m_pPayoff )
      {
        InstDataType* pInstDataTmp = (InstDataType*) m_path[nIdx].instdata.get();
        pInstDataTmp->SetOutsideInitialValue(m_pPayoff);       
      }

      m_path[nIdx].meshes->SetInitialState();
      m_path[nIdx].instdata->SetInitialValue();

    } // loop over paths
  }

  /**
      Sets an external initial condition/payoff.

      @param pPayoff The external initial condition
   */
  void SetInitialValue(const shared_ptr<Payoff>& pPayoff)
  {
    m_pPayoff = pPayoff;
  }

  /**
      Sets the current regime to which events will be applied.

      @param nCurrentRegime the current regime 
   */
  void SetCurrentRegime(size_t nCurrentRegime)
  {
    m_nCurrentRegime = nCurrentRegime;
  }

  /**
      Gets the number of regimes.

      @return the number of regimes
   */
  size_t GetNbRegimes()
  {
    return m_nNbRegimes;
  }


protected:

  /** 
      Constructs identical time meshes in all the paths. 

      Note: this function doesn't use anything specific to CB, so it can
      probably be implemented by using only base objets.

      @param ppMeshes (return) storage for the created mesh objects
      @param ppParams the param objects needed by mesh constructor
      @param pathDepEvents the path dependent events, needed for the times
      @param nNbPaths the number of paths
   */
  void ConstructIdenticalTimeMeshes( 
    shared_ptr< MeshManagerType >* ppMeshes,    
    std::vector< AutoPtr<ParamsType> >& ppParams, 
    const pricing::Model& model,  
    std::list< shared_ptr<pricing::PathDepEvent> >& pathDepEvents,
    size_t nNbPaths)
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
    numeric::mesh::SpecialTimes allSpecialTimes;
    size_t nIdx;
    for (nIdx = 0; nIdx < nNbPaths; nIdx++)
    {
      numeric::mesh::SpecialTimes specialTimes;
      numeric::mesh::SpecialTimes::const_iterator iter;

      ppParams[nIdx]->GetSpecialTimes(specialTimes);
      for (iter = specialTimes.begin(); iter != specialTimes.end(); ++iter)
        allSpecialTimes.push_back(*iter);
    }

    double dValuationTime = ppParams[0]->GetValuationTime();
    double dStoppingTime = ppParams[0]->GetStoppingTime();

    // Get the special times of model for the construction of the time mesh
    numeric::mesh::SpecialTimes specialTimesTmp;

    model.GetSpecialTimes(specialTimesTmp);

    // Add the useful special times in the model to the list
    numeric::mesh::SpecialTimes::const_iterator i;

    for (i = specialTimesTmp.begin(); i != specialTimesTmp.end(); ++i)
      if (i->GetTime() >= dValuationTime && i->GetTime() < dStoppingTime)
        allSpecialTimes.push_back(*i);

    // add path dependent event times
    std::list< shared_ptr<pricing::PathDepEvent> >::iterator iterEvents;
    for (iterEvents = pathDepEvents.begin();
         iterEvents != pathDepEvents.end();
         ++iterEvents)
    {
      allSpecialTimes.push_back( numeric::mesh::SpecialTime((*iterEvents)->GetTime()) );
    }

    allSpecialTimes.sort();

    numeric::mesh::SpecialTimes::iterator iter1 = allSpecialTimes.begin();
  
    // the final, well sorted, duplicated removed special times
    numeric::mesh::SpecialTimes specialTimes;
  
    specialTimes.push_back(*iter1++);

    for ( ; iter1 != allSpecialTimes.end(); ++iter1)
    {
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
    } // for loop over the times

    // Create a time mesh using the given special times for the master cb
    ppMeshes[0]->SetupTimeMesh(specialTimes);

    // Now create all the individual path meshes using the
    // master mesh time grid
    for (nIdx = 1; nIdx < nNbPaths; nIdx++)
      ppMeshes[nIdx]->SetupTimeMesh(ppMeshes[0]->GetTimes(),
                                    ppMeshes[0]->GetSchemeTypes(),
                                    ppMeshes[0]->GetNbTimes() );
  }

  /// The numerical output
  AutoPtr<NumOutputType> m_pNumOutput;

  /// instdata to update the output. Need the exact type, so can't use m_path. 
  InstDataType* m_pInstDataToSave;

  /// The computational flags (ie. what data to compute)
  const finance::ComputationalFlags& m_flags;

  /// The current regime.  Used in call to GetPriceData.
  size_t m_nCurrentRegime;

  /// The number of regimes
  size_t m_nNbRegimes;

  /// External payoff (eg. used for call notice)
  shared_ptr<Payoff> m_pPayoff;

private:

  NO_COPY_CLASS(PathDepStructure);

};

} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_PATHDEPSTRUCTURE_H_
