/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/pathdeppricer.cpp
// Purpose:     path dependent pricer class for ihg cb problems
// Author:      Yann and David
// Created:     8/03/2005
// RCS-ID:      $Id: cbpathdeppricer.cpp,v 1.4 2006/08/03 21:25:18 dave Exp $
// Copyright:   (c) 2004 - 2006  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/constants.h"

#include "ito33/pricing/pathdepstructure.h"

#include "ito33/numeric/predicatetime.h"

#include "ihg/instdatamulti.h"
#include "ihg/cbpathdeppricer.h"


namespace ito33
{

using numeric::AreTimesEqual;

namespace ihg
{

void CBPathDepPricer::AdvanceToTime(double dStopTime, 
                                    pricing::PathDepStructure& pathDepStruct)
{
  // Price until the stop time, which is typically an event time.
  // If there are no continuous events, the solution of all paths
  // can simply be advanced to the stop time. Otherwise, need to advance
  // all paths step by step while applying the continuous events
  size_t nNbPaths = pathDepStruct.m_path.size();
  if ( m_continuousEvents.size() == 0 )
  {
    // Should be more efficient to solve each path for as long as
    // possible, instead of proceeding in lock-step
    for (size_t nIdx = 0; nIdx < nNbPaths; nIdx++)
    {
      while 
      (
      !AreTimesEqual(pathDepStruct.m_path[nIdx].meshes->GetTime(), dStopTime)
      && pathDepStruct.m_path[nIdx].meshes->TryGoAhead() 
      )
      {
        pricing::BackwardMeshManagerMulti* mesh = 
          (pricing::BackwardMeshManagerMulti*) pathDepStruct.m_path[nIdx].meshes.get();

        if (mesh->IsEndOfGrid())
        {
          InstDataMulti *   
            instData = (InstDataMulti *) pathDepStruct.m_path[nIdx].instdata.get();

          instData->UpdateAtEndOfGrid();
          
          if (m_bSaveOutput)
            pathDepStruct.UpdateOutputEndOfGrid(nIdx);
          
        }
        else
        {
          AdvanceOneStep(nIdx, pathDepStruct);

          if ( m_bSaveOutput )
            pathDepStruct.UpdateOutput(nIdx);

        }
      }
    } // for loop over paths
  }
  else
  {
    // Go step by step
    double dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();
    while ( !AreTimesEqual(dCurrentTime, dStopTime) )
    {
      AdvanceOneStepForAll(pathDepStruct);
      ApplyContinuousEvents(pathDepStruct);
      dCurrentTime = pathDepStruct.m_path[0].meshes->GetTime();

      ASSERT_MSG( m_bSaveOutput == false, 
                  "Cannot save output while continuous events are active");

    } // while timestepping step by step

  } // if continuous events or not

}

} // namespace ihg

} // namespace ito33
