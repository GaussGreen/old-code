/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/mandatoryparams_pathdep.cpp
// Purpose:     PathDep stuff related to mandatory pricing
// Created:     2005/07/06
// RCS-ID:      $Id: mandatoryparams_pathdep.cpp,v 1.5 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/binarysearch.h"

#include "ito33/numeric/predicatedouble.h"
#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/mandatory.h"
#include "ito33/pricing/mandatoryparams.h"
#include "ito33/pricing/mandatoryconversion.h"
#include "ito33/pricing/pepsaveragingevent.h"
#include "ito33/pricing/constructeventtimes.h"

namespace ito33
{

namespace pricing
{


// Todo: should we disable explicitly n days call in case of averaging?
bool MandatoryParams::IsPathDependent() const
{
  return m_mandatory.HasAveragingPeriod() || CBLikeParams::IsPathDependent();
}

std::vector<double> 
MandatoryParams::ConstructPathDepGrid(Model& model) const
{
  if ( !m_mandatory.HasAveragingPeriod() )
    return CBLikeParams::ConstructPathDepGrid(model);

  double dCurrentAverage = m_mandatory.GetCurrentAverage();
  std::vector<double> pdGridY;

  // between avg end date and the maturity date
  // only one path need to be solved
  if ( numeric::IsEqualOrAfter( GetValuationTime(), 
       m_mandatory.GetAverageEndTime() ) )
  {  
    pdGridY.push_back(dCurrentAverage);
    return pdGridY;
  }

  // Get the mesh manager.  Work on clone to be safe 
  AutoPtr<MandatoryParams> params( Clone() );
  CBMeshManager cbmeshes(*params, model);

  params->Init();

  cbmeshes.SetupMe();
  cbmeshes.SetInitialState();

  size_t nGridSize     = cbmeshes.GetNbS();
  const double* pdGrid = cbmeshes.GetS();
 
 
  if ( m_mandatory.HasStockAveraging() )
  {  
    // Might be better to use a temporary vector and then copy back 
    // with current average inserted
    if ( pdGrid[0] > dCurrentAverage )
      pdGridY.push_back(dCurrentAverage);

    for (size_t nIdxS = 0; nIdxS < nGridSize; nIdxS++)
    {
      pdGridY.push_back(pdGrid[nIdxS]);
      
      // add current average if not in the grid
      if (   nIdxS < nGridSize - 1
          && pdGrid[nIdxS] < dCurrentAverage 
          && dCurrentAverage < pdGrid[nIdxS + 1]  
          && !numeric::IsEqual(pdGrid[nIdxS], dCurrentAverage)  
          && !numeric::IsEqual(pdGrid[nIdxS + 1], dCurrentAverage)
         )
        pdGridY.push_back(dCurrentAverage); 
    }

    if ( pdGrid[nGridSize - 1] < dCurrentAverage )
      pdGridY.push_back(dCurrentAverage);
  }
  else
  {
    double dFactor = dCurrentAverage / m_dSpot;

    for (size_t nIdS = 0; nIdS < nGridSize; nIdS++)  
      pdGridY.push_back( pdGrid[nIdS] * dFactor );
  }  

  return pdGridY;
}

std::vector< AutoPtr<CBLikeParams> >
MandatoryParams::ConstructPaths(const std::vector<double>& grid) const
{
  if ( !m_mandatory.HasAveragingPeriod() )
    return CBLikeParams::ConstructPaths(grid);

  size_t nNbPaths = grid.size();

  // Clone the params and create the paths
  std::vector< AutoPtr<CBLikeParams> > pMandatoryParams(nNbPaths);
  bool bHasStockAveraging = m_mandatory.HasStockAveraging(); 
  
  for (size_t nIdxPath = 0 ; nIdxPath < nNbPaths ; nIdxPath++)
  {
    AutoPtr<MandatoryParams> pMandTmp( Clone() );
  
    if ( bHasStockAveraging )
      pMandTmp->GetMandatory().SetStockAverage(grid[nIdxPath]);
    else
      pMandTmp->GetMandatory().SetConversionRatioAverage(grid[nIdxPath]);

    pMandatoryParams[nIdxPath] = AutoPtr<CBLikeParams>( pMandTmp.release() );  
  
  } //end loop constructing paths

  return pMandatoryParams;
}

std::list< shared_ptr<PathDepEvent> >
MandatoryParams::ConstructPathDepEvents
(const std::vector< AutoPtr<CBLikeParams> >& paths) const
{
  if ( !m_mandatory.HasAveragingPeriod() )
    return CBLikeParams::ConstructPathDepEvents(paths);

  std::list< shared_ptr<PathDepEvent> > pathDepEvents;

  double dStartSamplingTime = m_mandatory.GetAverageStartTime();
  double dStopSamplingTime  = m_mandatory.GetAverageEndTime();
  double dValuationTime     = GetValuationTime();
  size_t nNbSampling        = m_mandatory.GetNbSamplingAverages();
  size_t nNbSamplesUsed     = m_mandatory.GetNbSamplesUsed();

  std::vector<double> pdEventTimes;
  std::vector<size_t> pnSamplingNumbers;

  ConstructSamplingEventTimes(dValuationTime, dStartSamplingTime, dStopSamplingTime,
    nNbSampling, nNbSamplesUsed, pdEventTimes, pnSamplingNumbers);

  if ( pdEventTimes.empty() )
    return pathDepEvents;

  bool bHasStockAveraging = m_mandatory.HasStockAveraging(); 

  // First event
  bool bIsFirstEvent = true;
  shared_ptr<PEPSAveragingEvent> pEvent
       (
         new PEPSAveragingEvent(pdEventTimes[0], 
                                pnSamplingNumbers[0], 
                                bIsFirstEvent, 
                                m_mandatory.GetMandatoryConversion(), 
                                bHasStockAveraging)
        );

  pathDepEvents.push_back(pEvent);

  bIsFirstEvent = false;
  for (size_t nIdx = 1; nIdx < pdEventTimes.size() ; nIdx++)
  {
    shared_ptr<PEPSAveragingEvent> pEvent
       (
         new PEPSAveragingEvent(pdEventTimes[nIdx], 
                                pnSamplingNumbers[nIdx], 
                                bIsFirstEvent, 
                                m_mandatory.GetMandatoryConversion(), 
                                bHasStockAveraging)
        );
 
    pathDepEvents.push_back(pEvent);
  }

  return pathDepEvents;
}

size_t MandatoryParams::GetPathToSave(const std::vector<double>& pdGrid) const
{
  if ( !m_mandatory.HasAveragingPeriod() )
    return CBLikeParams::GetPathToSave( pdGrid);

  double dCurrentAverage = m_mandatory.GetCurrentAverage();

  if ( pdGrid.size() == 1 )
  {
    ASSERT_MSG( numeric::IsEqual(dCurrentAverage, pdGrid[0]),
      "The avg grid does not contain the current average.");
    
    return 0;
  }

  // the path to save is the one containing the current average
  size_t nPathToSave = BinSearch(& pdGrid[0], pdGrid.size(), dCurrentAverage);

  if ( !numeric::IsEqual( pdGrid[nPathToSave], dCurrentAverage) )
    nPathToSave--;

  return nPathToSave;
}

void MandatoryParams::InitPaths(PathDepStructure& pathDepStruct)
{
  if ( !m_mandatory.HasAveragingPeriod() )
    CBLikeParams::InitPaths(pathDepStruct);
}


} // namespace pricing

} // namespace ito33
