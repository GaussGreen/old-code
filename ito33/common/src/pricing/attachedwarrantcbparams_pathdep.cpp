/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/priicng/attachedwarrantcbparams_pathdep.cpp
// Purpose:     PathDep stuff related to attached warrant cb pricing 
// Created:     2005/07/01
// RCS-ID:      $Id: attachedwarrantcbparams_pathdep.cpp,v 1.4 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/binarysearch.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/attachedwarrantcbparams.h"
#include "ito33/pricing/sharedependentconversionevent.h"
#include "ito33/pricing/sharedependentconversion.h"

extern const double DEFAULT_BIG_RATIO;

namespace ito33
{

namespace pricing
{

  using namespace finance;
  using namespace numeric;


bool AttachedWarrantConvertibleBondParams::IsPathDependent() const
{
  
  if ( HasResetTime() )
  {
    if ( IsEqual( m_dStoppingTime, m_warrant.GetResetTime() ) )
      return CBLikeParams::IsPathDependent();
    else
      return !HasPutCallReset();
  }


  return CBLikeParams::IsPathDependent();

}

std::vector<double> 
AttachedWarrantConvertibleBondParams::ConstructPathDepGrid(Model& model) const
{
  if (    !HasResetTime() 
       || IsEqualOrBefore( m_dStoppingTime, m_warrant.GetResetTime() ) )
    return CBLikeParams::ConstructPathDepGrid(model);

  std::vector<double> pdGridY;

  // Basic idea: Use a cbmeshmanager to create a grid refined around the spot.
  // Shift this grid to be refined around the base conversion ratio
  // i.e the y grid corresponds to the conversion ratio
  // y=[cr_base , .... , cr_cap];
  
  // Get the mesh manager.  Work on clone to be safe 
  AutoPtr<CBLikeParams> cblikeparams( Clone() );
  CBMeshManager cbmeshes(*cblikeparams, model);

  cblikeparams->Init();

  cbmeshes.SetupMe();
  cbmeshes.SetInitialState();

  const size_t nGridSize = cbmeshes.GetNbS();
  const double* pdGrid = cbmeshes.GetS();

  size_t nIdxSpot = BinSearch(pdGrid, nGridSize, m_dSpot);

  // start at nIdxSpot - 1 since pdGrid[nIdxSpot] > dStrike
  // by default of the binarysearch
  nIdxSpot--;

  double dVal = pdGrid[nIdxSpot];
 
  double dBaseRatio = m_warrant.GetBaseRatio();
  double dCapRatio  = m_warrant.GetCapRatio();
  
  if ( IsEqual(dCapRatio, DEFAULT_BIG_RATIO) )
    dCapRatio = 20.* dBaseRatio;

  double dMult = dBaseRatio / dVal;

  // Make the Y grid of conversion prices. Limit the range due to caps/floors
  // starts only for the base conversion ratio
  // skip every_other point
  bool bSkip = false;
  for (size_t nIdS = nIdxSpot; nIdS < nGridSize; nIdS++)
  {
    if (bSkip)
    {
      bSkip = false;
      continue;
    }

    double dConvRatio = pdGrid[nIdS] * dMult;
  
    if ( dConvRatio < dCapRatio)
    {
      pdGridY.push_back(dConvRatio);
    }
    else 
    {
      //push in the cap ratio
      pdGridY.push_back(dCapRatio);
      break;
    }

    bSkip = true;
  } // end loop over grid

  return pdGridY;
} 

std::vector< AutoPtr<CBLikeParams> >
AttachedWarrantConvertibleBondParams::ConstructPaths
(const std::vector<double>& grid) const
{
  if (    !HasResetTime() 
       || IsEqualOrBefore( m_dStoppingTime, m_warrant.GetResetTime() ) )
    return CBLikeParams::ConstructPaths(grid);

  const size_t nNbPaths = grid.size();

  // Clone the params and create the paths
  std::vector< AutoPtr<CBLikeParams> > pCBLikeParams(nNbPaths);
   
  for ( size_t nIdPath = 0 ; nIdPath < nNbPaths ; nIdPath++)
  {
    pCBLikeParams[nIdPath] = AutoPtr<CBLikeParams>( Clone() );

    pCBLikeParams[nIdPath]->GetConversions()->SetRatios( grid[nIdPath] );

  } //end loop constructing paths

  return pCBLikeParams;
}

std::list< shared_ptr<PathDepEvent> >
AttachedWarrantConvertibleBondParams::ConstructPathDepEvents
(const std::vector< AutoPtr<CBLikeParams> >& paths) const
{
  if (    !HasResetTime() 
       || IsEqualOrBefore( m_dStoppingTime, m_warrant.GetResetTime() ) )
    return CBLikeParams::ConstructPathDepEvents(paths);

  // Create the list of events: there is only one single event
  std::list< shared_ptr<PathDepEvent> > pathDepEvents;

  ASSERT_MSG
  (
    dynamic_cast<ShareDependentConversion*>(paths[0]->GetConversions()), 
    "The type of Conversion in AWCB should be ShareDependentConversion."
  );

  shared_ptr<PathDepEvent>
    pEvent( new ShareDependentConversionEvent
                ( m_warrant.GetResetTime(),
                  static_cast<ShareDependentConversion*>(paths[0]->GetConversions())
                )
          );
   
  pathDepEvents.push_back(pEvent);

  return pathDepEvents;
}

size_t 
AttachedWarrantConvertibleBondParams::GetPathToSave
(const std::vector<double>& grid) const
{
  if (    !HasResetTime() 
       || IsEqualOrBefore( m_dStoppingTime, m_warrant.GetResetTime() ) )
    return CBLikeParams::GetPathToSave(grid);

  // NOTE: set the path to save index higher
  // than the actual total number of path
  // this way the surface is not saved since it does
  // not make sense.
  // The path to save is done in the function turn off path
  // in sharedependent conversion event
  return grid.size() + 1;
}

void 
AttachedWarrantConvertibleBondParams::InitPaths
(PathDepStructure& pathDepStruct)
{ 
  if (    !HasResetTime() 
       || IsEqualOrBefore( m_dStoppingTime, m_warrant.GetResetTime() ) )
    return CBLikeParams::InitPaths(pathDepStruct);
}


} // namespace pricing

} // namespace ito33
