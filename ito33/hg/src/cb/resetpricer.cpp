/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/resetpricer.cpp
// Purpose:     resettable cb pricer class
// Created:     2006/04/17
// RCS-ID:      $Id: resetpricer.cpp,v 1.5 2006/08/19 23:44:26 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/autoptr.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/pathdeppricer.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/resetparams.h"
#include "ito33/pricing/reset.h"

#include "hg/model.h"
#include "hg/cbinstdata.h"
#include "hg/cbnumoutput.h"
#include "hg/resetnumoutput.h"
#include "hg/cbstepper.h"
#include "hg/resetpricer.h"
#include "hg/cbpricer.h"
#include "hg/pathdepstructure.h"
#include "hg/pathdeppricer.h"


namespace ito33
{

namespace hg
{

using namespace pricing;

bool ResetPricer::IsPathDependent() const
{
  // The model is homogeneous, so only need to check params
  return m_resetparams.IsPathDependent();
}

AutoPtr<CBNumOutput> ResetPricer::Price()  
{
  // If there are no active reset times, then simply use the cb pricer
  if ( !m_resetparams.HasActiveResetDate() )
  {
    CBPricer cbpricer(m_resetparams, m_model, m_flags);

    if ( m_pPayoff )
      cbpricer.SetInitialValue(m_pPayoff);

    return cbpricer.Price();
  }

  // If path dependent, use path dependent pricer
  if ( IsPathDependent() )
    return PricePathDependent();

  // By default, use 1D pricing by setting y = kS
  return PriceOneD();
}


AutoPtr<CBNumOutput> ResetPricer::PricePathDependent()
{
  // Construct the relevent structures needed by the path dep pricer
  std::vector<double> pdGridY( m_resetparams.ConstructPathDepGrid(m_model) );

  std::vector< AutoPtr<CBLikeParams> >
    pParams( m_resetparams.ConstructPaths(pdGridY) );
  
  // Find the path to save
  size_t nIdxPathToSave = m_resetparams.GetPathToSave(pdGridY);

  // Create the event list for the path dep pricer
  std::list< shared_ptr<PathDepEvent> > 
    pathDepEvents( m_resetparams.ConstructPathDepEvents(pParams) );

  // Do the appropriate mesh construction
  size_t nNbPaths = pdGridY.size();
  std::vector< shared_ptr<CBMeshManager> > ppMeshes(nNbPaths);

  for (size_t nIdx = 0; nIdx < nNbPaths; nIdx++)
    ppMeshes[nIdx] = make_ptr( new CBMeshManager(*(pParams[nIdx]), m_model) );

  // Create the structure
  std::vector< std::vector<double> > ppdGrids(1);
  ppdGrids[0] = pdGridY;  

  PathDepStructure<CBInstData, 
                   CBMeshManager, 
                   CBLikeParams, 
                   CBNumOutput, 
                   CBStepper>
    pathDepStruct(ppMeshes, pParams, m_model, m_flags, 
                  pathDepEvents, ppdGrids, nIdxPathToSave);
  
  if ( m_pPayoff )
    pathDepStruct.SetInitialValue(m_pPayoff);

  pathDepStruct.PrepareForTimestepping();

  // Solve
  PathDepPricer<CBInstData, 
                CBMeshManager, 
                CBLikeParams, 
                CBNumOutput, 
                CBStepper>
    pathDepPricer(pathDepStruct);

  pathDepPricer.Price(pathDepEvents);

  return pathDepStruct.GetNumOutput();
} 

AutoPtr<CBNumOutput> ResetPricer::PriceOneD()
{
  // Since params are changed below, work on clone. This allows the same
  // params to be passed to the pricer multiple times (as with call notice)
  shared_ptr<ResetParams> pClonedParams( m_resetparams.Clone() );
  pClonedParams->SetAnalysisTime( m_resetparams.GetAnalysisTime() );

  // 1d pricing is based on the conversion ratio via y = kS. The caps
  // and floors are entered in terms of the conversion price.  Hence, need
  // to convert the caps and floors from price to ratio
  pClonedParams->GetReset().FlipFloorAndCap();

  // Conversion value is typically kS, where S is the grid values.  If we set 
  // k=1, then the conversion is based on the grid, which is correct for 1d
  pClonedParams->GetConversions()->SetRatios(1.0);

  // Shift the spot price so initial spot times initial conversion ratio is 
  // in the grid
  double dCurrentRatio =  m_resetparams.GetReset().GetNominal() 
   / m_resetparams.GetReset().GetCurrentConversionPrice() ;

  pClonedParams->SetSpotSharePrice
                  ( m_resetparams.GetSpotSharePrice() * dCurrentRatio );

  // Construct the objects needed for pricing
  CBMeshManager cbmeshes(*pClonedParams, m_model);

  CBInstData cbinstdata(*pClonedParams, m_model, cbmeshes);

  if ( m_pPayoff )
    cbinstdata.SetOutsideInitialValue(m_pPayoff);

  AutoPtr<ResetNumOutput> pNumOutput( new ResetNumOutput(*pClonedParams) );

  CBStepper cbstepper(cbinstdata, m_flags);

  ResetEngine 
    resetengine(*pClonedParams, cbmeshes, cbinstdata, cbstepper, *pNumOutput);

  cbinstdata.m_bComputeFugit = m_flags.GetComputeFugit();

  pClonedParams->Init();

  // add the reset events to the eventManager
  pClonedParams->ConstructResetEvents();

  cbmeshes.SetupMe();

  pNumOutput->SetFinalSave(true);

  pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

  pClonedParams->SetSpotSharePrice( m_resetparams.GetSpotSharePrice() );

  // Actually do the pricing and return the result
  resetengine.Run();

  AutoPtr<CBNumOutput> pCBNumOutput(pNumOutput.release());
  
  return pCBNumOutput;

} //AutoPtr<CBNumOutput> ResetPricer::PriceOneD()

} // namespace ihg

} // namespace ito33
