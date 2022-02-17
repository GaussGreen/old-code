/////////////////////////////////////////////////////////////////////////////
// Name:        cb/resetpricer.cpp
// Purpose:     resettable cb pricer class
// Author:      Yann and David
// Created:     2004/10/06
// RCS-ID:      $Id: resetpricer.cpp,v 1.21 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/vector.h"
#include "ito33/autoptr.h"

#include "ito33/pricing/engine.h"
#include "ito33/pricing/pathdeppricer.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/resetparams.h"

#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ihg/model.h"
#include "ihg/cbinstdata.h"
#include "ihg/cbnumoutput.h"
#include "ihg/resetnumoutput.h"
#include "ihg/cbstepper.h"
#include "ihg/resetpricer.h"
#include "ihg/cbpathdepstructure.h"
#include "ihg/cbpricer.h"
#include "ihg/cbpathdeppricer.h"


using namespace ito33::finance;
using namespace ito33::pricing;

namespace ito33
{

namespace ihg
{


bool ResetPricer::IsPathDependent() const
{
  // It is path dependent if at least one of the following is true
  // - has a cash or pseudocash dividend
  // - initial conversion ratio rule and more than one reset date
  // - stock dependent vol and hazard rate

  if ( m_resetparams.IsPathDependent() ) 
    return true;

  // hazard rate must be independent of S
  if ( !m_model.IsHazardRateTimeOnly() )
      return true;

  // volatility must be independent of S
  if ( !m_model.IsVolatilityTimeOnly() )
    return true;

  return false;
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

// todo: nearly the same as CBPricer::PricePathDep!
AutoPtr<CBNumOutput> ResetPricer::PricePathDependent()
{
  // Construct the relevent structures needed by the path dep pricer
  std::vector<double> pdGridY( m_resetparams.ConstructPathDepGrid(m_model) );

  std::vector< AutoPtr<pricing::CBLikeParams> >
    cblikeParams( m_resetparams.ConstructPaths(pdGridY) );
  
  // Create the event list for the path dep pricer
  std::list< shared_ptr<pricing::PathDepEvent> > 
    pathDepEventList( m_resetparams.ConstructPathDepEvents(cblikeParams) );
 
  size_t nPathtoSave = m_resetparams.GetPathToSave(pdGridY);

  // Create and initialize the path dependent structure
  CBPathDepStructure pathDepStruct(pdGridY, cblikeParams, m_model, m_flags,
                                   pathDepEventList, nPathtoSave);

  if ( m_pPayoff )
    pathDepStruct.SetInitialValue(m_pPayoff);

  pathDepStruct.PrepareForTimestepping(nPathtoSave);

  m_resetparams.InitPaths(pathDepStruct);

  // Finally, price
  CBPathDepPricer pricer;

  pricer.Price(pathDepStruct, pathDepEventList);

  AutoPtr<CBNumOutput> pNumOutput = pathDepStruct.GetOutput();

  return pNumOutput;
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
  pricing::CBMeshManager cbmeshes(*pClonedParams, m_model);

  CBInstData cbinstdata(*pClonedParams, m_model, cbmeshes);

  if ( m_pPayoff )
    cbinstdata.SetOutsideInitialValue(m_pPayoff);

  AutoPtr<ResetNumOutput> pNumOutput( new ResetNumOutput(*pClonedParams) );

  CBStepper cbstepper(cbinstdata, m_flags);

  ResetEngine 
    resetengine(*pClonedParams, cbmeshes, cbinstdata, cbstepper, *pNumOutput);

  // Now initialize the objects
  if ( !pClonedParams->GetCalls()->HasMakeWhole() && 
       !pClonedParams->GetCalls()->HasNoticePeriod() )
    cbinstdata.m_bComputeVega  = m_flags.GetComputeVega();

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
