/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/attachedwarrantcbpricer.cpp
// Purpose:     attached warrant cb pricer class
// Author:      Ito33
// Created:     2005/01/17
// RCS-ID:      $Id: attachedwarrantcbpricer.cpp,v 1.15 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/constants.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/finance/payoffdiscrete.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/attachedwarrantcbparams.h"
#include "ito33/pricing/sharedependentconversion.h"

#include "ito33/numeric/predicatedouble.h"

#include "ihg/model.h"
#include "ihg/cbinstdata.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbpathdepstructure.h"
#include "ihg/cbpathdeppricer.h"
#include "ihg/cbpricer.h"
#include "ihg/attachedwarrantcbpricer.h"

using namespace ito33::finance;
using namespace ito33::pricing;
using namespace ito33::numeric;

extern const double DEFAULT_BIG_RATIO;

namespace ito33
{

namespace ihg
{

 
bool AttachedWarrantConvertibleBondPricer::CanCalculatePriceAndGreeks()
{
  //When the current conversion ratio is stricly positive
  //this implies that the valuation time is after the reset
  //time.
  if ( !m_warrantparams.HasResetTime() 
      || m_warrantparams.HasPutCallReset()
      || m_warrantparams.GetCurrentConversionRatio() > 0. 
     )
    return true;

  return false;
}

AutoPtr<CBNumOutput> AttachedWarrantConvertibleBondPricer::Price()  
{

  // Price accordingly
  AutoPtr<CBNumOutput> pOutput;
 
  //reset the flags
  m_pFlagsPricer = make_ptr( new finance::ComputationalFlags(m_flags) );

  if ( CanCalculatePriceAndGreeks() )
  {
    //The only way to have the current conversion ratio strictly
    //positive is when the valuation date is after the reset time
    if ( !m_warrantparams.HasResetTime() 
      || m_warrantparams.GetCurrentConversionRatio() > 0. )
    {
      pOutput = Price1D();
    }
    else if ( m_warrantparams.HasPutCallReset() )
    {
      double dResetTime = m_warrantparams.GetResetTime();

      // Just build the mesh from the valuation time to the reset time
      m_warrantparams.SetStoppingTime(dResetTime);

      // Without coco from valuation to reset time, the payoff does
      // not matter.  However, if there is coco, there can be a trigger,
      // and the put/call constraints will not correctly set the values.
      // This function should take care of this
      m_warrantparams.GetCBLike().SetPayoff(ConstructPutCallPayoff()); 

      // Actually, just call the cb pricer.  The "1d" just means it is
      // not path dependent due to a share dependent feature

      // If there is new share, we can price this way, but not for the
      // greeks since the payoff will depend on the parameters
      // but I'm told that there is no new share with attached warrant cb

      pOutput = Price1D();
    }
   
  }
  else
    pOutput = PriceOnly();


  return pOutput;
}

AutoPtr<CBNumOutput> AttachedWarrantConvertibleBondPricer::PriceOnly()
{
  //Turn off all Greeks
  m_pFlagsPricer->TurnOffAllGreeks();

  // Price accordingly
  AutoPtr<CBNumOutput> pOutput;

  double dResetTime = m_warrantparams.GetResetTime();
  double dValuationTime = m_warrantparams.GetValuationTime();

  //solve from the maturity time to the reset time
  m_warrantparams.SetValuationTime(dResetTime);
  pOutput = Price2D();

  // A cloned params is used instead of the origin params, otherwise
  // the code below will change the payoff, and the payoff will be used
  // during the computation of greeks by perturbation.
  AutoPtr<CBLikeParams> clonedparams( m_warrantparams.Clone() );

  clonedparams->SetValuationTime(dValuationTime);

  // Analysis date is lost during the cloning
  clonedparams->SetAnalysisTime( m_warrantparams.GetAnalysisTime() );

  // Go just before reset time to avoid re-applying events which
  // should have been applied at the end of the previous solve
  // This should be considered as a hack and analysis time will
  // be ignored in this case.
  clonedparams->SetStoppingTime(dResetTime - TIMETOLERANCE * 2.);
  
  // Extract the final solution and use it to set the payoff for
  // a "1D" (normal cb) solve from the reset time to the valuation time
  const std::vector<double>& pdV = pOutput->GetFinalPrices();
  const std::vector<double>& pdS = pOutput->GetFinalMesh();    

  clonedparams->GetCBLike().SetPayoff(shared_ptr<finance::Payoff>(
      new finance::PayoffDiscrete(&pdS[0], &pdV[0], pdS.size())));

  // May not actually be 1D due to coco
  CBPricer cbPricer(*clonedparams, m_model, *m_pFlagsPricer);
  pOutput = cbPricer.Price();

  return pOutput;

}

AutoPtr<CBNumOutput> AttachedWarrantConvertibleBondPricer::Price1D()  
{
  // After the reset date (pricing backward), the problem is simply a 
  // normal cb with a share dependent conversion constraint. Simply
  // use the solver of a regular cb (which has the side benefit
  // of supporting coco)
  CBPricer cbPricer(m_warrantparams, m_model, *m_pFlagsPricer);
  return cbPricer.Price();
}

AutoPtr<CBNumOutput> AttachedWarrantConvertibleBondPricer::Price2D()  
{
  // Construct the relevent structures needed by the path dep pricer
  std::vector<double> pdGridY( m_warrantparams.ConstructPathDepGrid(m_model) );
  
  std::vector< AutoPtr<pricing::CBLikeParams> >
    cblikeParams( m_warrantparams.ConstructPaths(pdGridY) );
  
  // Create the event list for the path dep pricer
  std::list< shared_ptr<pricing::PathDepEvent> > 
    pathDepEventList( m_warrantparams.ConstructPathDepEvents(cblikeParams) );
 
  size_t nPathtoSave = m_warrantparams.GetPathToSave(pdGridY);

  // Create the path dependent structure
  CBPathDepStructure pathDepStruct(pdGridY, cblikeParams, m_model, m_flags,
                                   pathDepEventList, nPathtoSave);

  if ( m_pPayoff )
    pathDepStruct.SetInitialValue(m_pPayoff);

  pathDepStruct.PrepareForTimestepping(nPathtoSave);

  m_warrantparams.InitPaths(pathDepStruct);

  // Finally, price
  CBPathDepPricer pricer;

  pricer.Price(pathDepStruct, pathDepEventList);

  AutoPtr<CBNumOutput> pNumOutput = pathDepStruct.GetOutput();

  return pNumOutput;

} //  AutoPtr<CBNumOutput> AttachedWarrantConvertibleBondPricer::Price()

shared_ptr<finance::Payoff> 
AttachedWarrantConvertibleBondPricer::ConstructPutCallPayoff()
{
  // Get the mesh manager.  Work on clone to be safe 
  AutoPtr<CBLikeParams> clonedparams( m_warrantparams.Clone() );
  pricing::CBMeshManager cbmeshes(*clonedparams, m_model);

  // sets the params in the puts, calls, etc.  Sets up events, etc.
  clonedparams->Init();

  // make the mesh
  cbmeshes.SetupMe();

  // activate the appropriate sub-grid before getting the grid.
  // Also calls clonedparams->SetInitialState(resettime)
  cbmeshes.SetInitialState();

  size_t nNbS       = cbmeshes.GetNbS();
  const double* pdS = cbmeshes.GetS();

  Array<double> pdValues(nNbS);

  //New share treatment
  Array<double> pdNewSharePrices(nNbS);
  CBInstData cbinstdata(*clonedparams, m_model, cbmeshes);
  if( cbinstdata.HasNewShare() )
  {
    cbinstdata.ComputeNewSharePricesAtMaturity( pdS, nNbS, 
      pdNewSharePrices.Get() );
  }
  else
  {
    for(size_t nIdxS = 0; nIdxS < nNbS; nIdxS++)
      pdNewSharePrices[nIdxS] = pdS[nIdxS];
  }

  // If coco is set, the conversion may have a trigger that should not be 
  // applied at reset time. Note also that the coupon index has not been set
  clonedparams->GetConversions()->RemoveTriggers();
  //clonedparams->GetConversionConstraintValues(pdS, nNbS, pdValues.Get() );
  size_t nIdxConversion;
  clonedparams->GetCallConstraintValues(pdS, nNbS, pdValues.Get(), 
                                        nIdxConversion,
                                        pdNewSharePrices.Get(), false );

  return shared_ptr<finance::Payoff>(
            new finance::PayoffDiscrete(pdS, pdValues.Get(), nNbS));
}


} // namespace ihg

} // namespace ito33
