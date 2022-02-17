/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cb/cbpricer.cpp
// Purpose:     cb pricer class
// Author:      Nabil
// Created:     2004/04/14
// RCS-ID:      $Id: cbpricer.cpp,v 1.64 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/engine.h"
#include "ito33/pricing/cb.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/callprovisions.h"
#include "ito33/pricing/cbcalls.h"

#include "ito33/finance/computationalflags.h"

#include "ihg/model.h"
#include "ihg/cbinstdata.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbstepper.h"
#include "ihg/cbpricer.h"
#include "ihg/cbpathdepstructure.h"
#include "ihg/cbpathdeppricer.h"


using namespace ito33::finance;
using namespace ito33::pricing;

namespace ito33
{

namespace ihg
{


CBPricer::CBPricer(pricing::CBLikeParams& cbparams, 
                   Model& model, 
                   const finance::ComputationalFlags& flags)
                 : m_cbparams(cbparams), 
                   m_model(model),
                   m_flags(flags)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_cbparams.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}

AutoPtr<CBNumOutput> CBPricer::Price()
{    
  if ( m_cbparams.IsPathDependent() )
    return PricePathDep();
  else
    return PriceNormal();
}

AutoPtr<CBNumOutput> CBPricer::PriceNormal()
{
  // Construct the objects needed for pricing
  pricing::CBMeshManager cbmeshes(m_cbparams, m_model);

  CBInstData cbinstdata(m_cbparams, m_model, cbmeshes);
  
  cbinstdata.SetSolverType( m_flags.GetSolverType() );

  AutoPtr<CBNumOutput> pCBNumOutput( new CBNumOutput(m_cbparams) );

  CBStepper cbstepper(cbinstdata, m_flags);

  CBEngine cbengine(m_cbparams,cbmeshes,cbinstdata, cbstepper, *pCBNumOutput);
  
  // Now initialize the objects
  if (    !m_cbparams.GetCalls()->HasMakeWhole()  
       && !m_cbparams.GetCalls()->HasNoticePeriod()
       && !m_cbparams.HasNewShare() )
    cbinstdata.m_bComputeVega  = m_flags.GetComputeVega();

  cbinstdata.m_bComputeFugit = m_flags.GetComputeFugit();

  m_cbparams.Init();

  cbmeshes.SetupMe();

  pCBNumOutput->SetFinalSave(true);

  pCBNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( m_flags.GetComputeSurface() );

  if ( m_pPayoff )
    cbinstdata.SetOutsideInitialValue(m_pPayoff);

  // Actually do the pricing and return the result
  cbengine.Run();

  return pCBNumOutput;
}


AutoPtr<CBNumOutput> CBPricer::PricePathDep()
{
  ASSERT_MSG(m_cbparams.IsPathDependent(),
             "No path dependent clause found in PricePathDep");

  // Construct the relevent structures needed by the path dep pricer
  std::vector<double> pdGridY( m_cbparams.ConstructPathDepGrid(m_model) );
  
  std::vector< AutoPtr<pricing::CBLikeParams> >
    cblikeParams( m_cbparams.ConstructPaths(pdGridY) );
  
  // Create the event list for the path dep pricer
  std::list< shared_ptr<pricing::PathDepEvent> > 
    pathDepEventList( m_cbparams.ConstructPathDepEvents(cblikeParams) );
 
  size_t nPathtoSave = m_cbparams.GetPathToSave(pdGridY);

  // Create the path dependent structure
  CBPathDepStructure pathDepStruct(pdGridY, cblikeParams, m_model, m_flags,
                                   pathDepEventList, nPathtoSave);

  if ( m_pPayoff )
    pathDepStruct.SetInitialValue(m_pPayoff);

  pathDepStruct.PrepareForTimestepping(nPathtoSave);

  m_cbparams.InitPaths(pathDepStruct);

  // Finally, price
  CBPathDepPricer pricer;

  pricer.Price(pathDepStruct, pathDepEventList);

  AutoPtr<CBNumOutput> pNumOutput = pathDepStruct.GetOutput();

  return pNumOutput;
} 


} // namespace ihg

} // namespace ito33
