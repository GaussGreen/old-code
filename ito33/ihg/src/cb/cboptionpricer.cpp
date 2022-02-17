/////////////////////////////////////////////////////////////////////////////
// Name:        cb/cboptionpricer.cpp
// Purpose:     cb option pricer class
// Author:      Nabil
// Created:     2005/10/18
// RCS-ID:      $Id: cboptionpricer.cpp,v 1.4 2006/06/13 15:41:30 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cboptionparams.h"
#include "ito33/pricing/cboption.h"

#include "ito33/finance/computationalflags.h"

#include "ihg/cboptioninstdata.h"
#include "ihg/cboptionnumoutput.h"
#include "ihg/cboptionpricer.h"
#include "ihg/cboptionstepper.h"
#include "ihg/model.h"

using namespace ito33::finance;
using namespace ito33::pricing;

namespace ito33
{

namespace ihg
{


CBOptionPricer::CBOptionPricer(pricing::CBOptionParams& cboptionparams, 
                               Model& model, 
                               const finance::ComputationalFlags& flags)
                             : m_cboptionparams(cboptionparams), 
                               m_model(model),
                               m_flags(flags)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_cboptionparams.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );
}

AutoPtr<CBOptionNumOutput> CBOptionPricer::PriceNormal()
{
  // Construct the objects needed for pricing
  pricing::CBMeshManager cbmeshes(m_cboptionparams, m_model);

  CBOptionInstData cboptioninstdata(m_cboptionparams, m_model, cbmeshes);

  AutoPtr<CBOptionNumOutput> 
    pCBOptionNumOutput( new CBOptionNumOutput(m_cboptionparams) );

  CBOptionStepper cboptionstepper(cboptioninstdata, m_flags);

  CBOptionEngine 
    cboptionengine(m_cboptionparams, cbmeshes, cboptioninstdata,
                   cboptionstepper, *pCBOptionNumOutput);
  
  // Now initialize the objects
  if (    !m_cboptionparams.GetCalls()->HasMakeWhole()  
       && !m_cboptionparams.GetCalls()->HasNoticePeriod()
       && !m_cboptionparams.HasNewShare() )
  {
    cboptioninstdata.m_bComputeVega  = m_flags.GetComputeVega();
    // cbinstdata treatment
    cboptioninstdata.GetCBInstData()->m_bComputeVega = 
      m_flags.GetComputeVega();
  }

  cboptioninstdata.m_bComputeFugit = m_flags.GetComputeFugit();
  // cbinstdata treatment
  cboptioninstdata.GetCBInstData()->m_bComputeFugit = 
    m_flags.GetComputeFugit();

  m_cboptionparams.Init();

  cbmeshes.SetupMe();

  pCBOptionNumOutput->SetFinalSave(true);
  pCBOptionNumOutput->GetCBNumOutput()->SetFinalSave(true);

  pCBOptionNumOutput->GetComputationalFlags().SetComputeSurface
                                              ( m_flags.GetComputeSurface() );
  pCBOptionNumOutput->GetCBNumOutput()
    ->GetComputationalFlags().SetComputeSurface( m_flags.GetComputeSurface() );

  // Actually do the pricing and return the result
  cboptionengine.Run();
  
  return pCBOptionNumOutput;
}

// TODO: Path-dep ==> throw an exception for the moment.
AutoPtr<CBOptionNumOutput> CBOptionPricer::PricePathDep()
{
  ASSERT_MSG(!IsPathDependent(),
             "CB option with path-dependent clauses for "
             "the cb not treated yet");
 
  AutoPtr<CBOptionNumOutput> 
    pCBOptionNumOutput( new CBOptionNumOutput(m_cboptionparams) );

  return pCBOptionNumOutput;
} 

bool CBOptionPricer::IsPathDependent()
{
  return m_cboptionparams.HasPathDepCoCo() || 
         m_cboptionparams.HasPathDepCall();   
} 


} // namespace ihg

} // namespace ito33
