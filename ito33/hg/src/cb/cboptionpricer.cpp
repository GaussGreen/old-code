/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cboptionpricer.cpp
// Purpose:     cb option pricer class
// Created:     2006/01/19
// RCS-ID:      $Id: cboptionpricer.cpp,v 1.6 2006/08/02 20:57:22 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cboptionparams.h"
#include "ito33/pricing/cboption.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/hg/error.h"

#include "hg/cbnumoutput.h"
#include "hg/cboptioninstdata.h"
#include "hg/cboptionnumoutput.h"
#include "hg/cboptionpricer.h"
#include "hg/cboptionstepper.h"
#include "hg/model.h"

extern const ito33::hg::Error ITO33_HG_PATHDEPENDENT;

namespace ito33
{

namespace hg
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

bool CBOptionPricer::IsPathDependent()
{
  return m_cboptionparams.HasPathDepCoCo() || 
         m_cboptionparams.HasPathDepCall();   
} 

AutoPtr<CBOptionNumOutput> CBOptionPricer::PriceNormal()
{
  // Construct the objects needed for pricing
  pricing::CBMeshManager cbmeshes(m_cboptionparams, m_model);

  CBOptionInstData cboptioninstdata(m_cboptionparams, m_model, cbmeshes);

  cboptioninstdata.SetupFlags(m_flags);

  AutoPtr<CBOptionNumOutput> 
    pCBOptionNumOutput( new CBOptionNumOutput(m_cboptionparams) );

  CBOptionStepper cboptionstepper(cboptioninstdata, m_flags);

  CBOptionEngine 
    cboptionengine(m_cboptionparams, cbmeshes, cboptioninstdata,
                   cboptionstepper, *pCBOptionNumOutput);

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

AutoPtr<CBOptionNumOutput> CBOptionPricer::Price()
{    
  CHECK_COND(!IsPathDependent(), ITO33_HG_PATHDEPENDENT);
  
  return PriceNormal();
}

} // namespace hg

} // namespace ito33
