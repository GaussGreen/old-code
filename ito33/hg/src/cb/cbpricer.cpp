/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/cb/cbpricer.cpp
// Purpose:     cb pricer class
// Created:     2005/04/11
// RCS-ID:      $Id: cbpricer.cpp,v 1.12 2006/06/15 16:06:45 dave Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/engine.h"
#include "ito33/pricing/cb.h"
#include "ito33/pricing/cblikeparams.h"
#include "ito33/pricing/pathdeppricer.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/callprovisions.h"
#include "ito33/pricing/cbcalls.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/hg/error.h"

#include "hg/model.h"
#include "hg/cbinstdata.h"
#include "hg/cbnumoutput.h"
#include "hg/cbstepper.h"
#include "hg/cbpricer.h"

extern const ito33::hg::Error ITO33_HG_PATHDEPENDENT;

namespace ito33
{

namespace hg
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

bool CBPricer::IsPathDependent()
{
  // unlike reset, we just need to query the params 
  return m_cbparams.IsPathDependent();
}

AutoPtr<CBNumOutput> CBPricer::Price()
{    
  CHECK_COND(!IsPathDependent(), ITO33_HG_PATHDEPENDENT);
  
  return PriceNormal();
}

AutoPtr<CBNumOutput> CBPricer::PriceNormal()
{
  // Construct the objects needed for pricing
  pricing::CBMeshManager cbmeshes(m_cbparams, m_model);

  CBInstData cbinstdata(m_cbparams, m_model, cbmeshes);

  cbinstdata.SetupFlags(m_flags);
  
  AutoPtr<CBNumOutput> pCBNumOutput( new CBNumOutput(m_cbparams) );

  CBStepper cbstepper(cbinstdata, m_flags);

  CBEngine cbengine(m_cbparams,cbmeshes,cbinstdata, cbstepper, *pCBNumOutput);
 
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


} // namespace hg

} // namespace ito33
