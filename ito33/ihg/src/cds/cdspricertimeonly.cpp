/////////////////////////////////////////////////////////////////////////////
// Name:        cds/cdspricertimeonly.cpp
// Purpose:     time only cds pricer class
// Author:      Wang
// Created:     2004/03/26
// RCS-ID:      $Id: cdspricertimeonly.cpp,v 1.13 2006/08/21 14:54:32 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/cdsparams.h"

#include "ihg/model.h"
#include "ihg/cdsnumoutputtimeonly.h"
#include "ihg/cdspricertimeonly.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(ihg::CDSPricerTimeOnly);

} // namespace ito33


namespace ito33
{

namespace ihg
{


CDSPricerTimeOnly::CDSPricerTimeOnly(pricing::CDSParams& params, 
                                     Model& model,
                                     const finance::ComputationalFlags& flags)
  : m_params(params), 
    m_model(model),
    m_flags(flags),
    m_meshes(m_params, m_model),
    m_instdata(m_params, model, m_meshes),
    m_stepper(m_instdata),
    m_pNumOutput(new CDSNumOutputTimeOnly(params)),
    m_engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput)
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  // Instdata needs only to compute the price. The greeks can be computed
  // from the price

  // Vega even is trivial needs to be done in numoutput
  m_pNumOutput->GetComputationalFlags().SetComputeVega
                                        ( flags.GetComputeVega() );

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );
}

AutoPtr<CDSNumOutputTimeOnly> CDSPricerTimeOnly::Price()
{
  m_params.Init();

  m_meshes.SetupMeTimeOnly();

  m_engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
