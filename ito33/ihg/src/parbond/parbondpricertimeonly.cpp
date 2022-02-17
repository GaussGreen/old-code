/////////////////////////////////////////////////////////////////////////////
// Name:        parbond/parbondpricertimeonly.cpp
// Purpose:     time only parbond pricer class
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondpricertimeonly.cpp,v 1.3 2006/03/20 14:54:12 yann Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"

#include "ito33/pricing/parbondparams.h"

#include "ihg/model.h"
#include "ihg/parbondnumoutputtimeonly.h"
#include "ihg/parbondpricertimeonly.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(ihg::ParBondPricerTimeOnly);

} // namespace ito33


namespace ito33
{

namespace ihg
{


ParBondPricerTimeOnly::ParBondPricerTimeOnly(pricing::ParBondParams& params, 
                                     Model& model,
                                     const finance::ComputationalFlags& flags)
  : m_params(params), 
    m_model(model),
    m_flags(flags),
    m_meshes(m_params, m_model),
    m_instdata(m_params, model, m_meshes),
    m_stepper(m_instdata),
    m_pNumOutput(new ParBondNumOutputTimeOnly(params)),
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

AutoPtr<ParBondNumOutputTimeOnly> ParBondPricerTimeOnly::Price()
{
  m_params.Init();

  m_meshes.SetupMeTimeOnly();

  m_engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
