/////////////////////////////////////////////////////////////////////////////
// Name:        parbond/parbondpricer.cpp
// Purpose:     parbond pricer class
// Author:      Nabil
// Created:     2005/05/20
// RCS-ID:      $Id: parbondpricer.cpp,v 1.4 2006/06/13 15:43:51 wang Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/parbondparams.h"

#include "ihg/model.h"
#include "ihg/parbondnumoutput.h"
#include "ihg/parbondpricer.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
        <
          pricing::ParBondParams,
          pricing::ParBondMeshManager,
          ParBondInstData, 
          ParBondStepper, 
          ParBondNumOutput
        > ParBondEngine; 


ParBondPricer::ParBondPricer(pricing::ParBondParams& params, 
                     Model& model, 
                     const finance::ComputationalFlags& flags)
                   : m_params(params), 
                     m_model(model),
                     m_flags(flags),
                     m_meshes(m_params, m_model),
                     m_instdata(m_params, model, m_meshes),
                     m_stepper(m_instdata, m_flags),
                     m_pNumOutput(new ParBondNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom(flags.GetAnalysisDate()) );

  m_instdata.m_bComputeVega = flags.GetComputeVega();

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                          (flags.GetComputeSurface());
}

AutoPtr<ParBondNumOutput> ParBondPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMe();

  ParBondEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
