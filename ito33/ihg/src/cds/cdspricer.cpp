/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cds/cdspricer.cpp
// Purpose:     cds pricer class
// Author:      Nabil
// Created:     2003/10/29
// RCS-ID:      $Id: cdspricer.cpp,v 1.14 2006/08/21 14:50:08 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/computationalflags.h"

#include "ito33/pricing/cdsparams.h"

#include "ihg/model.h"
#include "ihg/cdsnumoutput.h"
#include "ihg/cdspricer.h"

namespace ito33
{

namespace ihg
{

typedef pricing::Engine
        <
          pricing::CDSParams,
          pricing::CDSMeshManager,
          CDSInstData, 
          CDSStepper, 
          CDSNumOutput
        > CDSEngine; 


CDSPricer::CDSPricer(pricing::CDSParams& params, 
                     Model& model, 
                     const finance::ComputationalFlags& flags)
                   : m_params(params), 
                     m_model(model),
                     m_flags(flags),
                     m_meshes(m_params, m_model),
                     m_instdata(m_params, model, m_meshes),
                     m_stepper(m_instdata, m_flags),
                     m_pNumOutput(new CDSNumOutput(params))
{
  if ( flags.GetAnalysisDate().IsValid() )
    m_params.SetAnalysisTime( GetDoubleFrom( flags.GetAnalysisDate() ) );

  m_instdata.m_bComputeVega = flags.GetComputeVega();

  m_pNumOutput->GetComputationalFlags().SetComputeSurface
                                        ( flags.GetComputeSurface() );
}

AutoPtr<CDSNumOutput> CDSPricer::Price()
{
  m_params.Init();

  m_meshes.SetupMe();

  CDSEngine
    engine(m_params, m_meshes, m_instdata, m_stepper, *m_pNumOutput);

  engine.Run();

  return m_pNumOutput;
}


} // namespace ihg

} // namespace ito33
