/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/transitionprobabilityparams.cpp
// Created:     2006/03/28
// RCS-ID:      $Id: transitionprobabilityparams.cpp,v 1.2 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/numeric/extrapolationmode.h"

#include "ito33/pricing/transitionprobabilityparams.h"

namespace ito33
{

namespace pricing
{

TransitionProbabilityParams::TransitionProbabilityParams
    (TransitionProbability& tp,    
     const shared_ptr<numeric::NumParams>& pNumParams,
     const shared_ptr<numeric::MeshParams>& pMeshParams)
   : Params(tp), m_tp(tp)
{
  m_pNumParams = pNumParams;
  m_pMeshParams = pMeshParams;

  m_pYieldCurve = make_ptr( new finance::YieldCurveFlat(0) );
  m_pYieldCurveForMesh = m_pYieldCurve;
  m_pForeignCurve = m_pYieldCurve;

  m_dValuationTime = 0;
  m_dStoppingTime = tp.GetMaturityTime();
  m_dAnalysisTime = -100.;

  m_dSpot = tp.GetStrike();
}

void TransitionProbabilityParams::Init()
{
  Params::Init();

  ConstructDividendEvents(numeric::ExtrapolationMode_Constant,
                          numeric::ExtrapolationMode_Linear,
                          numeric::InterpolationMethod_Quadratic);
}

} // namespace pricing

} // namespace ito33
