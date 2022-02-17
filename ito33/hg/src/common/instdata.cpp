/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/instdata.cpp
// Purpose:     implementation for HG base instdata class
// Created:     2005/01/13
// RCS-ID:      $Id: instdata.cpp,v 1.5 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/event.h"

#include "hg/model.h"
#include "hg/instdata.h"
#include "hg/payoff.h"

namespace ito33
{

namespace hg
{


InstData::InstData(pricing::Params& params,
                   Model& model,
                   pricing::MeshManager& meshes)
                 : pricing::InstData(params, meshes), m_model(model)
{
  m_nNbRegimes = m_model.GetNbRegimes();
}

void InstData::ApplyEvent(const ito33::pricing::Event *pEvent)
{
  double* pdPrices = m_pdPrices.Get();

  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
    pEvent->ApplyToPrice(m_pdS, pdPrices + nIdxR * m_nNbS, m_nNbS);
}

void InstData::ApplyBoundaryConditionToSensitivityRHS(double* pdRHS) const
{
  for (size_t nIdxR = 0; nIdxR < m_nNbRegimes; nIdxR++)
  {
    double* pdRHSTmp = pdRHS + nIdxR * m_nNbS;

    if (m_BoundaryCondition.GetLeftType() == numeric::BCType_Dirichlet )
      pdRHSTmp[0] = 0;

    if ( m_BoundaryCondition.GetRightType() == numeric::BCType_Dirichlet )
      pdRHSTmp[m_nNbS - 1] = 0;
  }
}

void InstData::UpdateBeforeStep()
{
  pricing::InstData::UpdateBeforeStep();

  // Copy m_pdOldprice inside m_pdPrice.
  // This is done in order to have consistent valid
  // initial values for the matrix solve
  // Note that for greeks the copy is not necessarly since
  // we are doing a direct solve
  memcpy(m_pdPrices.Get(), m_pdOldPrices.Get(), m_nNbX);
}

} // namespace hg

} // namespace ito33
