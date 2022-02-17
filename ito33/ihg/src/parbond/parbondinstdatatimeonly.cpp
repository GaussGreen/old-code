/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/parbond/parbondinstdatatimeonly.cpp
// Purpose:     time only parbond instdata class 
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondinstdatatimeonly.cpp,v 1.1 2005/06/08 16:00:08 zhang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
/**
    @file ihg/src/parbond/parbondinstdatatimeonly.cpp
    @brief Implementation of InstData class for time only parbond
 */

#include "ito33/useexception.h"

#include "ito33/numeric/numparams.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/parbondparams.h"
#include "ito33/pricing/parbondmeshmanager.h"

#include "ihg/model.h"
#include "ihg/parbondinstdatatimeonly.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace ihg
{


ParBondInstDataTimeOnly::ParBondInstDataTimeOnly
                                       ( pricing::ParBondParams& params,
                                         Model& model,
                                         pricing::ParBondMeshManager& meshes)
                                       : InstDataTimeOnly(params, meshes),  
                                         m_parbondParams(params), 
                                         m_model(model), 
                                         m_parbondMeshes(meshes)
{
}

void ParBondInstDataTimeOnly::Init()
{
}

void ParBondInstDataTimeOnly::SetInitialValue()
{
  m_dPrice = m_dOldPrice = 1.;

  DoEvents();
}

void ParBondInstDataTimeOnly::ApplyEvent(const pricing::Event* pEvent)
{
  double dS = 0.;

  if ( (pEvent->GetType()) == ito33::pricing::ET_Payment )
    pEvent->ApplyToPrice(&dS, &m_dPrice, 1);
}

void ParBondInstDataTimeOnly::UpdateBeforeStep()
{
  InstDataTimeOnly::UpdateBeforeStep();

  double dS = 0.;

  m_model.GetHazardRates(&dS, &m_dHazardRate, 1);

  m_dRecoveryValue = m_parbondMeshes.GetRecoveryValue();
  
  m_dCoeZero = m_dRate + m_dHazardRate;
  m_dCoeConst = m_dHazardRate * m_dRecoveryValue;
  
  // Swap the prices 
  m_dOldOldPrice = m_dOldPrice;
  m_dOldPrice = m_dPrice;
}

void ParBondInstDataTimeOnly::Run()
{
  double dRHS, dMatrix;

  switch(m_schemeType)
  {
  case ito33::numeric::SchemeType_Implicit:
    dRHS = m_dOldPrice * m_dOldTimeWeight + m_dCoeConst;
    dMatrix = m_dTimeWeight + m_dCoeZero;
    break;

  case ito33::numeric::SchemeType_ThreeLevel:
    dRHS = m_dOldPrice * m_dOldTimeWeight
         + m_dOldOldPrice * m_dOldOldTimeWeight 
         + m_dCoeConst;

    dMatrix = m_dTimeWeight + m_dCoeZero;

    break;
  
  case ito33::numeric::SchemeType_CrankNicolson:
    {
    double dTheta = 0.5;
    dRHS = m_dOldPrice * (m_dOldTimeWeight - dTheta * m_dOldCoeZero)
         + dTheta * m_dOldCoeConst + (1. - dTheta) * m_dCoeConst;

    dMatrix = m_dTimeWeight + (1. - dTheta) * m_dCoeZero;

    break;
    }

  default:
    throw EXCEPTION_MSG
         (
           ITO33_UNEXPECTED, 
           TRANS("Stepper doesn't handle the current scheme type yet!")
         );  
  }

  m_dPrice = dRHS / dMatrix;

  m_dOldCoeConst = m_dCoeConst;

  m_dOldCoeZero = m_dCoeZero;
}


} // namespace ihg

} // namespace ito33
