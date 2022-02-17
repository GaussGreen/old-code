/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/cds/cdsinstdatatimeonly.cpp
// Purpose:     time only cds instdata class 
// Author:      Wang
// Created:     2004/03/18
// RCS-ID:      $Id: cdsinstdatatimeonly.cpp,v 1.11 2005/02/03 14:02:12 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////
/**
    @file ihg/src/cds/cdsinstdatatimeonly.cpp
    @brief Implementation of InstData class for time only cds
 */

#include "ito33/useexception.h"

#include "ito33/numeric/numparams.h"

#include "ito33/pricing/event.h"
#include "ito33/pricing/cdsparams.h"
#include "ito33/pricing/cdsmeshmanager.h"

#include "ihg/model.h"
#include "ihg/cdsinstdatatimeonly.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace ihg
{


CDSInstDataTimeOnly::CDSInstDataTimeOnly(pricing::CDSParams& params,
                                         Model& model,
                                         pricing::CDSMeshManager& meshes)
                                       : InstDataTimeOnly(params, meshes),  
                                         m_cdsParams(params), 
                                         m_model(model), 
                                         m_cdsMeshes(meshes)
{
}

void CDSInstDataTimeOnly::Init()
{
}

void CDSInstDataTimeOnly::SetInitialValue()
{
  m_dPrice = m_dOldPrice = 0.;

  DoEvents();
}

void CDSInstDataTimeOnly::ApplyEvent(const pricing::Event* pEvent)
{
  double dS = 0.;

  if ( (pEvent->GetType()) == ito33::pricing::ET_Payment )
    pEvent->ApplyToPrice(&dS, &m_dPrice, 1);
}

void CDSInstDataTimeOnly::UpdateBeforeStep()
{
  InstDataTimeOnly::UpdateBeforeStep();

  double dS = 0.;

  m_model.GetHazardRates(&dS, &m_dHazardRate, 1);

  m_dRecoveryValue = m_cdsMeshes.GetRecoveryValue();
  
  m_dCoeZero = m_dRate + m_dHazardRate;
  m_dCoeConst = m_dHazardRate * m_dRecoveryValue;
  
  // Swap the prices 
  m_dOldOldPrice = m_dOldPrice;
  m_dOldPrice = m_dPrice;
}

void CDSInstDataTimeOnly::Run()
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
