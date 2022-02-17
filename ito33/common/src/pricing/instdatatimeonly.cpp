/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/instdatatimeonly.cpp
// Author:      Wang
// Created:     2004/02/10
// RCS-ID:      $Id: instdatatimeonly.cpp,v 1.12 2005/05/29 08:42:45 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/numeric/schemetype.h"

#include "ito33/pricing/params.h"
#include "ito33/pricing/meshmanager.h"
#include "ito33/pricing/instdatatimeonly.h"

extern const ito33::Error ITO33_UNEXPECTED;

using ito33::pricing::InstDataTimeOnly;

void InstDataTimeOnly::ComputeSchemeWeights()
{
  switch(m_schemeType)
  {
  case ito33::numeric::SchemeType_Implicit:
  
  case ito33::numeric::SchemeType_CrankNicolson:
  
  case ito33::numeric::SchemeType_Explicit:
  
    m_dTimeWeight = m_dInverseTimeStep;
    m_dOldTimeWeight = m_dInverseTimeStep;
    m_dOldOldTimeWeight = 0.;

    break;

  case ito33::numeric::SchemeType_ThreeLevel:

    m_dOldTimeWeight = m_dInverseTimeStep + m_dOldInverseTimeStep;
    m_dOldOldTimeWeight = - (m_dTimeStep * m_dOldInverseTimeStep)
                        / (m_dTimeStep + m_dOldTimeStep);
    m_dTimeWeight = m_dOldTimeWeight + m_dOldOldTimeWeight;
    
    break;

  default:

    throw EXCEPTION_MSG
         (
          ITO33_UNEXPECTED,
          TRANS("Pricing::InstData class doesn't handle current scheme type")
         );
  }
}

void InstDataTimeOnly::UpdateBeforeStep()
{
  m_dOldTimeStep = m_dTimeStep;

  m_meshes.GetInstValues( m_dTimeStep,
                          m_dRate, m_dDerivativeRate, m_dForeignRate,
                          m_schemeType );

  m_dOldInverseTimeStep = m_dInverseTimeStep;
  
  m_dInverseTimeStep = 1. / m_dTimeStep; 
  
  ComputeSchemeWeights();

  // By default, no event got applied
  m_bHasEvent = false;
}

void InstDataTimeOnly::ApplyEvents()
{
  const Event* pEvent;

  while ( ( pEvent = m_params.GetBasicEvent() ) != 0 )
  {
    m_bHasEvent = true;

    ApplyEvent(pEvent);
  }
}

void InstDataTimeOnly::DoEvents()
{ 
  ApplyEvents();  
}
