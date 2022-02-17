/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/derivatives.cpp
// Purpose:     do the necessary for Derivatives class
// Created:     2005/05/17
// RCS-ID:      $Id: derivatives.cpp,v 1.11 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"

extern const ito33::finance::Error
  ITO33_INVALID_DERIVATIVE,
  ITO33_INVALID_DERIVATIVEWEIGHT,
  ITO33_INCONSISTENT_SESSION_DATA;

namespace ito33
{

namespace finance
{

void 
Derivatives::CheckDerivative(const shared_ptr<finance::Derivative>& pDerivative)
{
  CHECK_COND(pDerivative, ITO33_INVALID_DERIVATIVE);
}

void 
Derivatives::CheckWeight(double dWeight)
{
  CHECK_COND(dWeight > 0.0, ITO33_INVALID_DERIVATIVEWEIGHT);
}

void Derivatives::SetSessionData(const shared_ptr<SessionData>& pSessionData)
{
  Derivatives::Elements::const_iterator i;
  for ( i = m_elements.begin(); i != m_elements.end(); ++i)
    i->first->SetSessionData(pSessionData);
}

void Derivatives::ValidateSessionDatas()
{
  // list empty or contains only one element
  if ( m_elements.size() <= 1 )
    return;

  Derivatives::Elements::const_iterator iter = m_elements.begin();
  shared_ptr<SessionData> pSessionData(iter->first->GetSessionData());

  ++iter; // increment pointer to start past first element

  for (iter; iter != m_elements.end(); ++iter)
    CHECK_COND( pSessionData == iter->first->GetSessionData(),
                ITO33_INCONSISTENT_SESSION_DATA );
}

shared_ptr<SessionData> Derivatives::GetSessionData() const
{
  ASSERT( !m_elements.empty() );

  return m_elements.front().first->GetSessionData();
}

} // namespace finance

} // namespace ito33
