/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/parbond.cpp
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbond.cpp,v 1.2 2006/04/04 16:29:48 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/parbond.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace pricing
{
  
ParBond::ParBond(const finance::ParBond &parbond)
    : m_dRecoveryRate( parbond.GetRecoveryRate() ),
      m_pSpreadStream( parbond.GetCashFlowStream() )
{  
  m_dMaturityTime = GetDoubleFrom( parbond.GetMaturityDate() );
  
  // TODO: The cross currency case.
  if( m_bIsCrossCurrency )
  {
    throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          TRANS("Cross currency parbond not treated yet")
        );      
  }
  else
    m_pDerivativeCurve = parbond.GetSessionData()->GetYieldCurve();
}

} // namespace pricing

} // namespace ito33
