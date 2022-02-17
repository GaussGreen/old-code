/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/onetouch.cpp
// Purpose:     OneTouch pricing level implementation
// Created:     2005/07/04
// RCS-ID:      $Id: onetouch.cpp,v 1.3 2006/04/04 16:29:48 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/exoticoption/onetouch.h"

#include "ito33/pricing/onetouch.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace pricing
{
  

OneTouch::OneTouch(const finance::OneTouch& oneTouch)
  : Contract(oneTouch),
    m_dBarrier( oneTouch.GetBarrier() ),
    m_barrierType( oneTouch.GetBarrierType() ),
    m_bImmmediate( oneTouch.GetRebateType() == finance::Rebate_Immediate )
{  
  // TODO: The cross currency case.
  if ( m_bIsCrossCurrency )
  {
    throw EXCEPTION_MSG
          (
            ITO33_UNEXPECTED,
            TRANS("Cross currency one touch not treated yet")
          );      
  }
  else
    m_pDerivativeCurve = oneTouch.GetSessionData()->GetYieldCurve();
}


} // namespace pricing

} // namespace ito33
