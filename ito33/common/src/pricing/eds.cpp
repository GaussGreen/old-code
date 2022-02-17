/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/eds.cpp
// Purpose:     EDS pricing level implementation
// Created:     2005/01/31
// RCS-ID:      $Id: eds.cpp,v 1.5 2006/04/04 16:29:48 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/eds.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace pricing
{
  
EDS::EDS(const finance::EDS& eds)
       : Contract(eds),
         m_dRecoveryRate( eds.GetRecoveryRate() ),
         m_spreadStream( *( eds.GetSpreadStream() ) ),
         m_dBarrier( eds.GetBarrier() )    
{  
  // TODO: The cross currency case.
  if( m_bIsCrossCurrency )
  {
    throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          TRANS("Cross currency eds not treated yet")
        );      
  }
  else
    m_pDerivativeCurve = eds.GetSessionData()->GetYieldCurve();
}

} // namespace pricing

} // namespace ito33
