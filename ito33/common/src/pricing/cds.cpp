/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cds.cpp
// Created:     2004/03/10
// RCS-ID:      $Id: cds.cpp,v 1.15 2006/08/21 14:43:14 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/cdslike.h"

#include "ito33/pricing/cds.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace pricing
{
  
CDS::CDS(const finance::CDSLike &cds)
       : m_dRecoveryRate( cds.GetRecoveryRate() ),
         m_pSpreadStream( cds.GetSpreadStream() )    
{  
  m_dMaturityTime = GetDoubleFrom( cds.GetMaturityDate() );
  
  // TODO: The cross currency case.
  if ( m_bIsCrossCurrency )
  {
    throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          TRANS("Cross currency cds not treated yet")
        );      
  }
  else
    m_pDerivativeCurve = cds.GetSessionData()->GetYieldCurve();
}

} // namespace pricing

} // namespace ito33
