/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/logcontract.cpp
// Purpose:     implementation of the pricing log contract class
// Created:     2006/07/18
// RCS-ID:      $Id: logcontract.cpp,v 1.1 2006/07/19 17:38:55 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/dateutils.h"

#include "ito33/finance/error.h"
#include "ito33/finance/logcontract.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/logcontract.h"

extern const ito33::finance::Error ITO33_START_SHARE_PRICE_UNDEFINED;

namespace ito33
{

namespace pricing
{

LogContract::LogContract(const finance::LogContract& logContract)
    : Contract(logContract),
      m_dT0( GetDoubleFrom(logContract.GetStartOfReturnPeriod()) ),
      m_dS0( logContract.GetStartSharePrice() )
{
  if ( m_dS0 < 0 )
  {
    CHECK_COND(   logContract.GetStartOfReturnPeriod()
                > logContract.GetSessionData()->GetValuationDate(),
               ITO33_START_SHARE_PRICE_UNDEFINED);

    m_dS0 = logContract.GetSessionData()->GetSpotSharePrice();
  }
}

} // namespace pricing

} // namespace ito33
