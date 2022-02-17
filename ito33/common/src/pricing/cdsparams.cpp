/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/cdsparams.cpp
// Author:      Wang
// Created:     2004/03/02
// RCS-ID:      $Id: cdsparams.cpp,v 1.19 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/problemtype.h"
#include "ito33/pricing/paymentevent.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/cdsparams.h"

using namespace ito33;
using namespace ito33::pricing;

using ito33::finance::CashFlowStream;

void CDSParams::Init()
{
  Params::Init();

  ConstructDividendEvents();

  // Add spread event
  CashFlowStream::const_iterator 
    pPayment = m_cds.GetSpreadStream().begin();

  for ( ; pPayment != m_cds.GetSpreadStream().end(); ++pPayment)
  {
    double dPaymentTime = GetDoubleFrom(pPayment->first);

    if ( !numeric::IsBefore(dPaymentTime, m_dValuationTime) )
    {     
      shared_ptr<Event>
        pEvent( new PaymentEvent(dPaymentTime, - pPayment->second) );

      m_eventManager.AddEvent(pEvent);
    }
  }
}


