/////////////////////////////////////////////////////////////////////////////
// Name:        pricing/parbondparams.cpp
// Author:      ZHANG
// Created:     2005/05/20
// RCS-ID:      $Id: parbondparams.cpp,v 1.3 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/pricing/problemtype.h"
#include "ito33/pricing/paymentevent.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/parbondparams.h"

using namespace ito33;
using namespace ito33::pricing;

using ito33::finance::CashFlowStream;

void ParBondParams::Init()
{
  Params::Init();

  ConstructDividendEvents();

  // Add spread event
  CashFlowStream::const_iterator 
    pPayment = m_parbond.GetSpreadStream()->begin();

  for ( ; pPayment != m_parbond.GetSpreadStream()->end(); ++pPayment)
  {
    double dPaymentTime = GetDoubleFrom(pPayment->first);

    if ( !numeric::IsBefore(dPaymentTime, m_dValuationTime) )
    {     
      shared_ptr<Event>
        pEvent( new PaymentEvent(dPaymentTime, pPayment->second) );

      m_eventManager.AddEvent(pEvent);
    }
  }
}


