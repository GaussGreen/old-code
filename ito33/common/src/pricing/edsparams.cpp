/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/edsparams.cpp
// Created:     2005/01/26
// RCS-ID:      $Id: edsparams.cpp,v 1.5 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/numeric/predicatetime.h"
#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/pricing/paymentevent.h"
#include "ito33/pricing/eventmanager.h"
#include "ito33/pricing/edsparams.h"

using namespace ito33;
using namespace ito33::pricing;

void EDSParams::Init()
{
  Params::Init();

  ConstructDividendEvents(numeric::ExtrapolationMode_Constant,
                          numeric::ExtrapolationMode_Linear,
                          numeric::InterpolationMethod_Quadratic);

  // Add spread event
  finance::CashFlowStreamUniform::const_iterator 
    pPayment = m_eds.GetSpreadStream().begin();

  for ( ; pPayment != m_eds.GetSpreadStream().end(); ++pPayment)
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
