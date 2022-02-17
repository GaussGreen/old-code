/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/paymentevent.cpp
// Purpose:     Implementation of the payment event
// Author:      Nabil
// Created:     2003/11/04
// RCS-ID:      $Id: paymentevent.cpp,v 1.7 2004/10/07 16:32:45 wang Exp $
// Copyright:   (c) 1999-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/pricing/paymentevent.h"

using ito33::pricing::PaymentEvent;

void PaymentEvent::ApplyToPrice(const double * /*pdS*/,
                                double *pdValues,
                                size_t nNbS) const
{
  for (size_t n = 0; n < nNbS; n++)
    pdValues[n] += m_dAmount;
}


void PaymentEvent::ApplyToGreek(const double * /*pdS*/,
                                double * /*pdValues*/,
                                size_t /*nNbS*/) const
{
  // do nothing
}
