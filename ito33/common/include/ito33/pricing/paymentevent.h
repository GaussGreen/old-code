/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/paymentevent.h
// Purpose:     The payment event 
// Author:      Nabil
// Created:     2003/11/04
// RCS-ID:      $Id: paymentevent.h,v 1.9 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/paymentevent.h
    @brief The declaration of payment event class.

    At the date of payment event, the holder of instrument is paid a amount.
    The payment can be negative, that is, it is the holder to pay back.
*/

#ifndef _ITO33_PRICING_PAYMENTEVENT_H_
#define _ITO33_PRICING_PAYMENTEVENT_H_

#include "ito33/pricing/event.h"
 
namespace ito33
{

namespace pricing
{


/** Payment event class
*/
class PaymentEvent : public Event
{

public:

  /** Constructor

      @param dTime time of the payment
      @param dAmount the payment amount
  */
  PaymentEvent(double dTime, double dAmount)
    : Event(dTime), m_dAmount(dAmount)
  {
    m_eventType = ET_Payment;
  }

  /** Apply the payment event

      @param pdS the spot prices
      @param pdValues  array of prices 
      @param nNbS the size of the array
   */
  void ApplyToPrice(const double *pdS, double* pdValues, size_t nNbS) const;
    
  /**
    Function to apply the event to Greek values pdValues on grid pdS

    @param pdS the spot prices
    @param pdValues array of Greek values
    @param nNbS the size of the array
   */
  void ApplyToGreek(const double *pdS, double* pdValues, size_t nNbS) const; 
  
  /**
    Function to get the amount of the payment.

   */
  double GetAmount(){ return m_dAmount; }

protected:

  /// The payment amount.
  double m_dAmount;

};

} // namespace pricing

} // namespace ito33 

#endif // #ifndef _ITO33_PRICING_PAYMENTEVENT_H_


