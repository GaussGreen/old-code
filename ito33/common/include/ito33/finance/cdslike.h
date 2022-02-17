/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/cdslike.h
// Purpose:     base financial cds-like class
// Created:     2006/05/19
// RCS-ID:      $Id: cdslike.h,v 1.5 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/cdslike.h
    @brief declaration of the base financial cds-like class.

    Normal CDS contracts, reference CDS contracts, and other cds-like
    contracts (if any) can be derived from this base class.
 */

#ifndef _ITO33_FINANCE_CDSLIKE_H_
#define _ITO33_FINANCE_CDSLIKE_H_

#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{

  class ITO33_DLLDECL CashFlowStreamUniform;

/**
   This class describes the financial aspects common to all
   Credit Default Swaps (CDSs).
 
   A Credit Default Swap (CDS) typically provides the buyer insurance against
   the risk of default of some company (usually the issuer of a bond).

   @nocreate
 */
class ITO33_DLLDECL CDSLike : public Derivative
{

public:
  // protected constructor

  /// Virtual destructor for base class
  virtual ~CDSLike() { }

  /**
     Gets the recovery rate.

     @return the recovery rate 
   */
  double GetRecoveryRate() const 
  { 
    return m_dRecoveryRate; 
  }  

  /**
     Gets the spread.

     @return the spread specified in the spread stream
   */
  virtual double GetSpread() const = 0;

  /**
     Gets the stream of spread payments.

     @return the stream of spread payments
   */
   virtual shared_ptr<CashFlowStreamUniform> GetSpreadStream() const = 0;

  /**
     Gets the premium.

     For a CDS-like contract, the premium is simply the market price.

     @return the premium (market price)
   */
  double GetPremium() const
  {
    return GetMarketPrice();
  }

  /**
     @internal
     @brief Gets the issue date.

     @return the issue date of the cds
     
     @noexport
   */
  Date GetIssueDate() const;

  /**
     Sets the market price for the CDS.

     Note: Re-implementation of virtual base class since CDS prices can be
     positive or negative

     @param dPrice market price
   */
  virtual void SetMarketPrice(double dPrice) 
  {
    DoSetMarketPrice(dPrice);
  }


protected:

  /**
     Constructor (only called by derived classes).
     
     All CDS contracts should specify a recovery rate. All derived
     classes should set the spread stream.

     @param dRecoveryRate the recovery rate of the CDS
   */
  CDSLike(double dRecoveryRate);

  /// The recovery value in case of default
  double m_dRecoveryRate;

  /// The spread stream specifying the CDS payments
  shared_ptr<CashFlowStreamUniform> m_pSpreadStream;

}; // class CDSLike


} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_CDSLIKE_H_
