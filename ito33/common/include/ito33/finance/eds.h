/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/eds.h
// Purpose:     financial EDS class 
// Created:     2005/01/26
// RCS-ID:      $Id: eds.h,v 1.8 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/eds.h
   @brief declaration of financial EDS class
 */

#ifndef _ITO33_FINANCE_EDS_H_
#define _ITO33_FINANCE_EDS_H_

#include "ito33/sharedptr.h"

#include "ito33/finance/derivative.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

  class ITO33_DLLDECL CashFlowStreamUniform;

  class DerivativeModifyingVisitor;

/**
   This class describes a Equity Default Swap (EDS).
 */
class ITO33_DLLDECL EDS : public Derivative
{

public:

  /**
     Ctor constructs a EDS from its recovery rate, spread stream, and the
     barrier.

     @param dRecoveryRate the recovery rate of the EDS
     @param pSpreadStream the spread stream of the EDS
     @param dBarrier The barrier of the EDS
   */
  EDS(double dRecoveryRate,
      const shared_ptr<CashFlowStreamUniform>& pSpreadStream,
      double dBarrier);

  // Default dtor is ok. CDS won't be derived

  /**
     Gets the recovery rate

     @return the recovery rate 
   */
  double GetRecoveryRate() const {  return m_dRecoveryRate; }  

  /**
     Gets the stream of spread payments.

     @return the stream of spreads
   */
  shared_ptr<finance::CashFlowStreamUniform> GetSpreadStream() const 
  { 
    return m_pSpreadStream; 
  }
  
  /**
     @internal
     @brief Replaces the existing spread stream

     @noexport
   */
  void SetSpreadStream(const shared_ptr<CashFlowStreamUniform>& pSpreadStream)
  {
    m_pSpreadStream = pSpreadStream;
  }

  /**
     Gets the barrier.

     @return The Barrier
   */
  double GetBarrier() const { return m_dBarrier; }

  /**
     Gets the maturity date.

     @return the maturity date 
   */
  Date GetMaturityDate() const;

  /**
     Sets the market price for the EDS.

     @param dPrice market price

     Note: Re-implementation of virtual base class since EDS prices can be
     positive or negative
   */
  void SetMarketPrice(double dPrice) 
  {
    DoSetMarketPrice(dPrice);
  }

  // implement base class pure virtuals
  void Visit(DerivativeVisitor& visitor) const;

  void Visit(DerivativeModifyingVisitor& visitor);
 
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// The recovery value in case of default
  double m_dRecoveryRate;

  /// The spread stream
  shared_ptr<CashFlowStreamUniform> m_pSpreadStream;

  /// The barrrier
  double m_dBarrier;

}; // class EDS


} //namespace finance

} //namespace ito33

#endif // #ifndef _ITO33_FINANCE_EDS_H_

