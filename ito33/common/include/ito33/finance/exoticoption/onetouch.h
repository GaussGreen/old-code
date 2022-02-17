///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/exoticoption/onetouch.h                                
// Purpose:          Description of One touch                                                 
// Created:          2005/07/05                                            
// RCS-ID:           $Id: onetouch.h,v 1.7 2006/01/31 10:26:06 wang Exp $
// Copyright         (c) 2005 -  Trilemma LLP                                        //
///////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/exoticoption/onetouch.h

    @todo Add financial description of one touch.
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCH_H_
#define _ITO33_FINANCE_EXOTICOPTION_ONETOUCH_H_

#include "ito33/finance/barriertype.h"
#include "ito33/finance/rebatetype.h"
#include "ito33/finance/derivative.h"

namespace ito33
{

namespace finance
{


/**
    One touch instrument.
 */
class ITO33_DLLDECL OneTouch : public Derivative
{
public:
  
  /** 
      Ctor creates an one touch instrument.

      @param maturityDate The maturity date of the one touch
      @param dBarrier The barrier
      @param barrierType The type of the barrier (up/down).
      @param rebateType If the rebate is paid immediately or at maturity
   */
  OneTouch(Date maturityDate, double dBarrier, BarrierType barrierType, 
           RebateType rebateType);
 
  /**
      Gets the maturity date of the one touch.

      @return maturity date of the option
   */
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }

  /**
      Gets the barrier of the one touch.
    
      @return dBarrier The barrier of the one touch.
   */
  virtual double GetBarrier() const { return m_dBarrier; }

  /** 
      Gets the type of the barrier (up/down).
     
      @return the type of the barrier
   */
  BarrierType GetBarrierType() const { return m_barrierType; }

  /**
      If the rebate will be paid immediately or at maturity.
     
      @return true if the rebate will be paid immediately,
             false if the rebate will be paid at maturity
   */
  RebateType GetRebateType() const { return m_rebateType; }

  /**
      @internal
      @brief Change the barrier level of the one touch.

      @noexport
   */
  virtual void SetBarrier(double dBarrier) { m_dBarrier = dBarrier; }

  // implement base class pure virtuals
  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

protected:
 
  /// Protected ctor called by derived class
  OneTouch(Date maturityDate) : m_maturityDate(maturityDate) { }

  /// The maturity date
  Date m_maturityDate;      
  
  /// The barrier level
  double m_dBarrier;

  /// The type of barrier (up-and-out, down-and-out, etc.)
  BarrierType m_barrierType;

  /// The type of rebate (immediate, at maturity, etc.)
  RebateType m_rebateType;

}; // class OneTouch

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_ONETOUCH_H_
