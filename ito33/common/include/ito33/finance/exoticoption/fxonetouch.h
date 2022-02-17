///////////////////////////////////////////////////////////////////////////////
// File:             ito33/finance/exoticoption/fxonetouch.h                                
// Purpose:          Description of One touch in FX                                              
// Created:          2006/01/31                                            
// RCS-ID:           $Id: fxonetouch.h,v 1.1 2006/01/31 10:26:06 wang Exp $
// Copyright         (c) 2006  Trilemma LLP                                        //
///////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/exoticoption/fxonetouch.h

    @todo Add financial description of one touch.
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_FXONETOUCH_H_
#define _ITO33_FINANCE_EXOTICOPTION_FXONETOUCH_H_

#include "ito33/finance/exoticoption/onetouch.h"

namespace ito33
{

namespace finance
{


/**
    One touch instrument in FX market, using market convention.
 */
class ITO33_DLLDECL FXOneTouch : public OneTouch
{
public:
  
  /** 
      Ctor creates an one touch instrument.

      @param maturityDate The maturity date of the one touch
      @param dBSBarrier Infers the barrier level via Black Scholes
      @param barrierType barrierType The type of the barrier (up/down).
      @param dVol The reference volatility     
   */
  FXOneTouch(Date maturityDate, double dBSBarrier, BarrierType barrierType,
             double dVol);
 
  /**
      Gets the barrier of the one touch from the Black-Scholes barrier.

      This function is costly so the computed barrier level is stored for
      subsequent calls.
    
      @return dBarrier The barrier of the one touch.
   */
  virtual double GetBarrier() const;

  /**
      The market quote.

      The market quote and Black-Scholes barrier price define the market price.

      @return The market quote
   */
  double GetMarketQuote() const;

  /**
      The market quote.

      The market quote and Black-Scholes barrier price define the market price.

      @param dMarketQuote The market quote
   */
  void SetMarketQuote(double dMarketQuote);


  /**
      Overload MarketPrice related functions so it can only be defined by the
      market quote and Black-Scholes barrier.
   */
  virtual void SetMarketPrice(double dPrice);

  // SetBarrier should not be called for FXOneTouch
  void SetBarrier(double dBarrier);

  // implement base class pure virtuals
  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// The Black-Scholes barrier price
  double m_dBSBarrier;

  /// The reference volatility used to compute the actual barrier
  double m_dReferenceVol;

  /// The market quote
  double m_dQuote;

}; // class FXOneTouch

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_FXONETOUCH_H_
