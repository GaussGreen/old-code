/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/varianceswapterms.h
// Purpose:     class for financial terms of variance and volatility swaps
// Created:     2006/07/20
// RCS-ID:      $Id: varianceswapterms.h,v 1.2 2006/08/10 14:47:30 dave Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/varianceswapterms.h
    @brief financial class for terms of variance and volatility swaps
 */

#ifndef _ITO33_FINANCE_VARIANCESWAPTERMS_H_
#define _ITO33_FINANCE_VARIANCESWAPTERMS_H_

#include "ito33/date.h"

#include "ito33/finance/swaptype.h"
#include "ito33/finance/returntype.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
    Class groups some financial terms of the VarianceSwap.

    Volatility strike is not included since it may not yet be determined for
    a variance swap to be issued.
 */
class ITO33_DLLDECL VarianceSwapTerms
{
public:
  /**
      Creates a variance swap terms. 
      
      Assumes a nominal of 1.  Log returns are used by default.  If the start
      of the sampling period is T_0, return calculations begin at T_1. The 
      return at T_1 will use the spot values S_1 and S_0.

      @param maturityDate The maturity date of the swap
      @param swapType The type of swap (variance or volatility)
      @param startOfSamplingPeriod The start of the sampling period
      @param nNbSamplingReturns The number of sampling days/returns
   */
  VarianceSwapTerms(Date maturityDate,
                    SwapType swapType,
                    Date startOfSamplingPeriod,
                    size_t nNbSamplingReturns);

  // Default dtor is ok

  /// @name Modifiers for variance swap terms.
  //@{

  /**
      The return calculation method (log or actual).

      @param returnType return calculation method
   */
  void SetReturnType(ReturnType returnType);  

  /**
      The cap multiplier. 

      It must be in the [1, 10] interval.

      @param dCapMultiplier the cap multiplier
   */
  void SetCapMultiplier(double dCapMultiplier);

  /**
      The number of sampling returns per year. 
      
      Set to 252 by default.

      @param nAnnualReturnFrequency the number of sampling returns per year
   */
  void SetAnnualReturnFrequency(size_t nAnnualReturnFrequency);

  /**
      The up corridor barrier.

      The return on a sampling date is only included if the previous share
      price is above the barrier.  Set to zero by default.  The value set
      must be greater than zero, and less than the down corridor barrier 
      if the down corridor barrier was previously set.

      @param dUpCorridorBarrier the barrier defining the up corridor
   */
  void SetUpCorridorBarrier(double dUpCorridorBarrier);

  /**
      The down corridor barrier.

      The return on a sampling date is only included if the previous share
      price is below the barrier.  Set to zero by default.  The value set
      must be greater than zero, and greater than the up corridor barrier
      if the up corridor barrier was previously set.

      @param dDownCorridorBarrier the barrier defining the down corridor
   */
  void SetDownCorridorBarrier(double dDownCorridorBarrier);

  //@}

  /// @name Accessors for variance swap.
  //@{

  /**
      Gets the maturity date of the swap.

      @return maturity date of the swap
   */
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }

  /**
      Gets the swap type (variance or volatility).

      @return the swap type
   */  
  SwapType GetSwapType() const
  {  
    return m_swapType;   
  }

  /**
      The return calculation method (log or actual). 

      @return the return calculation method
   */  
  ReturnType GetReturnType() const
  {  
    return m_returnType;   
  }

  /**
      Gets the start of the sampling period.

      @return the start of the sampling period
   */  
  Date GetStartOfSamplingPeriod() const
  {  
    return m_startOfSamplingPeriod;   
  }

  /**
      Gets the number of sampling days/returns.

      @return the number of sampling days/returns
   */  
  size_t GetNbSamplingReturns() const
  {  
    return m_nNbSamplingReturns;   
  }

  /**
      The cap multiplier.

      @return the cap multiplier
   */
  double GetCapMultiplier() const
  {
    return m_dCapMultiplier;
  }
  
  /**
      The number of sampling returns per year. 
      
      Set to 252 by default.

      @return the number of sampling returns per year
   */
  size_t GetAnnualReturnFrequency() const
  {
    return m_nAnnualReturnFrequency;
  }

  /**
      The up corridor barrier.

      @return the up corridor barrier
   */
  double GetUpCorridorBarrier() const
  {
    return m_dUpCorridorBarrier;
  }

  /**
      The down corridor barrier.

      @return the down corridor barrier
   */
  double GetDownCorridorBarrier() const
  {
    return m_dDownCorridorBarrier;
  }

  //@}
  
  /**
      Writes myself to tag parent which can be variance swap 
      or variance swaption or some other class tag using VarianceSwapTerms.

      @param tagParent tag of the class using VarianceSwapTerms

      @noexport
   */
  void DumpMe(XML::Tag& tagParent) const;


private:

  /// maturity date of the swap
  Date m_maturityDate;

  /// swap type (variance or volatility)
  SwapType m_swapType;

  /// return calculation method (actual or log)
  ReturnType m_returnType;

  /// start of the sampling period
  Date m_startOfSamplingPeriod;

  /// number of sampling days/returns
  size_t m_nNbSamplingReturns;

  /// cap multiplier
  double m_dCapMultiplier;

  /// number of sampling returns per year. 
  size_t m_nAnnualReturnFrequency;

  /// up corridor barrier
  double m_dUpCorridorBarrier;

  /// down corridor barrier
  double m_dDownCorridorBarrier;

}; // class VarianceSwapTerms


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_VARIANCESWAPTERMS_H_
