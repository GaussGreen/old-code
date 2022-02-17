/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/sharedependentconversion.h
// Purpose:     conversion class when the conversion depends on the
//              applicable underlying spot price
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: sharedependentconversion.h,v 1.8 2006/08/19 22:18:03 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_PRICING_SHAREDEPENDENTCONVERSION_H_
#define _ITO33_PRICING_SHAREDEPENDENTCONVERSION_H_

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/pricing/cbconversions.h"
#include "ito33/pricing/conversionprovisions.h"
#include "ito33/numeric/mesh/roots.h"

namespace ito33
{
 
namespace finance
{
  class ShareDependentConversion;
}

namespace pricing
{


class ShareDependentConversion : public CBConversions
{
public:

   /** 
       Constructor.

       @param pShareDeConv share dependent conversion parameters
    */
  ShareDependentConversion
    ( const shared_ptr<finance::ShareDependentConversion>& pShareDeConv,
      Date valuationDate );
  
  // Default dtor is ok

  /// @name implement virtual functions
  //@{

  virtual bool GetGrossParities(const double* pdS, size_t nNbS, 
    const double* pdNewShares, double* pdValues) const;
  //@}


  ///@name methods for accessing ShareDependentConversion
  //@{

 
  /**
      Gets the base conversion ratio.

      @return The base conversion ratio
   */
  double GetBaseRatio() const 
  { 
    return m_dBaseRatio; 
  }

  /**
      Gets the incremental share factor.

      @return the incremental share factor
   */
  double GetIncrementalShareFactor() const 
  { 
    return m_dIncrementalShareFactor; 
  }

  /**
      The strike.

      @return the incremental strike
   */
  double GetStrike() const; 

  /**
      Gets the maximum conversion ratio.

      @return the maximum conversion ratio
   */
  double GetCapRatio() const 
  { 
    return m_dCapRatio; 
  }


  /**
      Indicate whether or not a reset time has been specified.

      @return true/false if a reset date has been set/not set
   */
  bool HasResetTime() const
  {
    return m_bHasResetTime;
  }

  /**
      Gets the reset time if it has been set.
 
      @return reset time
   */
  double GetResetTime() const 
  { 
    ASSERT( HasResetTime() );

    return m_dResetTime; 
  }

 //@}


  double GetTriggerRatio(size_t /*nIdx*/) const
  {
    return m_dBaseRatio;
  }

private:

  /// Base conversion ratio
  double m_dBaseRatio;

  /// The incremental share factor
  double m_dIncrementalShareFactor;

  /// The fixed asset value above which the conversion ratio increases
  double m_dFixedStrike;

  /// The maximum conversion ratio
  double m_dCapRatio;

  /// store the reset time
  double m_dResetTime;

  /// indicate it there is a reset time
  bool m_bHasResetTime;

}; // class ShareDependentConversion

} // namespace pricing

} // namespace ito33

#endif  //  _ITO33_PRICING_SHAREDEPENDENTCONVERSION_H_
