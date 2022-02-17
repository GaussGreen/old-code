/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/mandatory.h
// Purpose:     contract class for mandatory (backward)
// Author:      Wang
// Created:     2004/08/16
// RCS-ID:      $Id: mandatory.h,v 1.16 2006/06/07 17:03:40 yann Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/mandatory.h
    @brief The declaration of the mandatory contract class.

    The base class for mandatory contracts.  
 */

#ifndef _ITO33_PRICING_MANDATORY_H_
#define _ITO33_PRICING_MANDATORY_H_

#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/mandatoryconversion.h"
#include "ito33/pricing/cbcalls.h"
#include "ito33/pricing/callfixedshare.h"
#include "ito33/pricing/callvariableshare.h"
#include "ito33/pricing/cblike.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL PERCSLike;
  class ITO33_DLLDECL GeneralizedPERCSLike;
  class ITO33_DLLDECL PEPSLike;
  class ITO33_DLLDECL ConvertibleLike;
}

namespace pricing
{

/// The declaration of the mandatory contract class.
class Mandatory : public CBLike 
{

public:

  /**
      Creates a mandatory by financial PERCS-like object.

      @param percs financial PERCS-like objet
   */
  Mandatory(const finance::PERCSLike& percs);

  /**
      Creates a mandatory by financial PEPS-like object.

      @param peps financial PEPS-like objet
   */
  Mandatory(const finance::PEPSLike& peps);

  /**
      Creates a mandatory by financial Generalized PEPS-like object.

      @param peps financial Generalized PEPS-like objet
   */
  Mandatory(const finance::GeneralizedPEPSLike& peps);

  /// default ctor calls the ctor of the base class
  Mandatory() : CBLike() { } 

  // default dtor is ok.
  
  /**
      Sets the calls.
   */
  void SetCalls(const CBCalls& calls) 
  { 
    m_calls = calls;  
  }

  // implement virtual function in base class
  CallProvisions* GetCalls()
  {
    if(m_callType == MandatoryCallType_FixedCash)
      return &m_calls;
    else if(m_callType == MandatoryCallType_VariableShare)
      return &m_callVariableShare;
    else
      return &m_callFixedShare;
  }

  ConversionProvisions* GetConversions() 
  { 
    return &m_conversions;
  }
 
  double GetClaim(double dTime, bool bPlus = true) const;

  /**
      Indicates if mandatory has an averaging period
      where the average is computed based on the stock
      price.

      @return true if mandatory has a stock based averaging period
   */
  bool HasStockAveraging() const
  {
    return m_bHasStockAveraging;
  }

  /**
      Indicates if mandatory has an averaging period.

      @return true if mandatory has averaging period
   */
  bool HasAveragingPeriod() const
  {
    return m_bHasAveragingPeriod;
  }

  /**
      Gets the average  start time.

      @return average start time
   */
  double GetAverageStartTime() const
  {
    return m_dAverageStartTime;
  }
  
  /**
      Gets the average  end time.

      @return observation end time
   */
  double GetAverageEndTime() const
  {
    return m_dAverageEndTime;
  }

  /**
      If the valuation date is in the middle of
      the av period, need to specify the path to save.

      @return current avg
   */
  double GetCurrentAverage() const
  {
    return m_dCurrentAverage;
  }
  
  /**
      The number of samples used to compute the current average.

      @return the number of samples used to compute the current average
   */
  size_t GetNbSamplesUsed() const
  {
    return m_nNbSamplesUsed;
  }

  
  /**
      Gets the number of sampling averages.

      @return the number of sampling averages
   */
  size_t GetNbSamplingAverages() const
  {
    return m_nNumberOfSamplingAverages;
  }

  /**
      Indicates that the Average of the stock should be used
      instead of the stock value for conversion.
    
      @param dStockAverage set the current stock average
   */
  void SetStockAverage(double dStockAverage)  
  {
    m_conversions.SetStockAverage(dStockAverage);
  }
   
  /**
      Indicates that the average conversion ratio must be used.

      @param dConversionRatio set the current conversion price
   */
  void SetConversionRatioAverage(double dConversionRatio)
  {
    m_conversions.SetConversionRatioAverage(dConversionRatio);
  }

  /**
      Helper function to get the mandatory conversion.

      @return pointer mandatory conversion
   */
  const MandatoryConversion *GetMandatoryConversion() const
  {
    return &m_conversions;
  }


private:

  enum MandatoryCallType
  {
    MandatoryCallType_FixedCash = 0,
    MandatoryCallType_FixedShare = 1,
    MandatoryCallType_VariableShare = 2,
    MandatoryCallType_FixedMax
  };

  /// call type
  MandatoryCallType m_callType;

  /// fixed share call provision
  CallFixedShare m_callFixedShare;

  /// variable share call provision
  CallVariableShare m_callVariableShare;
  
  /// the call provision for CallSchedule
  CBCalls m_calls;

  /// conversion provision
  MandatoryConversion m_conversions;

  /// Averaging observation start date
  double m_dAverageStartTime;

  /// Averaging observation end date
  double m_dAverageEndTime;

  /// Indicate if it is stock averaging
  bool m_bHasStockAveraging;

  /// Indicate that there is an averaging period
  bool m_bHasAveragingPeriod;

  /// current average
  double m_dCurrentAverage;

  /// The number of sampling averages
  size_t m_nNumberOfSamplingAverages;

  /// The number of samples used to compute the current average
  size_t m_nNbSamplesUsed;

}; // class Mandatory;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_MANDATORY_H_
