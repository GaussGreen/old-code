/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/pepsaveragingperiod.h
// Purpose:     define averaging period for peps
// Author:      ITO 33 Canada
// Created:     May 10, 2005
// RCS-ID:      $Id: pepsaveragingperiod.h,v 1.9 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/finance/bondlike/pepsaveragingperiod.h
    @brief declaration of the AveragingPeriod class    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_PEPS_AVERAGING_PERIOD_H_
#define _ITO33_FINANCE_BONDLIKE_PEPS_AVERAGING_PERIOD_H_

#include "ito33/dlldecl.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{


/**
   AveragingPeriod class for PEPS.

   @nocreate
 */
class ITO33_DLLDECL PEPSAveragingPeriod
{
public:

  /**
      Constructs an averaging period when the stock price is averaged.
   
      @param averageStartDate start of the average period
      @param averageEndDate   end of the average period
      @param nNumberOfSamplingAverages The number of sampling averages

      @noexport COM
   */
  static shared_ptr<PEPSAveragingPeriod> 
  CreateWithStock(Date averageStartDate, Date averageEndDate, 
    size_t nNumberOfSamplingAverages);

  /**
      Constructs an averaging period when the conversion ratio is averaged.
   
      @param averageStartDate start of the average period
      @param averageEndDate   end of the average period
      @param nNumberOfSamplingAverages The number of sampling averages

      @noexport COM
   */
  static shared_ptr<PEPSAveragingPeriod> 
  CreateWithConversionRatio(Date averageStartDate, Date AverageEndDate, 
    size_t nNumberOfSamplingAverages);


  // copy constructor is ok

  /// @name Modifiors
  //@{

  /**
      The current stock average.

      @param dCurrentStockAverage current average, can not be negative
      @param nNbSamplesUsed number of samples used in the current stock average
   */
  void SetCurrentStockAverage(double dCurrentStockAverage, 
    size_t nNbSamplesUsed);
  
  /**
      The current conversion ratio average.

      @param dCurrentConversionRatioAverage current average conversion ratio
      @param nNbSamplesUsed number of samples used in the current ratio average
   */
  void SetCurrentConversionRatioAverage(double dCurrentConversionRatioAverage, 
    size_t nNbSamplesUsed);

  //@}

  /// @name Acessors
  //@{

  /**
      Gets the average start date.

      @return the average start date
   */
  Date GetAverageStartDate() const
  {
    return m_averageStartDate;
  }

  /**
      Gets the averaging end date.

      @return the average end Date
   */
  Date GetAverageEndDate() const
  {
    return m_averageEndDate;
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
      The current running stock average.

      @return the current stock average
   */
  double GetCurrentStockAverage() const;

  /**
      The current running conversion ratio average.

      @return the current conversion average average
   */
  double GetCurrentConversionRatioAverage() const;

  /**
      Indicates if the averaging is based on the stock price.

      @return true if the averaging is based on the stock price
             false if the averaging is based on the conversion ratio
   */
  bool HasStockAveraging() const
  {
    return m_bHasStockAveraging;
  }

  //@}

  /**
      Validates the averaging period with the given valuation date.

      @param valuationDate the valuation date that the averaging period will
                           be validated against
   */
  void ValidateWith(Date valuationDate) const;

  /**
      @internal
     
      @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /**
      Constructs an Averaging Period.

      @param startDate start of the averaging period
      @param endDate   end of the averaging period
      @param nNumberOfSamplingAverages The number of sampling averages
   */
  PEPSAveragingPeriod(Date startDate, Date endDate, 
                  size_t nNumberOfSamplingAverages);

  /// start of averaging period
  Date m_averageStartDate;

  /// end of averaging period
  Date m_averageEndDate;

  /// Current average can be either stock average or conversion ratio
  double m_dCurrentAverage;

  /// indicates if averaging is based on the stock price
  bool m_bHasStockAveraging;

  /// The number of sampling averages
  size_t m_nNumberOfSamplingAverages;

  /// The number of samples used to compute the current average
  size_t m_nNbSamplesUsed;

}; // class PEPSAveragingPeriod


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_PEPS_AVERAGING_PERIOD_H_
