/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/exoticoption/asianoption.h
// Purpose:     financial for Asian option class
// Author:      ITO33 Canada
// Created:     March 29, 2005
// RCS-ID:      $Id: asianoption.h,v 1.8 2006/06/07 11:36:19 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/exoticoption/asianoption.h
    @brief declaration of the financial asian option class
 */

#ifndef _ITO33_FINANCE_EXOTICOPTION_ASIANOPTION_H_
#define _ITO33_FINANCE_EXOTICOPTION_ASIANOPTION_H_

#include "ito33/date.h"

#include "ito33/finance/optionlike.h"

namespace ito33
{

namespace finance
{


/**
    AsianOption represents the financial aspects of an Asian option.
 */
class  ITO33_DLLDECL AsianOption : public OptionLike
{
public:
  /**
      Creates a Fixed strike Asian Option object.
   
      @param dFixedStrike The strike value of the option which must be positive
      @param maturityDate The maturity date of the option
      @param optionType The type of Asian option
      @param exerciseType The exercise type of the option
      @param averageStartDate Date on which the averaging starts
      @param nNumberOfSamplingAverages The number of sampling averages
   */
  AsianOption(double dFixedStrike,
         Date maturityDate,
         OptionType optionType,
         ExerciseType exerciseType,
         Date averageStartDate,
         size_t nNumberOfSamplingAverages);


  /**
      Creates a Floating strike Asian Option object.
   
      @param maturityDate The maturity date of the option
      @param optionType The type of the option
      @param exerciseType The exercise type of the option
      @param averageStartDate Date on which the averaging starts
      @param nNumberOfSamplingAverages The number of sampling averages
   */
  AsianOption(Date maturityDate,
         OptionType optionType,
         ExerciseType exerciseType,
         Date averageStartDate,
         size_t nNumberOfSamplingAverages);

  /**
      Empty virtual destructor.
   */
  virtual ~AsianOption() { }


  /**
      @name Modifiors.
   */
  //@{

  /**
      Sets computed average at the valuation date.
     
      These values are only relevant if the valuation date is past the start 
      of the sampling period of the contract.  
      
      The current average defaults to zero. 

      The number of samples used should be less than the number of
      sampling returns specified in the constructor.

      @param dCurrentAverage current average
      @param nNbSamplesUsed number of samples used in the current average
   */
  void SetCurrentAverage(double dCurrentAverage, 
                         size_t nNbSamplesUsed);


  /**
      Date on which the averaging ends. 

      By default this is the maturity date.

      @param averageEndDate Date on which the averaging ends
   */
  void SetAverageEndDate(Date averageEndDate);

  //@}

  /**
      @name Accessors 
   */
  //@{

  /**
      Gets the number of sampling averages.

      @return the number of sampling averages
   */
  size_t GetNbSamplingAverages() const
  {
    return m_nNumberOfSamplingAverages;
  }
 
  /**
      Returns the strike for fixed strike Asian options, otherwise
      throws an exception.

      @return strike
   */
  double GetFixedStrike() const;

  /**
      The current average, relevant only if the valuation 
      date differs from the average start date.

      @return current average
   */
  double GetCurrentAverage() const
  {
    return m_dCurrentAverage;
  }

  /**
      Gets the average start date.

      @return the average start date
   */
  Date GetAverageStartDate() const
  {
    return m_averageStartDate;
  }

  /**
      Date on which the averaging ends.

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
      Indicates if this contract has a fixed strike. 	  	 
 	  	  	 
      @return true/false if it is/is not a fixed strike Asian option 	  	 
   */ 	  	 	
  bool HasFixedStrike() const 	  	 
  { 
    if ( m_dFixedStrike > 0.0 ) 
      return true; 
    
    return false; 	
  }

  //@}

  virtual void Visit(DerivativeVisitor& visitor) const;

  XML::Tag Dump(XML::Tag& tagParent) const;


protected:

  /// The number of sampling averages
  size_t m_nNumberOfSamplingAverages;

  /// The number of samples used to compute the current average
  size_t m_nNbSamplesUsed;

  /// Strike
  double m_dFixedStrike;

  /// Current average
  double m_dCurrentAverage;

  /// Average start date of the Asian Option
  Date m_averageStartDate;

  /// Average end date
  Date m_averageEndDate;

}; // class AsianOption


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_EXOTICOPTION_ASIANOPTION_H_
