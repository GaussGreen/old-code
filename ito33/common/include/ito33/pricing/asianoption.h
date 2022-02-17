/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/asianoption.h
// Purpose:     contracts class for asian option (backward)
// Author:      ITO33 Canada
// Created:     March 31, 2004
// RCS-ID:      $Id: asianoption.h,v 1.5 2006/06/02 18:08:16 yann Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/asianoption.h
    @brief The declaration of the contracts class for Asian options.  
 */

#ifndef _ITO33_PRICING_ASIANOPTION_H_
#define _ITO33_PRICING_ASIANOPTION_H_

#include "ito33/dateutils.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/pricing/option.h"


namespace ito33
{

namespace pricing
{


/// The declaration of the (backward) asian option contract class.
class AsianOption : public Option
{
public:

  /**
      The ctor.
    
      @param asianOption a reference to an Asian option
      @todo{Treat the cross currency case}
   */
  AsianOption(const finance::AsianOption & asianOption);

  virtual ~AsianOption() { }

  /**
      @name Accessors
   */
  //@{

  /**
      Gets end of the observation interval

      @return end of observation interval
   */
  double GetAverageEndTime() const
  {
    return m_dAverageEndTime;
  }

  /**
      Gets start of the observation interval

      @return start of observation interval
   */
  double GetAverageStartTime() const
  {
    return m_dAverageStartTime;
  }

  /**
      Helper function: Get the current average

      @return the current average
   */
  double GetCurrentAverage() const
  {
    return m_dCurrentAverage;
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
    if ( m_dStrike > 0.0 )
      return true;
    
    return false;
  }

  //@}


protected:

  /// The number of sampling averages
  size_t m_nNumberOfSamplingAverages;

  /// The number of samples used to compute the current average
  size_t m_nNbSamplesUsed;

  /// Begining of averaging period
  double  m_dAverageStartTime;
  
  /// End of averaging period
  double  m_dAverageEndTime;

  /// Current average
  double m_dCurrentAverage;

}; // class AsianOption;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ASIANOPTION_H_
