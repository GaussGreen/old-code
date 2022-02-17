/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/option.h
// Purpose:     contracts class for option (backward)
// Created:     2004/02/11
// RCS-ID:      $Id: option.h,v 1.18 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/option.h
    @brief The declaration of the contracts class.

    The base class for contracts.  
 */

#ifndef _ITO33_PRICING_OPTION_H_
#define _ITO33_PRICING_OPTION_H_

#include "ito33/finance/option.h"
#include "ito33/finance/optionlike.h"

#include "ito33/pricing/optionliketype.h" 
#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


/// The declaration of the (backward) option contract class.
class Option : public Contract 
{
public:

  /**
     The ctor.
    
     @param option a reference to an object of type finance::Option

     @todo{Treat the cross currency case}
   */
  Option(const finance::Option & option);

  /**
     The ctor.
    
     @param optionLike a reference to an object of type finance::OptionLike

     @todo{Treat the cross currency case}
   */
  Option(const finance::OptionLike & optionLike);

  virtual ~Option() { }
  
  /**
     Gets the pricing level option type information.

     @return pricing::OptionType
   */
  OptionLikeType GetOptionType() const 
  { 
    return m_optionType; 
  }

  /**
     Indicates if the optionlike contract is American or European.

     @return finance::ExerciseType
   */
  finance::ExerciseType GetExerciseType() const
  { 
    return m_exerciseType; 
  }

  /**
     Returns the strike value.

     @return strike
   */
  double GetStrike() const 
  { 
    return m_dStrike; 
  }


protected:

  /**
     Populates the information for common parameters for option.

     @param optionLike contract
   */
  void GetOptionLikeData(const finance::OptionLike &optionLike);

  /// Option type at the pricing level, different from the financial level
  OptionLikeType m_optionType;

  /// Exercise type European or American
  finance::ExerciseType m_exerciseType;

  /// strike
  double m_dStrike;

}; // class Option;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_OPTION_H_
