/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/optionlike.h
// Purpose:     financial option class
// Author:      ITO 33 Canada
// Created:     Apil 6, 2005
// RCS-ID:      $Id: optionlike.h,v 1.4 2005/06/29 10:23:24 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/optionlike.h
    @brief declaration of the financial optionlike class
 */

#ifndef _ITO33_FINANCE_OPTIONLIKE_H_
#define _ITO33_FINANCE_OPTIONLIKE_H_

#include "ito33/date.h"

#include "ito33/finance/exercisetype.h"
#include "ito33/finance/optiontype.h"

#include "ito33/finance/derivative.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{


/**
    OptionLike represents the financial aspects of an option.

    Classes corresponding to more complicated options (barrier and other
    exotics) derive from this one.

    @noexport COM
    @nocreate
 */
class ITO33_DLLDECL OptionLike : public Derivative
{
public:

  /// Empty virtual destructor for base class
  virtual ~OptionLike() { }

  /**
      @name Methods for accessing the OptionLike.
   */
  //@{

  /**
     Gets the type (Call/Put) of the option.

     @return the type of the option (Call/Put)
   */
  OptionType GetOptionType() const
  {
    return m_optionType;
  }

  /**
     Gets the maturity date of the option.

     @return maturity date of the option
   */
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }
 
  /**
     Gets the exercise type (American/European) of the option.

     @return the exercise type of the option (American/European)
   */
  ExerciseType GetExerciseType() const
  {
    return m_exerciseType;
  }

  //@}

  
protected:
  
  /**
     Creates an optionlike object.

     @param maturityDate the maturity date of the optionlike
     @param optionType the type of option put/call
     @param exerciseType the exercise type of the optionlike
   */
  OptionLike(Date maturityDate, OptionType optionType, 
             ExerciseType exerciseType);

  /**
     Writes myself to tag parent which can be derived options.

     @param tagParent tag of the derived class
   */
  void DumpMe(XML::Tag& tagParent) const;

  ///option type call or put
  OptionType m_optionType;

  ///maturity date of the option
  Date m_maturityDate;

  ///exercise type American or European
  ExerciseType m_exerciseType;


private:

  /// Don't allow generation of default ctor
  OptionLike();

}; // class OptionLike


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_OPTIONLIKE_H_
