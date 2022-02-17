/******************************************************************************
 * File.........: ito33/finance/optiontype.h
 * Purpose......: Definition of the enumeration of the type of option (put/call)
 * Author.......: laurence
 * Created......: 09/09/2003
 * RCS-ID.......: $Id: optiontype.h,v 1.9 2006/07/21 12:39:35 wang Exp $
 * Copyright....: (c) 2003 - 2006 Trilemma LLP
 ******************************************************************************/

/**
    @file ito33/finance/optiontype.h
    @brief Enumeration of the type of option (put/call) used by ITO33 projects
 */

#ifndef _ITO33_FINANCE_OPTIONTYPE_H_
#define _ITO33_FINANCE_OPTIONTYPE_H_

namespace ito33
{

namespace finance
{

/// type of the option: put or call
enum OptionType
{
  Option_Put,
  Option_Call

  #ifndef __CPP2ANY__
  , 
  Option_Digital,
  Option_Other,

  /// noexport
  OptionType_Max
  #endif

};

/**
    Checks if given option type is valid.

    @param optionType given value
    @return true if option type is valid, false if not 
    @noexport
 */
inline bool IsValidOptionType(OptionType optionType)
{
  return optionType == Option_Put ||
         optionType == Option_Call;
  
}

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_OPTIONTYPE_H_

