/******************************************************************************
 * File.........: ito33/pricing/optionliketype.h
 * Purpose......: Definition of the enumeration of the type of option 
 * Author.......: ITO 33 Canada
 * Created......: April 8, 2005
 * RCS-ID.......: $Id: optionliketype.h,v 1.3 2006/05/26 13:34:31 nabil Exp $
 * Copyright....: (c) 2005 Trilemma LLP
 ******************************************************************************/

/**
  @file ito33/pricing/optionliketype.h
  @brief Enumeration of the type of option used by ITO33 projects
 */

#ifndef _ITO33_PRICING_OPTIONLIKETYPE_H_
#define _ITO33_PRICING_OPTIONLIKETYPE_H_

namespace ito33
{

namespace pricing
{

/// type of the option: 
enum OptionLikeType
{
  Option_Put,
  Option_Call,
  Option_Digital,
  Option_Other,
  AsianOption_FixedStrikeCall,         
  AsianOption_FixedStrikePut,          
  AsianOption_FloatingStrikeCall,     
  AsianOption_FloatingStrikePut 

};



} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_OPTIONLIKETYPE_H_

