/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/mandatoryconversiontype.h
// Purpose:     enum of Mandatory conversion type
// Author:      Wang
// Created:     2004/08/23
// RCS-ID:      $Id: mandatoryconversiontype.h,v 1.7 2004/10/05 09:13:37 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/mandatoryconversiontype.h
    @brief enum of PEPS-like call type    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_MANDATORYCONVERSIONTYPE_H_
#define _ITO33_FINANCE_BONDLIKE_MANDATORYCONVERSIONTYPE_H_

namespace ito33
{

namespace finance
{


/// Conversion type of Mandatory
enum MandatoryConversionType
{
  /** 
  The issuer calls the bond by giving a fixed number of shares specified
  by the Call conversion ratio.
*/
  MandatoryConversionType_FixedShare,

/** 
  The issuer calls the bond for an amount obtained by the PEPS payoff
formula.
*/
  MandatoryConversionType_VariableShare
  
  #ifndef __CPP2ANY__
  ,
  // @noexport
  MandatoryConversionType_Max
  #endif

}; // enum MandatoryConversionType 


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_MANDATORYCONVERSIONTYPE_H_
