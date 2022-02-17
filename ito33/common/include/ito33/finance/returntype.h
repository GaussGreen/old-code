/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/returntype.h
// Purpose:     enum for different return calculation types
// Created:     2006/02/21
// RCS-ID:      $Id: returntype.h,v 1.2 2006/04/10 10:55:10 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/returntype.h
    @brief Enumeration of the type of return calculation (actual or log)
 */

#ifndef _ITO33_FINANCE_RETURNTYPE_H_
#define _ITO33_FINANCE_RETURNTYPE_H_

namespace ito33
{

namespace finance
{

/// type of the return calculation: actual or log
enum ReturnType
{
  Return_Actual,

  Return_Log

  #ifndef __CPP2ANY__
  , 

  /// @noexport
  Return_Max

  #endif

};

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_RETURNTYPE_H_
