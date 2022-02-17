/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/bonderror.h
// Purpose:     error codes for bondlike namespace
// Author:      ZHANG Yunzhi
// RCS-ID:      $Id: bonderror.h,v 1.7 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/finance/bondlike/bonderror.h
  @brief This file contains the Error class declaration in bondlike namespace.

  Bondkind error codes begins from 26000.
 */
#ifndef _ITO33_FINANCE_BONDLIKE_ERROR_H_
#define _ITO33_FINANCE_BONDLIKE_ERROR_H_

#include "ito33/finance/error.h"

namespace ito33
{

namespace finance
{


// ----------------------------------------------------------------------------
// Declares basic error codes, start from ITO33_BONDLIKE_ERROR_START
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_GIVEN_ERROR_CLASS( BondError, ITO33_BONDLIKE_ERROR_START );

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents error codes related to bond-like contracts

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
    
    The error codes starts at 26000.
  */
  class BondError;
#endif

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_BONDLIKE_ERROR_H_

