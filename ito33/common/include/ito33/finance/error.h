/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/error.h
// Purpose:     error codes of finance namespace
// Author:      ZHANG Yunzhi
// Created:     16.12.02
// RCS-ID:      $Id: error.h,v 1.1 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/error.h
    @brief This file contains the finance Error class declaration and also
           the definition of error ranges of sub-modules, such as option,
           bond-like etc..
 */
#ifndef _ITO33_FINANCE_ERROR_H_
#define _ITO33_FINANCE_ERROR_H_

#include "ito33/error_base.h"
#include "ito33/erroroffsets.h"

namespace ito33
{

namespace finance
{

// ----------------------------------------------------------------------------
// Declares basic error codes, start from ITO33_FINANCE_ERROR_START
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_ERROR_CLASS( ITO33_FINANCE_ERROR_START );

/// option error range goes from base+5000 to base+5999
#define ITO33_OPTION_ERROR_START (ITO33_FINANCE_ERROR_START + 5000)

/// bond-like error range goes from base+6000 to base+6999
#define ITO33_BONDLIKE_ERROR_START (ITO33_FINANCE_ERROR_START + 6000)

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents basic error codes which start from 20000.

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
    
    The finance error codes start at 20000.
  */
  class Error;
#endif

} // namespace finance

} // namespace ito33

#endif /* _ITO33_ERROR_H_ */

