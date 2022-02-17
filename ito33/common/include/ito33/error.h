/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/error.h
// Purpose:     error codes returned by ITO33 functions
// Author:      Vadim Zeitlin
// Created:     16.12.02
// RCS-ID:      $Id: error.h,v 1.30 2006/06/15 10:01:49 zhang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/error.h
  @brief This file contains the Error class declaration and some common error
         codes.

  We represent the error codes by class objects instead of the enum elements to
  make adding new error codes less painful (with an enum approach everything
  must be recompiled whenever a new error code is added even if is only used in
  one place).
 */
#ifndef _ITO33_ERROR_H_
#define _ITO33_ERROR_H_

#include "ito33/error_base.h"
#include "ito33/erroroffsets.h"

namespace ito33
{
// ----------------------------------------------------------------------------
// Declares basic error codes, start from 10000
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_ERROR_CLASS( ITO33_BASE_ERROR_START );

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents basic error codes which start from 10000.

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
    
    For compatibility with the old values for the other products the error codes
    in this enum start at 10000, not at 1 as might have been expected.
  */
  class Error;
#endif


} // namespace ito33

#endif /* _ITO33_ERROR_H_ */

