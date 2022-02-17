/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/ihg/error.h
// Purpose:     error codes for ihg namespace
// Author:      ZHANG Yunzhi
// RCS-ID:      $Id: error.h,v 1.4 2006/06/15 10:01:49 zhang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/ihg/error.h
    @brief This file contains the Error class declaration for ihg namespace
 */
#ifndef _ITO33_IHG_ERROR_H_
#define _ITO33_IHG_ERROR_H_

#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/error.h"
#include "ito33/erroroffsets.h"

namespace ito33
{

namespace ihg
{

// ----------------------------------------------------------------------------
// Declares IHG error codes, start from ITO33_IHG_ERROR_START
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_ERROR_CLASS( ITO33_IHG_ERROR_START );

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents error codes related to ihg namespace.

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
  */
  class Error;
#endif

}

} // namespace ito33


#endif /* _ITO33_IHG_ERROR_H_ */

