/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/hg/error.h
// Purpose:     HG error codes  
// RCS-ID:      $Id: error.h,v 1.3 2006/06/15 10:01:49 zhang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/hg/error.h
   @brief This file contains the Error class declaration for HG.
 */
#ifndef _ITO33_HG_ERROR_H_
#define _ITO33_HG_ERROR_H_

#include "ito33/error_base.h"
#include "ito33/erroroffsets.h"

namespace ito33
{

namespace hg
{

// ----------------------------------------------------------------------------
// Declares HG error codes, start from ITO33_HG_ERROR_START
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_ERROR_CLASS( ITO33_HG_ERROR_START );

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents error codes related to hg namespace.

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
  */
  class Error;
#endif



} // namespace hg

} // namespace ito33

#endif // #ifndef _ITO33_HG_ERROR_H_
