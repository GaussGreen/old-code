/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/optionerror.h
// Purpose:     error codes for option like
// Author:      ITO 33 Canada
// Created:     March 30, 2005
// RCS-ID:      $Id: optionerror.h,v 1.4 2006/06/15 18:57:04 zhang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/finance/optionerror.h
  @brief This file contains the Error class declaration for option like.

  Option error codes begins from 25000.
 */
#ifndef _ITO33_FINANCE_OPTION_ERROR_H_
#define _ITO33_FINANCE_OPTION_ERROR_H_


#include "ito33/finance/error.h"

namespace ito33
{

namespace finance
{


// ----------------------------------------------------------------------------
// Declares option error codes, start from ITO33_OPTION_ERROR_START
// ----------------------------------------------------------------------------
//
ITO33_DECLARE_GIVEN_ERROR_CLASS(OptionError, ITO33_OPTION_ERROR_START);

// ----------------------------------------------------------------------------
// To generate humain readable document
// ----------------------------------------------------------------------------
#ifdef DOXYGEN
  /**
    This class represents error codes related to option contracts

    An error code has numeric value and the associated error message. Error
    codes for different products belong to different, disjoint ranges.
    
    The error codes at 25000.
  */
  class OptionError;
#endif

} // namespace finance

} // namespace ito33

#endif // _ITO33_FINANCE_OPTION_ERROR_H_

