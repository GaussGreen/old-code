/////////////////////////////////////////////////////////////////////////////
// Name:        hg/cbutil.h
// Purpose:     Some helper functions for cb pricing with HG model
// Created:     2005/04/11
// RCS-ID:      $Id: cbutil.h,v 1.1 2005/04/12 13:43:40 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file hg/cbutil.h
   @brief Some helper functions for cb pricing with HG model
 */

#ifndef _HG_CBUTIL_H_
#define _HG_CBUTIL_H_

#include "ito33/hg/error.h"

// Exchangeable convertible is not yet supported. 
extern const ito33::hg::Error ITO33_HG_EXCHANGEABLE;

#endif // #ifndef _HG_CBUTIL_H_
