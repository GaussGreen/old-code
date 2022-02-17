/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/stepdebug.h
// Purpose:     macro for activate/desactive step by step debug
// Author:      Wang
// Created:     2004/09/23
// RCS-ID:      $Id: stepdebug.h,v 1.2 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/stepdebug.h
   @brief macro for activate/desactive step by step debug

   The macro shouldn't be defined when doing a release. So please don't commit
   your temporary change on the macro.

   Caution: Define the macro when you know what you are doing! It shouldn't be
   defined when doing testsuite or calibration. And undefine the macro when you
   have finished debug.

   This is probably a temporary solution. If we are not using a preprocessor 
   for the compiler, it's because that we need to change several projets with 
   the actual way the solution is organized.
 */

#ifndef _ITO33_STEPDEBUG_H_
#define _ITO33_STEPDEBUG_H_

// #define STEPBYSTEP_DEBUG


// When the macro is defined, we normaly want to dump data into files
#ifdef STEPBYSTEP_DEBUG

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#endif // #ifdef STEPBYSTEP_DEBUG


#endif // #ifdef _ITO33_STEPDEBUG_H_
