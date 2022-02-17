/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/beforestd.h
// Purpose:     Disable some warnings in standard headers.
// Author:      Vadim Zeitlin
// Created:     19.12.02
// RCS-ID:      $Id: beforestd.h,v 1.8 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/beforestd.h
    @brief  This file is used to disable the warnings in standard headers.

    Unfortunately, when compiling at maximum warning level, the standard
    headers themselves may generate warnings -- and really lots of them. So
    before including them, this header should be included to temporarily
    suppress the warnings and after this the header ito33/afterstd.h should be
    included to enable them back again.

    Note that there are intentionally no inclusion guards in this file, because
    it can be included several times.
 */
#ifdef _MSC_VER
  // these warning have to be disabled and not just temporarily disabled
  // because they will be given at the end of the compilation of the current
  // source -- and there is absolutely nothing we can do about them

  // 'foo': unreferenced inline function has been removed
  #pragma warning(disable:4514)

  // 'function': function not inlined
  #pragma warning(disable:4710)

  // 'id': identifier was truncated to 'num' characters in the debug info
  #pragma warning(disable:4786)

  // we have to disable (and reenable in afterstd.h) this one because, even
  // though it is of level 4, it is not disabled by warning(push, 1) below for
  // VC7.1!

  // unreachable code
  #pragma warning(disable:4702)

  #pragma warning(push, 1)
#endif


