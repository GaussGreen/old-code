/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/afterstd.h
// Purpose:     Enable warnings disabled by ito33/beforestl.h
// Author:      Vadim Zeitlin
// Created:     19.12.02
// RCS-ID:      $Id: afterstd.h,v 1.5 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/afterstd.h
    @brief  This header must be included after including ito33/beforestd.h.

    See the comments in beforestd.h.
 */

#ifdef _MSC_VER
  #pragma warning(pop)

  // don't restore this one for VC6, it gives it in each try/catch which is a
  // bit annoying to say the least
  #if _MSC_VER >= 0x1300
    // unreachable code
    #pragma warning(default:4702)
  #endif // VC++ >= 7
#endif


