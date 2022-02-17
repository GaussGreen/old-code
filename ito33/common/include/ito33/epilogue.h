/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/epilogue.h
// Purpose:     header to be included in the end of all public headers
// Author:      Vadim Zeitlin
// Created:     2005-09-29
// RCS-ID:      $Id: epilogue.h,v 1.2 2005/10/03 12:40:10 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/epilogue.h
    @brief Header to be included in the very end of all the public headers.

    This file restores the settings changed by ito33/prologue.h. It should only
    be included after the prologue header.
 */

#ifndef _ITO33_EPILOGUE_H_
#define _ITO33_EPILOGUE_H_

#ifdef _MSC_VER
  // disable the warning in case it was reenabled since ito33/prologue.h
  #pragma warning(disable: 4103)

  // restore default packing
  #pragma pack(pop)
#endif // _MSC_VER

#endif // _ITO33_EPILOGUE_H_


