//*************************-*-C++-*-***************************************************************
// file.............: ito33/cwrap.h
// purpose..........: definition of ITOCONVENTION macro
// author...........: ZHANG Yunzhi
// RCS-ID...........: $Id: cwrap.h,v 1.7 2004/10/05 09:13:35 pedro Exp $
//
// Copyright:   (c) 2003-2003 Trilemma LLP
//*************************************************************************************************

/**

  @file ito33/cwrap.h

  @brief definition of ITOCONVENTION macro
 */

#ifndef _ITO33_CWRAP_H_
#define _ITO33_CWRAP_H_

#ifdef WIN32
# ifdef ITO33_WIN32_DLL
#   define ITO33DLLCONV __declspec(dllexport)
# else
#   define ITO33DLLCONV
# endif

# define ITO33CALLCONV __stdcall
#else
# define ITO33DLLCONV
# define ITO33CALLCONV
#endif

/// Convention for all ITO33 public functions.
#ifndef ITOCONVENTION
# define ITOCONVENTION ITO33CALLCONV
#endif

/// Convention for C functions exported by the DLL
#ifndef ITO33_DLLEXP
# define ITO33_DLLEXP extern "C" ITO33DLLCONV
#endif

#endif // _ITO33_CWRAP_H_

