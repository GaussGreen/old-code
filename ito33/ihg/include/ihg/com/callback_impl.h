/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/com/callback_impl.h
// Purpose:     helper for calling ICallback::GetValues() from C++
// Author:      Vadim Zeitlin
// Created:     2004-06-29
// RCS-ID:      $Id: callback_impl.h,v 1.2 2004/10/04 18:04:04 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_COM_CALLBACK_H_
#define _IHG_COM_CALLBACK_H_

struct ICallback;

extern void
GetValuesFromCallback(ICallback *pCallback,
                      double dTime,
                      const double *pdS,
                      double *pdValues,
                      size_t nNbS);

#endif // _IHG_COM_CALLBACK_H_

