/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/init.h
// Purpose:     COM initialization
// Author:      Vadim Zeitlin
// Created:     2004-12-19
// RCS-ID:      $Id: init.h,v 1.1 2004/12/19 12:43:11 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/init.h
    @brief COM subsystem initialization for COM clients.
*/

#ifndef _ITO33_COM_INIT_H_
#define _ITO33_COM_INIT_H_

#include "ito33/com/exception.h"

namespace ito33
{

namespace COM
{

class Initialize
{
public:
    Initialize()
    {
      HRESULT hr = ::CoInitialize(NULL);
      if ( FAILED(hr) )
      {
        throw COM_EXCEPTION("CoInitialize", hr);
      }
    }

    ~Initialize() { ::CoUninitialize(); }
};

} // namespace ito33

} // namespace ito33

#endif // _ITO33_COM_INIT_H_
