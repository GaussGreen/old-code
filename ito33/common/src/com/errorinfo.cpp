/////////////////////////////////////////////////////////////////////////////
// Name:        com/errorinfo.cpp
// Purpose:     IErrorInfo-related stuff
// Author:      Vadim Zeitlin
// Created:     05.02.03
// RCS-ID:      $Id: errorinfo.cpp,v 1.11 2005/02/04 23:15:16 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/error.h"
#include "ito33/exception.h"
#include "ito33/win32/exception.h"
#include "ito33/com/exception.h"
#include "ito33/com/ptr.h"
#include "ito33/com/errorinfo.h"

extern const ito33::Error ITO33_BAD_PARAM, ITO33_NULL_PTR, ITO33_OUT_OF_MEMORY;

using namespace std;
using namespace ito33::COM;

// ============================================================================
// implementation
// ============================================================================

HRESULT
ito33::COM::DoSetErrorInfo(const IID& iid, const ito33::Exception& e)
{
  // as this may be called from a catch handler, ensure that we don't throw
  // anything
  try
  {
    ICreateErrorInfo *pCreateErrInfoRaw;
    if ( ::CreateErrorInfo(&pCreateErrInfoRaw) == S_OK )
    {
      Ptr<ICreateErrorInfo> pCreateErrInfo(pCreateErrInfoRaw);
      if ( SUCCEEDED(pCreateErrInfo->SetGUID(iid)) &&
        SUCCEEDED(pCreateErrInfo->
                SetDescription(
                  String::MB2WC(e.GetErrorMessage()))) )
      {
        Ptr<IErrorInfo> pErrInfo(pCreateErrInfo);
        if ( pErrInfo )
        {
          ::SetErrorInfo(0, pErrInfo.Get());

          // skip the assert at the bottom by returning from here

          // the exception object may carry the system error code
          // (this is the case of COM::Exception), ask it first
          long rc = e.GetSystemErrorCode();

          if ( !rc )
          {
            // it doesn't, [try to] find the standard error code
            // corresponding to our error code
            static const struct ErrorMap
            {
              long errorITO, errorCOM;
            } errorsMap[] =
            {
              { ITO33_BAD_PARAM,      E_INVALIDARG    },
              { ITO33_NULL_PTR,       E_POINTER       },
              { ITO33_OUT_OF_MEMORY,  E_OUTOFMEMORY   },
            };

            // fall back if we don't find anything in the map
            rc = E_UNEXPECTED;

            for ( size_t n = 0; n < SIZEOF(errorsMap); ++n )
            {
              if ( e.GetErrorCode() == errorsMap[n].errorITO )
              {
                rc = errorsMap[n].errorCOM;
                break;
              }
            }
          }

          return rc;
        }
      }
    }
  }
  catch ( ... )
  {
    // nothing to do, really
  }

  FAIL( "SetErrorInfo() failed" );

  return E_UNEXPECTED;
}


