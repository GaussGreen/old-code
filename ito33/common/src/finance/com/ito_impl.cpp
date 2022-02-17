/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/com/ito_impl.cpp
// Purpose:     implement main interface of the ito33 COM library
// Author:      Vadim Zeitlin
// Created:     2004-05-14
// RCS-ID:      $Id: ito_impl.cpp,v 1.92 2006/08/19 22:44:34 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/com/coclass.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/unknown_impl.h"
#include "ito33/com/c2a.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/showversion.h"

extern "C" const CLSID CLSID_ITO33;

#include "com/timefunction_impl.h"

#include "com/terms_impl.h"

#include "com/impliededsspreads_impl.h"
#include "com/impliedcdsspreads_impl.h"

#include "ito.h"

using namespace ito33;
using std::string;

DEFINE_COM_IFACE_FULL(ITO33, ITO33, IUnknown);

// ----------------------------------------------------------------------------
// ITO33Impl class declaration
// ----------------------------------------------------------------------------

class ITO33Impl : public COM::ImplementCoClass
                             <
                               IITO33,
                               ITO33Impl,
                               TYPE_LIST_2(ISupportErrorInfo, IDispatch)
                             >
{
public:
  // ctor and dtor
  ITO33Impl() { }
  virtual ~ITO33Impl() { }

  STDMETHODIMP get_Version(BSTR *version);

  // Methods in ITerms
#include "com/terms_impl_for_inclusion.h"

  STDMETHODIMP NewImpliedCDSSpreads(DATE contractingDate,
                                    DATE firstPaymentDate,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    double dRecoveryRate,
                                    IImpliedCDSSpreads **ppImpliedCDSSpreads);

  STDMETHODIMP NewImpliedEDSSpreads(DATE contractingDate,
                                    DATE firstPaymentDate,
                                    DayCountConvention dcc,
                                    Frequency freq,
                                    double dBarrier,
                                    double dRecoveryRate,
                                    IImpliedEDSSpreads **ppImpliedEDSSpreads);

private:
  
  ITO33Impl(const ITO33Impl&);
  ITO33Impl& operator=(const ITO33Impl&);
};

// ----------------------------------------------------------------------------
// definitions
// ----------------------------------------------------------------------------

DEFINE_COM_COCLASS_FULL(CLSID_ITO33, "ito33.ITO33", 1, ITO33Impl);

// ============================================================================
// Implementation for functions declared in IIerms interface
// ============================================================================

#define ClassImpl ITO33Impl 

#include "com/terms_impl_for_inclusion.cpp"

// ============================================================================
// Implementation for functions specific in IITO33 interface
// ============================================================================

STDMETHODIMP 
ITO33Impl::get_Version(BSTR *version)
{
  try
  {
    *version = C2A::COM::Translate<string>::To( ito33::finance::ShowVersion() );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ITO33Impl::NewImpliedCDSSpreads(DATE contractingDate,
                                DATE firstPaymentDate,
                                DayCountConvention dcc,
                                Frequency freq,
                                double dRecoveryRate,
                                IImpliedCDSSpreads **ppImpliedCDSSpreads)
{
  try
  {
    *ppImpliedCDSSpreads = new ImpliedCDSSpreadsImpl
                               (
                                 shared_ptr<finance::ImpliedCDSSpreads>
                                 (
                                   new finance::ImpliedCDSSpreads
                                       (
                                        C2A::COM::Translate<ito33::Date>::From(contractingDate),
                                        C2A::COM::Translate<ito33::Date>::From(firstPaymentDate),
                                        static_cast<Date::DayCountConvention>(dcc),
                                        static_cast<finance::Frequency>(freq),
                                        dRecoveryRate
                                       )
                                 )
                               );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}

STDMETHODIMP 
ITO33Impl::NewImpliedEDSSpreads(DATE contractingDate,
                                DATE firstPaymentDate,
                                DayCountConvention dcc,
                                Frequency freq,
                                double dBarrier,
                                double dRecoveryRate,
                                IImpliedEDSSpreads **ppImpliedEDSSpreads)
{
  try
  {
    *ppImpliedEDSSpreads = new ImpliedEDSSpreadsImpl
                               (
                                 shared_ptr<finance::ImpliedEDSSpreads>
                                 (
                                   new finance::ImpliedEDSSpreads
                                     (
                                       C2A::COM::Translate<ito33::Date>::From(contractingDate),
                                       C2A::COM::Translate<ito33::Date>::From(firstPaymentDate),
                                       static_cast<Date::DayCountConvention>(dcc),
                                       static_cast<finance::Frequency>(freq),
                                       dBarrier,
                                       dRecoveryRate
                                     )
                                 )
                               );

    return S_OK;
  }
  CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE
}
