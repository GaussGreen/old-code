/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/dispatch.h
// Purpose:     Everything related to IDispatch support
// Author:      Vadim Zeitlin
// Created:     25.03.03
// RCS-ID:      $Id: dispatch.h,v 1.14 2006/03/23 22:34:04 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/dispatch.h
    @brief Everything needed for transparent IDispatch support.

    To implement a dual interface, i.e. one which supports both direct calls
    to its methods via "normal" COM calls and also calls via IDispatch for
    late-binding languages (such as VBScript/JScript, for example) you simply
    have to:
        - use IDispatch instead of IUnknown as the base class in IDL file
        - add IDispatch to the list of extra interfaces to support in your
          COM::ImplementCoClass specialization

    Everything else is done automatically.
 */

#ifndef _ITO33_COM_DISPATCH_H_
#define _ITO33_COM_DISPATCH_H_

#include "ito33/com/unknown_impl.h"

/// Define the traits for IDispatch
DEFINE_COM_TRAITS(IDispatch, IUnknown);

#ifdef _MSC_VER
  // 'this' : used in base member initializer list
  #pragma warning(disable:4355)
#endif // _MSC_VER

namespace ito33
{

namespace COM
{

/**
    This class implements the meat of SupportIDispatch.

    It is used only to avoid the template bloat: the code we have here doesn't
    get instantiated multiple times but is only generated once.
 */
class SupportIDispatch
{
public:
  /**
      Ctor creates the object implementing IDispatch for the given interface.

      The interface is specified by its iid which we receive as parameter.

      @param pUnknown the pointer to the object implementing iid
      @param iid the interface we are supposed to implement
   */
  SupportIDispatch(IUnknown *pUnknown, const Uuid& iid)
    : m_iidThis(iid)
  {
    m_pThis = pUnknown;
    m_pDispatch = NULL;
  }

  /**
      Dtor frees the IDispatch object we create internally.

      IDispatch is created on demand -- only when it is first needed. But
      once it is created we keep it to avoid recreating it all the time and
      we only destroy in our destructor.
   */
  ~SupportIDispatch()
  {
    if ( m_pDispatch )
      m_pDispatch->Release();
  }

  /// helper for ImplementInterface<IDispatch>::RealImpl::QueryInterface()
  HRESULT DoQueryInterface(REFIID iid, void **ppObj);

  /**
      @name Forwarders for IDispatch methods.

      Note that we assume that DoQueryInterface() had been called before and
      so that m_pDispatch was initialized.
   */
  //@{

  HRESULT GetTypeInfoCount(UINT *pctinfo)
  {
    IDispatch *pDispatch = GetDispatch();

    return pDispatch ? pDispatch->GetTypeInfoCount(pctinfo) : E_FAIL;
  }

  HRESULT GetTypeInfo(UINT iTInfo, LCID lcid, ITypeInfo **ppTInfo)
  {
    IDispatch *pDispatch = GetDispatch();

    return pDispatch ? pDispatch->GetTypeInfo(iTInfo, lcid, ppTInfo)
            : E_FAIL;
  }

  HRESULT GetIDsOfNames(REFIID riid,
             LPOLESTR *rgszNames,
             UINT cNames,
             LCID lcid,
             DISPID *rgDispId)
  {
    IDispatch *pDispatch = GetDispatch();

    return pDispatch ? pDispatch->GetIDsOfNames(riid, rgszNames, cNames,
                          lcid, rgDispId)
            : E_FAIL;
  }

  HRESULT Invoke(DISPID dispIdMember,
         REFIID riid,
         LCID lcid,
         WORD wFlags,
         DISPPARAMS *pDispParams,
         VARIANT *pVarResult,
         EXCEPINFO *pExcepInfo,
         UINT *puArgErr)
  {
    IDispatch *pDispatch = GetDispatch();

    return pDispatch ? pDispatch->Invoke(dispIdMember, riid, lcid, wFlags,
                      pDispParams, pVarResult,
                      pExcepInfo, puArgErr)
            : E_FAIL;
  }

  //@}

private:
  /// get the IDispatch pointer creating it if necessary
  IDispatch *GetDispatch();


  /// IUnknown pointer we provide IDispatch for
  IUnknown *m_pThis;

  /// the IDispatch pointer is created on first use, never access it directly
  IDispatch *m_pDispatch;

  /// the main interface of m_pThis
  Uuid m_iidThis;
};

/**
    Specialization of ImplementInterface for IDispatch.

    The extra data for IDispatch implementation consists of an IDispatch
    pointer which is created on demand when it is first needed using
    CreateStdDispatch().
 */
template <>
struct ImplementInterface<IDispatch>
{
  /**
      The real implementation of ImplementInterface<ISupportErrorInfo>

      @sa ImplementInterface
   */
  template <class Base>
  class RealImpl : public Impl::DeriveFrom<Base>
  {
  public:
    RealImpl() : m_impl(this, Traits<Interface>::Uuid())
    {
    }

    /**
        Query the extra interfaces supported by this object.

        If we're asked for IDispatch, even not IUnknown), return
        m_pDispatch pointer instantiating it (for iid interface) if
        necessary.
     */
    STDMETHODIMP QueryInterface(REFIID iid, void **ppObj)
    {
      HRESULT hr = Base::QueryInterface(iid, ppObj);
      if ( hr == E_NOINTERFACE )
      {
        hr = m_impl.DoQueryInterface(iid, ppObj);
      }

      return hr;
    }

    // implement IDispatch methods
    // ---------------------------

    /// simply forward to the real implementation
    STDMETHODIMP GetTypeInfoCount(UINT *pctinfo)
      { return m_impl.GetTypeInfoCount(pctinfo); }

    /// simply forward to the real implementation
    STDMETHODIMP GetTypeInfo(UINT iTInfo, LCID lcid, ITypeInfo **ppTInfo)
      { return m_impl.GetTypeInfo(iTInfo, lcid, ppTInfo); }

    /// simply forward to the real implementation
    STDMETHODIMP GetIDsOfNames(REFIID riid,
                 LPOLESTR *rgszNames,
                 UINT cNames,
                 LCID lcid,
                 DISPID *rgDispId)
      { return m_impl.GetIDsOfNames(riid, rgszNames, cNames,
                     lcid, rgDispId); }

    /// simply forward to the real implementation
    STDMETHODIMP Invoke(DISPID dispIdMember,
              REFIID riid,
              LCID lcid,
              WORD wFlags,
              DISPPARAMS *pDispParams,
              VARIANT *pVarResult,
              EXCEPINFO *pExcepInfo,
              UINT *puArgErr)
      { return m_impl.Invoke(dispIdMember, riid, lcid, wFlags,
                 pDispParams, pVarResult, pExcepInfo,
                 puArgErr); }

  private:
    SupportIDispatch m_impl;
  };
};

} // namespace COM

} // namespace ito33

#ifdef _MSC_VER
  #pragma warning(default:4355)
#endif // _MSC_VER

#endif // _ITO33_COM_DISPATCH_H_

