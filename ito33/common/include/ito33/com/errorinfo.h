/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/errorinfo.h
// Purpose:     support for OLE Automation error handling
// Author:      Vadim Zeitlin
// Created:     05.02.03
// RCS-ID:      $Id: errorinfo.h,v 1.15 2006/06/07 23:48:59 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/errorinfo.h
    @brief OLE Automation error handling.

    All COM functions return an HRESULT which, in particular, is used as a
    success/failure indication. However this is not enough in most cases and so
    the callers may request more detailed error information (including at least
    a human readable error message) by asking for an IErrorInfo interface or,
    more precisely, for ISupportErrorInfo one (and only then for IErrorInfo
    itself).

    On the server side (the only one we're interested in so far) to support
    rich error information you need to derive your interface implementation
    class from ImplementCustomUnknown and include ISupportErrorInfo in the type
    list of extra interfaces instead of deriving ImplementUnknown, and call
    SetErrorInfo() whenever an error occurs.
 */

#ifndef _ITO33_COM_ERRORINFO_H_
#define _ITO33_COM_ERRORINFO_H_

#include "ito33/com/uuid.h"
#include "ito33/com/coclass.h"
#include "ito33/com/unknown_impl.h"

DEFINE_COM_TRAITS(IErrorInfo, IUnknown);
DEFINE_COM_TRAITS(ISupportErrorInfo, IUnknown);

namespace ito33
{

class Exception;

namespace COM
{

/**
    Sets the rich error information.

    For now we only support setting the error string but we should support
    the help topic and filename as well later.

    @todo add support for other ICreateErrorInfo methods

    @param iid the interface which generated the error
    @param e the exception object carrying information about the error
    @return S_OK if ok, otherwise an error occured and was ignored
 */
extern HRESULT DoSetErrorInfo(const IID& iid, const ito33::Exception& e);

/**
    As simple as it gets ISupportErrorInfo implementation.
 */
class ImplementSupportErrorInfo : public ImplementUnknown<ISupportErrorInfo>
{
public:
  /**
      Ctor creates an object which supports rich error info for the given
      interface.
   */
  ImplementSupportErrorInfo(const IID& iid) : m_iid(iid) { }

  /// return S_OK if the riid is the same as the IID of our interface
  STDMETHODIMP InterfaceSupportsErrorInfo(REFIID riid)
  {
    return riid == m_iid ? S_OK : S_FALSE;
  }

private:
  Uuid m_iid;
};

/**
    Specialization of ImplementInterface for ISupportErrorInfo.

    We add two things here: a pointer to an object supporting ISupportErrorInfo
    which we need internally and also a handy SetErrorInfo() function which
    should be caleld by the user code whenever an exception occurs,
 */
template <>
struct ImplementInterface<ISupportErrorInfo>
{
  /**
      The real implementation of ImplementInterface<ISupportErrorInfo>

      @sa ImplementInterface
   */
  template <class Base>
  class RealImpl : public Impl::DeriveFrom<Base>
  {
  public:
    /// default ctor
    RealImpl() { m_pSupportErrorInfo = NULL; }

    /// returns the interface pointer
    ImplementSupportErrorInfo *GetISupportErrorInfo(const IID& iid)
    {
      if ( !m_pSupportErrorInfo )
      {
        m_pSupportErrorInfo = new ImplementSupportErrorInfo(iid);
      }

      m_pSupportErrorInfo->AddRef();

      return m_pSupportErrorInfo;
    }

    /**
        Override QueryInterface() to indicate support for ISupportErrorInfo
     */
    STDMETHODIMP QueryInterface(REFIID iid, void **ppObj)
    {
      HRESULT hr = Base::QueryInterface(iid, ppObj);
      if ( hr == E_NOINTERFACE && iid == IID_ISupportErrorInfo )
      {
        if ( !ppObj )
        {
          hr = E_POINTER;
        }
        else // valid pointer provided, fill it
        {
          try
          {
            *ppObj = GetISupportErrorInfo(Traits<Interface>::Uuid());

            hr = S_OK;
          }
          catch ( std::bad_alloc& )
          {
            hr = E_OUTOFMEMORY;
          }
          catch ( ... )
          {
            hr = E_UNEXPECTED;
          }
        }
      }

      return hr;
    }

    /**
        Set COM error information from the given exception object.

        This method is only useful for the interfaces implementing IErrorInfo,
        but it doesn't hurt to have it for the others -- but calling it will
        result in compilation errors, of course.

        @param e the exception object containing the rich exception information
        @return COM error code
     */
    static HRESULT SetErrorInfo(const ito33::Exception& e)
    {
      return DoSetErrorInfo(Traits<Interface>::Uuid(), e);
    }

    /// dtor releases the interface pointer we hold
    ~RealImpl()
    {
      if ( m_pSupportErrorInfo )
        m_pSupportErrorInfo->Release();
    }


  private:
    /// the interface pointer (it is created on demand, @c NULL initially)
    ImplementSupportErrorInfo *m_pSupportErrorInfo;
  };
};

/**
    This class may be used as a base class instead of ImplementUnknown if
    support for IErrorInfo is desired.

    This is the simplest way to get support for rich error messages: simply
    inherit from ImplementUnknownAndErrorInfo<> instead of ImplementUnknown<>
    and call SetErrorInfo() when an error occurs.
 */
template <class Iface>
class ImplementUnknownAndErrorInfo : public ImplementCustomUnknown
                      <
                        Iface,
                        TYPE_LIST_1(ISupportErrorInfo)
                      >
{
};

/**
    The base class to be used for coclasses wishing to support IErrorInfo.

    @sa ImplementCoClass
 */
template <class Iface, class ImplClass>
class ImplementCoClassAndErrorInfo : public ImplementCoClass
                      <
                        Iface,
                        ImplClass,
                        TYPE_LIST_1(ISupportErrorInfo)
                      >
{
};

/// Should be used to handle C++ exceptions in COM methods
#define CATCH_EXCEPTION_AND_RETURN_COM_ERROR_CODE                             \
  catch ( const ::ito33::Exception& e )                                       \
  {                                                                           \
    return SetErrorInfo(e);                                                   \
  }                                                                           \
  catch ( std::bad_alloc& )                                                   \
  {                                                                           \
    return E_OUTOFMEMORY;                                                     \
  }                                                                           \
  catch ( ... )                                                               \
  {                                                                           \
    return E_UNEXPECTED;                                                      \
  }

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_ERRORINFO_H_
