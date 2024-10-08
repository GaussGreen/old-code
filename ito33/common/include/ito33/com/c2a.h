/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/c2a.h
// Purpose:     support file for C2A-generated code
// Author:      Vadim Zeitlin
// Created:     Nov 15, 2003
// RCS-ID:      $Id: c2a.h,v 1.33 2006/08/19 22:19:05 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/c2a.h
    @brief This file defines support functions used by C2A-generated code.

    Two functions To() and From() defined in C2A::COM namespace are used by the
    code generated by cpp2any to translate the COM types to C++ and from them,
    respectively. These functions should be defined for all types used in the
    exported API.
 */

#ifndef _ITO33_COM_C2A_H_
#define _ITO33_COM_C2A_H_

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
  // no kidding: VC++ 6 won't compile this file
  #error "Your compiler is broken, please update to at least 7.x"
#endif

#include "ito33/string.h"
#include "ito33/date.h"
#include "ito33/error.h"
#include "ito33/vector.h"
#include "ito33/autoptr.h"
#include "ito33/sharedptr.h"
#include "ito33/exception.h"

#include "ito33/com/bstr.h"
#include "ito33/com/traits.h"
#include "ito33/com/safearray.h"

// this definition is necessary in order to allow us to use the same type
// declarations for the type in both IDL and C++ sources
#define SAFEARRAY(type) SAFEARRAY *

extern const ito33::Error ITO33_NULL_PARAM;

namespace C2A
{

namespace COM
{


// ----------------------------------------------------------------------------
// default implementations
// ----------------------------------------------------------------------------

// these default implementations are necessary for ToVector<> and FromVector<>
// to work with all types, even primitive ones which don't need translation
// otherwise

/**
    Template used to translate COM type T to/from C++ type.
 */
template <typename T>
struct Translate
{
  typedef T CppType;
  typedef T COMType;

  /**
    Translate a COM object to C++ one.
    */
  static CppType From(COMType t) { return t; }

  /**
    Translate a C++ object to COM one.
    */
  static COMType To(CppType t) { return t; }
};

// ----------------------------------------------------------------------------
// pointers and objects
// ----------------------------------------------------------------------------

template <typename T>
struct Translate<T *>
{
  typedef T Iface;

  typedef Iface *COMType;

  // there is no CppType typedef because many C++ types map to COM interface
  // pointer type

  typedef typename ::ito33::COM::Traits<Iface> Traits;
  typedef typename Traits::Impl Impl;
  typedef typename Traits::Class Class;

  static ::ito33::shared_ptr<Class> From(Iface *p)
  {
    if ( !p )
    {
      using ito33::Exception;
      throw EXCEPTION_MSG(ITO33_NULL_PARAM, "Object is null.");
    }

    return static_cast<Impl *>(p)->GetImpl();
  }

  static Iface *To(const ::ito33::shared_ptr<Class>& ptr)
  {
    return ptr ? new Impl(ptr) : NULL;
  }

  static Iface *To(const ::ito33::AutoPtr<Class>& ptr)
  {
    return ptr ? To(ito33::shared_ptr<Class>(ptr.release())) : NULL;
  }

  static Iface *To(const ::std::auto_ptr<Class>& ptr)
  {
    return ptr.get() ? To(ito33::shared_ptr<Class>(
          const_cast<std::auto_ptr<Class>&>(ptr).release())
        )
      : NULL;
  }

  static Iface *To(const Class& obj)
  {
    return To(::ito33::shared_ptr<Class>(new Class(obj)));
  }

  static Iface *To(const Class *p)
  {
    return p ? To(*p) : NULL;
  }
};

// ----------------------------------------------------------------------------
// containers
// ----------------------------------------------------------------------------

template <typename T, typename U>
std::vector<T> ToVector(SAFEARRAY *pSA)
{
  using namespace ito33::COM;

  SafeArray<U> sa(pSA);
  SafeArrayAccessor<U> data(sa);

  const size_t count = data.GetCount();
  std::vector<T> vec(count);
  for ( size_t n = 0; n < count; ++n )
  {
    vec[n] = Translate<T>::From(data[n]);
  }

  return vec;
}

// this is for backwards compatibility with previous versions of cpp2any which
// generated an extra (and unnecessary) indirection for SAFEARRAY()s only,
// remove it as soon as the old version is not used any longer
template <typename T, typename U>
std::vector<T> ToVector(SAFEARRAY **ppSA)
{
  return ToVector<T,U>(*ppSA);
}

template <typename T>
SAFEARRAY *FromVector(const std::vector<T>& vec)
{
  using namespace ito33::COM;

  typedef typename Translate<T>::COMType U;

  const size_t count = vec.size();
  SafeArray<U> array(count);
  SafeArrayAccessor<U> arrayData(array);
  for ( size_t n = 0; n < count; ++n )
  {
    arrayData[n] = Translate<T>::To(vec[n]);
  }

  return array.Get();
}

// ----------------------------------------------------------------------------
// booleans
// ----------------------------------------------------------------------------

template <>
struct Translate<bool>
{
  typedef bool CppType;
  typedef VARIANT_BOOL COMType;

  static bool From(COMType b) { return b == VARIANT_TRUE; }
  static COMType To(bool b) { return b ? VARIANT_TRUE : VARIANT_FALSE; }
};

// ----------------------------------------------------------------------------
// dates
// ----------------------------------------------------------------------------

template <>
struct Translate<ito33::Date>
{
  typedef ito33::Date CppType;
  typedef ito33::Win32::OleDate COMType;

  static ito33::Date From(COMType d)
    { ito33::Date dt; dt.SetOleDate(d); return dt; }
  static COMType To(const ito33::Date& d)
    { return d.GetOleDate(); }
};

// ----------------------------------------------------------------------------
// string translations
// ----------------------------------------------------------------------------

template <>
struct Translate<std::string>
{
  static std::string From(BSTR bstr)
  {
    if ( !bstr )
    {
      using ito33::Exception;
      throw EXCEPTION_MSG(ITO33_NULL_PARAM, "String is null.");
    }

    return std::string(ito33::String::WC2MB(bstr));
  }

  static BSTR To(std::string str)
  {
    return ito33::COM::BasicString(str.c_str()).Detach();
  }
};

// TODO: support for C strings (char *)

} // namespace COM

} // namespace C2A

#endif // _ITO33_COM_C2A_H_
