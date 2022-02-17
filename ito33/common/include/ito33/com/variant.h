/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/variant.h
// Purpose:     OLE Variant related stuff
// Author:      Vadim Zeitlin
// Created:     26.01.03
// RCS-ID:      $Id: variant.h,v 1.23 2005/09/14 12:36:56 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/variant.h
    @brief Helpers for working with VARIANTs.

    For now there is not much here but later we're going to have all kinds of
    wrappers for OLE Automation functions working with VARIANTs.
 */

#ifndef _ITO33_COM_VARIANT_H_
#define _ITO33_COM_VARIANT_H_

#include "ito33/debug.h"
#include "ito33/string.h"

// have to include these headers to define IUnknown
#include "ito33/win32/winwrap.h"
#include <ole2.h>

#include "ito33/win32/oledate.h"

#include "ito33/com/bstr.h"

namespace ito33
{

namespace COM
{

/**
    VarType exists only to allow mapping the C++ types to OLE VT_XXX types.

    VarType<T>::Value yields the VT_XXX code for the given type T.
 */
template <typename T>
struct VarType
{
  /*
      Implementation details note: what we want to do here is to ensure that
      VarType<T>::Value for "T *" where T is any interface will be VT_UNKNOWN
      but without defining VarType<T>::Value for all types which would
      compromise the type safety.

      Unfortunately VC++ doesn't handle constructions such as

          template <class T> struct IsPointer<T *> { ... };

      at all which forces us to use an ugly workaround for it. Maybe even
      more unfortunately, g++ doesn't handle the MSVC workaround (because it
      considers conversion from any type to "void *" allowed!!) so we can't
      use the same code, however ugly, for both of them.

      Finally we use the workaround only for MSVC and the simpler version
      assuming standard compliant compiler for all the others (although
      surely not all of them are going to grok it...).
   */
#ifdef _MSC_VER
private:
  // we need 2 types with different sizeof()s:
  //
  // sizeof(Yes) == 1
  typedef char Yes;
  //
  // sizeof(No) > 1 unless compiler uses built in compression ;-)
  struct No { char dummy[3]; };

  // we also need a function which returns different results depending on
  // whether we feed it an "IUnknown *" or anything else -- this is achieved
  // by overloading with the first function chosen for such pointers and the
  // second one (the worst match for everybody except g++!) for everything
  // else
  //
  // NB: const volatile don't hurt and make this work with all pointers
  static Yes IsIUnknownPtr(const volatile IUnknown *);
  static No IsIUnknownPtr(...);

  // finally we also need a function being able to provide an object of class
  // T (can't use just T() because the class might not have a default ctor)
  static T MakeT();

  // and now we can finally know if T is a pointer to IUnknown or not
  enum { isIUnknownPtr = sizeof(IsIUnknownPtr(MakeT())) == sizeof(Yes) };

public:

  enum { Value = isIUnknownPtr ? VT_UNKNOWN : VT_EMPTY };

  /**
      As we can't avoid defining Value for all types, even those which are
      not pointers, catch using "undefined" VarType<T> specialization here
   */
  COMPILE_TIME_ASSERT( Value != VT_EMPTY, UndefinedVarType );
#endif // _MSC_VER
};

#ifndef _MSC_VER

// NB: this is not as good as it might be because we only want to do this for
//     T deriving from IUnknown, really, but I don't have time to spend on this
//     now -- to be looked at if/when we use any other compiler

/// specialization of VarType<> for the pointer types for standard compliant
/// compilers
template <typename T>
struct VarType<T *>
{
  enum { Value = VT_UNKNOWN };
};

#endif // !_MSC_VER

// ----------------------------------------------------------------------------
// specializations of VarType for the standard types
// ----------------------------------------------------------------------------

// as we use fully qualified name inside DEFINE_TYPE_VARTYPE macro (as we must
// because it may be used by the user code outside of any namespace), we can't
// specialize it inside the ito33::COM namespace -- close it and reopen it
// below
}

}

/**
    This macro is used to define an association between a C type and a VARTYPE
    constant.

    Note that the correspondence is not (and can not) be one to one,

    @param type is a (primitive) C type
    @param vt is a COM VT_XXX VARTYPE constant
 */
#define DEFINE_TYPE_VARTYPE(type, vt)                                         \
  template <>                                                                 \
  struct ::ito33::COM::VarType<type>                                          \
  {                                                                           \
    enum { Value = vt };                                                      \
  }

DEFINE_TYPE_VARTYPE(bool, VT_BOOL);
DEFINE_TYPE_VARTYPE(short, VT_I2);
DEFINE_TYPE_VARTYPE(int, VT_I4);
DEFINE_TYPE_VARTYPE(long, VT_I4);

DEFINE_TYPE_VARTYPE(unsigned short, VT_UI2);
DEFINE_TYPE_VARTYPE(unsigned int, VT_UI4);
DEFINE_TYPE_VARTYPE(unsigned long, VT_UI4);

DEFINE_TYPE_VARTYPE(float, VT_R4);
DEFINE_TYPE_VARTYPE(double, VT_R8);
DEFINE_TYPE_VARTYPE(::ito33::Win32::OleDate, VT_DATE);

DEFINE_TYPE_VARTYPE(BSTR, VT_BSTR);
DEFINE_TYPE_VARTYPE(VARIANT, VT_VARIANT);

namespace ito33
{

namespace COM
{

// ----------------------------------------------------------------------------
// VarType comparison
// ----------------------------------------------------------------------------

// the same type can have several corresponding VT_XXX constants which are all
// compatible between them, check for this using IsCompatible::With()
//
// NB: to care for broken VC++ 6 we need to use a class and not a simple
//     function here as it would be impossible to call it explicitely otherwise
//     with that compiler (it doesn't support it) and workaround with dummy
//     argument used below wouldn't allow us to define a generic implementation
//     as there is no such thing as partial function specialization
template <typename T>
struct IsCompatible
{
  static bool With(VARTYPE vt)
  {
    return vt == VarType<T>::Value;
  }
};

// under Win32 (but not Win64!), VT_INT is the same as VT_I4
#if !defined(_WIN64)

template <>
struct IsCompatible<short>
{
  static bool With(VARTYPE vt)
  {
    return vt == VT_I2 || vt == VT_BOOL;
  }
};

template <>
struct IsCompatible<int>
{
  static bool With(VARTYPE vt)
  {
    return vt == VT_I4 || vt == VT_INT;
  }
};

#endif // !Win64

// VT_DATE is the same as double (VT_R8) in C++
template <>
struct IsCompatible<double>
{
  static bool With(VARTYPE vt)
  {
    return vt == VT_R8 || vt == VT_DATE;
  }
};

// ----------------------------------------------------------------------------
// Variant class
// ----------------------------------------------------------------------------

#if 0

/**
    A very simple class allowing to store and extract values in/from VARIANT.
 */
template <typename T>
class Variant
{
public:
  /// The type of our value
  typedef T Type;

  /// The VT_XXX constant correspondingto our type
  enum { VT = VarType<T>::Value };

  /// @name Constructors
  //@{

  /**
      Initializes the object with an existing raw VARIANT.

      This Variant becomes associated with the provided VARIANT.
   */
  Variant(VARIANT& var);

  //@}

  /// @name Operations
  //@{

  /**
      Clears an uninitialized variant object.

      Clear() should be used if the object had been already initialized or
      memory leaks may occur.
   */
  void Init() { ::VariantInit(&m_var); }

  /**
      Clears an initialized variant object.

      Init() must be used if the object has never been initialized or the
      program may crash!
   */
  void Clear() { ::VariantClear(&m_var); }

  //@}

  /// @name Accessors
  //@{

  /// Get our value, throws if it is not of the correct type
  T Get() const
  {
    if ( m_var.vt != VT )
    {
    }
  }

  //@}

private:
  /// the VARIANT we're associated with
  VARIANT& m_var;
};

#endif // 0

// ----------------------------------------------------------------------------
// CopyTo/FromVariant() and related stuff
// ----------------------------------------------------------------------------

/**
    This function initializes a VARIANT object with the given value.

    The VARIANT object is supposed to be already initialized and empty.

    @param var the VARIANT which is going to be initialized with the value
    @param value the value to put in the VARIANT
 */
template <typename T>
void CopyToVariant(VARIANT& var, T value);

#if defined(_MSC_VER) && _MSC_VER <= 1200
  /// Helper macro for CopyToVariant()
  #define ITO33_COPY_FROM_VAR_EXTRA_PARAM(T) , T *
#else
  #define ITO33_COPY_FROM_VAR_EXTRA_PARAM(T)
#endif // VC++ <= 6.0

/**
    Extract the value of given type from the VARIANT.

    If the variant type is incorrect, an exception is thrown. Moreover, it must
    match exactly, i.e. even though it is possible to convert int to double in
    C an attempt to extract a double from a VARIANT of type VT_I4 will throw.

    Note that as VC++ 6.0 doesn't support explicit template function
    specification, we must have a dummy parameter to allow the compiler to
    select the function automatically but we hide this ugliness from the user
    by using a macro below.

    The template parameter is the type to be extracted.

    @param var the VARIANT to extract the value from
    @return the value of the given type
 */
template <typename T>
T DoCopyFromVariant(const VARIANT& var ITO33_COPY_FROM_VAR_EXTRA_PARAM(T));

/**
    Macro hiding the exact implementation of DoCopyFromVariant.

    @param T the type to be extracted
    @param var the VARIANT containing the value
 */
#if defined(_MSC_VER) && _MSC_VER <= 1200
  #define CopyFromVariant(T, var) DoCopyFromVariant(var, (T *)0)
#else
  #define CopyFromVariant(T, var) DoCopyFromVariant<T>(var)
#endif

/// specialization of CopyToVariant for VT_UNKNOWN
template <>
inline
void CopyToVariant(VARIANT& var, IUnknown *pUnknown)
{
  var.vt = VT_UNKNOWN;
  if ( (var.punkVal = pUnknown) != NULL )
  {
    // we must do this before giving the pointer to the VARIANT because it
    // is going to Release() it later
    pUnknown->AddRef();
  }
}

/// and another for VT_DISPATCH: this one must be used with VB(A)
template <>
inline
void CopyToVariant(VARIANT& var, IDispatch *pDispatch)
{
  var.vt = VT_DISPATCH;
  if ( (var.pdispVal = pDispatch) != NULL )
  {
    // we must do this before giving the pointer to the VARIANT because it
    // is going to Release() it later
    pDispatch->AddRef();
  }
}

/**
    Macro to declare specializations of CopyTo/FromVariant() for simple types.

    @param type to define the specialization for
    @param field the name of the VARIANT structure field to save the value in
 */
#define DEFINE_COPY_VARIANT(type, field)                                    \
  template <> inline                                                        \
  void CopyToVariant(VARIANT& var, type value)                              \
  {                                                                         \
    var.vt = VarType<type>::Value;                                          \
    var.field = value;                                                      \
  }                                                                         \
                                                                            \
  template <> inline                                                        \
  type DoCopyFromVariant<type>(const VARIANT& var                           \
                ITO33_COPY_FROM_VAR_EXTRA_PARAM(type))                      \
  {                                                                         \
    ASSERT_MSG( var.vt == VarType<type>::Value,                             \
            "incorrect VARIANT type in CopyFromVariant!" );                 \
                                                                            \
    return var.field;                                                       \
  }

DEFINE_COPY_VARIANT(short, iVal)
DEFINE_COPY_VARIANT(int, lVal)         // it's VT_I4 which counts, not int!
DEFINE_COPY_VARIANT(long, lVal)

DEFINE_COPY_VARIANT(unsigned short, uiVal)
DEFINE_COPY_VARIANT(unsigned int, ulVal)   // as above
DEFINE_COPY_VARIANT(unsigned long, ulVal)

DEFINE_COPY_VARIANT(float, fltVal)
DEFINE_COPY_VARIANT(double, dblVal)
DEFINE_COPY_VARIANT(::ito33::Win32::OleDate, date);

DEFINE_COPY_VARIANT(BSTR, bstrVal);

// special versions for some types
template <> inline
void CopyToVariant(VARIANT& var, bool b)
{
  var.vt = VT_BOOL;
  var.boolVal = b ? VARIANT_TRUE : VARIANT_FALSE;
}

template <> inline
bool DoCopyFromVariant(const VARIANT& var)
{
  ASSERT_MSG( var.vt == VT_BOOL,
              "incorrect VARIANT type in CopyFromVariant<bool>!" );

  return var.boolVal == VARIANT_TRUE;
}

template <> inline
void CopyToVariant(VARIANT& var, const char *str)
{
  var.vt = VT_BSTR;
  var.bstrVal = BasicString(str).Detach();
}

template <> inline
void CopyToVariant(VARIANT& var, const std::string& s)
{
  CopyToVariant(var, s.c_str());
}

#undef ITO33_COPY_FROM_VAR_EXTRA_PARAM

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_VARIANT_H_

