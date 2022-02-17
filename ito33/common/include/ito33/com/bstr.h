/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/bstr.h
// Purpose:     definition of COM string class
// Author:      Vadim Zeitlin
// Created:     02.04.03
// RCS-ID:      $Id: bstr.h,v 1.8 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/bstr.h
    @brief A simple COM string class.

    COM uses @c BSTR type which is simply a typedef for a pointer to @c
    wchar_t, i.e. a Unicode terminated string. However these strings are Pascal
    strings (length-prefixed) and they also must be allocated (and freed) by
    the system functions and so working with them directly is not as easy as it
    could be. So we have a small class making it easy.
 */

#ifndef _ITO33_COM_BSTR_H_
#define _ITO33_COM_BSTR_H_

#include "ito33/debug.h"
#include "ito33/string.h"

#include "ito33/win32/winwrap.h"

#include <ole2.h>

namespace ito33
{

namespace COM
{

/**
    A very simple wrapper around BSTR.

    @todo: - more ctors (create uninitialized string of given length...)
           - more assignment operators
           - comparisons
 */
class BasicString
{
public:
  /** @name Constructors and dtor */
  //@{

  /**
      Default ctor, the string is created uninitialized.

      You must assign something to it before using it.
   */
  BasicString();

  /**
      Creates a string containing the copy of the data of the given string.

      The BasicString object copies the data and so the string passed in may
      be deleted or modified without affecting this object.

      This ctor may throw an std::bad_alloc exception.

      @param sz a NUL-terminated, non NULL, ASCII string
   */
  BasicString(const char *sz);

  /**
      Copy ctor.
   */
  BasicString(const BasicString& str);

  /**
      Ctor from a raw BSTR.

      We don't take ownership of the string pointer but copy it. If you want
      to take ownership of it instead, use Attach().
   */
  BasicString(const wchar_t *bstr);

  /**
      Assignment operator from another BasicString object.
   */
  BasicString& operator=(const BasicString& str);

  /**
      Assignment operator from a raw BSTR.

      We don't take ownership of the string pointer but copy it. If you want
      to take ownership of it instead, use Attach().
   */
  BasicString& operator=(const wchar_t *bstr);

  /**
      Dtor frees the string.

      Note that the dtor is not virtual, this class shouldn't be used
      polymorphically (nor is it supposed to be derived from anyhow).

      If you don't want the string to be freed, use Detach().
   */
  ~BasicString();

  //@}

  /** @name Operations */
  //@{

  /**
      Attach to the given string, i.e. take ownership of it.

      The object takes ownership of the given string (freeing its old
      contents) and will free it in its dtor.
   */
  void Attach(BSTR bstr);

  /**
      Detach the string pointer from this object and return it.

      Detaching means that the pointer doesn't belong to this object any
      more, in particular it won't be freed in our dtor and so can be stored
      elsewhere safely.
   */
  BSTR Detach();

  /**
      Clear the string.

      This simple frees the current data, if any, and restes the string back
      to uninitialized state.
   */
  void Clear();

  //@}

  /** @name Accessors */
  //@{

  /// Returns true if the string is initialized, false otherwise
  bool IsOk() const { return m_bstr != NULL; }

  /// Return the length of the string, 0 if not initialized
  size_t Length() const { return ::SysStringLen(m_bstr); }

  /// Implicit conversion to BSTR
  operator BSTR() const { return m_bstr; }

  //@}

private:
  // init the string to an empty state
  void Init();

  // copy the contents of the specified NUL-terminated buffer to this string
  void Copy(const wchar_t *bstr);

  // free the string we're currently holding
  //
  // may be called even if m_bstr == NULL
  //
  // does *not* reset m_bstr to NULL, this is the callers responsability
  void Free();

  // the string itself (may be NULL)
  BSTR m_bstr;
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline
void BasicString::Init()
{
  m_bstr = NULL;
}

inline
BasicString::BasicString()
{
  Init();
}

inline
void BasicString::Free()
{
  // this is safe to call even if m_bstr == NULL
  ::SysFreeString(m_bstr);
}

inline
void BasicString::Attach(BSTR bstr)
{
  ASSERT_MSG( bstr, "NULL bstr in BasicString" );

  Free();

  m_bstr = bstr;
}

inline
void BasicString::Copy(const wchar_t *bstr)
{
  Attach(::SysAllocString(bstr));
}

inline
BasicString::BasicString(const char *sz)
{
  Init();

  Copy(String::MB2WC(sz));
}

inline
BasicString::BasicString(const BasicString& str)
{
  Init();

  Copy(str.m_bstr);
}

inline
BasicString::BasicString(const wchar_t *bstr)
{
  Init();

  Copy(bstr);
}

inline
BasicString& BasicString::operator=(const BasicString& str)
{
  Copy(str.m_bstr);

  return *this;
}

inline
BasicString& BasicString::operator=(const wchar_t *bstr)
{
  Copy(bstr);

  return *this;
}

inline
BasicString::~BasicString()
{
  Free();
}

inline
void BasicString::Clear()
{
  Free();

  Init();
}

inline
BSTR BasicString::Detach()
{
  BSTR bstr = m_bstr;

  Init();

  return bstr;
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_BSTR_H_

