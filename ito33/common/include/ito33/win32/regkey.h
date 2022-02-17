/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/regkey.h
// Purpose:     RegKey class encapsulating Win32 HKEY
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: regkey.h,v 1.19 2006/07/23 02:16:06 zeitlin Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/win32/regkey.h
    @brief  Class for working with the Win32 registry.
 */

#ifndef _ITO33_WIN32_REGKEY_H_
#define _ITO33_WIN32_REGKEY_H_

#include "ito33/debug.h"

#include "ito33/win32/winwrap.h"        // for HKEY &c
#include "ito33/win32/exception.h"
#include "ito33/error.h"
#include "ito33/gettext.h"

namespace ito33
{

namespace Win32
{

/**
    RegKey represents a registry key.

    RegKey will throw an exception if API call unexpectedly fails.

    @todo Note that it is a very simple class for now but will be extended
    later if needed.
 */
class RegKey
{
public:
  /**
      @name Constructors and dtor

      To construct a registry key you must normally specify its parent key
      and the path under parent. The parent may be either another RegKey
      object or one of the predefined keys (HKEY_XXX constants). You may also
      construct the key directly corresponding to one of the predefined keys.

      Note that there is (currently) no default ctor. Also, there is no (nor
      will be) copy ctor as it doesn't make sense to copy registry keys.
   */
  //@{

  /**
      Ctor creates the object corresonding to the key with the given name
      under the parent key.

      If the key doesn't exist, it is created. If this fails, an exception is
      thrown.

      @param hkeyParent the parent key, either RegKey or HKEY_XXX constant
      @param name the name under the parent key (without leading backslash)
   */
  RegKey(HKEY hkeyParent, const char *name);

  /// Same as above but taking a string
  RegKey(HKEY hkeyParent, const std::string& name);

  /**
      Dtor automatically closes the registry key.

      Notice that it isn't virtual and so this class shouldn't be used
      polymorphically.
   */
  ~RegKey();

  //@}

  /**
      @name Accessors

      @todo values and sub keys access
   */
  //@{

  /**
      Implicit conversion to HKEY.

      This makes it possible to use RegKey everywhere instead of HKEY,
      including as hkeyParent parameter of the other methods of this class.

      @return the HKEY corresponding to this key or @c NULL if it is invalid
   */
  operator HKEY() const;

  /**
      Check whether the given registry key exists.

      @param hkeyParent the parent HKEY
      @param name the name of the key to check for
      @return true if the key exists and may be opened, false otherwise
   */
  static bool Exists(HKEY hkeyParent, const char *name);

  /// Same as above but taking a string
  static bool Exists(HKEY hkeyParent, const std::string& name);

  /**
      Check whether this key has a value with this name.

      @param name the name of the value to check for
      @return true if the value exists or false otherwise
  */
  bool HasValue(const char *name) const;

  /// Same as above but taking a string
  bool HasValue(const std::string& name) const;

  /**
      Return the value of the given string key.

      If the value doesn't exist or is not of string type, an exception is
      thrown.

      @param name the name of the value, may be NULL or empty to get the value
                  of the unnamed (default) key value
   */
  std::string GetValue(const char *name) const;

  /// Same as above but taking a string
  std::string GetValue(const std::string& name) const;

  //@}

  /**
      @name Operations.

      @todo Set() for numeric values
   */
  //@{

  /**
      Set the value of the keys default value.

      This can only be used for the unnamed value of the key and only for
      string values. For anything else use Set().

      Note that unlike the usual assignment operator, this one doesn't return
      anything as it hardly makes sense to write chained assignments for this
      class.

      If an error occurs, an exception is thrown.

      @param value the value to be set
   */
  void operator=(const char *value);

  /// same as the previous version but taking a string
  void operator=(const std::string& value);

  /**
      Sets the given value to a string.

      If an error occurs, an exception is thrown.

      @param name the name of the value to be set, default value is set
                  if name is @c NULL or empty
      @param value the value to be set (cannot be NULL)
   */
  void Set(const char *name, const char *value);

  /// same as the previous version but taking a string
  void Set(const char *name, const std::string& value);

  /**
      Deletes the subkey of the given key.

      Note that the key must be empty, i.e. not have any subkeys or values,
      to be deleted. If this is not the case, the key is not deleted and
      false is returned.

      If another error occurs, an exception is thrown.

      @param hkeyParent the parent key, either RegKey or HKEY_XXX constant
      @param name the name under the parent key (without leading backslash)
      @return true of the key was deleted, false if it couldn't be deleted
              because it wasn't empty
   */
  static bool Delete(HKEY hkeyParent, const char *name);

  /// Same as above but takes a string
  static bool Delete(HKEY hkeyParent, const std::string& name);

  /**
      Deletes the value with the given name.

      If the value doesn't exist, an exception is thrown, use HasValue() first
      to check for this.

      @param name the name of the value, may be NULL or empty to reset the
                  default unnamed value (it can't be deleted as it always
                  exists)
   */
  void DeleteValue(const char *name);

  /// Same as above but takes a string
  void DeleteValue(const std::string& name);

  //@}

protected:
  /// Close the key checking whether there was an error (in debug build only)
  static void CloseKey(HKEY hkey);

private:
  // common part of both ctors
  void Init(HKEY hkeyParent, const char *name);

  // the key handle or 0 if none
  HKEY m_hkey;

  // registry keys can't be copied
  RegKey(const RegKey&);
  RegKey& operator=(const RegKey&);
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline void
RegKey::Init(HKEY hkeyParent, const char *name)
{
  LONG rc = ::RegCreateKey(hkeyParent, name, &m_hkey);
  if ( rc != 0 )
  {
    throw WIN32_EXCEPTION_ERR("RegCreateKey", rc);
  }
}

inline
RegKey::RegKey(HKEY hkeyParent, const char *name)
   : m_hkey(NULL)
{
  Init(hkeyParent, name);
}

inline
RegKey::RegKey(HKEY hkeyParent, const std::string& name)
   : m_hkey(NULL)
{
  Init(hkeyParent, name.c_str());
}

/* static */ inline
void
RegKey::CloseKey(HKEY hkey)
{
  long rc = ::RegCloseKey(hkey);
  if ( rc != 0 )
  {
#ifndef _NDEBUG
    // dtor shouldn't throw so we don't let the exception propagate but
    // we should still log it
    std::string
      msg = WIN32_EXCEPTION_ERR("RegCloseKey", rc).GetErrorMessage();
    ::OutputDebugString((msg + "\r\n").c_str());
#endif // _NDEBUG
  }
}

inline
RegKey::~RegKey()
{
  if ( m_hkey )
  {
    CloseKey(m_hkey);
  }
}

inline
RegKey::operator HKEY() const
{
  return m_hkey;
}

inline void
RegKey::Set(const char *name, const char *value)
{
  CHECK_VOID( value, "the value to set in RegKey::Set() can't be NULL" );
  ASSERT_MSG( m_hkey != NULL, "invalid HKEY in RegKey::Set()" );

  LONG rc = ::RegSetValueEx
        (
         m_hkey,                                 // key
         name,                                   // value
         0,                                      // [reserved]
         REG_SZ,                                 // value type
         reinterpret_cast<const BYTE *>(value),  // value itself
         static_cast<DWORD>(strlen(value))       // its size
        );

  if ( rc != 0 )
  {
    throw WIN32_EXCEPTION_ERR("RegSetValue", rc);
  }
}

inline void 
RegKey::Set(const char *name, const std::string& value)
{ 
  Set(name, value.c_str()); 
}

inline
void RegKey::operator=(const char *value)
{
  Set(NULL, value);
}

inline 
void RegKey::operator=(const std::string& value) 
{ 
  *this = value.c_str(); 
}

/* static */ inline
bool RegKey::Delete(HKEY hkeyParent, const char *name)
{
  CHECK( name, false, "the key to delete in RegKey::Delete() can't be NULL" );
  ASSERT_MSG( hkeyParent != NULL, "invalid HKEY in RegKey::Delete()" );

  LONG rc = ::RegDeleteKey(hkeyParent, name);
  if ( rc != 0 )
  {
    if ( rc == ERROR_ACCESS_DENIED )
    {
      FAIL( "Trying to delete a registry key which still has subkeys?" );
    }

    throw WIN32_EXCEPTION_ERR("RegDeleteKey", rc);
  }

  return true;
}

/* static */ inline
bool
RegKey::Delete(HKEY hkeyParent, const std::string& name)
{
  return Delete(hkeyParent, name.c_str());
}

/* static */ inline
bool RegKey::Exists(HKEY hkeyParent, const char *name)
{
  HKEY hkey;
  if ( ::RegOpenKey(hkeyParent, name, &hkey) != 0 )
  {
    return false;
  }

  CloseKey(hkey);

  return true;
}

/* static */ inline
bool
RegKey::Exists(HKEY hkeyParent, const std::string& name)
{
  return Exists(hkeyParent, name.c_str());
}

inline
bool RegKey::HasValue(const char *name) const
{
  return ::RegQueryValueExA(m_hkey, name, 0, NULL, NULL, NULL) == 0;
}

inline
bool
RegKey::HasValue(const std::string& name) const
{
  return HasValue(name.c_str());
}

inline
std::string RegKey::GetValue(const char *name) const
{
  // fixed size buffer for common case
  char buf[256];

  // but if necessary we allocate a larger one
  char *pBuf = buf;

  // try to read the value
  DWORD type;
  DWORD len = SIZEOF(buf);
  LONG rc = ::RegQueryValueExA(m_hkey, name, 0, &type,
                               reinterpret_cast<BYTE *>(pBuf), &len);

  // check what kind of value did we get
  std::string retval;
  switch ( type )
  {
    case REG_DWORD:
      if ( rc == ERROR_SUCCESS )
      {
        retval = String::Printf("%lu", *reinterpret_cast<DWORD *>(buf));
      }
      break;

    case REG_SZ:
      // check if the buffer was big enough
      if ( rc == ERROR_MORE_DATA )
      {
        pBuf = static_cast<char *>(malloc(len));
        rc = ::RegQueryValueExA(m_hkey, name, 0, NULL,
                                reinterpret_cast<BYTE *>(pBuf), &len);
      }

      if ( rc == ERROR_SUCCESS )
        retval = pBuf;

      if ( pBuf != buf )
        free(buf);
      break;

    default:
      extern const ito33::Error ITO33_BAD_PARAM;

      throw ito33::EXCEPTION_MSG
            (
              ITO33_BAD_PARAM,
              String::Printf
              (
                TRANS("Registry value \"%s\" is not of string type."),
                name ? name : ""
              )
            );
  }

  // check for errors
  if ( rc != ERROR_SUCCESS )
    throw WIN32_EXCEPTION_ERR("RegQueryValueExA", rc);

  return retval;
}

inline
std::string RegKey::GetValue(const std::string& name) const
{
  return GetValue(name.c_str());
}

inline
void RegKey::DeleteValue(const char *name)
{
  LONG rc = ::RegDeleteValue(m_hkey, name);
  if ( rc != ERROR_SUCCESS )
    throw WIN32_EXCEPTION_ERR("RegDeleteValue", rc);
}

inline
void RegKey::DeleteValue(const std::string& name)
{
  return DeleteValue(name.c_str());
}

} // namespace Win32

} // namespace ito33

#endif // _ITO33_WIN32_REGKEY_H_

