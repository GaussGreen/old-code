/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/dynlib.h
// Purpose:     class working with dynamic libraries (DLLs)
// Author:      Vadim Zeitlin
// Created:     21.03.03
// RCS-ID:      $Id: dynlib.h,v 1.12 2006/02/24 14:45:32 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/dynlib.h
    @brief  Encapsulates DLL functions support on different platforms.

    This class uses either Win32 or dlopen() Unix functions to work with the
    shared libraries. It is a low level class and as such makes no attempt to
    hide the semantic differences between the platforms (e.g. under Unix you
    can dlopen() the program itself and backlink to it, under Win32 this is
    impossible), only syntaxic ones.

    Also, it doesn't support the old HP-UX shl_open() way although adding
    support for it would be trivial but then HP-UX 11 (and maybe 10?) provides
    dlopen() anyhow so this doesn't seem very useful. Remember that HAVE_DLOPEN
    must be defined for this class to compile under Unix (normally done by
    autoconf) and the programs using it must link with libdl (use -ldl switch)
 */

#ifndef _ITO33_DYNLIB_H_
#define _ITO33_DYNLIB_H_

#include "ito33/debug.h"
#include "ito33/string.h"

#ifdef _WIN32
  #include "ito33/win32/winwrap.h"

  /// type of the DLL handle
  typedef HMODULE DllHandle_t;
#elif defined(HAVE_DLOPEN)
  #include <dlfcn.h>

  typedef void *DllHandle_t;
#else // unknown platform
  #error "DynamicLibrary can't be compiled on this platform, sorry.\n"\
     "(maybe you forgot to define HAVE_DLOPEN?)"
#endif

namespace ito33
{

/**
    DynamicLibrary encapsulates a low level DLL handle.

    The library is typically opened either in the ctor or when Load() is called
    later on and is unloaded automatically in its dtor (although this can also
    be forced by calling Unload() explicitly or prevented from happening by
    using Detach()).

    After the library is loaded successfull, GetSymbol() may be called to
    dynamically resolve the symbols in the shared library as many times as
    needed.
 */
class DynamicLibrary
{
public:
  /// flags for the ctor and Load()
  enum
  {
    /// force resolving all symbols on load (default, no effect under Windows)
    Load_Now = 0,

    /// resolve symbols only when they're used, not immediately when the
    /// library is loaded (faster to load; no effect under Windows)
    Load_Lazy = 1,

    /// make all loaded symbols available to the subsequently loaded
    /// libraries and not only to the caller (no effect under Windows)
    Load_Global = 4,

    /// use the provided library name as is, don't append extension to it
    Load_Verbatim = 8,

    /// the default value of the flags
    Load_Default = Load_Now
  };


  /**
      @name Constructors and destructor

      This class does provide a default ctor but the object is in an unusable
      state when it is used and Load() must be called to change this.

      As semantics of copying dynamic libraries is not really well-defined,
      this class does not provide neither copy ctor nor the assignment
      operator.
   */
  //@{

  /**
      Default ctor, Load() must be called before the library may be used.
   */
  DynamicLibrary();

  /**
      Ctor which loads the given library.

      This ctor will throw if the library can't be loaded.

      @param libname the base or full library name
      @param flags combination of Load_XXX constants defined above
   */
  DynamicLibrary(const char *libname, int flags = Load_Default);

  /**
      Dtor unloads the library automatically if it had been loaded.

      If you want to prevent this from happening (because the library has to
      stay in memory but for whatever reason the object which loaded it must
      be destroyed) use Detach().

      Note that the dtor is not virtual as this class is not supposed to be
      used as a base class.
   */
  ~DynamicLibrary();

  //@}


  /// @name Operations
  //@{

  /**
      Load the given library.

      If another library had been previously loaded by this object it is
      unloaded first.

      Note that the name is interpreted differently depending on whether
      flags parameter contains Load_Verbatim flag or not. If it doesn't
      (default), the system-dependent extension is added to the library name
      if it doesn't already end in one and, if the library path is relative
      and not absolute, "lib" prefix is prepended under Linux. This allows to
      use the same library name in portable code under all systems.

      @param libname the base or full library name
      @param flags combination of Load_XXX constants defined above
      @return true if the library could be loaded, false otherwise
   */
  bool Load(const char *libname, int flags = Load_Default);

  /**
      Detach this object from the library it represents.

      If for some reason you don't want the object to automatically unload
      the library when it is destroyed, this method may be used to take the
      library handle away from it. In this case you, the caller, become
      responsable for unloading it later.

      @return the handle to the dynamic library
    */
  DllHandle_t Detach();

  /**
      Unload the library now.

      Normally there is no need to call this function as it is going to be
      done from dtor anyhow. Also note that you must @b not call it as long
      as any object or functions from the library are still used as it is
      going to result in mysterious crashes during the next access to them.

      After a call to Unload() the object transits to an uninitialized state
      and Load() must be called before ti can be used again.

      As this function is most often called from dtor, it doesn't throw even
      if unloading the library fails -- but just returns false.

      @return true if successfully unloaded, false if unloading failed
   */
  bool Unload();

  /**
      Unload the specified library handle.

      This function should only be used if you had used Detach() before.

      As this function is often called from dtors, it doesn't throw even if
      unloading the library fails.

      @param handle the DLL handle previously returned by Detach()
      @return true if successfully unloaded, false if unloading failed
   */
  static bool Unload(DllHandle_t handle);

  //@}


  /// @name Accessors
  //@{

  /// Returns true if the library had been successfully loaded
  bool IsOk() const;

  /**
      Returns true if the given handle is valid.

      This is only useul for testing the return value of Detach().

      @param handle the DLL handle previously returned by Detach()
   */
  static bool IsOk(DllHandle_t handle);

  /**
      Checks for the symbol existence in the library.

      If the given symbol (function or variable) exists, it is loaded and
      a (non @c NULL) pointer to it is returned, otherwise @c NULL is
      returned. This function may be used instead of GetSymbol() if it is
      expected that the symbol might not be there and you don't want to catch
      the exception which would be thrown by GetSymbol() in this case.

      If the library hadn't been loaded, this method will throw, otherwise it
      won't.

      @param name the name of the symbol (system-dependent)
   */
  void *HasSymbol(const char *name) const;

  /**
      Loads the symbol from the library.

      If the given symbol doesn't exist in the library an exception is
      thrown, use HasSymbol() if you want to avoid this.

      @param name the name of the symbol (system-dependent)
   */
  void *GetSymbol(const char *name) const;

  //@}


  /**
      Returns the full path to the file this DLL was loaded from.

      Under Unix we don't have a portable way to determine the full shared
      library path so sometimes this function just returns the base name of the
      DLL. Under Windows the full path is always returned.
   */
  std::string GetFullPath() const;


  /**
      @name Error handling helpers.

      If you want a function which normally just returns an error code to
      throw an exception you may do something like
      @code
          // suppose we can't use non default ctor for some reason
          DynamicLibrary dll;
          ...
          if ( !dll.Load(name) )
              DynamicLibrary::ThrowOnLoad();
      @endcode
   */
  //@{

  // this function throws an exception if Load() fails
  static void ThrowOnLoad();

  // this function throws an exception if GetSymbol() fails
  static void ThrowOnGetSymbol();

  //@}

private:
  // initialize m_handle to an initial, invalid value
  void Init();

  // this function throws the exception using the last error
  static void Throw(const char *funcname);

  // our library handle
  DllHandle_t m_handle;

  // the name of the library as passed to ctor
  std::string m_name;


  NO_COPY_CLASS(DynamicLibrary);
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline
void DynamicLibrary::Init()
{
  m_handle = 0;
}

inline
DynamicLibrary::DynamicLibrary()
{
  Init();
}

inline
DynamicLibrary::DynamicLibrary(const char *libname, int flags)
{
  // call Init() so that Load() doesn't try to Unload(some junk) first
  Init();

  if ( !Load(libname, flags) )
  {
    ThrowOnLoad();
  }
}

/* static */ inline
bool DynamicLibrary::IsOk(DllHandle_t handle)
{
  return handle != 0;
}

inline
bool DynamicLibrary::IsOk() const
{
  return IsOk(m_handle);
}

inline
bool DynamicLibrary::Unload()
{
  m_name.clear();

  return Unload(m_handle);
}

inline
DllHandle_t DynamicLibrary::Detach()
{
  DllHandle_t handle = m_handle;
  m_handle = 0;
  return handle;
}

inline
DynamicLibrary::~DynamicLibrary()
{
  Unload();
}

inline
void *DynamicLibrary::GetSymbol(const char *name) const
{
  void *symbol = HasSymbol(name);
  if ( !symbol )
    ThrowOnGetSymbol();

  return symbol;
}

} // namespace ito33

#endif // _ITO33_DYNLIB_H_


