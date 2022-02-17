/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/kernobj.h
// Purpose:     Win32 KernelObject class
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: kernobj.h,v 1.4 2004/10/05 09:13:40 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_WIN32_KERNOBJ_H_
#define _ITO33_WIN32_KERNOBJ_H_

/**
    @file   ito33/win32/kernobj.h
    @brief  Class encapsulating a Win32 kernel object.

    Win32 kernel objects are identified by an opaque HANDLE. The KernelObject
    class hides this handle and provides generic functions for working with it,
    including Wait()ing on it and, especially, destroying it automatically in
    the objects dtor.
 */

#include "ito33/win32/winwrap.h"

#include "ito33/win32/exception.h"

namespace ito33
{

namespace Win32
{

/**
    KernelObject: any Win32 kernel object, i.e. any HANDLE.

    You may construct a KernelObject from an existing HANDLE or using more
    specific Win32 functions in the derived classes.
 */
class KernelObject
{
public:
  /// Any kernel object may be in either signaled or non signaled state.
  enum State
  {
    // don't change the values of the enum elements, they correspond to
    // Win32 convention: TRUE == signaled, FALSE == not

    /// Non signalled state
    NonSignaled,

    /// Signaled (active) state
    Signaled
  };

  /**
      Default ctor, creates an invalid kernel object.

      You must use Attach() later if you want to do anything useful with the
      object.
   */
  KernelObject();

  /**
      Constructor from an existing kernel object.

      We take ownership of the HANDLE and will close it in our destructor.

      @param handle a valid Win32 HANDLE, must be non @c NULL
   */
  KernelObject(HANDLE handle);

  /**
      Destructor closes the kernel object.

      If the object had never been initialized, nothing is done.
   */
  ~KernelObject();

  /**
      Implicit conversion to a HANDLE.

      This allows passing KernelObjects to any Win32 function expecting a
      HANDLE.

      @return the handle for our object
   */
  operator HANDLE() const;

  /**
      Take ownership of the given handle.

      We free our existing handle, if any, and take ownership (i.e. we're
      going to close it later) of the given one.

      @param handle a valid Win32 HANDLE, may be @c NULL
   */
  void Attach(HANDLE handle);

  /**
      Yield ownership of our handle.

      This function may be used if the handle has to be stored beyond life
      time of this object. After calling this method the caller is
      responsible for closing the handle.
   */
  HANDLE Detach();

  /**
      Wait for the object.

      Wait until this kernel object becomes signaled, returns true if this
      happened or false if timeout expired.

      Throws an exception if an error occurs.
   */
  bool Wait(DWORD timeout = INFINITE) const;

protected:
  /// reset the handle
  void Init();

  /// frees the kernel object if necessary
  void Free();

  /// the kernel object handle or NULL if invalid
  HANDLE m_handle;

//private: -- avoid warnings for all derived classes
  /// kernel objects can't be copied
  KernelObject(const KernelObject&);

  /// kernel objects can't be copied
  KernelObject& operator=(const KernelObject&);
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline
void KernelObject::Init()
{
  m_handle = NULL;
}

inline
KernelObject::KernelObject()
{
  Init();
}

inline
KernelObject::KernelObject(HANDLE handle)
      : m_handle(handle)
{
  ASSERT_MSG( handle, "KernelObject: invalid HANDLE" );
}

inline
void KernelObject::Free()
{
  if ( m_handle )
  {
    if ( !::CloseHandle(m_handle) )
    {
      FAIL( "KernelObject: CloseHandle() failed" );
    }

    Init();
  }
}

inline
KernelObject::~KernelObject()
{
  Free();
}

inline
KernelObject::operator HANDLE() const
{
  return m_handle;
}

inline
void KernelObject::Attach(HANDLE handle)
{
  Free();

  m_handle = handle;
}

inline
HANDLE KernelObject::Detach()
{
  HANDLE handle = m_handle;

  Free();

  return handle;
}

inline
bool KernelObject::Wait(DWORD timeout) const
{
  switch ( ::WaitForSingleObject(m_handle, timeout) )
  {
    case WAIT_OBJECT_0:
      return true;

    case WAIT_TIMEOUT:
      return false;

    default:
      throw WIN32_EXCEPTION("WaitForSingleObject");
  }
}

} // namespace Win32

} // namespace ito33

#endif // _ITO33_WIN32_KERNOBJ_H_

