/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/atomic.h
// Purpose:     functions for atomically inc/decrementing a value
// Author:      Vadim Zeitlin
// Created:     27.01.04
// RCS-ID:      $Id: atomic.h,v 1.7 2006/02/18 19:28:13 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_THREAD_INCDEC_H_
#define _ITO33_THREAD_INCDEC_H_

/**
    @file   ito33/thread/atomic.h
    @brief  Functions for atomically incrementing and decrementing integers.

    Although atomically changing a value can be done using a critical section
    or a mutex, it is usually much faster to use the native functions (which
    often just expand into inline assembler statements) to do it. Win32 and
    many other systems have such functions, this header defines a type and two
    functions which allow to use them in a portable way.
 */

namespace ito33
{

/**
    Contains everything related to elementary atomic operations on integers.
 */
namespace Atomic
{

#ifdef _WIN32
  #include "ito33/win32/winwrap.h"

  /**
      The integer type on which atomic operations may be performed

      Note that currently we provide a single type for atomic operations only
      but we could do better and have a Atomic::Type<T> template which would
      yield the smallest type on which atomic operations can be performed and
      which can contain all values of the type T (e.g. long for short, int,
      long but int64 for int64 under Windows). We don't need this for now
      however.
   */
  typedef long IntType;

  /**
      Increment the given value.

      There is no return value as it would be unsafe to use it anyhow.
   */
  inline void Inc(IntType *pValue)
  {
    InterlockedIncrement(pValue);
  }

  /**
      Decrement the given value.

      @return false if the value became 0 after the decrement, true otherwise
              (even if it is negative)
   */
  inline bool Dec(IntType *pValue)
  {
    return InterlockedDecrement(pValue) != 0;
  }

#else // !_WIN32
  // TODO: implement faster atomic operations using CPU-specific code, see
  //       http://bugzilla.gnome.org/show_bug.cgi?id=63621 for a long
  //       discussion about how to do it in the best way
  class IntType
  {
  public:
    IntType() { /* intentionally don't initialize, as with primitive type */ }

    IntType(int n) : m_n(n) { }

    IntType& operator=(int n)
    {
      Lock<CriticalSection> lock(m_cs);
      m_n = n;
      return *this;
    }

    bool operator==(int n) const
    {
      Lock<CriticalSection> lock(m_cs);
      return m_n == n;
    }

    bool operator!=(int n) const
    {
      return !(*this == n);
    }

    void Inc()
    {
      Lock<CriticalSection> lock(m_cs);
      m_n++;
    }

    bool Dec()
    {
      Lock<CriticalSection> lock(m_cs);
      return --m_n != 0;
    }

  private:
    int m_n;

    mutable CriticalSection m_cs;
  };

  inline void Inc(IntType *pValue) { pValue->Inc(); }
  inline bool Dec(IntType *pValue) { return pValue->Dec(); }
#endif // OS

} // namespace Atomic

} // namespace ito33

#endif // _ITO33_THREAD_INCDEC_H_
