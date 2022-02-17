/////////////////////////////////////////////////////////////////////////////
// Name:        src/win32/thread.cpp
// Purpose:     implementation of thread-related classes using Win32 API
// Author:      Vadim Zeitlin
// Created:     18.02.03
// RCS-ID:      $Id: thread.cpp,v 1.9 2005/05/02 12:10:03 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/string.h"

#include "ito33/thread.h"
#include "ito33/win32/exception.h"

#ifdef HAVE__BEGINTHREADEX
  #include <process.h>
#endif

// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

// TLSIndexUntyped is one slot in the TLS (thread locale storage)
//
// this class is not well documented but if you refer to the MSDN documentation
// about TLS its purpose should be quite clear: basicly it is just an
// application of RIAA idea to the TLS slots
//
// the only "subtlety" is that TLSIndexUntyped is not used at all -- only
// TLSIndex defined below is
class TLSIndexUntyped
{
public:
  TLSIndexUntyped()
  {
    m_tls = ::TlsAlloc();
    if ( m_tls == TLS_OUT_OF_INDEXES )
    {
      throw WIN32_EXCEPTION("TlsAlloc");
    }
  }

  ~TLSIndexUntyped()
  {
    if ( !::TlsFree(m_tls) )
    {
#ifndef _NDEBUG
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      std::string
        msg = WIN32_EXCEPTION("TlsFree").GetFullMessage();
      ::OutputDebugString((msg + "\r\n").c_str());
#endif // _NDEBUG
    }
  }

  // these functions are only to be called by TLSIndexUntyped<>
protected:
  void Set(void *value)
  {
    if ( !::TlsSetValue(m_tls, value) )
    {
      throw WIN32_EXCEPTION("TlsSetValue");
    }
  }

  void *Get() const
  {
    void *value = ::TlsGetValue(m_tls);
    if ( !value && ::GetLastError() != NO_ERROR )
    {
      throw WIN32_EXCEPTION("TlsGetValue");
    }

    return value;
  }

private:
  // the TLS index
  DWORD m_tls;
};

// a type safe wrapper around TLSIndexUntyped
template <typename T>
class TLSIndex : public TLSIndexUntyped
{
public:
  // chained assignments don't make sense with TLS slots so operator=()
  // exceptionally returns nothing in this class
  void operator=(T value) { Set(reinterpret_cast<void *>(value)); }
  T operator()() const { return reinterpret_cast<T>(Get()); }
};

// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

// TLS index of the slot where we store the pointer to the current thread
//
// we use a function returning a reference to the static object to create it on
// first use only and not on the program startup
static inline TLSIndex<ito33::Thread *>& TlsThisThread()
{
  static TLSIndex<ito33::Thread *> s_tlsThisThread;

  return s_tlsThisThread;
}

// ============================================================================
// Thread implementation
// ============================================================================

namespace ito33
{

// ----------------------------------------------------------------------------
// Thread
// ----------------------------------------------------------------------------

/* static */
THREAD_RETTYPE THREAD_CALLCONV
Thread::ThreadStart(void *data)
{
  Thread * const thread = reinterpret_cast<Thread *>(data);

  // remember the thread pointer in the TLS to be able to retrieve it later
  TlsThisThread() = thread;

  // do run the thread now
  thread->m_run->Run();

  // the detached threads should be deleted when they terminate because
  // there is no other place to do it -- joinable threads, however, are
  // deleted by the user code (presumably after it Join()s them)
  if ( !thread->m_isJoinable )
  {
    delete thread;
  }

  return 0;
}

/* static */
Thread *Thread::GetCurrent()
{
  return TlsThisThread()();
}

void Thread::Init(bool isJoinable)
{
  m_isJoinable = isJoinable;

  // we need to create the thread in the suspended state because otherwise
  // it could start running and try to use m_hThread before CreateThread()
  // returns and so before we write m_hThread to memory
  m_hThread =
#ifdef HAVE__BEGINTHREADEX
        (HANDLE)_beginthreadex
#else // compiler doesn't have _beginthreadex
        ::CreateThread
#endif // _beginthreadex/CreateThread
        (
          NULL,                     // default security
          0,                        // default stack size
          Thread::ThreadStart,      // entry point function
          this,                     // parameter to pass it
          CREATE_SUSPENDED,         // don't run it just yet
          &m_id                     // [out] thread id
        );

  if ( !m_hThread )
  {
    throw WIN32_EXCEPTION("CreateThread");
  }

  // now that we have m_hThread, we can run it
  if ( ::ResumeThread(m_hThread) == (DWORD)-1 )
  {
    throw WIN32_EXCEPTION("ResumeThread");
  }
}

Thread::~Thread()
{
  Close();

  m_run->OnThreadExit();
}

void Thread::Close()
{
  if ( m_hThread )
  {
    if ( !::CloseHandle(m_hThread) )
    {
      FAIL( "CloseHandle(thread) failed" );
    }

    m_hThread = NULL;
  }
}

// for other compilers this function is defined (as empty) in the header
#ifdef _MSC_VER

void Thread::SetName(const char *name)
{
  if ( !name )
  {
    // not an error: simply no name specified in thread ctor
    return;
  }

  // the code is copy-and-pasted from MSDN (topic "Setting a thread name")
  typedef struct tagTHREADNAME_INFO
  {
    DWORD dwType; // must be 0x1000
    LPCSTR szName; // pointer to name (in user addr space)
    DWORD dwThreadID; // thread ID (-1=caller thread)
    DWORD dwFlags; // reserved for future use, must be zero
  } THREADNAME_INFO;

  THREADNAME_INFO info;
  info.dwType = 0x1000;
  info.szName = name;
  info.dwThreadID = m_id;
  info.dwFlags = 0;

  __try
  {
    ::RaiseException(0x406D1388, 0, sizeof(info)/sizeof(DWORD), (DWORD *)&info);
  }
  __except( EXCEPTION_CONTINUE_EXECUTION )
  {
  }
}

#endif // VC++

// ----------------------------------------------------------------------------
// JoinableThread
// ----------------------------------------------------------------------------

JoinableThread::~JoinableThread()
{
  if ( m_hThread )
  {
    try
    {
      Join();
    }
#ifndef _NDEBUG
    catch ( Win32::Exception& e )
    {
      // dtor shouldn't throw so we don't let the exception propagate but
      // we should still log it
      ::OutputDebugString((e.GetFullMessage() + "\r\n").c_str());
    }
#endif // _NDEBUG
    catch ( ... )
    {
      // do nothing -- dtor must not throw
    }
  }
}

void JoinableThread::Join()
{
  CHECK_VOID( m_hThread, "Join() had been already called or invalid thread" );

  if ( ::WaitForSingleObject(m_hThread, INFINITE) != WAIT_OBJECT_0 )
  {
    throw WIN32_EXCEPTION("WaitForSingleObject(thread)");
  }

  // we ignore the exit code for now anyhow
  DWORD exitCode;
  ::GetExitCodeThread(m_hThread, &exitCode);

  ASSERT_MSG( exitCode != STILL_ACTIVE, "Thread not terminated in Join()?" );

  Close();
}

} // namespace ito33

