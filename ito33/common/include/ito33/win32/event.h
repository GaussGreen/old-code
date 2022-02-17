/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/win32/event.h
// Purpose:     Win32 Event class
// Author:      Vadim Zeitlin
// Created:     25.06.03
// RCS-ID:      $Id: event.h,v 1.3 2004/10/05 09:13:40 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_WIN32_EVENT_H_
#define _ITO33_WIN32_EVENT_H_

/**
    @file   ito33/win32/event.h
    @brief  Class encapsulating a Win32 event object.

    Events are low level Win32 synchronization objects, we need them in our
    higher level thread code but they may also be used directly if needed.
 */

#include "ito33/win32/kernobj.h"

namespace ito33
{

namespace Win32
{

/**
    Event encapsulates Win32 kernel event object.

    Win32 events may be automatic or manual: automatic ones are reset by the
    system immediately after a (single!) thread has successfully Wait()ed on
    them while the manual events stay set until they are manually Reset().
 */
class Event : public KernelObject
{
public:
  /// Kind of the event.
  enum Kind
  {
    // don't change the values of the enum elements, they correspond to
    // Win32 convention for CreateEvent() parameters

    /// Automatic Win32 event, reset by the system
    Automatic,

    /// Manual Win32 event, must be Reset()
    Manual
  };

  /**
      Creates a new Win32 event object.

      Throws if event can't be created.

      @param kind kind of the event to create, automatic or manual
      @param state the initial event state, signaled or not
   */
  Event(Kind kind = Automatic, State state = NonSignaled);

  /**
      Signal the event.

      For a manual event, all threads waiting on it are woken up. For an
      automatic event only a single ("random") thread is woken up and the
      others continue to wait.

      Throws if the event couldn't be signalled.
   */
  void Set();

  /**
      Resets the event back to non signaled state.

      This function mostly makes sense for manual events, automatic ones are
      always reset automatically if another thread waited for them.

      Throws if the event couldn't be reset.
   */
  void Reset();
};

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

inline
Event::Event(Kind kind, State state)
{
  m_handle = ::CreateEvent
        (
          NULL,   // security descriptor: none (default)
          kind,   // auto reset or manual?
          state,  // initial state: signaled or nonsignaled?
          NULL    // unnamed event
        );

  if ( !m_handle )
  {
    throw WIN32_EXCEPTION("CreateEvent");
  }
}

inline
void Event::Set()
{
  if ( !::SetEvent(m_handle) )
  {
    throw WIN32_EXCEPTION("SetEvent");
  }
}

inline
void Event::Reset()
{
  if ( !::ResetEvent(m_handle) )
  {
    throw WIN32_EXCEPTION("ResetEvent");
  }
}

} // namespace Win32

} // namespace ito33

#endif //  _ITO33_WIN32_EVENT_H_

