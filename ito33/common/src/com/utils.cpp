/////////////////////////////////////////////////////////////////////////////
// Name:        com/utils.cpp
// Purpose:     implementation of various COM-related utility functions
// Author:      Vadim Zeitlin
// Created:     04.01.03
// RCS-ID:      $Id: utils.cpp,v 1.9 2006/03/30 13:16:06 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/com/unknown_impl.h"
#include "ito33/com/server.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// global data defined here
// ----------------------------------------------------------------------------

#ifndef NDEBUG
  IUnknown *ito33::COM::g_objectDebugRef = NULL;
#endif // Debug

// ----------------------------------------------------------------------------
// private data
// ----------------------------------------------------------------------------

// this struct holds the number of currently active COM objects and the
// critical section protecting it
static struct ObjectCount
{
  // the critical section protecting count
  CriticalSection cs;

  // the count of COM objects still alive, always >= 0
  unsigned long count;

  // ctor sets the number of objects to 0
  ObjectCount()
  {
    // it's safe to not use cs here, we're called before any threads are
    // created
    count = 0;
  }

private:
  // we can't be copied
  ObjectCount(const ObjectCount&);
  ObjectCount& operator=(const ObjectCount&);
} gs_objects;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// active object counting stuff
// ----------------------------------------------------------------------------

void ito33::COM::IncNumberOfActiveObjects()
{
  Lock<CriticalSection> lock(gs_objects.cs);

  ++gs_objects.count;
}

void ito33::COM::DecNumberOfActiveObjects()
{
  {
    Lock<CriticalSection> lock(gs_objects.cs);

    ASSERT_MSG( gs_objects.count > 0, "error in active object counting" );

    if ( --gs_objects.count )
      return;
  } // unlock gs_objects.cs as OnNoMoreObjects() could want to lock it too

  // no more active objects left, notify the server about it
  COM::Server *server = Server::Get();
  if ( server )
    server->OnNoMoreObjects();
}

unsigned long ito33::COM::GetNumberOfActiveObjects()
{
  Lock<CriticalSection> lock(gs_objects.cs);
  return gs_objects.count;
}

