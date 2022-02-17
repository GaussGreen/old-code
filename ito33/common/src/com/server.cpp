/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/server.cpp
// Purpose:     supporting code for COM servers
// Author:      Vadim Zeitlin
// Created:     Jan 21, 2004
// RCS-ID:      $Id: server.cpp,v 1.3 2004/10/05 09:13:43 pedro Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declaration
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"

#include "ito33/com/server.h"
#include "ito33/com/unknown_impl.h"   // for GetNumberOfActiveObjects()

#include "ito33/thread/atomic.h"

using namespace ito33;
using namespace ito33::COM;

// ----------------------------------------------------------------------------
// variables
// ----------------------------------------------------------------------------

Server *Server::ms_server = NULL;

// the lock count of the global server
static Atomic::IntType gs_nLocks = 0;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Server
// ----------------------------------------------------------------------------

Server::~Server()
{
  // nothing to do here
}

/* static */
void Server::Set(Server *server)
{
  ASSERT_MSG( !ms_server, "Calling Server::Set() twice?" );

  ms_server = server;
}

/* static */
Server *Server::Reset()
{
  ASSERT_MSG( ms_server, "Calling Server::Reset() without previous Set()?" );

  Server *server = ms_server;
  ms_server = NULL;

  return server;
}

void Server::Lock()
{
  Atomic::Inc(&gs_nLocks);
}

bool Server::Unlock()
{
  return Atomic::Dec(&gs_nLocks);
}

// ----------------------------------------------------------------------------
// EXEServer
// ----------------------------------------------------------------------------

bool EXEServer::Unlock()
{
  const bool rc = Server::Unlock();
  if ( !rc )
  {
    // no more locks, check if there are any objects alive
    if ( !COM::GetNumberOfActiveObjects() )
    {
      // and exit if there are none
      Exit();
    }
  }

  return rc;
}

void EXEServer::OnNoMoreObjects()
{
  // if we're not locked, we should exit immediately, otherwise we'll have to
  // wait until the last lock is removed
  //
  // to avoid race condition here (a lock could appear while we're deciding
  // whether to exit...), provoke a check in Unlock()
  Lock();
  Unlock();
}

int EXEServer::Run()
{
  if ( !OnInit() )
    return -1;

  int rc = MainLoop();

  OnExit();

  return rc;
}

