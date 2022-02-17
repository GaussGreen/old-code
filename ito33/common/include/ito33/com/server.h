/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/server.h
// Purpose:     Base class for out of process COM servers
// Author:      Vadim Zeitlin
// Created:     Jan 20, 2004
// RCS-ID:      $Id: server.h,v 1.7 2005/03/01 13:24:34 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/server.h
    @brief Base class for out of process COM servers.


 */

#ifndef _ITO33_COM_SERVER_H_
#define _ITO33_COM_SERVER_H_

#include "ito33/common.h"

namespace ito33
{

namespace COM
{

/**
    Base class for all COM servers.
 */
class Server
{
public:
  /**
      Base class ctor is trivial.
   */
  Server() {} 

  /**
      Virtual dtor for the base class.
   */
  virtual ~Server();

  /**
      @name Server object,

      There is always exactly one server object, these functions provide access
      to it.
   */
  //@{

  /**
      Set the server object, for internal use only.

      This function is only called from the server startup code and must not be
      called from anywhere else.

      We don't take ownership of the provided pointer, the caller is still
      responsible for deleting it if needed.

      @param server the server object
   */
  static void Set(Server *server);

  /**
      Reset the server object, for internal use only.

      This function is only called from the server shutdown code and must not
      be called from anywhere else.

      It does not delete the server object but simply returns it to the caller
      who can delete it if necessary as it is not stored inside this class any
      longer.
   */
  static Server *Reset();

  /**
      Get the current server object.

      This function can only ever return @c NULL during server startup or
      shutdown.

      The returned pointer must not be deleted.

      @return pointer to the global server object
   */
  static Server *Get();

  //@}

  /**
      @name Server callbacks.

      The server receives notifications when it is about to be started and
      before it is exited as well as in some other situations. An in process
      server doesn't have to do anything in these callbacks but an out of
      process has to implement at least OnNoMoreObjects().
   */
  //@{

  /**
      Initialize the server.

      This method is called as soon as possible after the server startup,
      i.e. when COM initialization has already been done but not much else.

      It should only return false if a catastrophic error occured as the server
      is going to shut down then.

      @return true to continue or false to abandon
   */
  virtual bool OnInit() { return true; }

  /**
      Notify the server that there are no more active objects alive.

      Implementation of this callback depends on the server kind, please see
      DLLServer and EXEServer.
   */
  virtual void OnNoMoreObjects() = 0;

  /**
      Cleanup the server.

      This method is called when the server shutdown is initiated but before
      cleaning up COM.
   */
  virtual void OnExit() { }

  //@}

  /**
      @name Locking the server.

      If the server is locked, it is not destroyed, even if there are no more
      active objects alive in it. For the in process servers there is normally
      no reason to lock them from inside but for the out of process servers
      with GUI, the GUI should lock the server for as long as there are any
      open windows as we definitely don't want the server to disappear while
      the user is working with it.
   */
  //@{

  /**
      Increment server lock count.
   */
  virtual void Lock();

  /**
      Decrement the server lock count.

      This may result in the EXE server termination (or DLL unload) if it was
      the last lock.

      @return true if this wasn't the last lock or false if it was
   */
  virtual bool Unlock();

  //@}

private:
  // the one and only server
  static Server *ms_server;

  NO_COPY_CLASS(Server);
};

/**
    Base class for in process COM servers.
 */
class DLLServer : public Server
{
public:
  /**
      Notify the server that there are no more active objects alive.

      There is no need to do anything here for an in process server, we're
      polled using DllCanUnloadNow() anyhow and don't have to do anything
      actively.
   */
  virtual void OnNoMoreObjects() { }


  NO_COPY_CLASS(DLLServer);
};

/**
    Base class for out of process COM servers.
 */
class EXEServer : public Server
{
public:
  /// Trivial default constructor.
  EXEServer() { }

  /**
      Enter the main server loop.

      When this function returns, the server exits and its return value becomes
      the exit code of the process.
   */
  virtual int MainLoop() = 0;

  /**
      Exit the server.

      Calling this function terminates the main loop and exits the server
      unconditionally. This shouldn't be done as long as we have any objects or
      locks.
   */
  virtual void Exit() = 0;

  /**
      Notify the server that there are no more active objects alive.

      By default, an out of process server terminates when there are no more
      active objects and no more locks on the server.
   */
  virtual void OnNoMoreObjects();

  /**
      Convenience function for running a server.

      Run the server: this calls OnInit(), MainLoop() and OnExit() in
      succession.

      May throw an exception if either of the methods above throws.

      @return the value returned by MainLoop() or -1 if OnInit() failed
   */
  int Run();

  /**
      Override base class Unlock() to exit the server when last lock is
      removed.
   */
  virtual bool Unlock();


  NO_COPY_CLASS(EXEServer);
};

/**
    This macro must be used exactly once in all out of process COM servers to
    allow the library code to create a EXEServer-derived object.

    This macro should be used at a global scope, i.e. not inside a namespace.

    @param server the name of a EXEServer-derived class
 */
#define COM_DEFINE_EXE_SERVER(server)                                         \
  extern ito33::COM::EXEServer *CreateEXEServer() { return new server; }      \
  struct DummyCreateStruct ## server /* just to force semicolon after macro */

// ----------------------------------------------------------------------------
// Server inline functions implementation
// ----------------------------------------------------------------------------

inline
Server *Server::Get()
{
  return ms_server;
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_SERVER_H_

