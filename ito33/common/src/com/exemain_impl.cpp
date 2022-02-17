/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/exemain_impl.cpp
// Purpose:     supporting code for COM EXE (out of process) servers
// Author:      Vadim Zeitlin
// Created:     Jan 20, 2004
// RCS-ID:      $Id: exemain_impl.cpp,v 1.7 2006/08/19 23:18:54 wang Exp $
// Copyright:   (c) 2003-2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declaration
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"
#include "ito33/gettext.h"
#include "ito33/vector.h"
#include "ito33/sharedptr.h"
#include "ito33/error.h"

#include "ito33/win32/winwrap.h"
#include "ito33/win32/module.h"

#include "ito33/com/coclass.h"
#include "ito33/com/register.h"
#include "ito33/com/server.h"
#include "ito33/com/coinit.h"

#include <stdlib.h>         // for EXIT_SUCCESS/EXIT_FAILURE

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// for ParseServerCmdLine() below
static const int EXIT_DONT = -1;

// ----------------------------------------------------------------------------
// function prototypes
// ----------------------------------------------------------------------------

// this function is defined using COM_DEFINE_SERVER in the user code
extern COM::EXEServer *CreateEXEServer();

// register/unregister the server if needed, return the exit code or -1 if we
// shouldn't exit
static int ParseServerCmdLine(LPSTR lpCmdLine);

// show the error message to the user
static void ShowError(const std::string& msg);

// ----------------------------------------------------------------------------
// local classes
// ----------------------------------------------------------------------------

// set up the server as the global one during the life time of this object and
// delete it when we're deleted
class ScopedServer
{
public:
  ScopedServer(COM::EXEServer *server)
  {
    if ( !server )
    {
      throw EXCEPTION_MSG(ITO33_UNEXPECTED,
                          TRANS("Server initialization failed."));
    }

    COM::Server::Set(server);
  }

  COM::EXEServer *operator->() const
  {
    return static_cast<COM::EXEServer *>(COM::Server::Get());
  }

  ~ScopedServer()
  {
    delete COM::Server::Reset();
  }

private:
  NO_COPY_CLASS(ScopedServer);
};

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// helpers
// ----------------------------------------------------------------------------

static int ParseServerCmdLine(LPSTR lpCmdLine)
{
  if ( lpCmdLine && (lpCmdLine[0] == '/' || lpCmdLine[0] == '-') )
  {
    if ( String::Strnicmp(lpCmdLine + 1, "REGSERVER", 9) == 0 )
    {
      COM::RegisterAll(COM::Server_OutOfProcess);

      // if it didn't throw, it was successful
      return EXIT_SUCCESS;
    }

    if ( String::Strnicmp(lpCmdLine + 1, "UNREGSERVER", 11) == 0 )
    {
      return COM::UnregisterAll() ? EXIT_SUCCESS : EXIT_FAILURE;
    }
  }

  return EXIT_DONT;
}

static void ShowError(const std::string& msg)
{
  ::MessageBox(NULL,
               msg.c_str(),
               TRANS("Fatal COM Server Error"),
               MB_ICONERROR | MB_OK);
}

// ----------------------------------------------------------------------------
// entry point
// ----------------------------------------------------------------------------

int APIENTRY
WinMain(HINSTANCE hInstance,
        HINSTANCE /* hPrevInstance */,    // unused in Win32
        LPSTR lpCmdLine,
        int /* nCmdShow */)
{
  try
  {
    // save our module handle, it is needed in many places
    Win32::Module::SetHandle(hInstance);

    // check for [un]registration command line switch
    int rc = ParseServerCmdLine(lpCmdLine);
    if ( rc != EXIT_DONT )
      return rc;

    // initialize COM
    COM::CoInitializer coInit;

    // make all our coclasses known
    std::vector< shared_ptr<COM::ClassObjectRegistrar> > registeredObjects;
    for ( const COM::CoClassInfo *coclass = COM::CoClassInfo::GetFirst();
          coclass;
          coclass = coclass->GetNext() )
    {
      COM::Ptr<IUnknown> pObject(coclass->CreateInstance());

      // the class object registrat takes ownership of the object pointer and
      // will release it in its own dtor
      shared_ptr<COM::ClassObjectRegistrar> ptr(new
          COM::ClassObjectRegistrar(coclass->GetClsid(), pObject));

      registeredObjects.push_back(ptr);
    }

    // create the server object
    ScopedServer server(CreateEXEServer());

    // enter the main server loop
    return server->Run();
  }
  catch ( const std::bad_alloc& )
  {
    ShowError(TRANS("Out of memory."));
  }
  catch ( const ito33::Exception& e )
  {
    ShowError(e.GetFullMessage());
  }
  catch ( ... )
  {
    ShowError(TRANS("Unexpected error -- sorry."));
  }

  return EXIT_FAILURE;
}

