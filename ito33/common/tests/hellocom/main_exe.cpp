/////////////////////////////////////////////////////////////////////////////
// Name:        hellocom/main_exe.cpp
// Purpose:     COM support for the out process (DLL) servers
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: main_exe.cpp,v 1.3 2004/10/05 09:13:50 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/com/server.h"

#include "ito33/win32/winwrap.h"

using namespace ito33;

class SimpleEXEServer : public COM::EXEServer
{
public:
  SimpleEXEServer() { }

  virtual bool OnInit()
  {
    ::MessageBox(NULL, "Starting...", "Out of Process Hello COM Server",
                 MB_ICONINFORMATION);

    return true;
  }

  virtual void OnExit()
  {
    ::MessageBox(NULL, "Exiting...", "Out of Process Hello COM Server",
                 MB_ICONINFORMATION);
  }

  virtual int MainLoop()
  {
    MSG msg;
    while ( ::GetMessage(&msg, 0, 0, 0) > 0 )
    {
      ::TranslateMessage(&msg);
      ::DispatchMessage(&msg);
    }

    return msg.wParam;
  }

  virtual void Exit()
  {
    ::PostMessage(NULL, WM_QUIT, 0, 0);
  }
};

COM_DEFINE_EXE_SERVER(SimpleEXEServer);

