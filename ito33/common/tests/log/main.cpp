/////////////////////////////////////////////////////////////////////////////
// Name:        test/log/main.cpp
// Purpose:     test for logging macros
// Author:      Vadim Zeitlin
// Created:     2004-12-11
// RCS-ID:      $Id: main.cpp,v 1.2 2006/03/07 15:26:21 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/log.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

ITO33_DEFINE_LOG_CATEGORY(LogTest, "logtest");

#define TEST_TRACE    ITO33_TRACE_CATEGORY(LogTest)

int main()
{
#ifdef _WIN32
  log::Win32DebugSink sink;
#else
  log::StdioSink sink;
#endif // _WIN32
  sink.SubscribeAll();

  ITO33_TRACE(LogTest, "This message uses %s", "ITO33_TRACE");
  ITO33_TRACE_CATEGORY(LogTest)("This one uses %s", "ITO33_TRACE_CATEGORY");
  TEST_TRACE("And the last one is output by %s", "TEST_TRACE");

  return 0;
}

