/////////////////////////////////////////////////////////////////////////////
// Name:        test/debug/main.cpp
// Purpose:     main file of debug test program
// Author:      Vadim Zeitlin
// Created:     29.04.03
// RCS-ID:      $Id: main.cpp,v 1.6 2004/10/05 09:13:49 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"

#include <stdio.h>
#include <string.h>

#ifdef _MSC_VER
  #include <crtdbg.h>
#endif

// ----------------------------------------------------------------------------
// local functions
// ----------------------------------------------------------------------------

// this function shows how FAIL/ASSERT behave
static void Fail()
{
  FAIL( "Don't call Fail(), it always fails." );
}

// this function leaks memory, you will see the summary of the leaked memory
// in VC++ debug window on program termination
static void LeakMemory()
{
  char *buf = strdup("Hello, world!");
  buf[10] = '?';

  // don't leak all memory to be "more realistic"
  delete [] new double[10];

  int *arr = new int[10];
  arr[0] = 0xAAAAAAAA;
  arr[1] = 0xBBBBBBBB;
  arr[2] = 0xCCCCCCCC;
  arr[3] = 0xDDDDDDDD;
  arr[4] = 0xEEEEEEEE;
  arr[5] = 0xFFFFFFFF;
}

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main()
{
  puts("Welcome to a buggy program!\n");

  printf("Calling function which fails an assert: ");
  Fail();
  puts("done.\n");

  puts("Now we're going to allocate some memory:");

  ito33::MemoryUsageSnapshot ms1;
  printf("\tMax consumption  = %lu bytes\n", ms1.GetMaxUsage());

  LeakMemory();

  ito33::MemoryUsageSnapshot ms2;

  printf("After allocating memory:\n"
     "\tTotal bytes diff = %ld bytes\n"
     "\tAllocations diff = %ld blocks\n"
     "\tMax consumption  = %ld bytes\n\n",
     ms2.GetUsage() - ms1.GetUsage(),
     ms2.GetAllocCount() - ms1.GetAllocCount(),
     ms2.GetMaxUsage());

  printf("Another assertion failure: ");
  ASSERT_MSG( 2 + 2 == 5, "Unexpected addition result" );
  puts("done.\n");

  return 0;
}

