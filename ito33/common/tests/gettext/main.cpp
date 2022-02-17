/////////////////////////////////////////////////////////////////////////////
// Name:        test/gettext/main.cpp
// Purpose:     main file of gettext test program
// Author:      Vadim Zeitlin
// Created:     18.12.02
// RCS-ID:      $Id: main.cpp,v 1.3 2004/10/05 09:13:50 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/gettext.h"

#include <stdio.h>

// ----------------------------------------------------------------------------
// program entry point
// ----------------------------------------------------------------------------

int main(int argc, char **argv)
{
  static const char *testStrings[] =
  {
    TRANS_NOOP("an example string inside TRANS_NOOP"),
    TRANS_NOOP("second one"),
    "untranslated string",      // this 1 won't get into msg catalog
  };

  puts(TRANS("Hello, world!"));

  unsigned int n;
  if ( argc == 2 )
  {
    if ( sscanf(argv[1], "%u", &n) != 1 || n >= SIZEOF(testStrings) )
    {
      printf(TRANS("Usage: %s [n]\nwhere n must be less than %u.\n"),
         argv[0], SIZEOF(testStrings));
      return 1;
    }
  }
  else
  {
    n = 0;
  }

  printf(TRANS("Test string %u is \"%s\"\n"), n, TRANS(testStrings[n]));

  return 0;
}

