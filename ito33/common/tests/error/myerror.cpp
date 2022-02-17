/////////////////////////////////////////////////////////////////////////////
// Name:        test/error/myerror.cpp
// Purpose:     example of file with error codes definitions
// Author:      Vadim Zeitlin
// Created:     12.03.04
// RCS-ID:      $Id: myerror.cpp,v 1.3 2006/06/13 14:46:33 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/gettext.h"

#include "myerror.h"

namespace MyNameSpace
{
  ITO33_IMPLEMENT_ERROR_CLASS;
}

// ----------------------------------------------------------------------------
// the error messages table
// ----------------------------------------------------------------------------

extern const MyNameSpace::Error ITO33_MYERROR_FIRST
("First error.");

extern const MyNameSpace::Error ITO33_MYERROR_SECOND
("Second error.");

extern const MyNameSpace::Error ITO33_MYERROR_FOO
("Foo error.");

extern const MyNameSpace::Error ITO33_MYERROR_BAR
("Bar error.");

extern const MyNameSpace::Error ITO33_MYERROR_BAZ
("Baz error.");

extern const MyNameSpace::Error ITO33_MYERROR_COMPUTER_POSSESSED
("Demogorgon gazes into your eyes, you are turning to stone.");

extern const MyNameSpace::Error ITO33_MYERROR_FROBNICATION_INCOMPLETE
("The error FROBNICATION_INCOMPLETE is self explanatory.");
