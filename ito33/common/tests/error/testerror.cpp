/////////////////////////////////////////////////////////////////////////////
// Name:        test/error/main.cpp
// Purpose:     test file of error example/test program
// Author:      Vadim Zeitlin
// Created:     12.03.04
// RCS-ID:      $Id: testerror.cpp,v 1.3 2006/06/13 14:46:33 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/cppunit.h"

#include "myerror.h"

#include "ito33/tests/testerror.h"

#include <stdio.h>

typedef MyNameSpace::Error MyError;

// declare some error codes we're going to use
extern const MyError ITO33_MYERROR_FIRST;
extern const MyError ITO33_MYERROR_SECOND;


void ErrorCodeTestCase::GetCode()
{
  // we don't know what their error codes are, but they _should_ be
  // consecutive
  CPPUNIT_ASSERT( ITO33_MYERROR_SECOND == ITO33_MYERROR_FIRST + 1 );
}

void ErrorCodeTestCase::OperatorInt()
{
  extern const MyError ITO33_MYERROR_FROBNICATION_INCOMPLETE;

  CPPUNIT_ASSERT( ITO33_MYERROR_FROBNICATION_INCOMPLETE ==
                  ITO33_MYERROR_FROBNICATION_INCOMPLETE );
}

void ErrorCodeTestCase::GetMessage()
{
  CPPUNIT_ASSERT( !strcmp(ITO33_MYERROR_FIRST.GetMessage(), "First error.") );
  CPPUNIT_ASSERT( !strcmp(ITO33_MYERROR_SECOND.GetMessage(), "Second error.") );
}
