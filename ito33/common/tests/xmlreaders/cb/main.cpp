/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/main.cpp
// Purpose:     main file of bondlike test
// Author:      Zhang
// Created:     24.06.04
// RCS-ID:      $Id: main.cpp,v 1.5 2004/10/05 09:13:53 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "testxmlreader_bondlike.h"


int main()
{
  
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(XMLReaderBondLikeTest::suite());

    return runner.run("") ? 0 : 1;
}

