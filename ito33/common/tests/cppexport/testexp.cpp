/////////////////////////////////////////////////////////////////////////////
// Name:        tests/cppexport/testexp.cpp
// Purpose:     small test to validate the macros for class exporting
// Author:      Pedro Ferreira
// Created:     17.02.05
// RCS-ID:      $Id: testexp.cpp,v 1.1 2005/02/17 18:00:16 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"
#include "ito33/dlldecl.h"

class ITO33_DLLDECL commonTest
{
public:
  void hello()
  {
    std::cout << "Hello (common)" << std::endl;
  }
};

