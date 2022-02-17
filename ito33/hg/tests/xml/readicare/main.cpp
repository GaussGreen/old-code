/////////////////////////////////////////////////////////////////////////////
// Name:        hg/tests/xml/readicare/main.cpp
// Purpose:     Run HG calibration by using icare output file
// Created:     2005/06/21
// RCS-ID:      $Id: main.cpp,v 1.2 2006/08/19 23:47:18 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "ito33/useexception.h"

#include "ito33/finance/derivatives.h"

#include "ito33/hg/parametrization.h"
#include "ito33/hg/theoreticalmodel.h"
#include "hg/tests/readicare.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceOption);
ITO33_FORCE_LINK_MODULE(HGPriceCDS);

using namespace ito33;
using namespace ito33::finance;

int main()
{
  try
  {
    std::cout.precision(16);

    std::ifstream sIn("C:\\ito33\\output\\icare.txt");

    shared_ptr<hg::Parametrization> pParametrization;
    shared_ptr<finance::Derivatives> pDerivatives;

    hg::ReadICARE(sIn, pParametrization, pDerivatives);

    pParametrization->EnableDebugOutput();

    pParametrization->Calibrate(*pDerivatives);
  }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    printf("std exception: %s\n", e.what());
  }
  catch ( ... )
  {
    puts("unknown exception!");

  }

  return 0;
}
