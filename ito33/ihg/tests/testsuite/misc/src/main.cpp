#include <iostream>

#include "ito33/link.h"
#include "ito33/exception.h"

#include "testcdp.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);


int main()
{
  try
  {

    // For acceptance testing
    CppUnit::TextUi::TestRunner runner;

    // Cumulative default probability tests
    runner.addTest(ito33::ihg::CDPTEST::CDPTest::suite());

    // Run the tests
    runner.run("");

    return 0;
  }
  catch (const ito33::Exception& e)
  {
    std::cerr << "Exception caught:" << std::endl;
    std::cerr << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught." << std::endl;

    return 2;
  }

}