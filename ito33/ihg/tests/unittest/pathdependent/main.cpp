
#include "ito33/cppunit.h"

#include "ihg/tests/testpathdep.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceCB);

int main()
{

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( PathDepTest::suite() );
  return runner.run("");

}
