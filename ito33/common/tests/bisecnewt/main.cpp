
#include "ito33/cppunit.h"

#include "ito33/tests/testbisecnewt.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(BisecNewtTest::suite());

  return runner.run("") ? 0 : 1;
}
