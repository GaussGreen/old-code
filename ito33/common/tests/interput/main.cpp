
#include "ito33/cppunit.h"
#include "ito33/tests/testinterput.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(InterpUtTest::suite());

  return runner.run("") ? 0 : 1;
}
