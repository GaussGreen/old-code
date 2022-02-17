
#include "ito33/cppunit.h"
#include "ito33/tests/testeventmanager.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(EventManagerTest::suite());

  return runner.run("") ? 0 : 1;
}
