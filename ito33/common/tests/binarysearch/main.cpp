
#include "ito33/cppunit.h"

#include "ito33/tests/testbinarysearch.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(BinarySearchTest::suite());

  return runner.run("") ? 0 : 1;
}
