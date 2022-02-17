
#include "ito33/tests/testonetouch.h"

int main()
{  
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(OneTouchTest::suite());

  return runner.run("");
}
