

#include "ito33/tests/testequity.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(EquityTest::suite());

   return runner.run("");
 
}
