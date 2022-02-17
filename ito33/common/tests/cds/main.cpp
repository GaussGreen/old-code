

#include "ito33/tests/testcds.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(CDSTest::suite());

   return runner.run("");
 
}
