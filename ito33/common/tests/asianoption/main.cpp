

#include "ito33/tests/testasianoption.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(AsianOptionTest::suite());
   runner.addTest(CurranTest::suite());

   return runner.run("");
 
}
