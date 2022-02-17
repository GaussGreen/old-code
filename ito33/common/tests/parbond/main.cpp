

#include "ito33/tests/testparbond.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(ParBondTest::suite());

   return runner.run("");
 
}
