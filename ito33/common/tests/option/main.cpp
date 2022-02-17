

#include "ito33/tests/testoption.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(OptionTest::suite());

   return runner.run("");
 
}
