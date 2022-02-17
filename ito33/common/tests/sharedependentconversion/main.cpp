
#include "ito33/tests/testsharedependentconversion.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(ShareDependentConversionTest::suite());

   return runner.run("");
 
}
