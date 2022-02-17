
#include "ito33/tests/testattachedwarrantconvertiblebond.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(AttachedWarrantConvertibleBondTest::suite());

   return runner.run("");
 
}
