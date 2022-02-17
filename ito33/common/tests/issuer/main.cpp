

#include "ito33/tests/testissuer.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(IssuerTest::suite());

   return runner.run("");
 
}
