

#include "ito33/tests/testmoneymarket.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(MoneyMarketTest::suite());

   return runner.run("");
 
}
