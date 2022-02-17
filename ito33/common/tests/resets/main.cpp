

#include "ito33/tests/testconversionpricereset.h"
#include "ito33/tests/testresetconversionschedule.h"
#include "ito33/tests/testresets.h"


int main()
{
  
   CppUnit::TextUi::TestRunner runner;

   runner.addTest(ConversionPriceResetTest::suite());
   runner.addTest(ConversionScheduleResetTest::suite());
   runner.addTest(ResetTest::suite());

   return runner.run("");
 
}
