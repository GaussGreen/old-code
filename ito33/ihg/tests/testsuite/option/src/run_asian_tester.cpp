#include "ito33/beforestd.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/link.h"

#include "ito33/cppunit.h"

#include "utiltest.h"
#include "testparam.h"

#include "cppunit_convergence_test.h"
#include "cppunit_asianoption_test.h"

using namespace ito33;
using namespace ito33::ihg::test;

extern ito33::ihg::test::TestParam testParam;

ITO33_FORCE_LINK_MODULE(IHGPriceAsianOption);

void RunAsianTester(bool bIsAcceptance)
{

  try
  {
  
    //Run the files that are permanently stored
    //inside the cvs directory
    size_t nIdx = 0;
    size_t nNbFiles = 3;
    for ( nIdx = 0 ; nIdx < nNbFiles ; nIdx++)
    {  
      char fileName[1024];
      sprintf(fileName,"xmlfiles/ihgasian%d.xml",nIdx);  
      testParam.SetFileName(fileName);

      std::cout <<"Testing: " << testParam.GetFileName() << std::endl;

      CppUnit::TextUi::TestRunner runner; 
      runner.addTest( CppUnitAsianOptionTest::suite() );

      if ( !bIsAcceptance ) //convergence testing could be slow
       runner.addTest( CppUnitConvergenceTest::suite() );

      runner.run("");
    }
    

 
    /**************************************************************************
    *
    * Here are stored the specific tests. By specific we mean that
    * the developpers is looking at a particular things he really wants
    * to check for.
    * e.g. interpolation is fine for dividends. forwardoptions, etc ...
    *
    **************************************************************************/
    
    //Reproduce the numbers in Curran's paper using the
    //pde code. The analytical curran formula implementation
    //already does it
    // see common/tests/asianoption
    //Note to reproduce the exact numbers of the paper
    //it is easier to set ONEWEEK to be 1./52. instead
    //of 7./365.
    {
      CppUnit::TextUi::TestRunner runner;
      std::cout <<"Testing: curran formula's." << std::endl;
      runner.addTest( CppUnitCurranAsianOptionTest::suite() );
      runner.run("");
    } 
    


    //Test different properties of Asian option
    {
      CppUnit::TextUi::TestRunner runner;
      std::cout <<"Testing: Asian option properties." << std::endl;
      runner.addTest( CppUnitAsianOptionSpecificTest::suite() );
      runner.run("");
    }


   }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    printf("std exception: %s\n", e.what());
  }
  catch ( ... )
  {
    puts("unknown exception!");
  }

}
