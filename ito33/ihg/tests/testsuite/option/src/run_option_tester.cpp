#include "ito33/beforestd.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/link.h"

#include "ito33/cppunit.h"

#include "utiltest.h"
#include "testparam.h"

#include "cppunit_option_test.h"
#include "cppunit_convergence_test.h"
#include "cppunit_forwardoption_test.h"
#include "testputcallparity.h"

using namespace ito33;
using namespace ito33::ihg::test;

extern ito33::ihg::test::TestParam testParam;

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

void RunOptionTester(bool bIsAcceptance)
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
      sprintf(fileName,"xmlfiles/ihg_option%d.xml",nIdx);  
      testParam.SetFileName(fileName);

      std::cout <<"Testing: " << testParam.GetFileName() << std::endl;

      CppUnit::TextUi::TestRunner runner; 
      runner.addTest( CppUnitOptionTest::suite() );
      runner.addTest( CppUnitConvergenceTest::suite() );

      runner.run("");
    }

  
    if ( !bIsAcceptance ) 
    {    
      size_t nIdx   = 0;
      size_t nFileCounter = 0; 
      GenerateXMLTestingFiles(nFileCounter);
     
      //Loop over each test file and do cppUnit
      for ( nIdx = 0; nIdx < nFileCounter; nIdx++)
      {
        char fileName[1024];
        sprintf(fileName,"xmlfiles/ihg_option_temp%d.xml",nIdx);  
        testParam.SetFileName(fileName);

        std::cout <<"Testing: " << testParam.GetFileName() << std::endl;
 
        CppUnit::TextUi::TestRunner runner; 
        runner.addTest( CppUnitOptionTest::suite() );
        runner.addTest( CppUnitConvergenceTest::suite() );

        runner.run("");
      }

      std::cout <<"Remove temporarely created files in xmlfiles directory " << std::endl;
      system(ito33::String::Printf("del %s\\ihg_option_temp*.xml", "xmlfiles").c_str());

    }//end if

    /**************************************************************************
    *
    * Here are stored the specific tests. By specific we mean that
    * the developpers is looking at a particular things he really wants
    * to check for.
    * e.g. interpolation is fine for dividends. forwardoptions, etc ...
    *
    **************************************************************************/
    

    //forward option
    {
      CppUnit::TextUi::TestRunner runner;
      std::cout <<"Testing: forward option." << std::endl;
      runner.addTest( CppUnitForwardOptionTest::suite() );
      runner.run("");
    }   
    
    //regular option tests
    {
      CppUnit::TextUi::TestRunner runner;
      std::cout <<"Testing: Legacy option test code." << std::endl;
      runner.addTest(  CppUnitOptionLegacyTests::suite() );
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
