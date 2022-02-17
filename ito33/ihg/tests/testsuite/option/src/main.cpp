#include "ito33/beforestd.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/link.h"
#include "ito33/string.h"
#include "ito33/exception.h"

#include "ito33/xml/write.h"

#include "optionlike_test_runner.h"
#include "testparam.h"

std::ofstream MainReportFile("reports/main.xml");

ito33::XML::RootTag root("root", MainReportFile, "../../test_style_sheet.xsl");  

std::ofstream ConvergenceReportFile("reports/convergence.xml");

ito33::XML::RootTag 
  convergenceRootTag( "root",
                      ConvergenceReportFile,
                      "../../convergence_test_style_sheet.xsl");  

ito33::ihg::test::TestParam testParam;

int main(int ,char **)
{
 try
  {

    /*******************************************************************
    * bIsAcceptance is flag to activate/desactivate acceptance test.
    * acceptance test is a small sub-set of full test cases. It is a
    * quick test that Developement team must run before pass the codes 
    * to Testing team.
    *******************************************************************/
    //bool bIsAcceptance = true;
    bool bIsAcceptance = false;

    std::cout <<"Remove old files in failedcases directory " << std::endl;
    system(ito33::String::Printf("del %s\\*.xml", "failedcases").c_str());

    RunOptionTester(bIsAcceptance);
    RunAsianTester(bIsAcceptance);


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

  return 0;

}
