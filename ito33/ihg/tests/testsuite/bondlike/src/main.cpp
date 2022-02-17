#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/cppunit.h"

#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/convertiblebond.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"

#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

// local files
#include "bondlike_reader.h"
#include "basic_features_tester.h"
#include "conversion_tester.h"
#include "reset_tester.h"
#include "newshare_tester.h"
#include "bondlike_test_runners.h"
#include "call_tester.h"
#include "cboption_tester.h"

using namespace ito33;
using namespace ito33::finance;


extern void GlobalInitialization(size_t);

namespace ito33
{

extern size_t g_nNbFailedSelfRegressionTests;

}

ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceBond);
ITO33_FORCE_LINK_MODULE(IHGPriceAttachedWarrantConvertibleBond);
ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceReset);
ITO33_FORCE_LINK_MODULE(IHGPricePEPSLike);
ITO33_FORCE_LINK_MODULE(IHGPriceGeneralizedPEPSLike);
ITO33_FORCE_LINK_MODULE(IHGPricePERCSLike);
ITO33_FORCE_LINK_MODULE(IHGPriceCBOption);

int main(int ,char **)
{
/*******************************************************************
 * nConvergenceNbTests is number of levels for convergence test run
 * for each test case. the recommended number is 3. 
 * Note that when bIsAcceptance is activated, this number becomes 
 * irrelevant.
 *******************************************************************/
  size_t nConvergenceNbTests = 3;

/*******************************************************************
 * bIsAcceptance is flag to activate/desactivate acceptance test.
 * acceptance test is a small sub-set of full test cases. It is a
 * quick test that Developement team must run before pass the codes 
 * to Testing team.
 *******************************************************************/
  //bool bIsAcceptance = true;
  bool bIsAcceptance = false;

//-------------------------------------------------------------------
  if(bIsAcceptance)
    nConvergenceNbTests = 1;

  try
  {  
    GlobalInitialization(nConvergenceNbTests);

    std::cout <<
      "Remove old files in \"" << ITO33_BONDLIKE_FAILEDCASES_DIR << "\"" << std::endl;
    system(String::Printf("del %s\\*.xml", ITO33_BONDLIKE_FAILEDCASES_DIR).c_str());

    RunBondTester(bIsAcceptance);
    RunBasicFeatureTester(bIsAcceptance);
    RunConversionTester(bIsAcceptance);
    RunResetTester(bIsAcceptance);    
    RunPEPSTester(bIsAcceptance);
    RunPERCSTester(bIsAcceptance);  
    RunGeneralizedPEPSTester(bIsAcceptance);    
    RunWarrantTester(bIsAcceptance);
    RunCallTester(bIsAcceptance);
    RunNewShareTester(bIsAcceptance);
    RunCBOptionTester(bIsAcceptance);

    std::cout << "---------------------------------------------------------\n";
    std::cout << "     REPORT\n";
    std::cout << "---------------------------------------------------------\n";

    std::cout << "\n";

    std::cout << "with stricter rule\n";

    std::cout << "\n";

    if(g_nNbWrongConvergenceTests > 0)
      std::cout << g_nNbWrongConvergenceTests
                << " convergence tests have failed!!!!!!" << "\n";
    if(g_nNbWarningConvergenceTests)
      std::cout << g_nNbWarningConvergenceTests
                << "  convergence tests have given warning." << "\n";
    if(g_nNbWrongConvergenceTests == 0 && g_nNbWarningConvergenceTests == 0)
      std::cout << "All convergence tests have passed.\n";

    std::cout << "\n"
              << "Please check " << ITO33_BONDLIKE_REPORTS_DIR << "/" 
              << ITO33_BONDLIKE_CONVERGENCE_REPORT << " for more details."
              << std::endl;
    std::cout << "---------------------------------------------------------\n";

    std::cout << "\n";

    if( g_nNbFailedSelfRegressionTests > 0)
      std::cout << g_nNbFailedSelfRegressionTests
                << " regression tests have failed!!!!!!" << "\n";
    else
      std::cout << "All regression tests have passed.\n";

    std::cout << "\n"
              << "Please check " << ITO33_BONDLIKE_REPORTS_DIR << "/" 
              << ITO33_BONDLIKE_REGRESSION_SELF_REPORT << " for more details."
              << std::endl;
    std::cout << "---------------------------------------------------------\n";

    return 0; 

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
