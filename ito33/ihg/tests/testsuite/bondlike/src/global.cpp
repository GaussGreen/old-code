
#include <fstream>

#include "ito33/string.h"
#include "ihg/tests/testdata.h"

// local files
#include "bondlike_reader.h"
#include "bond_tester.h"
#include "basic_features_tester.h"
#include "conversion_tester.h"
#include "reset_tester.h"
#include "peps_tester.h"
#include "generalized_peps_tester.h"
#include "percs_tester.h"

#include "warrant_tester.h"
#include "call_tester.h"

#include "newshare_tester.h"
#include "cboption_tester.h"

namespace ito33
{

// Although there is a value given here, it will probably be overwritten 
// by the call to GlobalInitialization, so will not be used.
size_t g_nConvergenceNbTests = 5;

size_t g_nNbWrongConvergenceTests = 0;
size_t g_nNbWarningConvergenceTests = 0;

size_t g_nNbFailedSelfRegressionTests = 0;
}

using namespace std;
using namespace ito33;

std::ofstream fileConvergence
      (
      String::Printf("%s/%s", ITO33_BONDLIKE_REPORTS_DIR,
                              ITO33_BONDLIKE_CONVERGENCE_REPORT).c_str()
      );

ito33::XML::RootTag tagConvergenceTest
                        (
                          "root",
                          fileConvergence,
                          "../../convergence_test_style_sheet.xsl"
                        );

std::ofstream fileRegre
      (
      String::Printf("%s/%s", ITO33_BONDLIKE_REPORTS_DIR,
                              ITO33_BONDLIKE_REGRESSION_SELF_REPORT).c_str()
      );
ito33::XML::RootTag tagRegressionTest
                        (
                          "root",
                          fileRegre,
                          "../../comparison_test_style_sheet.xsl"
                        );

TagsForIhgBondLikeTest testTags(tagConvergenceTest, tagRegressionTest);

BondTester bondTester(testTags);
BasicFeaturesTester basicFeatureTester(testTags);
ConversionTester conversionTester(testTags);
ResetTester resetTester(testTags);
CBUsingResetPricerTester resetTester100(testTags);
Reset1DModelTester reset1DTester(testTags);
PEPSTester pepsTester(testTags);
GeneralizedPEPSTester generalized_pepsTester(testTags);
PERCSTester percsTester(testTags);
WarrantTester warrantTester(testTags);
CallTester callTester(testTags);
NewShareTester newshareTester(testTags);
CBOptionTester cboptionTester(testTags);

ito33::RegressionOutputSet<finance::ModelOutput> g_outputsCB;

void GlobalInitialization(size_t nNbConvergenceTests = 1)
{  
  if(nNbConvergenceTests != 0)
    g_nConvergenceNbTests = nNbConvergenceTests;
  // if nNbConvergenceTests = 0, just use initial one

  tagConvergenceTest.precision(10);
  tagRegressionTest.precision(10);
}
