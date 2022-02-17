
#include "ito33/hg/bondlikeoutput.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/xml/reader.h"

#include "hg/tests/generic.h"
#include "hg/tests/sensitivity.h"
#include "hg/tests/sparseconstraintsolver.h"
#include "hg/tests/hedging.h"
#include "hg/tests/hero.h"

#include "hg/optionstepper.h"

#include <iostream>
#include <list>

using namespace ito33;

int RunGenericTests();
int RunSensitivityTests();
int RunSparseConstraintSolverTests();
int RunHedgingTests();
int RunHeroTests();


int main()
{
  try
  {    
    int iFailed = 0;
    iFailed += RunGenericTests();
    iFailed += RunSensitivityTests();
    iFailed += RunSparseConstraintSolverTests();
    iFailed += RunHedgingTests();
    iFailed += RunHeroTests();

    std::cout << "Total number of failed tests: " << iFailed << std::endl;
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


int RunGenericTests()
{
  std::cout << "Generic tests" << std::endl;

  hg::XML::Reader reader("xmlfiles/option01.xml");
    
  shared_ptr<finance::Derivative> pDerivative;

  reader.ReadDerivative(pDerivative);

  CppUnit::TextUi::TestRunner runner;
    
  // Setup the generic test on models
  hg::GenericTest::Setup(pDerivative); 
  runner.addTest(hg::GenericTest::GenericTest::suite());

  return runner.run("") ? 0 : 1;
}


int RunSensitivityTests()
{

  std::cout << "Sensitivity tests" << std::endl;

  // Make a list of xml files to test
  std::list<std::string> fileList;

  fileList.push_back("xmlfiles/option01.xml");
  fileList.push_back("xmlfiles/option02.xml");
  fileList.push_back("xmlfiles/option03.xml");
  fileList.push_back("xmlfiles/option04.xml");
  fileList.push_back("xmlfiles/option05.xml");

  fileList.push_back("xmlfiles/eds01.xml");
  fileList.push_back("xmlfiles/eds02.xml");

  fileList.push_back("xmlfiles/cb01.xml");

  fileList.push_back("xmlfiles/cds01.xml");  

  std::list<std::string>::const_iterator iterFileList;

  // Run the tests
  int iFailed = 0;
  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {
    hg::XML::Reader reader( (*iterFileList).c_str() );

    shared_ptr<finance::Derivative> pDerivative;

    reader.ReadDerivative(pDerivative);

    shared_ptr<hg::TheoreticalModel> pModel;
    
    reader.ReadTheoreticalModel(pModel);

    hg::SensitivityTest::Setup(pDerivative, pModel);

    CppUnit::TextUi::TestRunner runner;

    runner.addTest(hg::SensitivityTest::SensitivityTest::suite());

    iFailed += runner.run("") ? 0 : 1;

  }

  // Run partial sensitivity tests on each contract type.  No need
  // to run on all the files.
  fileList.clear();

  fileList.push_back("xmlfiles/option04.xml");
  fileList.push_back("xmlfiles/eds02.xml");
  fileList.push_back("xmlfiles/cb01.xml");
  fileList.push_back("xmlfiles/cds01.xml");  

  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {
    hg::XML::Reader reader( (*iterFileList).c_str() );

    shared_ptr<finance::Derivative> pDerivative;

    reader.ReadDerivative(pDerivative);

    shared_ptr<hg::TheoreticalModel> pModel;
    
    reader.ReadTheoreticalModel(pModel);

    hg::SensitivityTest::Setup(pDerivative, pModel);

    CppUnit::TextUi::TestRunner runner;

    runner.addTest(hg::SensitivityTest::PartialSensitivityTest::suite());

    iFailed += runner.run("") ? 0 : 1;

  }


  return iFailed;
}


int RunSparseConstraintSolverTests()
{
  std::cout << "Sparse constraint solver tests" << std::endl;

  // Make a list of xml files to test
  std::list<std::string> fileList;

  // Standard American put option
  fileList.push_back("xmlfiles/option03.xml");

  // American put option with large dividend
  fileList.push_back("xmlfiles/option04.xml");

  // American put option with large jump values
  fileList.push_back("xmlfiles/option05.xml");


  std::list<std::string>::const_iterator iterFileList;

  // Run the tests
  int iFailed = 0;
  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {
    hg::XML::Reader reader( (*iterFileList).c_str() );

    shared_ptr<finance::Derivative> pDerivative;

    reader.ReadDerivative(pDerivative);

    shared_ptr<hg::TheoreticalModel> pModel;
    
    reader.ReadTheoreticalModel(pModel);

    hg::SparseConstraintSolverTest::Setup(pDerivative, pModel);

    CppUnit::TextUi::TestRunner runner;

    runner.addTest(hg::SparseConstraintSolverTest::SparseConstraintSolverTest::suite());

    iFailed += runner.run("") ? 0 : 1;

  }

  return iFailed;

}


int RunHedgingTests()
{
  std::cout << "Hedging tests" << std::endl;

  CppUnit::TextUi::TestRunner runner;

  runner.addTest(hg::HedgingTest::HedgingTest::suite());

  int iFailed = runner.run("") ? 0 : 1;

  return iFailed;

}


int RunHeroTests()
{
  std::cout << "Hero tests" << std::endl;

  CppUnit::TextUi::TestRunner runner;

  runner.addTest(hg::HeroTest::HeroTest::suite());

  int iFailed = runner.run("") ? 0 : 1;

  return iFailed;

}
