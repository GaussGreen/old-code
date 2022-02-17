/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testdata.h
// Author:      ZHANG Yunzhi
// Created:     16/08/2004
// RCS-ID:      $Id: testdata.h,v 1.18 2006/08/20 09:36:16 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _IHG_TESTS_TESTDATA_H_
#define _IHG_TESTS_TESTDATA_H_

#include "ito33/beforestd.h"
#include <iostream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/debug.h"

#include "ihg/tests/showconvergence.h"
#include "ihg/tests/regression_outputs_set.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ito33/xml/write.h"

#include "ito33/tests/error_quality.h"


///////////////////////////////////////////////////////////////////////////////
#define ITO33_TESTDATA_ADD_TO_PARAMVALUE_IFPOSSIBLE(name)                     \
{                                                                             \
  if( outputOld.Has ##name() && outputNew.Has ##name())                       \
    ITO33_TESTDATA_ADD_TO_PARAMVALUE(name);                                   \
} 

#define ITO33_TESTDATA_ADD_TO_PARAMVALUE(name)                                \
{                                                                             \
  paramvalue.m_names.push_back( #name );                                      \
  paramvalue.m_firsts.push_back(outputOld.Get ##name ());                     \
  paramvalue.m_seconds.push_back(outputNew.Get ##name ());                    \
  paramvalue.m_diffs.push_back                                                \
      (outputNew.Get ##name () - outputOld.Get ##name ());                    \
} 

namespace ito33
{

extern size_t g_nNbWrongConvergenceTests;
extern size_t g_nNbWarningConvergenceTests;

namespace finance
{
  class SessionData;
}

namespace ihg
{
  class TheoreticalModel;
}

namespace XML
{
  class Tag;
}

extern size_t g_nConvergenceNbTests;


static const char STR_FAIL[] = "fail";
static const char STR_PASS[] = "pass";
static const char STR_WARN[] = "warn";

/// temporary output xml file
const std::string ITO33_TEST_OUTPUT_FILE = "./output___.xml";

/**
  TestData is somehow a test interface for ihg model

  It contains ihg model data, derviative data and session data.
  It also has other information such as test file name, test name etc.
  */
template<class T>
class TestData
{
public:

  /// creates a default TestData
  TestData() :
      m_strFailedCasesDir("failedcases"),
      m_nNbErrors(0),
      m_bOutputForRegressionTestSet(false)
  {
  }

  /// type of derivative
  typedef T DerivativeType;

  /// session data
  shared_ptr<finance::SessionData> m_pSessionData;

  /// derivative data
  shared_ptr<DerivativeType> m_pDeriv;

  /// model data
  shared_ptr<ihg::TheoreticalModel> m_pModel;

  /// name of the test
  std::string m_strTestName;

  /// name of the file, if any, where the test data is stored
  std::string m_strFileName;

  /// directory where xml file should be saved when the test fails.
  std::string m_strFailedCasesDir;

  /**
    runs the computation

    @return theorectical price value
    */
  double Price()
  {
    try
    {
      return m_pModel->Compute(*m_pDeriv)->GetPrice();
    }
    catch(const ito33::Exception& e)
    {
      std::cout << e.GetFullMessage() << std::endl;
      HandleException();
    }
    catch(...) 
    {
      HandleException();
    }

    return 0;
  }

  /**
    Gets name of failed test

    @return name of failed test
    */
  std::string GetFailedTestNames()
  {
    return m_strFailedTest;
  }

  /**
    Gets the number of errors

    @return number of errors.
    */
  size_t GetNumberErrors()
  {
    return m_nNbErrors;
  }


  

  /**
    Create output file name
    
    @param sTestName name of the test

    return output file name of the form "failedcases/<testname>_<filename>.xml"
  */
  std::string GetXMLOutputFile(std::string sTestName)
  {
    return  String::Printf("%s/%s_%s.xml", m_strFailedCasesDir.c_str(),
                         sTestName.c_str(), GetTestFileName().c_str() );
  }

  protected:

   /**
    Get the file name that is currently being tested.
    strip .xml and the xmfiles/ information

    return test file name
  */
  std::string GetTestFileName()
  {
    std::string::size_type nBegin = m_strFileName.find_last_of('/');
  
    if(nBegin == std::string::npos)
      nBegin = 1;
    else
      nBegin++;

    std::string::size_type nEnd = m_strFileName.find_last_of('.');

    return m_strFileName.substr(nBegin, nEnd - nBegin);
  }

  /**
    This function should be called before running all tests in
    implementing Test() function for derived class
    */
  void BeforeAllTests()
  {
    // do nothing for undefine data
    CHECK_VOID( m_strTestName != "" || m_strFileName != "",
        "Undefined test data.");
    
    std::string::size_type nBegin = m_strFileName.find_last_of('/');
    if(nBegin == std::string::npos)
      nBegin = 1;
    else
      nBegin++;
    std::string::size_type nEnd = m_strFileName.find_last_of('.');

    m_strOutputName
        = "_" + m_strFileName.substr(nBegin, nEnd - nBegin) + 
          "___" +  m_strTestName;

    m_pDeriv->SetSessionData(m_pSessionData);

    m_pModel->SetDebugOutputFile(ITO33_TEST_OUTPUT_FILE);
  }

  /**
    This function should be called after running all tests in
    implementing Test() function for derived class
    */
  void AfterAllTests()
  {
    // save the case when something is wrong
    if(m_strFailedTest != "")
    {
      rename
        (
          ITO33_TEST_OUTPUT_FILE.c_str(),
          (String::Printf("%s/%s_%s.xml",
                          m_strFailedCasesDir.c_str(),
                          m_strFailedTest.c_str(),
                          m_strOutputName.c_str())).c_str()
        );
    }
    else
      unlink(ITO33_TEST_OUTPUT_FILE.c_str());
  }

  /**
    This function should be called when an exception is thrown
    through the tests.
    */
  void HandleException()
  {
    // we catch exception here, as we want to save failed cases.
    // In fact, if we let CppUnit handle exceptions, we lose the xml file.
    m_nNbErrors++;

    m_strFailedTest += "ExceptionThrown_";
    
    rename
      (
        ITO33_TEST_OUTPUT_FILE.c_str(),
        (String::Printf("%s/%s_%s.xml",
                        m_strFailedCasesDir.c_str(),
                        m_strFailedTest.c_str(),
                        m_strOutputName.c_str())).c_str()
      );

    // the file has been saved, throw exception again to CppUnit
    throw "Exception caught!";
  }

  /**
    Sets output for regression test

    @param po shared pointer to current model output 
    */
  void SetOutputForRegressionTest(shared_ptr<finance::ModelOutput> po)
  {
    m_pOutput = po;
    m_bOutputForRegressionTestSet = true;
  }

  /**
    Gets the reference output result. It can be used for regression
    test or other comparison tests.

    @return reference output result
    */
  shared_ptr<finance::ModelOutput> GetOutputForRegressionTest()
  {
    ASSERT_MSG
        (
          m_bOutputForRegressionTestSet,
          "You have to run the pricer for once "\
          "and set the output for regression test."
        );
    return m_pOutput;
  }

  /// name of failed tests
  std::string m_strFailedTest;

  /// number of failed tests
  size_t m_nNbErrors;

private:

  shared_ptr<finance::ModelOutput> m_pOutput;

  bool m_bOutputForRegressionTestSet;

  std::string m_strOutputName;

protected:

  /**
    Prints out results of regression test

    @param tag xml tag
    @param outputOld output result from an older version
    @param outputNew result obtained by current version
    @param errorQuality comparison result between old and new outputs.
    */
  void WriteRegressionResults(ito33::XML::RootTag& tag,
                              const finance::ModelOutput& outputOld,
                              const finance::ModelOutput& outputNew,
                              ErrorQuality errorQuality)
  {
    XML::Tag localtest("test", tag);
    localtest.Element("name")(m_strTestName.c_str());
    localtest.Element("file")(m_strFileName.c_str());
    if(errorQuality == ErrorQuality_pass)
      localtest.Element("result")(STR_PASS);
    else
      localtest.Element("result")(STR_FAIL);
    ito33::XML::ComparisonParameterValue paramvalue;

    paramvalue.Init("name", "old version", "current", "difference");

    ITO33_TESTDATA_ADD_TO_PARAMVALUE(Price);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE(Delta);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE(Gamma);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE(Theta);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE_IFPOSSIBLE(Rho);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE_IFPOSSIBLE(Vega);
    ITO33_TESTDATA_ADD_TO_PARAMVALUE_IFPOSSIBLE(Fugit);

    //output details
    localtest.Element("summary",paramvalue);
  }

  /**
    Runs Regression test

    @param tagRegressionTest xml tag that regression test report should
                            be stored
    @param outputs set of old outputs that new output should be compared to
    */
  void RegressionTest
        (
          ito33::XML::RootTag& tagRegressionTest,
          ito33::RegressionOutputSet<finance::ModelOutput>& outputs
        )
  {
    //------ self regression test ------------------
    shared_ptr<finance::ModelOutput> pOutputOld;
    finance::ModelOutput* po = GetOutputForRegressionTest().get();
    ErrorQuality errorRegression 
      = outputs.CheckModelOutput(m_strOutputName, *po, pOutputOld);

    // if old one exists, we print the comparison result.
    if(pOutputOld)
      WriteRegressionResults
          (
            tagRegressionTest,
            *pOutputOld,
            *po,
            errorRegression
          );

    if(errorRegression == ErrorQuality_fail)
    {
      m_strFailedTest += "SELFREGRESSION_";
      m_nNbErrors ++;
    }
  }

  /**
    Runs convergence test

    @param tagConvergenceTest xml tag will report should be saved
    @return convergence quality
    */
  DifferenceQuality ConvergenceTest
                      (
                        ito33::XML::RootTag& tagConvergenceTest
                      )
  { 

    DifferenceQuality result;
    shared_ptr<finance::ModelOutput> pOutput;

    try
    {

    size_t nNbTests = g_nConvergenceNbTests;
    std::vector< shared_ptr<finance::ModelOutput> > pResults(nNbTests);
    std::vector< shared_ptr<finance::ModelOutput> > pDiffs;
    std::vector< shared_ptr<finance::ModelOutput> > pRatings;

    ihg::ShowConvergence(*m_pModel, *m_pDeriv, pResults,
                    pDiffs, pRatings);

    pOutput = pResults[0];

    std::vector<double> pdPrices(nNbTests);
    std::vector<double> pdDiffs(nNbTests);
    std::vector<double> pdRatings(nNbTests);

    size_t n;

    for(n = 0; n < nNbTests; n++)
    {
      pdPrices[n] = pResults[n]->GetPrice();
      if(n > 0)
        pdDiffs[n] = pDiffs[n]->GetPrice();
      if(n > 1)
        pdRatings[n] = pRatings[n]->GetPrice();
    }        
    
    double dMin = 1.e10, dMax = 0;

    for(n = 0; n < nNbTests; n++)
    {
      if(dMin > pdPrices[n])
        dMin = pdPrices[n];
      if(dMax < pdPrices[n])
        dMax = pdPrices[n];
    }
    result = ComparePrice(dMin, dMax);
    //=======================================================================
    // if the prices are closed between each other, we will let it pass the 
    // test at cppunit level no matter whether there is an oscillation or not
    if(result.lowLevel == ErrorQuality_fail) 
      m_nNbErrors++;
    //=======================================================================

    // now apply strict rule
    if(result.lowLevel == ErrorQuality_fail)
      result.strictLevel = ErrorQuality_fail;
    else
    {
      for(n = 2; n < nNbTests; n++)
      {
        if(fabs(pdDiffs[n]) > 1.1e-3) // let it pass if we have 3 correct digit
        {
          if(pdRatings[n] < 0)
          {
            result.strictLevel = ErrorQuality_fail;
            break;
          }
          else if(pdRatings[n] < 1.1)
          {
            result.strictLevel = ErrorQuality_warn;
            break;
          }
        }
      }
      if(n == nNbTests)
        result.strictLevel = ErrorQuality_pass;
    }

    XML::Tag localtest("test", tagConvergenceTest);
    localtest.Element("name")(m_strTestName.c_str());
    localtest.Element("file")(m_strFileName.c_str());
    switch (result.strictLevel)
    {
    case ErrorQuality_pass:
      localtest.Element("result")(STR_PASS);
      break;
    case ErrorQuality_warn:
      localtest.Element("result")(STR_WARN);
      break;
    case ErrorQuality_fail:
      localtest.Element("result")(STR_FAIL);
      break;
    default:
      localtest.Element("result")("");
      break;
    }
    ito33::XML::ConvergenceParameterValue paramvalue;

    for( n = 0; n < nNbTests; n++)
    {
      paramvalue.m_prices.push_back(  pdPrices[n] );
      paramvalue.m_diffs.push_back( pdDiffs[n] ) ;
      paramvalue.m_ratings.push_back( pdRatings[n] ) ;
    }

    //output details
    localtest.Element("summary",paramvalue);

    SetOutputForRegressionTest(pOutput);
    

    if(result.strictLevel == ErrorQuality_fail)
    {
      m_strFailedTest += "CONVERGENCE_";
      g_nNbWrongConvergenceTests++;
    }
    else if(result.strictLevel == ErrorQuality_warn)
      g_nNbWarningConvergenceTests++;


    } // end try
    catch(const ito33::Exception& e)
    {
      std::cout << e.GetFullMessage() << std::endl;
      HandleException();
    }
    catch(...) 
    {
      HandleException();
    }

    return result;
  }
};

///////////////////////////////////////////////////////////////////////////////



}

#endif
