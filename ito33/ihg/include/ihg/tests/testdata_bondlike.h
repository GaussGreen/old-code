
#ifndef _IHG_TEST_DATA_BONDLIKE_H_
#define _IHG_TEST_DATA_BONDLIKE_H_

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/useexception.h"
#include "ito33/tests/comparison_parameter_value.h"

#include "ihg/tests/testdata.h"

extern ito33::RegressionOutputSet<ito33::finance::ModelOutput> g_outputsCB;

namespace ito33
{

namespace finance
{
  class Bond;
  class ConvertibleBond;
  class Reset;
  class PEPSLike;
  class GeneralizedPEPSLike;
  class PERCSLike;
  class ModelOutput;
  class AttachedWarrantConvertibleBond;
  class CBOption;
}


class TagsForIhgBondLikeTest
{
public:
  TagsForIhgBondLikeTest
    (
      ito33::XML::RootTag& tagConvergence,
      ito33::XML::RootTag& tagRegressionTest
    )
    : m_tagConvergence(tagConvergence),
      m_tagRegression(tagRegressionTest)
  {
  }

  ito33::XML::RootTag& GetTagConvergence()
  {
    return m_tagConvergence;
  }
  
  ito33::XML::RootTag& GetTagRegression()
  {
    return m_tagRegression;
  }

private:

  ito33::XML::RootTag& m_tagConvergence;

  ito33::XML::RootTag& m_tagRegression;

private:

  NO_COPY_CLASS(TagsForIhgBondLikeTest);
};

/**
  TestData for bondlike instruments, such as bond, convertible bond, reset
  or PEPS etc.
  */
template <class T>
class TestDataBondLike : public ito33::TestData<T>
{
public:

  /// creates a default TestDataBondLike
  TestDataBondLike() : TestData<T>()
  {
  }

  /**
    Run all possible tests

    @param testTags set of xml tags for reporting different kind of tests
    @return 0 if all tests pass.
    */
  size_t Test(TagsForIhgBondLikeTest& testTags)
  {
    BeforeAllTests();

    //------ convergence test ------------------
    DifferenceQuality 
        qualityConvergence = ConvergenceTest(testTags.GetTagConvergence());

    //----- regression test ---------------------
    RegressionTest(testTags.GetTagRegression(), g_outputsCB);

    AfterAllTests();

    return GetNumberErrors();
  }

};

/// Bond test data
typedef TestDataBondLike<finance::Bond> BondData;

/// Convertible bond test data
typedef TestDataBondLike<finance::ConvertibleBond> CBData;

/// Reset test data
typedef TestDataBondLike<finance::Reset> ResetData;

/// PEPS test data
typedef TestDataBondLike<finance::GeneralizedPEPSLike> GeneralizedPEPSData;

/// PEPS test data
typedef TestDataBondLike<finance::PEPSLike> PEPSData;

/// PERCS test data
typedef TestDataBondLike<finance::PERCSLike> PERCSData;

/// Warrant test data
typedef TestDataBondLike<finance::AttachedWarrantConvertibleBond> WarrantData;

/// Cb option test data
typedef TestDataBondLike<finance::CBOption> CBOptionData;
}

#endif // #define _IHG_TEST_DATA_BONDLIKE_H_
