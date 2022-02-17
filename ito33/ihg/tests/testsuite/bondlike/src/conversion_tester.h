#ifndef _ITO33_CB_CONVERSION_TESTER_H_
#define _ITO33_CB_CONVERSION_TESTER_H_


#include "ito33/common.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"
#include "basic_features_tester.h"

#define XML_IHGTEST_CB_DEFAULT_CONVERSION         "default_conversion_schedule"
#define XML_IHGTEST_CB_ADD_ANYTIME_CHECKDATE      "add_anytime_checkdate"
#define XML_IHGTEST_CB_ADD_ANYTIME_FOREVER        "add_anytime_forever"
#define XML_IHGTEST_CB_ADD_QUARTERLY_FOREVER      "add_quarterly_forever"
#define XML_IHGTEST_CB_ADD_QUARTERLY_NEXT_QUARTER "add_quarterly_next_quarter"

namespace ito33
{


class ConversionTester
{
public:

  ConversionTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  size_t Test()
  {
    size_t nNbErrors = 0;

    nNbErrors += AddAnytimeForever().Test(m_tags);

    return nNbErrors;
  }


  void RestoreInitData(CBData &data);

  CBData AddAnytimeCheckDate()
  {
    return ReadTestFromXML(XML_IHGTEST_CB_ADD_ANYTIME_CHECKDATE);
  }

  CBData AddAnytimeForever()
  {
    return ReadTestFromXML(XML_IHGTEST_CB_ADD_ANYTIME_FOREVER);
  }

  CBData AddQuarterlyForever()
  {
    return ReadTestFromXML(XML_IHGTEST_CB_ADD_QUARTERLY_FOREVER);
  }

  CBData AddQuarterlyNextQuarter()
  {
    return ReadTestFromXML(XML_IHGTEST_CB_ADD_QUARTERLY_NEXT_QUARTER);
  }

  void TestForcedConversion();

  void TestTriggerContinuity();

  void TestComparisons();

  void TestKeepAccrued();

  void TestAnytimeForever();

protected:

  CBData ReadTestFromXML(const char* sTestName);

private:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::ConvertibleBond> m_pCBInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;
  

  NO_COPY_CLASS(ConversionTester);
};

}

#endif // #define _ITO33_CB_CONVERSION_TESTER_H_
