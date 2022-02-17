#ifndef _ITO33_CB_CALL_TESTER_H_
#define _ITO33_CB_CALL_TESTER_H_


#include "ito33/common.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"
#include "basic_features_tester.h"

#define XML_IHGTEST_CB_CALL_BASIC_DATA         "call_basic_data"


namespace ito33
{


class CallTester
{
public:

  CallTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  size_t Test()
  {
    size_t nNbErrors = 0;

    nNbErrors += GetBasicData().Test(m_tags);

    return nNbErrors;
  }


  void RestoreInitData(CBData &data);

  CBData GetBasicData();

  void PriceDecreaseAsTriggerHistoryIncreases();

  void PriceIncreaseAsCallPeriodIncreases();

private:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::ConvertibleBond> m_pCBInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;
  

  NO_COPY_CLASS(CallTester);
};

}

#endif // #define _ITO33_CB_CALL_TESTER_H_
