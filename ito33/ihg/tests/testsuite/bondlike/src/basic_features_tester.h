
#ifndef _ITO33_CB_BASIC_FEATURE_TESTER_H_
#define _ITO33_CB_BASIC_FEATURE_TESTER_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class BasicFeaturesTester
{
public:

  BasicFeaturesTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  size_t Test()
  {
    size_t nNbErrors = 0;
    nNbErrors += GetBasicData().Test(m_tags);

    nNbErrors += AddSoftCall().Test(m_tags);

    nNbErrors += AddPVCouponMakeWhole().Test(m_tags);
    
    nNbErrors += AddNoPVCouponMakeWhole().Test(m_tags);
    nNbErrors += AddPremiumMakeWhole().Test(m_tags);
    nNbErrors += AddCallNoticeToBasicData().Test(m_tags);

    nNbErrors += AddCallNoticeToAddSoftCallData().Test(m_tags);

    return nNbErrors;

  }


  void RestoreInitData(CBData &data);

  CBData GetBasicData();
  
  CBData AddSoftCall();

  CBData AddPVCouponMakeWhole();
  
  CBData AddNoPVCouponMakeWhole();
  
  CBData AddPremiumMakeWhole();

  CBData AddCallNoticeToBasicData();

  CBData AddCallNoticeToAddSoftCallData();

  CBData AddCallClaimTriggerAsPercentageOfToBasicData();

  CBData AddCallIssuePriceTriggerAsPercentageOfToBasicData();

  CBData AddCallPrincipalTriggerAsPercentageOfToBasicData();

  void TestPriceIncreasesAsYTCIncreases();

  void TestPriceIncreasesAsYTPIncreases();

  void TestPriceIncreasesAsCallNoticePeriodIncreases();

private:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::ConvertibleBond> m_pCBInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;
  
  NO_COPY_CLASS(BasicFeaturesTester);
};

}

#endif // #define _ITO33_CB_BASIC_FEATURE_TESTER_H_
