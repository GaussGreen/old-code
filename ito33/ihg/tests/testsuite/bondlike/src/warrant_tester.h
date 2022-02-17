
#ifndef _ITO33_WARRANT_TESTER_H_
#define _ITO33_WARRANT_TESTER_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class WarrantTester
{
public:

  WarrantTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(WarrantData &data);

  WarrantData GetBasicData();
  WarrantData AddCallNoticeToBasicData();
  WarrantData DividendOnResetDate();
  WarrantData IsLastTriggerConditionMet();
      
  void TestCapRatioIncreases();
  void TestStrike();
  void TestShareFactor();
  void TestBaseRatio();
  void TestResetDate();
  void TestOneDvsTwoD();
  void TestLastTriggerConditionMetProperties();

protected:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::AttachedWarrantConvertibleBond> m_pWarrantInit;

  shared_ptr<ihg::TheoreticalModel> m_pModelInit;
  
  //store the initial call schedule to avoid
  //problems when the call notice flag is added
  shared_ptr<finance::CallSchedule> m_pCallScheduleInit;

  //helper function
  double Price(const shared_ptr<finance::ShareDependentConversion> &pShareDeConv,
               const shared_ptr<finance::BondTerms> &pBondTerms,              
               std::string &xmlOutputFile, 
               std::string sTestName);

  void CompleteShareDependentConversion(
    shared_ptr<finance::ShareDependentConversion> &pShareDeConv);
  
  NO_COPY_CLASS(WarrantTester);
};



}

#endif // #define _ITO33_WARRANT_TESTER_H_
