
#ifndef _ITO33_GENERALIZEDPEPS_TESTER_H_
#define _ITO33_GENERALIZEDPEPS_TESTER_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/ihg/theoreticalmodel.h"


#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class GeneralizedPEPSTester
{
public:

  GeneralizedPEPSTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(GeneralizedPEPSData &data);

  GeneralizedPEPSData GetInitialData();
  
  GeneralizedPEPSData GetDataWithOptionalConversion();

  GeneralizedPEPSData AddCallExchange();

  GeneralizedPEPSData AddCallFixedCash();

  GeneralizedPEPSData AddDividends();

  GeneralizedPEPSData AddNewShare();

  GeneralizedPEPSData AddNewShareAndCrossCurrency();

  GeneralizedPEPSData AddCrossCurrency();

  GeneralizedPEPSData AddCallFixedCashWithNotice();
  
  GeneralizedPEPSData AddCallFixedCashWithTriggerPeriod();

protected:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::GeneralizedPEPSLike> m_pGeneralizedPEPSInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;
  
  shared_ptr<finance::CallSchedule> m_pCallFixedCash;

  shared_ptr<finance::GeneralizedPEPSLikeCall> m_pCallExchange;

  shared_ptr<finance::MoneyMarket> m_pMoneyMarket;

  shared_ptr<finance::SpotFXRates> m_pSpotFXRates;

  shared_ptr<finance::Dividends> m_pDividends;

  NO_COPY_CLASS(GeneralizedPEPSTester);
};

}
#endif // #DEFINE _ITO33_GENERALIZEDPEPS_TESTER_H_
