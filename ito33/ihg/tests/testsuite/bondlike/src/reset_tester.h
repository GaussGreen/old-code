
#ifndef _ITO33_RESET_TESTER_H_
#define _ITO33_RESET_TESTER_H_


#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/moneymarket.h"

#include "ito33/finance/bondlike/reset.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class ResetTester
{
public:

  ResetTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(ResetData &data);

  ResetData GetBasicData();

  ResetData AddDividends();

  ResetData AddNewShare();

  ResetData AddCrossCurrency();

  ResetData AddNewShareAndCrossCurrency();
  
  ResetData AddSoftCall();

  ResetData AddPVCouponMakeWhole();
  
  ResetData AddNoPVCouponMakeWhole();
  
  ResetData AddPremiumMakeWhole();

  ResetData AddCallNoticeToBasicData();

  ResetData AddCallNoticeToAddSoftCallData();

  void PriceDecreaseCapWhenNumberOfResetIncreases();
  void PriceIncreaseFloorWhenNumberOfResetIncreases();
  void PriceStaySameWhenFloorAndCapAreOne();
  void PriceStaySameWithResetDateBeforeValuationDate();

  CBData GetEquivalentCB(Date valuationDate);

protected:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::Reset> m_pResetInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;

  shared_ptr<finance::MoneyMarket> m_pMoneyMarket;

  shared_ptr<finance::Dividends> m_pDividends;
  
  /**
    helper function to price a reset

    @param pBondTerms bond terms
    @param pResetSchedule conversion reset schedule
    @param xmlOutputFile name of the output file created
    @param sTestName name of the test that we want a price for
    @param valuationDate valuation date

    @return price price of the resets
    @return xmlOutpuFile name of the outputfile in case of failure
  */
  double Price(const shared_ptr<finance::BondTerms> &pBondTerms,
             const shared_ptr<finance::ResetConversionSchedule> &pResetSchedule,
             std::string &xmlOutputFile, 
             std::string sTestName, Date valuationDate);

  /**
    helper function to ensure that cash, forfeit coupon and keep
    accrued flags are the same as the xml file read

    @param pResetSchedule reset conversion schedule
  */
  void CompleteConversionSchedule(
    shared_ptr<finance::ResetConversionSchedule> &pResetSchedule);


    NO_COPY_CLASS(ResetTester);
};

/**
  This is to verify that with floor rate and cap floor both 100%,
  reset pricer gives same result as that of cb pricer.
  */
class CBUsingResetPricerTester : public ResetTester
{
public:
  CBUsingResetPricerTester(TagsForIhgBondLikeTest& testTags)
    : ResetTester(testTags)
  {
  }

  bool Compare();

  

private:
  

  NO_COPY_CLASS(CBUsingResetPricerTester);
};


/**
  This is to verify that 1D model and 2D model give same result for
  resets to which 1D model can be applied.
  */
class Reset1DModelTester : public ResetTester
{
public:
  Reset1DModelTester(TagsForIhgBondLikeTest& testTags) : ResetTester(testTags)
  {
  }

  bool Compare();

private:

  NO_COPY_CLASS(Reset1DModelTester);
};

}

#endif // #define _ITO33_RESET_TESTER_H_
