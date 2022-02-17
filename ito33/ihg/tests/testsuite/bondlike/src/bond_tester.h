
#ifndef _ITO33_BOND_TESTER_H_
#define _ITO33_BOND_TESTER_H_

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/bondlike/bond.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/bondlikeoutput.h"

#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class BondTester
{
public:

  BondTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(BondData &data);

  BondData GetBasicData();

  void TestPriceIncreasesAsYTCIncreases();

  void TestPriceIncreasesAsYTPIncreases();

  // (FDF for Forward Discount Factor)
  void TestPriceWithCouponsDefinedByFDF();

  /*
    For a bond without constraints and without default (ie: HR null), we have:
    F(t) = T-t where F is the fugit, T the maturity and t the valuation time. 
    T - t is the time to maturity and must be in years.
  */
  void TestFugitForBondWithNoConstraintsAndNoDefault();

  /*
    For a bond without constraints and with a flat hazard rate p != 0, we have:
    F(t) = (1/p)[ 1 - exp( -p(T-t) ) ] where F is the fugit, T the maturity and
    t the valuation time. T - t must be in years.

    Remark: Test is done with low, medium and high hazard rate.
  */
  void TestFugitForBondWithNoConstraintsAndNonNullFlatHR();

private:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::Bond> m_pBondInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;
  
  NO_COPY_CLASS(BondTester);
};

}

#endif // #define _ITO33_BOND_TESTER_H_
