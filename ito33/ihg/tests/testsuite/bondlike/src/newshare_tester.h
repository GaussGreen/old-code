/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/newshare_tester.h
// Purpose:     Tests for new share feature
// Author:      Nabil
// Created:     2005/04/12
// RCS-ID:      $Id: newshare_tester.h,v 1.6 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NEWSHARE_TESTER_H_
#define _ITO33_NEWSHARE_TESTER_H_


#include "ito33/beforestd.h"
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

class NewShareTester
{
public:

  NewShareTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(CBData &data);

  CBData GetBasicData();

  //void PriceStaySameWhenNewShareButNoDividends();
  void PriceDecreaseWhenNewShare();

  // The same tests must be done for Rest, Peps, Percs...
  // Tests with cross-currency too must be done for all these instruments 

protected:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::ConvertibleBond> m_pCBInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;


  NO_COPY_CLASS(NewShareTester);
};

}

#endif // #define _ITO33_NEWSHARE_TESTER_H_
