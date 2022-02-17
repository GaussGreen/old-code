/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/cboption_tester.h
// Purpose:     Tests for cb option instrument
// Author:      Nabil
// Created:     2005/09/20
// RCS-ID:      $Id: cboption_tester.h,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_CBOPTION_TESTER_H_
#define _ITO33_CBOPTION_TESTER_H_

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include "ito33/finance/bondlike/cboption.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/cboptionoutput.h"

#include "ihg/tests/testdata.h"
#include "ihg/tests/testdata_bondlike.h"
#include "ihg/tests/regression_outputs_set.h"

#include "test_bondlike_common.h"

namespace ito33
{

class CBOptionTester
{
public:

  CBOptionTester(TagsForIhgBondLikeTest& testTags) : m_tags(testTags)
  {
  }

  void Setup(std::string strInputFilename);

  void RestoreInitData(CBOptionData &data);

  CBOptionData GetBasicData();

  void PriceIncreaseWithRecallSpread();

protected:
  
  TagsForIhgBondLikeTest& m_tags;

  std::string m_strInputFile;

  shared_ptr<finance::SessionData> m_pSessionDataInit;

  shared_ptr<finance::CBOption> m_pCBOptionInit;

  shared_ptr<ito33::ihg::TheoreticalModel> m_pModelInit;

  NO_COPY_CLASS(CBOptionTester);
};

}

#endif // #define _ITO33_CBOPTION_TESTER_H_
