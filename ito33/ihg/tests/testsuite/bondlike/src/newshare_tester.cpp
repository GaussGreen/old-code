/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/newshare_tester.cpp
// Purpose:     Tests for new share feature
// Author:      Nabil
// Created:     2005/04/12
// RCS-ID:      $Id: newshare_tester.cpp,v 1.6 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/sharedptr.h"

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"

#include "ito33/xml/finance/bondlike/callschedule.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "newshare_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void NewShareTester::Setup(std::string strInputFilename)
{
  try
   {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pCBInit = visitor.GetConvertibleBond();

    if(!m_pCBInit)
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no cb in input xml file");

    m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);
  }
  catch (const ito33::Exception& e)
  {
    std::cout << e.GetErrorMessage() << "\n";
    throw e;
  }
}

void NewShareTester::PriceDecreaseWhenNewShare()
{
  CBData cb = GetBasicData();
  
  double dPriceInit = cb.Price();

  // Deactivate new share feature
  cb.m_pDeriv->SetConvertIntoNewShare(false);

  double dPriceWithoutNewShare = cb.Price();

  RestoreInitData(cb);

  CPPUNIT_ASSERT( numeric::IsEqualOrLess(dPriceInit,dPriceWithoutNewShare) );

}

void NewShareTester::RestoreInitData(CBData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = m_pCBInit;

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

CBData NewShareTester::GetBasicData()
{ 
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();
  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_BASICDDATA);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

