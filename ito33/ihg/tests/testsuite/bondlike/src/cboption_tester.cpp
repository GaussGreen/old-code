/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testsuite/bondlike/src/cboption_tester.cpp
// Purpose:     Tests for cb option instument
// Author:      Nabil
// Created:     2005/09/20
// RCS-ID:      $Id: cboption_tester.cpp,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/numeric/predicatedouble.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "cboption_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void CBOptionTester::Setup(std::string strInputFilename)
{
  try
   {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pCBOptionInit = visitor.GetCBOption();

    if(!m_pCBOptionInit)
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no cb option in input xml file");

    m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);
  }
  catch (const ito33::Exception& e)
  {
    std::cout << e.GetErrorMessage() << "\n";
    throw e;
  }
}

void CBOptionTester::RestoreInitData(CBOptionData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = m_pCBOptionInit;

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

CBOptionData CBOptionTester::GetBasicData()
{ 
  CBOptionData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();
  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_BASICDDATA);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

void CBOptionTester::PriceIncreaseWithRecallSpread()
{
  CBOptionData cboption = GetBasicData();
  
  double dPriceInit = cboption.Price();
  
  CBOptionData newcboption = cboption;

  shared_ptr<finance::FloatingRates>
    pFloatingRatesInit = newcboption.m_pDeriv->GetASWFloatingRates();
  
  double
    dMarginInit = pFloatingRatesInit->GetMargin(),
    dMargin = dMarginInit * 2.;

  shared_ptr<finance::FloatingRates> 
    pASWFloatingRates (pFloatingRatesInit);
  
  pASWFloatingRates->SetMargin(dMargin);
  newcboption.m_pDeriv->SetASWFloatingRates(pASWFloatingRates);

  double dPriceWithHigherRecallSpread = newcboption.Price();
  
  //Restore the initial margin
  pASWFloatingRates->SetMargin(dMarginInit);

  CPPUNIT_ASSERT( 
    numeric::IsEqualOrLess(dPriceInit, dPriceWithHigherRecallSpread) );
}
