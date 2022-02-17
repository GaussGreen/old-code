#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/sharedptr.h"
#include "ito33/debug.h"

#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/pepslike.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/finance/moneymarket.h"
#include "ito33/xml/finance/spotfxrates.h"
#include "ito33/xml/finance/bondlike/callfixedshare.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "peps_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void PEPSTester::Setup(std::string strInputFilename)
{
  try
  {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pPEPSInit = visitor.GetPEPSLike();

    if(!m_pPEPSInit)
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no PEPS in input xml file");

    m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);

    // fixed share
    Restore(reader.GetRootNode(), m_pCallFixedShare);

    // fixed cash
    Restore(reader.GetRootNode(), m_pCallFixedCash);

    // cross currency
    const xml::node& nodeRoot = reader.GetRootNode();

    xml::node::const_iterator
      pNode = nodeRoot.find(XML_TAG_MONEYMARKET_ROOT);
    
    if( pNode != nodeRoot.end() )
      m_pMoneyMarket = GetMoneyMarketFromNode(*pNode);

    pNode = nodeRoot.find(XML_TAG_SPOT_FX_RATES_ROOT);
  
    if( pNode != nodeRoot.end() )
      m_pSpotFXRates = GetSpotFXRatesFromNode(*pNode);

    // Dividends
    const xml::node::const_iterator 
      div_node = reader.GetRootNode().find(XML_TAG_FINANCE_DIVIDENDS);
    if ( div_node != reader.GetRootNode().end() )
    {
      m_pDividends = ReadDividends( *div_node );
    }
    else
      m_pDividends.reset();
  }
  catch (const ito33::Exception& e)
  {
    std::cout << e.GetErrorMessage() << "\n";
    throw e;
  }
}


void PEPSTester::RestoreInitData(PEPSData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = make_ptr( new PEPSLike(*m_pPEPSInit) );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

PEPSData PEPSTester::GetInitialData()
{ 
  PEPSData data;

  RestoreInitData(data);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

PEPSData PEPSTester::GetDataWithOptionalConversion()
{ 
  PEPSData data(GetInitialData());

  data.m_pDeriv->EnableOptionalConversion();

  data.m_strTestName = std::string("add_optional_conversion");

  return data;
}


PEPSData PEPSTester::AddCallFixedShare()
{ 
  ASSERT_MSG(m_pCallFixedShare, "m_pCallFixedShare must exist");

  PEPSData data(GetDataWithOptionalConversion());

  data.m_pDeriv->SetCallFixedShare(m_pCallFixedShare);

  data.m_strTestName = std::string("add_call_fixed_share");

  return data;
}


PEPSData PEPSTester::AddCallFixedCash()
{ 
  ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");

  PEPSData data(GetDataWithOptionalConversion());

  data.m_pDeriv->SetCallFixedCash(m_pCallFixedCash);

  data.m_strTestName = std::string("add_call_fixed_cash");

  return data;
}

PEPSData PEPSTester::AddCallFixedCashWithNotice()
{ 
  ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");

  PEPSData data(GetDataWithOptionalConversion());

  shared_ptr<CallSchedule> pCalls(new CallSchedule(*m_pCallFixedCash));
  pCalls->SetNoticePeriod(21);
  data.m_pDeriv->SetCallFixedCash(pCalls);

  data.m_strTestName = std::string("add_call_fixed_cash_with_notice");

  return data;
}

PEPSData PEPSTester::AddCallFixedCashWithTriggerPeriod()
{ 
  ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");

  PEPSData data(GetDataWithOptionalConversion());

  shared_ptr<CallSchedule> pCalls(new CallSchedule(*m_pCallFixedCash));
  pCalls->SetTriggerCheckPeriod(21, 0);
  data.m_pDeriv->SetCallFixedCash(pCalls);

  data.m_strTestName = std::string("add_call_fixed_cash_with_trigger_period");

  return data;
}

PEPSData PEPSTester::AddCrossCurrency()
{ 
  ASSERT_MSG(m_pSpotFXRates && m_pMoneyMarket, 
             "cross currency property undefined");

  PEPSData data(GetDataWithOptionalConversion());

  data.m_pDeriv->SetNumeraire( m_pMoneyMarket->GetNumeraire() );

  data.m_pSessionData->GetRateData()->SetYieldCurve(
    m_pMoneyMarket->GetNumeraire(), m_pMoneyMarket->GetYieldCurve() );

  data.m_pSessionData->GetRateData()->SetSpotFXRates(m_pSpotFXRates);

  data.m_strTestName = std::string("add_cross_currency");

  return data;
}

PEPSData PEPSTester::AddDividends()
{ 
  ASSERT_MSG(m_pDividends, "dividends must exist");

  PEPSData data(GetDataWithOptionalConversion());

  // Dividends added
  shared_ptr<Dividends> pDividends(new Dividends(*m_pDividends));

  data.m_pSessionData->GetEquity()->SetDividends(pDividends);

  data.m_strTestName = std::string("add_dividends");

  return data;
}

PEPSData PEPSTester::AddNewShare()
{ 
  PEPSData data(GetDataWithOptionalConversion());

  // Dividends added
  if( m_pDividends )
  {
    shared_ptr<Dividends> pDividends(new Dividends(*m_pDividends));
    data.m_pSessionData->GetEquity()->SetDividends(pDividends);
  }

  data.m_pDeriv->SetConvertIntoNewShare(true);

  data.m_strTestName = std::string("add_newshare");

  return data;
}

PEPSData PEPSTester::AddNewShareAndCrossCurrency()
{ 
  ASSERT_MSG(m_pSpotFXRates && m_pMoneyMarket, 
             "cross currency property undefined");

  PEPSData data(GetDataWithOptionalConversion());

  // Dividends added
  if( m_pDividends )
  {
    shared_ptr<Dividends> pDividends(new Dividends(*m_pDividends));
    data.m_pSessionData->GetEquity()->SetDividends(pDividends);
  }

  data.m_pDeriv->SetConvertIntoNewShare(true);

  data.m_pDeriv->SetNumeraire( m_pMoneyMarket->GetNumeraire() );

  data.m_pSessionData->GetRateData()->SetYieldCurve(
    m_pMoneyMarket->GetNumeraire(), m_pMoneyMarket->GetYieldCurve() );

  data.m_pSessionData->GetRateData()->SetSpotFXRates(m_pSpotFXRates);

  data.m_strTestName = std::string("add_newshare_and_cross_currency");

  return data;
}

