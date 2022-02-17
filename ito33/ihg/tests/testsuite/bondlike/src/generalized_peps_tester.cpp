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
#include "ito33/finance/bondlike/generalizedpepslike.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/finance/moneymarket.h"
#include "ito33/xml/finance/spotfxrates.h"
#include "ito33/xml/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "generalized_peps_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void GeneralizedPEPSTester::Setup(std::string strInputFilename)
{
  try
  {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pGeneralizedPEPSInit = visitor.GetGeneralizedPEPSLike();

    if(!m_pGeneralizedPEPSInit)
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no GeneralizedPEPS in input xml file");

    m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);

    // conversion at issuer's option
    Restore(reader.GetRootNode(), m_pCallExchange);

    //m_bIsCallExchange = false;

    //if (m_pCallExchange)
    //  m_bIsCallExchange = true;

    // fixed share
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


void GeneralizedPEPSTester::RestoreInitData(GeneralizedPEPSData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = make_ptr( new GeneralizedPEPSLike(*m_pGeneralizedPEPSInit) );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

GeneralizedPEPSData GeneralizedPEPSTester::GetInitialData()
{ 
  GeneralizedPEPSData data;

  RestoreInitData(data);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::GetDataWithOptionalConversion()
{ 
  GeneralizedPEPSData data(GetInitialData());

  data.m_pDeriv->EnableOptionalConversion();

  data.m_strTestName = std::string("add_optional_conversion");

  return data;
}


GeneralizedPEPSData GeneralizedPEPSTester::AddCallExchange()
{ 
  GeneralizedPEPSData data(GetDataWithOptionalConversion());
  data.m_strTestName = std::string("add_call_exchange");

  if ( !data.m_pDeriv->HasAveragingPeriod() )
  {
    ASSERT_MSG(m_pCallExchange, "m_pCallFixedShare must exist");
    data.m_pDeriv->SetGeneralizedPEPSLikeCall(m_pCallExchange);
  }
 

  return data;
}


GeneralizedPEPSData GeneralizedPEPSTester::AddCallFixedCash()
{ 
  GeneralizedPEPSData data(GetDataWithOptionalConversion());
  data.m_strTestName = std::string("add_call_fixed_cash");

  if ( !data.m_pDeriv->HasAveragingPeriod() )
  {
    ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");
    data.m_pDeriv->SetCallFixedCash(m_pCallFixedCash);
  }

  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::AddCallFixedCashWithNotice()
{ 
  GeneralizedPEPSData data(GetDataWithOptionalConversion());
  
  data.m_strTestName = std::string("add_call_fixed_cash_with_notice");

  if ( !data.m_pDeriv->HasAveragingPeriod() )
  {
    ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");

    shared_ptr<CallSchedule> pCalls(new CallSchedule(*m_pCallFixedCash));
    pCalls->SetNoticePeriod(21);
    data.m_pDeriv->SetCallFixedCash(pCalls);
  }


  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::AddCallFixedCashWithTriggerPeriod()
{ 

  GeneralizedPEPSData data(GetDataWithOptionalConversion());

  data.m_strTestName = std::string("add_call_fixed_cash_with_trigger_period");

  if ( !data.m_pDeriv->HasAveragingPeriod())
  {
    ASSERT_MSG(m_pCallFixedCash, "m_pCallFixedCash must exist");
    shared_ptr<CallSchedule> pCalls(new CallSchedule(*m_pCallFixedCash));
    pCalls->SetTriggerCheckPeriod(21, 0);
    data.m_pDeriv->SetCallFixedCash(pCalls);
  }  

  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::AddCrossCurrency()
{  
  GeneralizedPEPSData data(GetDataWithOptionalConversion());

  data.m_strTestName = std::string("add_cross_currency");

  if (  !data.m_pDeriv->HasAveragingPeriod() )
  {
    ASSERT_MSG(m_pSpotFXRates && m_pMoneyMarket, 
               "cross currency property undefined");

    data.m_pDeriv->SetNumeraire( m_pMoneyMarket->GetNumeraire() );

    data.m_pSessionData->GetRateData()->SetYieldCurve(
      m_pMoneyMarket->GetNumeraire(), m_pMoneyMarket->GetYieldCurve() );
  }

  data.m_pDeriv->SetConvertIntoNewShare(false);

  data.m_pSessionData->GetRateData()->SetSpotFXRates(m_pSpotFXRates);

  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::AddDividends()
{ 
  ASSERT_MSG(m_pDividends, "dividends must exist");

  GeneralizedPEPSData data(GetDataWithOptionalConversion());

  // Dividends added
  shared_ptr<Dividends> pDividends(new Dividends(*m_pDividends));

  data.m_pSessionData->GetEquity()->SetDividends(pDividends);

  data.m_strTestName = std::string("add_dividends");

  return data;
}

GeneralizedPEPSData GeneralizedPEPSTester::AddNewShare()
{
  GeneralizedPEPSData data(GetDataWithOptionalConversion());

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

GeneralizedPEPSData GeneralizedPEPSTester::AddNewShareAndCrossCurrency()
{ 
  ASSERT_MSG(m_pSpotFXRates && m_pMoneyMarket, 
             "cross currency property undefined");

  GeneralizedPEPSData data(GetDataWithOptionalConversion());

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
