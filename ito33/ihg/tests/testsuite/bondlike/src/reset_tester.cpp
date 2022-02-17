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
#include "ito33/finance/issuer.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/spotfxrates.h"

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

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/moneymarket.h"
#include "ito33/xml/finance/spotfxrates.h"
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
#include "reset_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void ResetTester::Setup(std::string strInputFilename)
{
  try
   {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pResetInit = visitor.GetReset();

    if(!m_pResetInit)
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no reset in input xml file");

    m_pModelInit = make_ptr( new ito33::ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);
    
    // cross currency
    const xml::node& nodeRoot = reader.GetRootNode();

    xml::node::const_iterator
      pNode = nodeRoot.find(XML_TAG_MONEYMARKET_ROOT);
    
    if( pNode != nodeRoot.end() )
    {
      m_pMoneyMarket = GetMoneyMarketFromNode(*pNode);
    }

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

double ResetTester::Price(const shared_ptr<finance::BondTerms> &pBondTerms,
              const shared_ptr<finance::ResetConversionSchedule> &pResetSchedule,
              std::string &xmlOutputFile, 
              std::string sTestName, Date valuationDate)
{
  ResetData data;
  RestoreInitData(data);

  data.m_pDeriv = make_ptr( new finance::Reset(pBondTerms, pResetSchedule) );

  if ( m_pResetInit->GetCallSchedule() )
    data.m_pDeriv->SetCallSchedule(m_pResetInit->GetCallSchedule());

  if ( m_pResetInit->GetPutSchedule() )
    data.m_pDeriv->SetPutSchedule(m_pResetInit->GetPutSchedule());

  if ( m_pResetInit->IsCrossCurrency() )
  {
    data.m_pDeriv->SetTriggerInCurrencyOf( 
      m_pResetInit->GetTriggerInCurrencyOf() );
    
    if ( m_pResetInit->GetFixedFXRate() > 0. )
      data.m_pDeriv->SetFixedFXRate(m_pResetInit->GetFixedFXRate());
  }

  data.m_pDeriv->SetConversionTriggerAsPercentageOf
      ( m_pResetInit->GetConversionTriggerAsPercentageOf() );

  data.m_pSessionData->SetValuationDate( valuationDate );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  if ( m_pResetInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pResetInit->GetNumeraire() );

  xmlOutputFile = data.GetXMLOutputFile(sTestName);

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);

  try
  {
    return data.Price();
  }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }

  return 0;
}

void ResetTester::CompleteConversionSchedule(
    shared_ptr<finance::ResetConversionSchedule> &pResetSchedule)
{

   pResetSchedule->SetForfeitCoupon
    (
      m_pResetInit->GetResetConversionSchedule()->GetForfeitCoupon()
    );

   pResetSchedule->SetKeepAccrued
     (
      m_pResetInit->GetResetConversionSchedule()->GetKeepAccrued()
     );

  double dCash = m_pResetInit->GetResetConversionSchedule()->GetCash();
  
  if ( dCash != 0. )
    pResetSchedule->SetCash(dCash);

}


void ResetTester::PriceDecreaseCapWhenNumberOfResetIncreases()
{
 
  ResetFlooredBy resetFlooredBy = ResetFlooredBy_PrevailingConversionPrice;

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms( *m_pResetInit->GetBondTerms() ) );
 
  Date maturityDate  = pBondTerms->GetMaturityDate();
  Date issueDate     = pBondTerms->GetIssueDate();
  Date valuationDate = m_pSessionDataInit->GetValuationDate();

  double dLocalCapRate   = 1.2;
  double dLocalFloorRate = 1.0;
  double dInitialConvPrice = pBondTerms->GetNominal();
  double dCurrentConvPrice = pBondTerms->GetNominal();
  double dMultiplier     = 1.0;

  long lDays = Date::DaysDiff(valuationDate,maturityDate);

  // Setup the reset data
  shared_ptr<ResetConversionSchedule> 
    pResetConversionSchedule 
    (
      new ResetConversionSchedule(issueDate, maturityDate, dInitialConvPrice, 
                                  dCurrentConvPrice, resetFlooredBy)
    );

  CompleteConversionSchedule(pResetConversionSchedule);

  Date resetDate1 = valuationDate;
  resetDate1.AddDays( int(lDays/3) );

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset 
    ( 
      new ConversionPriceReset(resetDate1,  dLocalFloorRate)
    );

  pConversionPriceReset->SetMultiplier(dMultiplier);

  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
 
  std::string xmlOutputFileOneResetDate;

  double dPriceOneReset = Price(pBondTerms, 
                     pResetConversionSchedule, 
                     xmlOutputFileOneResetDate,
                     "PriceDecreaseCapWhenNumberOfResetIncreasesOneResetDate",
                     valuationDate);

  //two resets date
  Date resetDate2 = resetDate1;
  resetDate2.AddDays( int(lDays/3) );
  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate2, dLocalFloorRate ) );

  pConversionPriceReset->SetMultiplier(dMultiplier);

  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  std::string xmlOutputFileTwoResetDate;

  double dPriceTwoReset = Price(pBondTerms, 
                     pResetConversionSchedule, 
                     xmlOutputFileTwoResetDate,
                     "PriceDecreaseCapWhenNumberOfResetIncreasesTwoResetDate",
                     valuationDate);

  CPPUNIT_ASSERT( numeric::IsEqualOrLess(dPriceTwoReset,dPriceOneReset) );
 
  //threee reset dates
  Date resetDate3= resetDate2;
  resetDate3.AddDays( int(lDays/4) );

  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        (resetDate3,dLocalFloorRate) );
  pConversionPriceReset->SetMultiplier(dMultiplier);
  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
  
  std::string xmlOutputFileThreeResetDate;

  double dPriceThreeReset = Price(pBondTerms, 
                   pResetConversionSchedule, 
                   xmlOutputFileThreeResetDate,
                   "PriceDecreaseCapWhenNumberOfResetIncreasesThreeResetDate",
                   valuationDate);

  CPPUNIT_ASSERT( 
                numeric::IsEqualOrLess(dPriceThreeReset, dPriceOneReset)
                );
  CPPUNIT_ASSERT( 
                numeric::IsEqualOrLess(dPriceThreeReset, dPriceTwoReset)
                );
                
  //test is passed so removed xml output
  //otherwise, cpp unit will fail and go to the next
  //test inside the testsuite, that is why we only need
  //to remove the file at the end
  CPPUNIT_ASSERT( unlink( xmlOutputFileOneResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileTwoResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileThreeResetDate.c_str() ) == 0 );

}//ResetAcceptanceTest::PriceIncreaseCapWhenNumberOfResetIncreases()

void ResetTester::PriceIncreaseFloorWhenNumberOfResetIncreases()
{
 
   finance::ResetFlooredBy resetFlooredBy =
     finance::ResetFlooredBy_InitialConversionPrice;

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pResetInit->GetBondTerms()) );

  Date maturityDate  = pBondTerms->GetMaturityDate();
  Date issueDate     = pBondTerms->GetIssueDate();
  Date valuationDate = m_pSessionDataInit->GetValuationDate();

  double dLocalCapRate   = 1.2;
  double dLocalFloorRate = .1;
  double dInitialConvPrice = pBondTerms->GetNominal();
  double dCurrentConvPrice = pBondTerms->GetNominal();
  double dMultiplier     = 1.;

  long lDays = Date::DaysDiff(valuationDate,maturityDate);

  // Setup the reset data
  shared_ptr<ResetConversionSchedule> 
    pResetConversionSchedule 
    (
      new ResetConversionSchedule(issueDate, maturityDate, dInitialConvPrice, 
                                  dCurrentConvPrice, resetFlooredBy)
    );

  CompleteConversionSchedule(pResetConversionSchedule);

  Date resetDate1 = valuationDate;
  resetDate1.AddDays( int(lDays/3) );

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset 
    ( 
      new ConversionPriceReset(resetDate1, dLocalFloorRate)
    );

  pConversionPriceReset->SetMultiplier(dMultiplier);
  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
   
  std::string xmlOutputFileOneResetDate;

  double dPriceOneReset = Price(pBondTerms, 
                   pResetConversionSchedule, 
                   xmlOutputFileOneResetDate,
                   "PriceIncreaseFloorWhenNumberOfResetIncreasesOneResetDate",
                   valuationDate);

  //two resets date
  Date resetDate2 = resetDate1;
  resetDate2.AddDays( int(lDays/3) );

  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate2, dLocalFloorRate ) );
  pConversionPriceReset->SetMultiplier(dMultiplier);
  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  std::string xmlOutputFileTwoResetDate;

  double dPriceTwoReset = Price(pBondTerms, 
                   pResetConversionSchedule, 
                   xmlOutputFileTwoResetDate,
                   "PriceIncreaseFloorWhenNumberOfResetIncreasesTwoResetDate",
                   valuationDate);
 
  CPPUNIT_ASSERT( numeric::IsEqualOrGreater(dPriceTwoReset, dPriceOneReset) );

  //three reset dates
  Date resetDate3 = resetDate2;
  resetDate3.AddDays( int(lDays/4) );

  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate3, dLocalFloorRate ) );
  pConversionPriceReset->SetMultiplier(dMultiplier);
  pConversionPriceReset->SetCap(dLocalCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  std::string xmlOutputFileThreeResetDate;

  double dPriceThreeReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileThreeResetDate,
                 "PriceIncreaseFloorWhenNumberOfResetIncreasesThreeResetDate",
                 valuationDate);

  CPPUNIT_ASSERT( 
                numeric::IsEqualOrGreater(dPriceThreeReset, dPriceOneReset)
                );
  CPPUNIT_ASSERT( 
                numeric::IsEqualOrGreater(dPriceThreeReset, dPriceTwoReset)
                );

  //test is passed so removed xml output
  //otherwise, cpp unit will fail and go to the next
  //test inside the testsuite, that is why we only need
  //to remove the file at the end
  CPPUNIT_ASSERT( unlink( xmlOutputFileOneResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileTwoResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileThreeResetDate.c_str() ) == 0 );

} //PriceIncreaseFloorWhenNumberOfResetIncreases()

void ResetTester::PriceStaySameWithResetDateBeforeValuationDate()
{

  Date newValuationDate =  m_pSessionDataInit->GetValuationDate();
  newValuationDate.AddMonths(1);
 
  finance::ResetFlooredBy resetFlooredBy 
    = finance::ResetFlooredBy_InitialConversionPrice;

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pResetInit->GetBondTerms()) );
 
  Date maturityDate  = pBondTerms->GetMaturityDate();
  Date issueDate     = pBondTerms->GetIssueDate();
  //Date valuationDate = m_pSessionDataInit->GetValuationDate();

  double dLocalFloorRate    = 1.;
  double dInitialConvPrice  = pBondTerms->GetNominal();
  double dCurrentConvPrice  = pBondTerms->GetNominal();
  double dMultiplier        = 1.;
  
  // Setup the reset data
  shared_ptr<ResetConversionSchedule> 
    pResetConversionSchedule 
    ( 
      new ResetConversionSchedule(issueDate, maturityDate, dInitialConvPrice, 
                                  dCurrentConvPrice, resetFlooredBy)
    );

  CompleteConversionSchedule(pResetConversionSchedule);
 
  //Reset date is before valuation Date going forward
  Date resetDate = m_pSessionDataInit->GetValuationDate();
  resetDate.AddDays(1);
 
  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset
    ( 
      new ConversionPriceReset(resetDate, dLocalFloorRate)
    );

  pConversionPriceReset->SetMultiplier(dMultiplier);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
 
  std::string xmlOutputFileOneResetDate;

  double dPriceOneReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileOneResetDate,
                 "PriceStaySameWhenValuationDateBeforeResetDate",
                 newValuationDate);

  std::string xmlOutputFileTwoResetDate;

  
 //two resets date
  Date resetDate2 = resetDate;
  resetDate2.AddDays( 2 );

  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate2,dLocalFloorRate ) );
  pConversionPriceReset->SetMultiplier(dMultiplier);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

 double dPriceTwoReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileTwoResetDate,
                 "PriceStaySameWhenValuationDateBeforeTwoResetDate",
                 newValuationDate);

  CPPUNIT_ASSERT( fabs(dPriceOneReset - dPriceTwoReset) < DOUBLETOLERANCE );

  // Test is passed so removed xml output
  // otherwise, cpp unit will fail and go to the next
  // test inside the testsuite, that is why we only need
  // to remove the file at the end
  CPPUNIT_ASSERT( unlink( xmlOutputFileOneResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileTwoResetDate.c_str() ) == 0 );

}

void ResetTester::PriceStaySameWhenFloorAndCapAreOne()
{

  finance::ResetFlooredBy resetFlooredBy 
    = finance::ResetFlooredBy_InitialConversionPrice;

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pResetInit->GetBondTerms()) );
 
  Date maturityDate  = pBondTerms->GetMaturityDate();
  Date issueDate     = pBondTerms->GetIssueDate();
  Date valuationDate = m_pSessionDataInit->GetValuationDate();

  double dLocalFloorRate = 1.;
  double dInitialConvPrice = pBondTerms->GetNominal();
  double dCurrentConvPrice  = pBondTerms->GetNominal();

  double dMultiplier     = 1.;
  
  long lDays = Date::DaysDiff(valuationDate,maturityDate);
  
  // Setup the reset data
  shared_ptr<ResetConversionSchedule> 
    pResetConversionSchedule 
    ( 
      new ResetConversionSchedule(issueDate, maturityDate, dInitialConvPrice, 
                                  dCurrentConvPrice, resetFlooredBy)
    );

  CompleteConversionSchedule(pResetConversionSchedule);

  Date resetDate1 = valuationDate;
  resetDate1.AddDays( int(lDays/3) );

  shared_ptr<ConversionPriceReset> 
    pConversionPriceReset
    ( 
      new ConversionPriceReset(resetDate1, dLocalFloorRate)
    );

  pConversionPriceReset->SetMultiplier(dMultiplier);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
 
  std::string xmlOutputFileOneResetDate;

  double dPriceOneReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileOneResetDate,
                 "PriceStaySameWhenFloorAndCapAreOneOneResetDate",
                 valuationDate);

  //two resets date
  Date resetDate2 = resetDate1;
  resetDate2.AddDays( int( lDays/3) );

  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate2, dLocalFloorRate ) );
  pConversionPriceReset->SetMultiplier(dMultiplier);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  std::string xmlOutputFileTwoResetDate;

  double dPriceTwoReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileTwoResetDate,
                 "PriceStaySameWhenFloorAndCapAreOneTwoResetDate",
                 valuationDate);
    
  CPPUNIT_ASSERT( 
    fabs(dPriceTwoReset-dPriceOneReset)/dPriceOneReset < 1.e-2 );

  //threee reset dates
  Date resetDate3 = resetDate2;
  resetDate3.AddDays( int(lDays/4) );
  pConversionPriceReset = make_ptr( new ConversionPriceReset
                                        ( resetDate3, dLocalFloorRate ) );
  pConversionPriceReset->SetMultiplier(dMultiplier);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);
    
  std::string xmlOutputFileThreeResetDate;

  double dPriceThreeReset = Price(pBondTerms, 
                 pResetConversionSchedule, 
                 xmlOutputFileThreeResetDate,
                 "PriceStaySameWhenFloorAndCapAreOneThreeResetDate",
                 valuationDate);

  CPPUNIT_ASSERT( 
    fabs(dPriceThreeReset-dPriceTwoReset)/dPriceTwoReset < 1.e-2 );

  CPPUNIT_ASSERT( 
    fabs(dPriceThreeReset-dPriceOneReset)/dPriceOneReset < 1.e-2 );

  //test is passed so removed xml output
  //otherwise, cpp unit will fail and go to the next
  //test inside the testsuite, that is why we only need
  //to remove the file at the end
  CPPUNIT_ASSERT( unlink( xmlOutputFileOneResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileTwoResetDate.c_str() ) == 0 );
  CPPUNIT_ASSERT( unlink( xmlOutputFileThreeResetDate.c_str() ) == 0 );

}//PriceStaySameWhenFloorAndCapAreOne

void ResetTester::RestoreInitData(ResetData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = m_pResetInit;

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pDeriv->SetConvertIntoNewShare(false);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
  /*
  data.m_pSessionData = m_pSessionDataInit;

  // data.m_pDeriv = m_pResetInit;
  shared_ptr<BondTerms> bondTerms = new BondTerms(*m_pResetInit->GetBondTerms());
  shared_ptr<CallSchedule> calls = new CallSchedule(*m_pResetInit->GetCallSchedule());
  shared_ptr<ConversionSchedule> 
  conversions = new ConversionSchedule(*m_pResetInit->GetConversionSchedule());
  
  data.m_pDeriv = new finance::ConvertibleBond(bondTerms, conversions);
  data.m_pDeriv->SetCallSchedule(calls);

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
  */
}

ResetData ResetTester::GetBasicData()
{ 
  ResetData data;

  RestoreInitData(data);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

ResetData ResetTester::AddDividends()
{ 
  ASSERT_MSG(m_pDividends, "dividends must exist");

  ResetData data(GetBasicData());

  // Dividends added
  shared_ptr<Dividends> pDividends(new Dividends(*m_pDividends));

  data.m_pSessionData->GetEquity()->SetDividends(pDividends);

  data.m_strTestName = std::string("add_dividends");

  return data;
}

ResetData ResetTester::AddNewShare()
{ 
  ResetData data(GetBasicData());

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

ResetData ResetTester::AddCrossCurrency()
{ 
  ASSERT_MSG(m_pMoneyMarket, "cross currency property undefined");

  ResetData data(GetBasicData());

  data.m_pDeriv->SetNumeraire( m_pMoneyMarket->GetNumeraire() );

  data.m_pSessionData->GetRateData()->SetYieldCurve( 
    m_pMoneyMarket->GetNumeraire(), m_pMoneyMarket->GetYieldCurve() );

  double dFixedRate = .8126;

  data.m_pDeriv->SetTriggerInCurrencyOf( TriggerInCurrencyOf_Underlying );
      
  data.m_pDeriv->SetFixedFXRate(dFixedRate);

  data.m_strTestName = std::string("add_cross_currency");

  return data;
}

ResetData ResetTester::AddNewShareAndCrossCurrency()
{ 
  ASSERT_MSG(m_pMoneyMarket, "cross currency property undefined");

  ResetData data(GetBasicData());

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
  
  double dFixedRate = 1.1;

  data.m_pDeriv->SetTriggerInCurrencyOf( TriggerInCurrencyOf_Underlying );
      
  data.m_pDeriv->SetFixedFXRate(dFixedRate);

  data.m_strTestName = std::string("add_newshare_and_cross_currency");

  return data;
}

bool CBUsingResetPricerTester::Compare()
{
  

  ResetData reset = GetBasicData();
  CBData cb = GetEquivalentCB(reset.m_pSessionData->GetValuationDate());

  if( fabs(reset.Price() - cb.Price() ) > 0.03)
    return false;

  double dSpot = reset.m_pSessionData->GetSpotSharePrice();

  reset.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot * 0.5);
  cb.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot * 0.5);
  
  if( fabs(reset.Price() - cb.Price() ) > 0.03)
    return false;

  reset.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot * 1.5);
  cb.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot * 1.5);
  
  if( fabs(reset.Price() - cb.Price() ) > 0.03)
    return false;

  return true;
}


CBData ResetTester::GetEquivalentCB(Date valuationDate)
{
  CBData data;

  data.m_pSessionData = make_ptr( new SessionData(*m_pSessionDataInit) );
  data.m_pSessionData->SetValuationDate( valuationDate );

  shared_ptr<ResetConversionSchedule>
    rcs = m_pResetInit->GetResetConversionSchedule();

  shared_ptr<ConversionSchedule> pConv ( new ConversionSchedule() );

  double dCurrentConversionPrice = rcs->GetCurrentConversionPrice();

  if ( m_pResetInit->IsCrossCurrency() )
  {
    double dFixedFXRate = m_pResetInit->GetFixedFXRate();

    dCurrentConversionPrice *= dFixedFXRate;
  }

  double dCurrentRatio = m_pResetInit->GetBondTerms()->GetNominal()
    / dCurrentConversionPrice;

  pConv->AddConversionPeriod
            ( shared_ptr<ConversionPeriod> 
                (new ConversionPeriod
                    ( rcs->GetStartDate(), rcs->GetEndDate(), 
                      dCurrentRatio ))
            );

  data.m_pDeriv = make_ptr( new finance::ConvertibleBond
                                ( m_pResetInit->GetBondTerms(), pConv ) );

  if(m_pResetInit->GetCallSchedule())
    data.m_pDeriv->SetCallSchedule(m_pResetInit->GetCallSchedule());

  if (m_pResetInit->GetPutSchedule())
    data.m_pDeriv->SetPutSchedule(m_pResetInit->GetPutSchedule());

  if ( m_pResetInit->IsCrossCurrency() )
  {
    data.m_pDeriv->SetTriggerInCurrencyOf(
      m_pResetInit->GetTriggerInCurrencyOf() );
    
    if ( m_pResetInit->GetFixedFXRate() > 0. )
      data.m_pDeriv->SetFixedFXRate(m_pResetInit->GetFixedFXRate());
  }
    
  data.m_pDeriv->SetConversionTriggerAsPercentageOf
      ( m_pResetInit->GetConversionTriggerAsPercentageOf() );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  if ( m_pResetInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pResetInit->GetNumeraire() );

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;

  return data;
}

//-------------------------------------------------------------------------
// Reset1DModelTester
//-------------------------------------------------------------------------

bool Reset1DModelTester::Compare()
{
  std::cout.precision(13);

  ResetData reset = GetBasicData();

  reset.m_pModel->SetDebugOutputFile("ihg.xml");

  double dPrice1D = reset.Price();
  std::cout << "\n" << dPrice1D << "\n";

  if(reset.m_pSessionData->GetDividends())
    reset.m_pSessionData->GetDividends()->Add
                      (Dividend::Cash, Date(1900, Date::Jan, 1), 1);
  else
  {
    shared_ptr<Dividends> pDivs ( new Dividends() );
    pDivs->Add(Dividend::Cash, Date(1900, Date::Jan, 1), 1);
    reset.m_pSessionData->GetEquity()->SetDividends(pDivs);
  }

  reset.m_pModel->SetDebugOutputFile("ihg2.xml");

  double dPrice2D = reset.Price();

  std::cout << dPrice2D << std::endl;

  if( fabs(dPrice1D - dPrice2D ) > 0.03)
    return false;

  return true;
}

