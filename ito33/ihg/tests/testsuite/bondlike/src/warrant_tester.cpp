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
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"
#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/cocotype.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/dividends.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#include "ito33/numeric/predicatedouble.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "warrant_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void WarrantTester::Setup(std::string strInputFilename)
{
  try
   {
    m_strInputFile = strInputFilename;

    ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

    m_pSessionDataInit = reader.ReadSessionData();

    finance::BondLikeVisitor visitor;
    reader.ReadDerivatives(visitor);
    m_pWarrantInit = visitor.GetAttachedWarrantConvertibleBond();

    m_pCallScheduleInit = make_ptr
        ( new finance::CallSchedule(*m_pWarrantInit->GetCallSchedule()) );

    if ( !m_pWarrantInit )
      throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no warrant in input xml file");

    m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
    reader.ReadTheoreticalModel(m_pModelInit);
  }
  catch (const ito33::Exception& e)
  {
    std::cout << e.GetErrorMessage() << "\n";
    throw e;
  }
}

void WarrantTester::RestoreInitData(WarrantData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  data.m_pDeriv = m_pWarrantInit;

  shared_ptr<finance::CallSchedule> pCallTmp = 
    shared_ptr<finance::CallSchedule> 
    ( new finance::CallSchedule(*m_pCallScheduleInit) );

  data.m_pDeriv->SetCallSchedule( pCallTmp );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

WarrantData WarrantTester::GetBasicData()
{ 
  WarrantData data;

  RestoreInitData(data);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}


WarrantData WarrantTester::AddCallNoticeToBasicData()
{
  WarrantData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  data.m_strTestName = std::string(XML_IHGTEST_CB_CALLNOTICE_TO_BASIC);

  //do not run the call notice test
  //when there is a reset data
 if ( data.m_pDeriv->GetShareDependentConversion()->HasResetDate() )
    return data;

  xml::node::const_iterator
    pNode = nodeRoot.find(XML_IHGTEST_CB_CALLNOTICE_PERIOD);

  if ( pNode != nodeRoot.end() )
  {
    shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

    pCallSchedule->SetNoticePeriod(GetLongFromNode(*pNode));
  }

  return data;
}

WarrantData WarrantTester::DividendOnResetDate()
{
  WarrantData data;

  RestoreInitData(data);

  shared_ptr<ShareDependentConversion> pConv = 
  data.m_pDeriv->GetShareDependentConversion();

  data.m_strTestName = std::string(XML_IHGTEST_WARRANT_DIVIDENDS);

  if ( !pConv->HasResetDate() )
    return data;

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  double dSpot = data.m_pSessionData->GetSpotSharePrice();
  shared_ptr<Numeraire> 
    pCurrency = data.m_pSessionData->GetEquity()->GetNumeraire();

  shared_ptr<Equity> pEquity ( new Equity(dSpot, pCurrency) );
  
  pEquity->SetBorrowCurve(data.m_pSessionData->GetEquity()->GetBorrowCurve());

  xml::node::const_iterator pNode;
  if ( (pNode = nodeRoot.find(XML_IHGTEST_WARRANT_DIVIDENDS )) != nodeRoot.end() )
    pEquity->SetDividends(ReadDividends(*pNode));
  else
    return data;

  shared_ptr<Issuer> pIssuer = data.m_pSessionData->GetEquity()->GetIssuer();
  pEquity->SetIssuer(pIssuer);

  Date valuationDate = data.m_pSessionData->GetValuationDate();

  shared_ptr<RateData> pRateData = data.m_pSessionData->GetRateData();

  data.m_pSessionData = make_ptr
      ( new SessionData(pRateData, pEquity, valuationDate) );

  return data;
  
}

WarrantData WarrantTester::IsLastTriggerConditionMet()
{
  WarrantData data;

  RestoreInitData(data);

  data.m_strTestName = std::string(XML_IHGTESTS_WARRANT_LASTTRIGGERCONDMET);

  shared_ptr<ShareDependentConversion> 
    pConv ( new ShareDependentConversion
                ( *m_pWarrantInit->GetShareDependentConversion() ) 
          );

  //no coco so the code does not do anything
  if ( !pConv->HasCoCo() )
    return data;

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  bool bIsLastTriggerConditionMet  = false;

  if ( nodeRoot.find(XML_IHGTESTS_WARRANT_LASTTRIGGERCONDMET) != nodeRoot.end() ) 
    bIsLastTriggerConditionMet = GetBoolFromName
        ( nodeRoot, XML_IHGTESTS_WARRANT_LASTTRIGGERCONDMET );
  else
    return data;

  //somehow the tests is doing nothing
  //so make the code fail
  if ( pConv->GetIsLastTriggerConditionMet() == bIsLastTriggerConditionMet )
    CPPUNIT_FAIL("The last trigger condition met is the same as the one" \
    " already in the input xml file. Please rectify.");

  pConv->SetCoCo(pConv->GetTrigger(), pConv->GetCoCoType(),
                 pConv->GetChangeRate(), pConv->GetExtremeTrigger(),
                 bIsLastTriggerConditionMet);
  
  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  data.m_pDeriv = make_ptr
      ( new AttachedWarrantConvertibleBond(pBondTerms,pConv) );

  if ( m_pWarrantInit->GetCallSchedule() )
    data.m_pDeriv->SetCallSchedule(m_pWarrantInit->GetCallSchedule());

  if ( m_pWarrantInit->GetPutSchedule() )
    data.m_pDeriv->SetPutSchedule(m_pWarrantInit->GetPutSchedule());

  if ( m_pWarrantInit->IsCrossCurrency() )
  {
    data.m_pDeriv->SetTriggerInCurrencyOf(
      m_pWarrantInit->GetTriggerInCurrencyOf() );
    
    if ( m_pWarrantInit->GetFixedFXRate() > 0. )
      data.m_pDeriv->SetFixedFXRate(m_pWarrantInit->GetFixedFXRate());
  }

  data.m_pDeriv->SetConversionTriggerAsPercentageOf(m_pWarrantInit->GetConversionTriggerAsPercentageOf());

  return data;

}

double WarrantTester::Price(
              const shared_ptr<finance::ShareDependentConversion> &pShareDeConv,
              const shared_ptr<finance::BondTerms> &pBondTerms,
              std::string &xmlOutputFile, 
              std::string sTestName)
{
  WarrantData data;
  RestoreInitData(data);

  data.m_pDeriv = make_ptr
     ( new AttachedWarrantConvertibleBond(pBondTerms, pShareDeConv) );

  if ( m_pWarrantInit->GetCallSchedule() )
    data.m_pDeriv->SetCallSchedule(m_pWarrantInit->GetCallSchedule());

  if ( m_pWarrantInit->GetPutSchedule() )
    data.m_pDeriv->SetPutSchedule(m_pWarrantInit->GetPutSchedule());

  if ( m_pWarrantInit->IsCrossCurrency() )
  {
    data.m_pDeriv->SetTriggerInCurrencyOf(
      m_pWarrantInit->GetTriggerInCurrencyOf());

    if ( m_pWarrantInit->GetFixedFXRate() > 0. )
      data.m_pDeriv->SetFixedFXRate(m_pWarrantInit->GetFixedFXRate());
  }

  data.m_pDeriv->SetConversionTriggerAsPercentageOf
      ( m_pWarrantInit->GetConversionTriggerAsPercentageOf() );

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

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

void WarrantTester::CompleteShareDependentConversion(
    shared_ptr<finance::ShareDependentConversion> &pShareDeConv)
{

  if ( m_pWarrantInit->GetShareDependentConversion()->HasResetDate() )
    pShareDeConv->SetResetDate
      ( 
       m_pWarrantInit->GetShareDependentConversion()->GetResetDate() 
      );

  if ( m_pWarrantInit->GetShareDependentConversion()->GetCurrentRatio() > 0)
   pShareDeConv->SetCurrentRatio(
    m_pWarrantInit->GetShareDependentConversion()->GetCurrentRatio() 
    );

  if ( m_pWarrantInit->GetShareDependentConversion()->GetFixedStrike() >= 0. )
    pShareDeConv->SetFixedStrike
      ( 
        m_pWarrantInit->GetShareDependentConversion()->GetFixedStrike() 
      );
      
  pShareDeConv->SetKeepAccrued(
    m_pWarrantInit->GetShareDependentConversion()->GetKeepAccrued() );
 
  pShareDeConv->SetForfeitCoupon(
    m_pWarrantInit->GetShareDependentConversion()->GetForfeitCoupon() );
}


void WarrantTester::TestCapRatioIncreases()
{
 /*
   As cap ratio increases the price should not decrease
 */

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  double dCapRatio    = 
    m_pWarrantInit->GetShareDependentConversion()->GetCapRatio();
  
  double dCapRatioMin = 
    m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio();
  
  //avoid too big loop when there is not cap ratio (e.g. cap ratio = 1.e8)
  double dCapRatioMax  = std::min(100.*dCapRatioMin,dCapRatio);
  double dCapRatioStep = .15*dCapRatioMin;
  dCapRatio            = dCapRatioMin;
 
  shared_ptr<ShareDependentConversion> 
    pShareDeConv
    (
      new ShareDependentConversion
      ( 
        m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
        m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
        m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio(), 
        m_pWarrantInit->GetShareDependentConversion()
        ->GetIncrementalShareFactor()  
      )
    );

  double dCurrentRatio =  
    m_pWarrantInit->GetShareDependentConversion()->GetCurrentRatio();
     
  if ( dCurrentRatio > 0 )
    pShareDeConv->SetCurrentRatio( dCurrentRatio );
  
  pShareDeConv->SetCapRatio(dCapRatio);

  CompleteShareDependentConversion(pShareDeConv);

  std::string xmlOutputFile;

  double dOldPrice = Price( pShareDeConv, pBondTerms, 
                            xmlOutputFile, "TestCapRatioIncreases");

  while ( dCapRatio < dCapRatioMax + dCapRatioStep)
  {
    dCapRatio += dCapRatioStep;

    shared_ptr<ShareDependentConversion> 
      pShareDeConv
      (
        new ShareDependentConversion
        ( 
          m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
          m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
          m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio(), 
          m_pWarrantInit->GetShareDependentConversion()
          ->GetIncrementalShareFactor()
        )
      );
      
    if ( dCurrentRatio > 0 )
      pShareDeConv->SetCurrentRatio( dCurrentRatio );

    pShareDeConv->SetCapRatio(dCapRatio);
     
    CompleteShareDependentConversion(pShareDeConv);

    double dPrice = Price(pShareDeConv, pBondTerms, 
                          xmlOutputFile, "TestCapRatioIncreases");

    CPPUNIT_ASSERT( numeric::IsEqualOrGreater(dPrice, dOldPrice) );

    //if old price is equal to new price no point
    //in continuing as the price will not increase any longer
    if ( numeric::IsEqual(dPrice, dOldPrice) )
      break;

    dOldPrice = dPrice;
  }

  unlink( xmlOutputFile.c_str() );
}

void WarrantTester::TestStrike()
{

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  shared_ptr<ShareDependentConversion> 
    pConv( 
            new ShareDependentConversion
                ( *m_pWarrantInit->GetShareDependentConversion() )
         );

  double dFixedStrike     = pConv->GetFixedStrike();
  
  if (dFixedStrike < 0 )
    return;

  double dFixedStrikeMin  = .1*dFixedStrike;
  double dFixedStrikeMax  = 3.*dFixedStrike;
  double dFixedStrikeStep = .2*dFixedStrike;
  dFixedStrike = dFixedStrikeMin;

  pConv->SetFixedStrike(dFixedStrike);
  
  std::string xmlOutputFile;

  double dOldPrice = Price( pConv, pBondTerms, xmlOutputFile, "TestStrike");

  while ( dFixedStrike < dFixedStrikeMax + dFixedStrikeStep )
  {
    dFixedStrike += dFixedStrikeStep;
    
    pConv->SetFixedStrike(dFixedStrike);

    double dPrice = Price( pConv, pBondTerms, xmlOutputFile, "TestStrike");

    //as strike increases price should not increase
    CPPUNIT_ASSERT( numeric::IsEqualOrGreater(dOldPrice, dPrice) );

    dOldPrice = dPrice;
  }

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 );
}
 

void WarrantTester::TestShareFactor()
{
  /*
   As incremental share factor increases price should increase
 */

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  double dParam     = 
    m_pWarrantInit->GetShareDependentConversion()->GetIncrementalShareFactor();
  double dParamMin  = .1*dParam;
  double dParamMax  = 5.*dParam;
  double dParamStep = .2*dParam;
  dParam            = dParamMin;

  shared_ptr<ShareDependentConversion> 
    pShareDeConv
    (
      new ShareDependentConversion 
      (
      m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
      m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
      m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio(), 
      dParam
      )
    );

  double dCurrentRatio = 
    m_pWarrantInit->GetShareDependentConversion()->GetCurrentRatio();
  
  if ( dCurrentRatio > 0 )
    pShareDeConv->SetCurrentRatio( dCurrentRatio );

  pShareDeConv->SetCapRatio
    ( m_pWarrantInit->GetShareDependentConversion()->GetCapRatio() );

  CompleteShareDependentConversion(pShareDeConv);

  std::string xmlOutputFile;

  double dOldPrice = Price( pShareDeConv, pBondTerms , 
    xmlOutputFile, "TestShareFactor");

  while ( dParam < dParamMax + dParamStep)
  {
    dParam += dParamStep;

    shared_ptr<ShareDependentConversion> 
      pShareDeConv
      (
        new ShareDependentConversion 
        (
        m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
        m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
        m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio(), 
        dParam
        )
      );
            
    if ( dCurrentRatio > 0 )
      pShareDeConv->SetCurrentRatio( dCurrentRatio );

    pShareDeConv->SetCapRatio
      ( m_pWarrantInit->GetShareDependentConversion()->GetCapRatio() );

    CompleteShareDependentConversion(pShareDeConv);

    double dPrice = Price( pShareDeConv, pBondTerms,
      xmlOutputFile, "TestShareFactor");

    
    //as param increases price should increase
    CPPUNIT_ASSERT( numeric::IsEqualOrGreater(dPrice, dOldPrice) );

    dOldPrice = dPrice;
  }

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 );

}
  
void WarrantTester::TestBaseRatio()
{

  /*
   As the bas ratio increase price should increase
 */

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  double dParam     = m_pWarrantInit->GetShareDependentConversion()->GetBaseRatio();
  double dParamMin  = .1*dParam;
  double dParamMax  = std::min
     (
     m_pWarrantInit->GetShareDependentConversion()->GetCapRatio(),
     2.*dParam
     );
  double dParamStep = .2*dParam;
  dParam            = dParamMin;
  
  shared_ptr<ShareDependentConversion> 
    pShareDeConv
    (
      new ShareDependentConversion 
      (
        m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
        m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
        dParam,
        m_pWarrantInit->GetShareDependentConversion()
        ->GetIncrementalShareFactor()
      )
    );

  double dCurrentRatio = 
    m_pWarrantInit->GetShareDependentConversion()->GetCurrentRatio();

  if ( dCurrentRatio > 0 )
    pShareDeConv->SetCurrentRatio( dCurrentRatio );

  pShareDeConv->SetCapRatio
    ( m_pWarrantInit->GetShareDependentConversion()->GetCapRatio() );

  CompleteShareDependentConversion(pShareDeConv);

  std::string xmlOutputFile;

  double dOldPrice = Price(pShareDeConv, pBondTerms, 
    xmlOutputFile, "TestBaseRatio");

  while ( dParam < dParamMax )
  {
   
    dParam += dParamStep;
    dParam = std::min(dParam, dParamMax);

    shared_ptr<ShareDependentConversion> 
      pShareDeConv
      (
        new ShareDependentConversion 
        (
          m_pWarrantInit->GetShareDependentConversion()->GetStartDate(), 
          m_pWarrantInit->GetShareDependentConversion()->GetEndDate(),
          dParam,
          m_pWarrantInit->GetShareDependentConversion()
          ->GetIncrementalShareFactor()
        )
      );
       
    if ( dCurrentRatio > 0 )
      pShareDeConv->SetCurrentRatio( dCurrentRatio );

    pShareDeConv->SetCapRatio
      ( m_pWarrantInit->GetShareDependentConversion()->GetCapRatio() );

    CompleteShareDependentConversion(pShareDeConv);

    double dPrice = Price( pShareDeConv, pBondTerms,
    xmlOutputFile, "TestBaseRatio");

    //as param increases price should increase
    CPPUNIT_ASSERT( numeric::IsEqualOrGreater(dPrice, dOldPrice) );

    dOldPrice = dPrice;
  }

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0);

}
  
void WarrantTester::TestResetDate()
{
/*
   As the reset date is increased
   from the valuation date to the original
   reset date price should increase
 */

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  shared_ptr<ShareDependentConversion> 
    pConv ( new ShareDependentConversion
                ( *m_pWarrantInit->GetShareDependentConversion() )
          );

  //if no reset date is specified no point in doing the test
  if ( !pConv->HasResetDate() )
    return;

  Date valuationDate     = m_pSessionDataInit->GetValuationDate();
  Date originalResetDate = pConv->GetResetDate();
  Date resetDate         = valuationDate;
  int resetDateStepper   = 6;
  resetDate.AddMonths(resetDateStepper);

  pConv->SetResetDate(resetDate);

  std::string xmlOutputFile;

  double dOldPrice = Price( pConv, pBondTerms, xmlOutputFile, "TestResetDate");

  while ( resetDate < originalResetDate )
  {
    if ( resetDate.AddMonths(resetDateStepper) > originalResetDate )
      resetDate = originalResetDate;
    else
      resetDate.AddMonths(resetDateStepper);
    
    pConv->SetResetDate( resetDate );

    double dPrice = Price( pConv, pBondTerms, xmlOutputFile, "TestResetDate");

    //as param increases price should increase

    //std::cout << "dPrice: " << dPrice << " dOldPrice: " << dOldPrice << std::endl;
    CPPUNIT_ASSERT( dPrice > dOldPrice );

    dOldPrice = dPrice;
  }

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 );
}
  
void WarrantTester::TestLastTriggerConditionMetProperties()
{

  shared_ptr<BondTerms> 
    pBondTerms ( new BondTerms(*m_pWarrantInit->GetBondTerms()) );

  shared_ptr<ShareDependentConversion> 
    pConv ( new ShareDependentConversion
                ( *m_pWarrantInit->GetShareDependentConversion() )
          );

  if ( !pConv->HasCoCo() )
    return;


  shared_ptr<ShareDependentConversion> 
    pNewShare
    (
      new ShareDependentConversion(
      pConv->GetStartDate(), pConv->GetEndDate(), pConv->GetBaseRatio(),
      pConv->GetIncrementalShareFactor() )
    );

  pNewShare->SetCapRatio( pConv->GetCapRatio() );

  CompleteShareDependentConversion(pNewShare);

}
  
void WarrantTester::TestOneDvsTwoD()
{

}

