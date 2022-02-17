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
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"

#include "ito33/xml/finance/bondlike/callschedule.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"


#include "ito33/finance/theoreticalmodel.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/tests/convergence_parameter_value.h"

#include "ihg/tests/testdata.h"

#include "test_bondlike_common.h"

// local files
#include "bondlike_reader.h"
#include "basic_features_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void BasicFeaturesTester::Setup(std::string strInputFilename)
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


void BasicFeaturesTester::RestoreInitData(CBData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  // data.m_pDeriv = m_pCBInit;
  shared_ptr<BondTerms> bondTerms ( new BondTerms(*m_pCBInit->GetBondTerms()) );
  
  shared_ptr<CallSchedule> 
    calls ( new CallSchedule(*m_pCBInit->GetCallSchedule()) );
  
  shared_ptr<ConversionSchedule> 
    conversions ( new ConversionSchedule(*m_pCBInit->GetConversionSchedule()) );
  
  data.m_pDeriv = make_ptr( new finance::ConvertibleBond
                                (bondTerms, conversions) );
  data.m_pDeriv->SetCallSchedule(calls);

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  if ( m_pCBInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pCBInit->GetNumeraire() );

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

CBData BasicFeaturesTester::GetBasicData()
{ 
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();
  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_BASICDDATA);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

CBData BasicFeaturesTester::AddSoftCall()
{ 
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_CB_TRIGGER);

  xml::node::const_iterator pNodePeriod;

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();
  for(pNodePeriod = pNode->begin(); pNodePeriod != pNode->end(); pNodePeriod++)
  {
    shared_ptr<finance::CallPeriod> pCallPeriod;
    if(XML::Restore(*pNodePeriod, pCallPeriod))
        pCallSchedule->AddCallPeriod(pCallPeriod);
  }

  data.m_strTestName = std::string(XML_IHGTEST_CB_TRIGGER);


  return data;
}

CBData BasicFeaturesTester::AddCallClaimTriggerAsPercentageOfToBasicData()
{
   CBData data;

  RestoreInitData(data);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();
 
  if ( pCallSchedule )
  {
   pCallSchedule->SetTriggerAsPercentageOf(finance::TriggerAsPercentageOf_Claim);

  CallSchedule::Elements::iterator iter;
  CallSchedule::Elements element = pCallSchedule->GetAll();

  for ( iter = element.begin(); iter != element.end(); iter++)
  {
    (*iter)->SetTrigger(1.2);
  }

  }

  data.m_strTestName = std::string(XML_IHGTEST_CB_TRIGGERASPERCENTAGEOF_CLAIM);


  return data;

}

CBData BasicFeaturesTester::AddCallIssuePriceTriggerAsPercentageOfToBasicData()
{
   CBData data;

  RestoreInitData(data);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  if ( pCallSchedule )
  {
    pCallSchedule->SetTriggerAsPercentageOf(finance::TriggerAsPercentageOf_IssuePrice);
   
    CallSchedule::Elements::iterator iter;
    CallSchedule::Elements element = pCallSchedule->GetAll();

    for ( iter = element.begin(); iter != element.end(); iter++)
    {
      (*iter)->SetTrigger(1.1);
    }
  }

  data.m_strTestName = 
    std::string(XML_IHGTEST_CB_TRIGGERASPERCENTAGEOF_ISSUEPRICE );


  return data;

}

CBData BasicFeaturesTester::AddCallPrincipalTriggerAsPercentageOfToBasicData()
{
   CBData data;

  RestoreInitData(data);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  if ( pCallSchedule )
  {
    pCallSchedule->SetTriggerAsPercentageOf(finance::TriggerAsPercentageOf_Principal);
   
    CallSchedule::Elements::iterator iter;
    CallSchedule::Elements element = pCallSchedule->GetAll();

    for ( iter = element.begin(); iter != element.end(); iter++)
    {
      (*iter)->SetTrigger(.9);
    }
  }

  data.m_strTestName = 
    std::string(XML_IHGTEST_CB_TRIGGERASPERCENTAGEOF_PRINCIPAL);


  return data;

}


CBData BasicFeaturesTester::AddPVCouponMakeWhole()
{
  CBData data(AddSoftCall());

  data.m_pDeriv->GetCallSchedule()->SetCouponMakeWhole(true);

  data.m_strTestName = std::string(XML_IHGTEST_CB_PVCOUPON_MAKEWHOLE);

  return data;

}

CBData BasicFeaturesTester::AddNoPVCouponMakeWhole()
{
  CBData data(AddSoftCall());

  data.m_pDeriv->GetCallSchedule()->SetCouponMakeWhole(false);

  data.m_strTestName = std::string(XML_IHGTEST_CB_NOPVCOUPON_MAKEWHOLE);

  return data;
}

CBData BasicFeaturesTester::AddPremiumMakeWhole()
{
  CBData data(AddSoftCall());

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  xml::node::const_iterator
    pNode = nodeRoot.find(XML_IHGTEST_CB_PREMIUM_MAKEWHOLE);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  pCallSchedule->SetPremiumMakeWhole(GetDoubleFromNode(*pNode));

  data.m_strTestName = std::string(XML_IHGTEST_CB_PREMIUM_MAKEWHOLE);

  return data;
}

CBData BasicFeaturesTester::AddCallNoticeToBasicData()
{
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  xml::node::const_iterator
    pNode = nodeRoot.find(XML_IHGTEST_CB_CALLNOTICE_PERIOD);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  pCallSchedule->SetNoticePeriod(GetLongFromNode(*pNode));

  data.m_strTestName = std::string(XML_IHGTEST_CB_CALLNOTICE_TO_BASIC);

  return data;
}

CBData BasicFeaturesTester::AddCallNoticeToAddSoftCallData()
{
  CBData data(AddSoftCall());

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  xml::node::const_iterator
    pNode = nodeRoot.find(XML_IHGTEST_CB_CALLNOTICE_PERIOD);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  pCallSchedule->SetNoticePeriod(GetLongFromNode(*pNode));

  data.m_strTestName = std::string(XML_IHGTEST_CB_CALLNOTICE_TO_SOFT);

  return data;
}


void BasicFeaturesTester::TestPriceIncreasesAsYTCIncreases()
{
 
  double dYieldToCallMax   = .18;
  double dYieldToCallStart = 1.e-8;
  double dYieldToCallStep  = .01;
  double dYieldToCall      = dYieldToCallStart;
  double dPriceOld         = 0.0;

  Date maturityDate  = m_pCBInit->GetMaturityDate();
  Date issueDate     = m_pCBInit->GetBondTerms()->GetIssueDate();

  CBData data;

  data.m_pSessionData = m_pSessionDataInit;
  data.m_pModel       = m_pModelInit;
  data.m_strFileName  = m_strInputFile;

  // data.m_pDeriv = m_pCBInit;
  shared_ptr<BondTerms> bondTerms ( new BondTerms(*m_pCBInit->GetBondTerms()) );
  
  shared_ptr<ConversionSchedule> 
    conversions 
    ( 
      new ConversionSchedule( *m_pCBInit->GetConversionSchedule() )
    );
  
  data.m_pDeriv = make_ptr( new finance::ConvertibleBond
                                (bondTerms, conversions) );
 
  data.m_pDeriv->SetSessionData(data.m_pSessionData);
  
  if ( m_pCBInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pCBInit->GetNumeraire() );

  std::string xmlOutputFile = 
    data.GetXMLOutputFile("TestPriceIncreasesAsYTCIncreases");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);

  while ( dYieldToCall < dYieldToCallMax + dYieldToCallStep )
  { 
    shared_ptr<CallPeriod> 
      pCallPeriod 
      (     
        CallPeriod::CreateWithYield(issueDate, maturityDate, dYieldToCall)
      );

    shared_ptr<CallSchedule> pCalls ( new CallSchedule() );
    pCalls->AddCallPeriod(pCallPeriod);

    data.m_pDeriv->SetCallSchedule(pCalls);

    double dPrice = data.Price();     

    CPPUNIT_ASSERT( dPrice >= dPriceOld );
    
    dPriceOld = dPrice;
    dYieldToCall += dYieldToCallStep;
  } //end while

  //test is passed so removed xml output
  //otherwise, cpp unit will fail and go to the next
  //test inside the testsuite, that is why we only need
  //to remove the file at the end
  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 );

}//BasicFeaturesTester::TestPriceIncreasesAsYTCIncreases()

void BasicFeaturesTester::TestPriceIncreasesAsYTPIncreases()
{

  double dYieldToPutMax   = .18;
  double dYieldToPutStart = 1.e-8;
  double dYieldToPutStep  = .01;
  double dYieldToPut      = dYieldToPutStart;
  double dPriceOld        = 0.0;

  Date maturityDate  = m_pCBInit->GetMaturityDate();
  Date issueDate     = m_pCBInit->GetBondTerms()->GetIssueDate();
  Date putDate       = issueDate;
  long lDays = Date::DaysDiff(issueDate, maturityDate);
  putDate.AddDays(lDays/2);

  //create the cb
  CBData data;
  RestoreInitData(data);
  
  std::string xmlOutputFile = 
    data.GetXMLOutputFile("TestPriceIncreasesAsYTPIncreases");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);

  while ( dYieldToPut < dYieldToPutMax + dYieldToPutStep )
  {  
    shared_ptr<finance::PutSchedule> pPutSchedule ( new PutSchedule() );  
    pPutSchedule->AddPutWithYield(putDate,dYieldToPut);
  
    data.m_pDeriv->SetPutSchedule(pPutSchedule);

    double dPrice = data.Price();

    CPPUNIT_ASSERT( dPrice >= dPriceOld );

    dPriceOld    = dPrice;
    dYieldToPut += dYieldToPutStep;
  } //end while

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 ) ;

}//BasicFeaturesTester::TestPriceIncreasesAsYTPIncreases()

void BasicFeaturesTester::TestPriceIncreasesAsCallNoticePeriodIncreases()
{
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  xml::node::const_iterator
    pNode = nodeRoot.find(XML_IHGTEST_CB_CALLNOTICE_PERIOD);

  shared_ptr<finance::CallSchedule> pCallSchedule = data.m_pDeriv->GetCallSchedule();

  std::string xmlOutputFile = 
  data.GetXMLOutputFile("TestPriceIncreasesAsCallNoticePeriodIncreases");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);

  double dPriceOld = data.Price();

  size_t nIdx; //start at 1, call notice period of 0 does not make sense
  for ( nIdx = 1 ; nIdx < 20 ; nIdx += 2 )
  {
    pCallSchedule->SetNoticePeriod( nIdx );

    data.m_pDeriv->SetCallSchedule(pCallSchedule);

    double dPrice = data.Price();

    CPPUNIT_ASSERT( dPrice > dPriceOld );

    dPriceOld = dPrice;
  }

  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0);

}//BasicFeaturesTester::TestPriceIncreasesAsYTPIncreases()
