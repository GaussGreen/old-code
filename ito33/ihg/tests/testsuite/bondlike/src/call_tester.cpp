#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"


#include "ito33/sharedptr.h"
#include "ito33/cppunit.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/xml/finance/bondlike/conversionschedule.h"

#include "ihg/tests/testdata.h"

// local files
#include "bondlike_reader.h"
#include "call_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;


void CallTester::Setup(std::string strInputFilename)
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


void CallTester::RestoreInitData(CBData &data)
{ 
  // Need to read the input file again.
  Setup(m_strInputFile);

  // Does not work unless the input file is re-read
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

CBData CallTester::GetBasicData()
{ 
  CBData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();
  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_BASICDDATA);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}


void CallTester::PriceDecreaseAsTriggerHistoryIncreases()
{
  
 CBData data;

  RestoreInitData(data);

  shared_ptr<CallSchedule> 
    calls ( new CallSchedule(*m_pCBInit->GetCallSchedule()) );

  /*
    Generate new call schedule where the trigger
    has been lowered othwerwise price does not change.
  */
  shared_ptr<CallSchedule> pCallSchedule ( new CallSchedule() );

  CallSchedule::Elements pList = calls->GetAll();

  CallSchedule::Elements::iterator iter;

  double dStrike = 1.0;
  double dTrigger = .5;

  for ( iter = pList.begin() ;  iter != pList.end(); iter++)
  {
    shared_ptr<CallPeriod> pCallPeriod = *iter;

    Date startDate = pCallPeriod->GetStartDate();
    Date endDate   = pCallPeriod->GetEndDate();

    shared_ptr<CallPeriod> 
      pCallPeriodNew
      (
        CallPeriod::CreateWithStrike(startDate, endDate, dStrike)
      );

    pCallPeriodNew->SetTrigger( dTrigger );

    pCallSchedule->AddCallPeriod(pCallPeriodNew);
  }

  
  size_t nIdx;
  double dPriceOld = 1.e8;
  
  for ( nIdx = 2; nIdx <= 30; nIdx += 2)
  {
    pCallSchedule->SetTriggerCheckPeriod(30, nIdx);
    data.m_pDeriv->SetCallSchedule(pCallSchedule);

    double dPrice = data.Price();

    CPPUNIT_ASSERT( dPriceOld >= dPrice );
    dPriceOld = dPrice;

  } //end loop
}


void CallTester::PriceIncreaseAsCallPeriodIncreases()
{

  CBData data;

  RestoreInitData(data);

  shared_ptr<CallSchedule> 
    calls ( new CallSchedule(*m_pCBInit->GetCallSchedule()) );


  size_t nIdx;
  double dPriceOld = 0.;
  size_t nTriggerHistory = 0;
  
  for ( nIdx = 5; nIdx <= 30; nIdx += 5)
  {
    calls->SetTriggerCheckPeriod(nIdx, nTriggerHistory);
    data.m_pDeriv->SetCallSchedule(calls);

    double dPrice = data.Price();

    CPPUNIT_ASSERT( dPriceOld < dPrice );
    
    dPriceOld = dPrice;

  } //end loop

}