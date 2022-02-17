#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/cppunit.h"
#include "ito33/dateutils.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivativevisitors/bondlikevisitor.h"

#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"

#include "ito33/xml/read.h"
#include "ito33/xml/write.h"

#include "ito33/xml/finance/bondlike/callschedule.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/hazardrateflat.h"

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
#include "bond_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::XML;

void BondTester::Setup(std::string strInputFilename)
{
  m_strInputFile = strInputFilename;

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());

  m_pSessionDataInit = reader.ReadSessionData();

  finance::BondLikeVisitor visitor;
  reader.ReadDerivatives(visitor);
  m_pBondInit = visitor.GetBond();

  if ( !m_pBondInit )
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "no bond in input xml file");

  m_pModelInit = make_ptr( new ihg::TheoreticalModel() );
  reader.ReadTheoreticalModel(m_pModelInit);
}

void BondTester::RestoreInitData(BondData &data)
{ 
  data.m_pSessionData = m_pSessionDataInit;

  shared_ptr<BondTerms> 
    bondTerms ( new BondTerms(*m_pBondInit->GetBondTerms()) );  
  
  data.m_pDeriv = make_ptr( new finance::Bond(bondTerms) );
  
  if ( m_pBondInit->GetCallSchedule() )
  {
    shared_ptr<CallSchedule> 
      calls ( new CallSchedule(*m_pBondInit->GetCallSchedule()) );
    
    data.m_pDeriv->SetCallSchedule(calls);
  }
  
  if ( m_pBondInit->GetPutSchedule() )
  {
    shared_ptr<PutSchedule> 
      puts ( new PutSchedule(*m_pBondInit->GetPutSchedule()) );
    
    data.m_pDeriv->SetPutSchedule(puts);
  }

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  if ( m_pBondInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pBondInit->GetNumeraire() );

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}

BondData BondTester::GetBasicData()
{ 
  BondData data;

  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();
  xml::node::const_iterator pNode = nodeRoot.find(XML_IHGTEST_BASICDDATA);

  data.m_strTestName = std::string(XML_IHGTEST_BASICDDATA);

  return data;
}

void BondTester::TestPriceIncreasesAsYTCIncreases()
{
 
  double dYieldToCallMax   = .18;
  double dYieldToCallStart = 1.e-8;
  double dYieldToCallStep  = .01;
  double dYieldToCall      = dYieldToCallStart;
  double dPriceOld         = 0.0;

  Date maturityDate  = m_pBondInit->GetMaturityDate();
  Date issueDate     = m_pBondInit->GetBondTerms()->GetIssueDate();

  BondData data;

  data.m_pSessionData = m_pSessionDataInit;
  data.m_pModel       = m_pModelInit;
  data.m_strFileName  = m_strInputFile;

  shared_ptr<BondTerms> 
    bondTerms ( new BondTerms(*m_pBondInit->GetBondTerms()) );

  data.m_pDeriv = make_ptr( new finance::Bond(bondTerms) );
 
  data.m_pDeriv->SetSessionData(data.m_pSessionData);
  
  if ( m_pBondInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pBondInit->GetNumeraire() );

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

}//BondTester::TestPriceIncreasesAsYTCIncreases()

void BondTester::TestPriceIncreasesAsYTPIncreases()
{

  double dYieldToPutMax   = .18;
  double dYieldToPutStart = 1.e-8;
  double dYieldToPutStep  = .01;
  double dYieldToPut      = dYieldToPutStart;
  double dPriceOld        = 0.0;

  Date maturityDate  = m_pBondInit->GetMaturityDate();
  Date issueDate     = m_pBondInit->GetBondTerms()->GetIssueDate();
  Date putDate       = issueDate;
  long lDays = Date::DaysDiff(issueDate, maturityDate);
  putDate.AddDays(lDays/2);

  //create the bond
  BondData data;
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

}//BondTester::TestPriceIncreasesAsYTPIncreases()

void BondTester::TestPriceWithCouponsDefinedByFDF()
{
  //create the bond
  BondData data;
  RestoreInitData(data);
  
  // No put
  shared_ptr<finance::PutSchedule> pPutSchedule ( new PutSchedule() ); 
  data.m_pDeriv->SetPutSchedule(pPutSchedule);
  
  // No call
  shared_ptr<finance::CallSchedule> pCallSchedule ( new CallSchedule() ); 
  data.m_pDeriv->SetCallSchedule(pCallSchedule);

  // No default
  shared_ptr<ihg::HazardRate> pHRNull (new ihg::HazardRateFlat(0.));  
  data.m_pModel->SetHazardRate( pHRNull );
  
  std::string xmlOutputFile = 
    data.GetXMLOutputFile("TestPriceWithCouponsDefinedByFDF");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);
  
  shared_ptr<CashFlowStream> 
    pCFSInit = data.m_pDeriv->GetBondTerms()->GetCashDistribution();
  
  // If the bond has no non floating coupons, test will stop here.
  if (! pCFSInit )
    return;
  
  shared_ptr<YieldCurve> pYC = m_pBondInit->GetSessionData()->GetYieldCurve();
  
  std::vector<CashFlow> pCashFlows = pCFSInit->GetAll();
  
  size_t
    nIdx,
    nNbCF = pCashFlows.size();

  Array<double> pdMaturities(nNbCF + 1);
  Array<double> pdDFF(nNbCF);
  
  pdMaturities[0] = GetDoubleFrom( pCFSInit->GetContractingDate() );
  for(nIdx = 0; nIdx < nNbCF; ++nIdx)
    pdMaturities[nIdx+1] = GetDoubleFrom( pCashFlows[nIdx].GetDate() );

  pYC->GetForwardDiscountFactor(pdMaturities.Get(), pdDFF.Get(), nNbCF + 1);

  std::vector<Date> paymentDates;
  std::vector<double> paymentAmounts;

  for( nIdx = 0; nIdx < nNbCF; ++nIdx )
  {
    paymentDates.push_back( pCashFlows[nIdx].GetDate() );
    paymentAmounts.push_back( 1 / pdDFF[nIdx] - 1 );
  }

  shared_ptr<CashFlowStream> 
    pCFS (new CashFlowStreamGeneral( pCFSInit->GetContractingDate(),
                paymentDates, paymentAmounts, 
                pCFSInit->GetDayCountConvention(),                         
                pCFSInit->GetPaymentFrequency() ));

  data.m_pDeriv->GetBondTerms()->SetCashDistribution( pCFS );

  double 
    dPrice,
    dExpectedPrice,
    dNominal = data.m_pDeriv->GetBondTerms()->GetNominal();

  for( nIdx = 0; nIdx < nNbCF - 1; ++nIdx )
  {
    data.m_pDeriv->GetSessionData()->SetValuationDate( paymentDates[nIdx] );
    
    dExpectedPrice = dNominal * (1. + paymentAmounts[nIdx]);

    dPrice = data.Price();
  
    CPPUNIT_ASSERT_DOUBLES_EQUAL( dExpectedPrice, dPrice, 1.e-6);
  }
  
  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 ) ;

}//BondTester::TestPriceWithCouponsDefinedByFDF()

void BondTester::TestFugitForBondWithNoConstraintsAndNoDefault()
{
  //create the bond
  BondData data;
  RestoreInitData(data);
  
  // No put
  shared_ptr<finance::PutSchedule> pPutSchedule ( new PutSchedule() ); 
  data.m_pDeriv->SetPutSchedule(pPutSchedule);
  
  // No call
  shared_ptr<finance::CallSchedule> pCallSchedule ( new CallSchedule() ); 
  data.m_pDeriv->SetCallSchedule(pCallSchedule);

  // No default
  shared_ptr<ihg::HazardRate> pHRNull (new ihg::HazardRateFlat(0.));  
  data.m_pModel->SetHazardRate( pHRNull );

  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);

  pFlags->SetComputeFugit(true);

  data.m_pDeriv->SetComputationalFlags(pFlags);
  
  std::string xmlOutputFile = 
    data.GetXMLOutputFile("TestFugitForBondWithNoConstraintsAndNoDefault");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);

  shared_ptr<finance::ModelOutput>
    pMO = data.m_pModel->Compute(*data.m_pDeriv.get());
  
  double
    dFugit,
    dExpectedFugit;

  dFugit = pMO->GetFugit();

  // dExpectedFugit = T - t
  dExpectedFugit = GetDoubleFrom(data.m_pDeriv->GetMaturityDate()) 
    - GetDoubleFrom( data.m_pDeriv->GetSessionData()->GetValuationDate() );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( dExpectedFugit, dFugit, 1.e-6);
  
  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 ) ;

}//BondTester::TestFugitForBondWithNoConstraintsAndNoDefault()

void BondTester::TestFugitForBondWithNoConstraintsAndNonNullFlatHR()
{
  //create the bond
  BondData data;
  RestoreInitData(data);
  
  // No put
  shared_ptr<finance::PutSchedule> pPutSchedule ( new PutSchedule() ); 
  data.m_pDeriv->SetPutSchedule(pPutSchedule);
  
  // No call
  shared_ptr<finance::CallSchedule> pCallSchedule ( new CallSchedule() ); 
  data.m_pDeriv->SetCallSchedule(pCallSchedule);

  shared_ptr<finance::ComputationalFlags> 
    pFlags(new finance::ComputationalFlags);

  pFlags->SetComputeFugit(true);

  data.m_pDeriv->SetComputationalFlags(pFlags);
  
  std::string xmlOutputFile = 
    data.GetXMLOutputFile("TestFugitForBondWithNoConstraintsAndNoDefault");

  data.m_pModel->SetDebugOutputFile(xmlOutputFile);
  
  double pdDefaultIntensities[3] = {0.02, 0.2, 0.9};
  
  double
    dTimeToMaturity,
    dFugit,
    dExpectedFugit;

  shared_ptr<ihg::HazardRate> pHR;

  shared_ptr<finance::ModelOutput> pMO;

  for( size_t nIdx = 0; nIdx < 3; ++nIdx)
  {    
    // Flat hazard rate
    pHR = shared_ptr<ihg::HazardRate> 
          (new ihg::HazardRateFlat(pdDefaultIntensities[nIdx]) );
    data.m_pModel->SetHazardRate( pHR );

    pMO = data.m_pModel->Compute(*data.m_pDeriv.get());
    
    dFugit = pMO->GetFugit();

    ASSERT_MSG( pdDefaultIntensities[nIdx] > 0., 
      "Default intensity must be strictly positif for this test.");

    dTimeToMaturity = GetDoubleFrom(data.m_pDeriv->GetMaturityDate()) 
      - GetDoubleFrom( data.m_pDeriv->GetSessionData()->GetValuationDate() );
    
    // dExpectedFugit = (1/p)[ 1 - exp( -p(T-t) ) ]  
    dExpectedFugit = ( 1. / pdDefaultIntensities[nIdx] ) * 
            ( 1. - exp( - pdDefaultIntensities[nIdx] * dTimeToMaturity) );

    CPPUNIT_ASSERT_DOUBLES_EQUAL( dExpectedFugit, dFugit, 1.e-2);
  }
  
  CPPUNIT_ASSERT( unlink( xmlOutputFile.c_str() ) == 0 ) ;

}//BondTester::TestFugitForBondWithNoConstraintsAndNonNullFlatHR()
