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
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/xml/finance/bondlike/conversionschedule.h"

#include "ihg/tests/testdata.h"

// local files
#include "bondlike_reader.h"
#include "conversion_tester.h"

extern const ito33::Error ITO33_UNEXPECTED;

using namespace ito33;
using namespace ito33::finance;


void ConversionTester::Setup(std::string strInputFilename)
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


void ConversionTester::RestoreInitData(CBData &data)
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
    conversions
    ( 
      new ConversionSchedule(*m_pCBInit->GetConversionSchedule())
    );
  
  data.m_pDeriv = make_ptr( new finance::ConvertibleBond
                                (bondTerms, conversions) );
  data.m_pDeriv->SetCallSchedule(calls);

  data.m_pDeriv->SetSessionData(data.m_pSessionData);

  if ( m_pCBInit->GetNumeraire() )
    data.m_pDeriv->SetNumeraire( m_pCBInit->GetNumeraire() );

  data.m_pModel = m_pModelInit;

  data.m_strFileName = m_strInputFile;
}


CBData ConversionTester::ReadTestFromXML(const char* sTestName)
{
  // Restore the original data
  CBData data;
  RestoreInitData(data);

  ihg::XML::BondLikeReader reader(m_strInputFile.c_str());
  const xml::node& nodeRoot = reader.GetRootNode();

  // Find the appropriate section in the xml file
  xml::node::const_iterator pNode = nodeRoot.find(sTestName);

  // Get and set the new conversion period
  shared_ptr<finance::ConversionSchedule> pConversionSchedule;  
  XML::Restore(*pNode, pConversionSchedule);
  data.m_pDeriv->SetConversionSchedule(pConversionSchedule);

  // Set the test name
  data.m_strTestName = std::string(sTestName);

  return data;

}

void ConversionTester::TestForcedConversion()
{

  // Use the default conversion schedule: no coco, ratio of 1, full call
  CBData data = ReadTestFromXML(XML_IHGTEST_CB_DEFAULT_CONVERSION);

  // Force conversion by setting large spot, and ensuring a full
  // call period
  double dNominal = data.m_pDeriv->GetBondTerms()->GetNominal();
  double dSpot = dNominal * 2.0;
  data.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot);

  Date issueDate = data.m_pDeriv->GetBondTerms()->GetIssueDate();
  Date maturityDate = data.m_pDeriv->GetMaturityDate();
  shared_ptr<CallSchedule> pCallSchedule ( new CallSchedule() );
  pCallSchedule->SetKeepAccrued(false);
  pCallSchedule->SetForfeitCoupon(false);
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                             (issueDate, maturityDate, 1.0) );
  pCallSchedule->AddCallPeriod(pCallPeriod);

  data.m_pDeriv->SetCallSchedule(pCallSchedule);

  // The price should be conversion ratio times the spot
  double dPrice = data.Price();
  double dExpectedPrice = dSpot * 1.0;

  ITO33_ASSERT_DOUBLES_EQUAL(dPrice, dExpectedPrice);

}


void ConversionTester::TestTriggerContinuity()
{

  // Use the quarterly conversion schedule. Should be
  // discontinuous at a discrete observation date, and
  // continuous at all other dates
  CBData data = ReadTestFromXML(XML_IHGTEST_CB_ADD_QUARTERLY_FOREVER);

  // Remove the call window, if any
  shared_ptr<CallSchedule> pCallSchedule ( new CallSchedule() );
  data.m_pDeriv->SetCallSchedule(pCallSchedule);

  // Move to the last quarterly observation date. Continuity effects
  // will be largest at this time.  If the issue date is used, the solution
  // will be smooth due to the length of the contract
  shared_ptr<ConversionSchedule> pConversionSchedule =
    data.m_pDeriv->GetConversionSchedule();

  ConversionSchedule::Elements periods = pConversionSchedule->GetAll();
  ConversionSchedule::Elements::iterator iterPeriods;

  CPPUNIT_ASSERT( periods.size() == 1);

  iterPeriods = periods.begin();
  Date observationDate = (*iterPeriods)->GetStartDate();  
  Date maturityDate = data.m_pDeriv->GetMaturityDate();
  while (observationDate < maturityDate)
    observationDate.AddMonths(3);

  observationDate.AddMonths(-3);

  // Set the valuation date to be 1 day before an observation, at the
  // observation, and 1 day after the observation.
  // Get prices at the trigger, and to either side. 
  double dNominal = data.m_pDeriv->GetBondTerms()->GetNominal();
  double dTriggerRate = (*iterPeriods)->GetTrigger();
  double dTriggerLevel = dNominal * dTriggerRate;

  std::vector<double> pdPrices(9);

  size_t nIdxPrices = 0;
  for (int iDateShift = -1; iDateShift <= 1 ; iDateShift++)
  {
    // Set the new valuation date
    Date valuationDate = observationDate;
    valuationDate.AddDays(iDateShift);
    data.m_pSessionData->SetValuationDate(valuationDate);

    for (double dSpotShift = -1.0; dSpotShift < 1.5; dSpotShift += 1.0)
    {
      // Shift the spot and save the price
      double dSpot = dTriggerLevel + dSpotShift;
      data.m_pSessionData->GetEquity()->SetSpotSharePrice(dSpot);

      double dPrice = data.Price();
      pdPrices[nIdxPrices] = dPrice;
      nIdxPrices++;

    } // shifting the spot

  } // shifting valuation date

  /*
  for (size_t nIdx = 0; nIdx < 9; nIdx++)
    std::cout << "price " << nIdx << " = " << pdPrices[nIdx] << std::endl;
  */

  // Test continuity using left and right-sided differences
  // (equivalent to finite differences)
 
  // Just before the observation date, the solution should have
  // been smoothed by one timestep
  double dLeftDifference = pdPrices[1] - pdPrices[0];
  double dRightDifference = pdPrices[2] - pdPrices[1];

  // If this generates divide by zero, then there is another problem...
  double dRelDiff = fabs(dLeftDifference - dRightDifference)/dLeftDifference; 
  CPPUNIT_ASSERT( dRelDiff < 0.3 );

  // At the observation date, should have a discontinuity
  dLeftDifference = pdPrices[4] - pdPrices[3];
  dRightDifference = pdPrices[5] - pdPrices[4];

  CPPUNIT_ASSERT(dLeftDifference/dRightDifference > 5.0);

  // After the observation, nothing has happened.  No calls, no puts,
  // no conversions, so we just have a bond. The prices should
  // be the same (if we went past a previous observation date,
  // the solution should be smooth)
  dLeftDifference = pdPrices[7] - pdPrices[6];
  dRightDifference = pdPrices[8] - pdPrices[7];

  ITO33_ASSERT_DOUBLES_EQUAL(dLeftDifference, 0.0);
  ITO33_ASSERT_DOUBLES_EQUAL(dRightDifference, 0.0);

}

void ConversionTester::TestComparisons()
{

  // Get prices for the base cases, as well as the different coco clauses

  CBData data = ReadTestFromXML(XML_IHGTEST_CB_DEFAULT_CONVERSION);
  double dPriceNoTrigger = data.Price();

  data = ReadTestFromXML(XML_IHGTEST_CB_ADD_QUARTERLY_FOREVER);
  double dPriceQuarterlyForever = data.Price();

  data = ReadTestFromXML(XML_IHGTEST_CB_ADD_QUARTERLY_NEXT_QUARTER);
  double dPriceQuarterlyNextQuarter = data.Price();

  data = ReadTestFromXML(XML_IHGTEST_CB_ADD_ANYTIME_FOREVER);
  double dPriceAnytimeForever = data.Price();

  data = ReadTestFromXML(XML_IHGTEST_CB_ADD_ANYTIME_CHECKDATE);
  double dPriceAnytimeCheckDate = data.Price();

  // Compare the prices. Use a tolerance to avoid numerical errors
  // Must have:
  //   NoTrigger >= AnytimeForever >= QuarterlyForever >= QuarterlyNextQuarter
  //             >= NoConversion
  //   AnytimeForever >= AnytimeCheckDate  
  // This assumes that all other conversion properties (trigger, lengths, etc)
  // are the same
  double dTol = dPriceNoTrigger * 1.e-6;

  CPPUNIT_ASSERT( dPriceNoTrigger > dPriceAnytimeForever - dTol );

  CPPUNIT_ASSERT( dPriceAnytimeForever > dPriceQuarterlyForever - dTol );

  CPPUNIT_ASSERT( dPriceQuarterlyForever > dPriceQuarterlyNextQuarter - dTol );

  CPPUNIT_ASSERT( dPriceAnytimeForever > dPriceAnytimeCheckDate - dTol );

  /*
  std::cout << "no trigger = " << dPriceNoTrigger << std::endl;
  std::cout << "anytime forever = " << dPriceAnytimeForever << std::endl;
  std::cout << "quarterly forever = " << dPriceQuarterlyForever << std::endl;
  std::cout << "quarterly next quarter = " << dPriceQuarterlyNextQuarter << std::endl;
  std::cout << "anytime checkdate = " << dPriceAnytimeCheckDate << std::endl;
  */
}

void ConversionTester::TestKeepAccrued()
{
  // Price with keep accrued should be greater than without

  CBData data = ReadTestFromXML(XML_IHGTEST_CB_DEFAULT_CONVERSION);

  data.m_pDeriv->GetConversionSchedule()->SetKeepAccrued(true);
  double dPriceKAtrue = data.Price();

  data.m_pDeriv->GetConversionSchedule()->SetKeepAccrued(false);
  double dPriceKAfalse = data.Price();


  double dTol = 1.e-6;
  CPPUNIT_ASSERT( dPriceKAtrue > dPriceKAfalse - dTol );

  /*
  std::cout << "Price with keep accrued = " << dPriceKAtrue << std::endl; 
  std::cout << "Price without keep accrued = " << dPriceKAfalse << std::endl; 
  */
}

void ConversionTester::TestAnytimeForever()
{
  // Add anytime forever periods of different lengths, and different rates
  // Compare with cb with and without conversion

  // Full conversion, no triggers.  Should be max price
  CBData data = ReadTestFromXML(XML_IHGTEST_CB_DEFAULT_CONVERSION);
  double dPriceNoTrigger = data.Price();

  // The CoCo anytime forever should vary between the above prices
  // Restore the original data
  RestoreInitData(data);

  Date periodStart = data.m_pDeriv->GetBondTerms()->GetIssueDate();
  Date periodEnd = data.m_pDeriv->GetMaturityDate();

  // No trigger, full period
  shared_ptr<ConversionSchedule> 
    pConversionSchedule ( new ConversionSchedule() );

  shared_ptr<ConversionPeriod> 
    pConversionPeriod ( new ConversionPeriod(periodStart, periodEnd, 1.0) );
  pConversionPeriod->SetCoCo(1.e-6, 
    CoCoType_CheckAnyTimeAndConvertAsOfCheckDate, 0.0, 1.e-6); 
  pConversionSchedule->AddConversionPeriod(pConversionPeriod);
  data.m_pDeriv->SetConversionSchedule(pConversionSchedule);

  double dPriceCoCoNoTrigger = data.Price();

  // Trigger rate of 1.0, full period
  pConversionSchedule = make_ptr( new ConversionSchedule() );
  pConversionPeriod = make_ptr( new ConversionPeriod(periodStart, periodEnd, 1.0) );
  pConversionPeriod->SetCoCo(1.0, 
    CoCoType_CheckAnyTimeAndConvertAsOfCheckDate, 0.0, 1.0); 
  pConversionSchedule->AddConversionPeriod(pConversionPeriod);
  data.m_pDeriv->SetConversionSchedule(pConversionSchedule);

  double dPriceCoCoWithTrigger = data.Price();

  // Start the ordering checks
  double dTol = 1.e-6;
  ITO33_ASSERT_DOUBLES_EQUAL(dPriceNoTrigger, dPriceCoCoNoTrigger);
  CPPUNIT_ASSERT(dPriceCoCoNoTrigger > dPriceCoCoWithTrigger - dTol); 

  // Save old price for later ordering test
  double dPreviousPrice = dPriceCoCoWithTrigger;

  // Reduce the conversion period
  periodStart.AddYears(1);
  periodEnd.AddYears(-1);
  while (periodStart < periodEnd)
  {
    // Trigger rate of 1.0, reduced period
    pConversionSchedule = make_ptr( new ConversionSchedule() );
    pConversionPeriod = make_ptr( new ConversionPeriod
                                      (periodStart, periodEnd, 1.0) );
    pConversionPeriod->SetCoCo(1.0, 
      CoCoType_CheckAnyTimeAndConvertAsOfCheckDate, 0.0, 1.0); 
    pConversionSchedule->AddConversionPeriod(pConversionPeriod);
    data.m_pDeriv->SetConversionSchedule(pConversionSchedule);

    double dPriceCoCoReduced = data.Price();

    periodStart.AddYears(1);
    periodEnd.AddYears(-1);

    CPPUNIT_ASSERT(dPreviousPrice > dPriceCoCoReduced - dTol);
    dPreviousPrice = dPriceCoCoReduced;

    //std::cout << "anytime forever, with trigger: " << dPriceCoCoReduced << std::endl;
  }

  /*
  std::cout << "normal, no triggers: " << dPriceNoTrigger << std::endl;
  std::cout << "normal, no conversion: " << dPriceNoConversion << std::endl;
  std::cout << "anytime forever, no trigger: " << dPriceCoCoNoTrigger << std::endl;
  std::cout << "anytime forever, with trigger: " << dPriceCoCoWithTrigger << std::endl;
  */

}
