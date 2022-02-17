#include "ito33/beforestd.h"
#include <iostream>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/domain.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"

#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"
#include "ito33/finance/bondlike/callschedule.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"

#include "hg/model.h"

using namespace ito33;
using namespace ito33::finance;

ITO33_FORCE_LINK_MODULE(HGPriceReset);

shared_ptr<hg::TheoreticalModel> InitModel(size_t nNbRegimes);
shared_ptr<SessionData> InitSessionData(bool bHasDividend);
shared_ptr<Reset> InitReset();

int main()
{
  try
  {
    bool bFakeDividend = true;
    shared_ptr<SessionData> pSessionData = InitSessionData(bFakeDividend);
   
    shared_ptr<Reset> pReset = InitReset();
    pReset->SetSessionData(pSessionData);

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetAnalysisDate( pSessionData->GetValuationDate() );
    flags->SetComputeSurface(true);

    pReset->SetComputationalFlags(flags);

    shared_ptr<hg::TheoreticalModel> pModel = InitModel(1);

    // Actually price
    pModel->SetDebugOutputFile("hg_reset.xml");
    shared_ptr<ModelOutput> output = pModel->Compute(*pReset);

    // Report
    std::cout.precision(15);
    std::cout << "HG reset test" << std::endl;
    std::cout << "  Issue date = " << pReset->GetBondLikeTerms()->GetIssueDate() << std::endl;
    std::cout << "  Maturity date = " << pReset->GetMaturityDate() << std::endl;
    std::cout << "  Reset dates = ";
    ResetConversionSchedule::Elements resetDates =
      pReset->GetResetConversionSchedule()->GetAll();

    ResetConversionSchedule::Elements::const_iterator iterDates;
    iterDates = resetDates.begin();
    std::cout << (*iterDates)->GetDate();
    ++iterDates;
    for ( ; iterDates != resetDates.end(); ++iterDates)
      std::cout << ", " << (*iterDates)->GetDate();
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "  Price = " << output->GetPrice() << std::endl;
    std::cout << std::endl;

    if ( output->HasPriceSurface() )
    {
      std::cout << "Surface computed: ";

      finance::SharedSurface pSurface = output->GetPriceSurface();

      shared_ptr<finance::Domain> pDomain = pSurface->GetDomain();

      finance::Domain::Spots pdSpots(20);
      pdSpots[0] = 0.8 * pSessionData->GetSpotSharePrice();
      double dStep = 0.1 * (pSessionData->GetSpotSharePrice() - pdSpots[0]);
      for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
        pdSpots[nIdxS] = pdSpots[nIdxS-1] + dStep;

      pDomain->SetUnderlyingSharePrices(pdSpots);

      const finance::Domain::Dates& pdDates = pDomain->GetDates();
      size_t nNbDates = pdDates.size();
      finance::SurfaceDouble::Doubles pValues 
        = pSurface->GetValuesAt(nNbDates-1);
      finance::SurfaceDouble::Doubles pValues2 
        = pSurface->GetValuesAt(0);

      std::cout << nNbDates << " dates saved" << std::endl;
      std::cout << "First date = " << pdDates[0] << std::endl;
      std::cout << "Last date = " << pdDates[nNbDates-1] << std::endl;

      for (size_t nIdx = 0; nIdx < pValues.size(); nIdx++)
        std::cout << pdSpots[nIdx] << ", value = " << pValues[nIdx]
                  << ", " << pValues2[nIdx]
                  << std::endl;
      
      std::cout << std::endl;
    }

    if ( output->HasPriceAtAnalysisDate() )
    {

      std::cout << "Analysis date values computed:" << std::endl;
      std::vector<double> pdSpots = output->GetSpotsAtAnalysisDate();
      std::vector<double> pdPrices = output->GetPricesAtAnalysisDate();

      for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
      {
        std::cout << pdSpots[nIdx] << " "
                  << pdPrices[nIdx] << " "
                  << std::endl;
      }
      std::cout << std::endl;
    }

    //ShowConvergence(*pModel, *pReset);

    return 0;
    
  }

  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
      << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught.\n";

    return 2;
  }
 
}


shared_ptr<SessionData> InitSessionData(bool bHasDividend)
{
  Date valuationDate(2002, Date::Jan, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  shared_ptr<RateData> pRateData(new RateData);

  double dContinuousRate = 0.05;
  double dAnnualRate = exp(dContinuousRate) - 1.0;
  shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dAnnualRate));

  pRateData->SetYieldCurve(pCurrency, pyc);

  double dS0 = 100.0;
  shared_ptr<Equity> pEquity(new Equity(dS0, pCurrency));

  if ( bHasDividend )
  {
    shared_ptr<Dividends> pDividends(new Dividends());
    pDividends->AddCash(Date(2002, Date::Jul, 1), 2.0e-8);
    pEquity->SetDividends(pDividends);
  }

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


shared_ptr<hg::TheoreticalModel> InitModel(size_t nNbRegimes)
{
  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.05);

  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  hg::Jumps jumps;

  //jumps.push_back(hg::Jump(0.1, -0.2));
  //pUnderlyingProcess->SetJumps(0, 0, jumps);
    
  if (nNbRegimes > 1)
  {    
    jumps.clear();
    jumps.push_back(hg::Jump(0.15, -0.25));
    pUnderlyingProcess->SetJumps(0, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.2, -0.5));
    pUnderlyingProcess->SetJumps(1, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.11, 0.35));
    pUnderlyingProcess->SetJumps(1, 0, jumps);    
  }

  if (nNbRegimes > 2)
  {
    jumps.clear();
    jumps.push_back(hg::Jump(0.1, -0.2));
    pUnderlyingProcess->SetJumps(0, 2, jumps); 
 
    jumps.clear();
    jumps.push_back(hg::Jump(0.1, -0.3));
    pUnderlyingProcess->SetJumps(1, 2, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.3, -0.2));
    pUnderlyingProcess->SetJumps(2, 0, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.2, -0.4));
    pUnderlyingProcess->SetJumps(2, 1, jumps); 

    jumps.clear();
    jumps.push_back(hg::Jump(0.05, -0.25));
    pUnderlyingProcess->SetJumps(2, 2, jumps); 
  }

  shared_ptr<hg::TheoreticalModel> 
    pModel(new hg::TheoreticalModel(pUnderlyingProcess));

  return pModel;
}


shared_ptr<finance::Reset> InitReset()
{

  //--------------------------------------------------------------------------
  // Create the different windows of constraints
  //--------------------------------------------------------------------------
  Date issueDate     = Date(2003,Date::Jan, 1);
  Date maturityDate  = Date(2004,Date::Jan, 1);

  Date startCallDate = Date(2003,Date::Feb, 1);
  Date endCallDate   = Date(2004,Date::May, 1);

  Date startConvDate = Date(2003,Date::Feb, 1);
  Date endConvDate   = Date(2004,Date::May, 1);


  //--------------------------------------------------------------------------
  // Reset parameters
  //--------------------------------------------------------------------------
  double dCapRate      = 1.;
  double dFloorRate    = 0.6;
  double dMultiplier   = 1.0;
  double dParValue     = 110;
  double dInitialConvPrice = dParValue;
  double dCurrentConvPrice = dParValue;

  finance::ResetFlooredBy resetFlooredBy 
    = finance::ResetFlooredBy_PrevailingConversionPrice;

  // Setup the reset data
  std::list<Date> resetDates;
  shared_ptr<ResetConversionSchedule> pResetConversionSchedule( new 
    ResetConversionSchedule(issueDate, maturityDate, dInitialConvPrice, 
                dCurrentConvPrice, resetFlooredBy) );

  Date resetDate1(2003, Date::Mar, 1);
  resetDates.push_back(resetDate1);
  shared_ptr<ConversionPriceReset> pConversionPriceReset 
    ( new ConversionPriceReset( resetDate1,  dFloorRate) );

  pConversionPriceReset->SetMultiplier(dMultiplier);

  pConversionPriceReset->SetCap(dCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  Date resetDate2(2003, Date::Aug, 1);
  resetDates.push_back(resetDate2);
  pConversionPriceReset = make_ptr( new ConversionPriceReset(resetDate2, dFloorRate) );
  
  pConversionPriceReset->SetMultiplier(dMultiplier);

  pConversionPriceReset->SetCap(dCapRate);

  pResetConversionSchedule->AddConversionPriceReset(pConversionPriceReset);

  shared_ptr<finance::CallSchedule> pCall( new CallSchedule() );
  
  shared_ptr<finance::CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                              (issueDate,maturityDate,1.7) );

  pCallPeriod->SetTrigger(1.3);

  pCall->AddCallPeriod(  pCallPeriod );

  
  double dIssuePrice = 1;
  double dRedemptionRate = 1;
  double dRecoveryRate = 0.;  

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                       maturityDate, dParValue, dRedemptionRate,
                       dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream>
   pInterests( new CashFlowStreamGeneral
                   (issueDate, paymentDates, paymentRates,
                   Date::DayCountConvention_Act365,finance::Frequency_Annual)
             );

  bc->SetCashDistribution(pInterests);

  // finally create and return
  shared_ptr<finance::Reset> 
    pReset( new finance::Reset(bc, pResetConversionSchedule) );

  pReset->SetCallSchedule(pCall);

  return pReset;
}

