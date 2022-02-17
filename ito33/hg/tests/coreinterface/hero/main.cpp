#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/timer.h"

#include "ito33/finance/domain.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/option.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/hedgeratiodata.h"
#include "ito33/hg/heroflags.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_MODULE(HGPriceOption);
ITO33_FORCE_LINK_MODULE(HGPriceCDS);
ITO33_FORCE_LINK_MODULE(HGPriceEDS);
ITO33_FORCE_LINK_MODULE(HGPriceCB);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

shared_ptr<SessionData> InitSessionData()
{
  Date valuationDate(2003, Date::Feb, 1);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  shared_ptr<Equity> pEquity(new Equity(60.0, pCurrency));
  //shared_ptr<Equity> pEquity(new Equity(40.0, pCurrency));

  shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.01));
  pEquity->SetBorrowCurve(pyf);

  shared_ptr<Dividends> pDividends( new Dividends );
  pDividends->AddCash(Date(2003, Date::Jul, 1), 1.0);
  pEquity->SetDividends(pDividends);
    
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.05) );

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pyc);
    
  shared_ptr<SessionData> pSessionData( 
    new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}


shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date issueDate = pSessionData->GetValuationDate();
  Date maturityDate = issueDate;
  maturityDate.AddYears(1);

  Date conversionStart = issueDate;
  conversionStart.AddMonths(3);
  Date conversionEnd = maturityDate;

  double
    dIssuePrice = 1,
    dParValue = pSessionData->GetSpotSharePrice() * 1.1,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                      maturityDate, dParValue, dRedemptionRate,
                      dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date payDate = issueDate;
  payDate.AddMonths(4);
  paymentDates.push_back(payDate);
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral(issueDate, paymentDates, paymentRates,
    Date::DayCountConvention_Act365,Frequency_SemiAnnual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  shared_ptr<ConversionPeriod> 
    pConvPeriod( new ConversionPeriod(conversionStart, conversionEnd, 1) );
  conv->AddConversionPeriod(pConvPeriod);

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  
//  size_t nNoticePeriod = 20;
//    pCallSchedule->SetNoticePeriod(nNoticePeriod);
//    startCallDate.AddMonths(4);
//    endCallDate.AddMonths(-4);
        
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                             (startCallDate, endCallDate, 1.) );

  pCallSchedule->AddCallPeriod( pCallPeriod);

        
  pConvertibleBond->SetCallSchedule(pCallSchedule);// 

  pConvertibleBond->SetSessionData(pSessionData);
  
  return pConvertibleBond;
}


void InitHedgeContracts(Derivatives& hedgeContracts, 
                        shared_ptr<SessionData> pSessionData)
{

  double dSpot = pSessionData->GetSpotSharePrice();

  Date valuationDate = pSessionData->GetValuationDate();
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);

  Date maturityDate2 = maturityDate;
  maturityDate2.AddMonths(-4);

  shared_ptr<Option> opt;

  opt = make_ptr( new Option(dSpot - 4.0, maturityDate2, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  hedgeContracts.Add(opt);


  opt = make_ptr( new Option(dSpot + 8.0, maturityDate, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  hedgeContracts.Add(opt);


  opt = make_ptr( new Option(dSpot + 8.0, maturityDate, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  hedgeContracts.Add(opt);


  opt = make_ptr( new Option(dSpot + 2.0, maturityDate2, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  //hedgeContracts.Add(opt);


  opt = make_ptr( new Option(dSpot - 4.0, maturityDate, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  //hedgeContracts.Add(opt);

  opt = make_ptr( new Option(dSpot + 16.0, maturityDate, 
                             Option_Call, ExerciseType_European) );
  opt->SetSessionData(pSessionData);

  //hedgeContracts.Add(opt);


  Date issueDate = valuationDate;
  Date firstDate = valuationDate;
  firstDate.AddMonths(3);
  Date lastDate = maturityDate;

  double dSpread = 0.25;
   
  shared_ptr<CashFlowStreamUniform>
     pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_Quarterly
                        )
                  );

  double dRecoveryRate = 0.5;
  shared_ptr<finance::CDS>
     pCDS( new ito33::finance::CDS(dRecoveryRate, pSpreadStream) );

  pCDS->SetSessionData(pSessionData);

  //hedgeContracts.Add(pCDS);
  
}


shared_ptr<Derivative> InitTargetContract(shared_ptr<SessionData> pSessionData)
{
/*
  double dSpot = pSessionData->GetSpotSharePrice();
  Date valuationDate = pSessionData->GetValuationDate();
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);
  shared_ptr<Option> opt(new Option(dSpot, maturityDate, Option_Put,
                                   ExerciseType_European));

  opt->SetSessionData(pSessionData);

  return opt;
*/
/*
  double dSpot = pSessionData->GetSpotSharePrice();
  Date valuationDate = pSessionData->GetValuationDate();
  Date maturityDate = valuationDate;
  maturityDate.AddYears(1);

  Date issueDate = valuationDate;
  Date firstDate = valuationDate;
  firstDate.AddMonths(3);
  Date lastDate = maturityDate;

  double dSpread = 0.25;
   
  shared_ptr<CashFlowStreamUniform>
     pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_Quarterly
                        )
                  );

  double dRecoveryRate = 0.2;
  double dBarrier = dSpot * 0.3;
  shared_ptr<finance::EDS>
     pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);

  return pEDS;
*/

  shared_ptr<finance::ConvertibleBond> pBond = InitCB(pSessionData);
  return pBond;
}

shared_ptr<hg::TheoreticalModel> InitModel()
{
  size_t nNbRegimes = 3;

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.1);
  //pdDefaultIntensities.resize(nNbRegimes, 0.0);

  if ( nNbRegimes > 1 )
  {
    pdVols[1] = 0.35;
    pdDefaultIntensities[1] = 0.15;
  }

  if ( nNbRegimes > 2 )
  {
    pdVols[2] = 0.45;
    pdDefaultIntensities[2] = 0.25;
  }



  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess( new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

  hg::Jumps jumps;

  jumps.push_back(hg::Jump(0.1, -0.2));
  pUnderlyingProcess->SetJumps(0, 0, jumps);
  
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

  pModel->SetSharpeRatio(0.5);

  return pModel;
}



void ShowHeroConvergence(shared_ptr<finance::Derivative> pTarget,
                         Derivatives& hedgeDerivatives,
                         shared_ptr<hg::TheoreticalModel> pModel,
                         size_t nNbTests)
{
  std::cout << "Convergence testing" << std::endl;
  std::cout << std::endl;

  bool bQuiet = false;

  size_t nNbTimes = numeric::NumParamsReference::GetNbTimeStepsFor5Years();
  size_t nNumberSpots = numeric::NumParamsReference::GetMinNbSpaceSteps();
  double dX = numeric::NumParamsReference::GetMaxDeltaLogS() * nNumberSpots;

  std::vector<double> pdDiffs;
  std::vector<double> pdRatings;
  std::vector<double> pdPrices;

  pdPrices.resize(nNbTests);
  pdDiffs.resize(nNbTests);
  pdRatings.resize(nNbTests);

  size_t n;
  for (n = 0; n < nNbTests; n++, nNbTimes *=2 , nNumberSpots*=2 )
  {
    numeric::NumParamsReference::SetNbTimeStepsFor5Years(nNbTimes);
    numeric::NumParamsReference::SetMinNbSpaceSteps(nNumberSpots);
    numeric::NumParamsReference::SetMaxDeltaLogS(dX / nNumberSpots);

    shared_ptr<HedgeOutput> 
      pHedgeOutput( pModel->ComputeHERO(*pTarget, hedgeDerivatives) );

    double dHero = pHedgeOutput->GetHERO();
    
    pdPrices[n] = dHero;

    if(n > 0)
      pdDiffs[n] = pdPrices[n] - pdPrices[n-1];
    if(n > 1)
      pdRatings[n] = pdDiffs[n - 1] / pdDiffs[n];
    
    if(!bQuiet)
    {
      if (n == 0)
      std::cout << pdPrices[n] << std::endl;
      else if (n == 1)
        std::cout << pdPrices[n] << "  " << pdDiffs[n] << std::endl;
      else
        std::cout << pdPrices[n] << "  " << pdDiffs[n] << "  "
                  << pdRatings[n] << std::endl;
    }
  } // refinement loop

}



int main()
{
  try
  {
    StopWatch sw;

    bool bComputeSurface = true;

    shared_ptr<SessionData> pSessionData = InitSessionData();

    Derivatives hedgeDerivatives;
    InitHedgeContracts(hedgeDerivatives, pSessionData);

    shared_ptr<Derivative> pTarget = InitTargetContract(pSessionData);

    shared_ptr<hg::TheoreticalModel> pModel = InitModel();

    //pModel->SetDebugOutputFile("hero_test.xml");

    shared_ptr<hg::HEROFlags> pHEROFlags(new hg::HEROFlags);
    pHEROFlags->SetComputeSurface(bComputeSurface);

    shared_ptr<HedgeOutput> 
      pHedgeOutput( pModel->ComputeHERO(*pTarget, hedgeDerivatives, pHEROFlags) );

    std::cout.precision(10);

    // Output HERO
    if ( pHedgeOutput->HasHERO() )
      std::cout << "HERO = " << pHedgeOutput->GetHERO()
                << std::endl 
                << std::endl;

    // Output hedge ratios
    double dUnderlyingRatio = pHedgeOutput->GetUnderlyingHedgeRatio();
    std::cout << "Underlying ratio = " << dUnderlyingRatio << std::endl;

    std::vector< shared_ptr<hg::HedgeRatioData> > ppHedgeRatioData 
      = pHedgeOutput->GetHedgeRatioData();

    for (size_t nIdx = 0; nIdx < ppHedgeRatioData.size(); nIdx++)
    {
      double dRatio = ppHedgeRatioData[nIdx]->GetRatio();
      std::cout << " Ratio " << nIdx << " = " << dRatio << std::endl;
    }
    std::cout << std::endl;

    // Output prices
    shared_ptr<finance::ModelOutput> pModelOutput;
    pModelOutput = pHedgeOutput->GetTargetModelOutput();
    std::cout << "Target price = " << pModelOutput->GetPrice() << std::endl;
    std::cout << "Target delta = " << pModelOutput->GetDelta() << std::endl;

    for (size_t nIdx = 0; nIdx < ppHedgeRatioData.size(); nIdx++)
    {
      double dPrice = ppHedgeRatioData[nIdx]->GetModelOutput()->GetPrice();
      std::cout << " Price " << nIdx << " = " << dPrice << std::endl;
    }
    std::cout << std::endl;

    // Output analysis date values
    if ( pHedgeOutput->HasHEROAtAnalysisDate() )
    {
      std::cout << "Analysis date values:" << std::endl;

      const finance::Values& pdValues 
        = pHedgeOutput->GetHEROValuesAtAnalysisDate();

      const finance::Values& pdSpots 
        = pHedgeOutput->GetHEROSpotsAtAnalysisDate();

      for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
        std::cout << pdSpots[nIdx] << " " << pdValues[nIdx] << std::endl;

    } // if analysis date 
    std::cout << std::endl;
    
    // Output part of hero surface
    if ( pHedgeOutput->HasHEROSurface() )
    {
      std::cout << "surface computed" << std::endl;
      
      finance::SharedSurface pSurface = pHedgeOutput->GetHEROSurface();

      shared_ptr<finance::Domain> pDomain = pSurface->GetDomain();

      finance::Domain::Spots pdSpots(20);
      pdSpots[0] = 0.8 * pSessionData->GetSpotSharePrice();
      double dStep = 0.1 * (pSessionData->GetSpotSharePrice() - pdSpots[0]);
      for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
        pdSpots[nIdxS] = pdSpots[nIdxS-1] + dStep;

      pDomain->SetUnderlyingSharePrices(pdSpots);

      size_t nNbDates = pDomain->GetDates().size();
      finance::SurfaceDouble::Doubles pValues 
        = pSurface->GetValuesAt(nNbDates-1);

      for (size_t nIdx = 0; nIdx < pValues.size(); nIdx++)
        std::cout << pdSpots[nIdx] << ", value = " << pValues[nIdx]
                  << std::endl;

    } // if surface
    std::cout << std::endl;

    //ShowHeroConvergence(pTarget, hedgeDerivatives, pModel, 4);

    std::cout << "Total time = " << sw() << std::endl;

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
