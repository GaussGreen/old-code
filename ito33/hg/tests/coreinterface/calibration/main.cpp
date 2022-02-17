#include <iostream>
#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/finance/derivatives.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/option.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/parametrization.h"
#include "ito33/hg/multioutput.h"

#include "hg/priceforwardoption.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceEDS);
ITO33_FORCE_LINK_MODULE(HGPriceCDS);
ITO33_FORCE_LINK_MODULE(HGPriceOption);
ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwap);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

shared_ptr<SessionData> InitSessionDataGenerated()
{
  Date valuationDate(2002, Date::Jul, 1);

  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR") );

  double dSpot = 45;
  shared_ptr<Equity> pEquity(new Equity(45, pCurrency));
  
  double dPreviousSpot = dSpot;
  pEquity->SetPreviousSharePrice(dPreviousSpot);

  shared_ptr<Dividends> pDividends( new Dividends() );
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
  pEquity->SetBorrowCurve(pyf);
    
  double dRate = 0.05;
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

  shared_ptr<RateData> pRateData( new RateData());
  pRateData->SetYieldCurve(pCurrency, pyc);

  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

shared_ptr<Option> InitOption(double dStrike, 
                             Date maturityDate, 
                             OptionType optionType)
{
  shared_ptr<Option> pOption(new 
    Option(dStrike, maturityDate, optionType, ExerciseType_European) );

  return pOption;
}

shared_ptr<EDS> InitEDS(Date issueDate,
                       Date firstDate,
                       Date lastDate,
                       double dSpread,
                       double dBarrier)
{
  double dRecoveryRate = 0.5;

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                        )
                  );

  shared_ptr<EDS>
    pEDS( new EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  return pEDS;
}

shared_ptr<CDS> InitCDS(Date issueDate,
                       Date firstDate,
                       Date lastDate,
                       double dSpread)
{
  double dRecoveryRate = 0.0;

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                        )
                  );

  shared_ptr<CDS>
    pCDS( new CDS(dRecoveryRate, pSpreadStream) );

  return pCDS;
}

std::list< shared_ptr<Option> > 
ConstructOptionList(shared_ptr<hg::TheoreticalModel> pModel, 
                     shared_ptr<SessionData> pSessionData )
{

  // the return derivatives
  std::list< shared_ptr<Option> > optionList;

  // Get session data
  double dSpot = pSessionData->GetSpotSharePrice();
  Date valuationDate = pSessionData->GetValuationDate();

  // Construct some options
  double dStrikeMin = dSpot / 2.0;
  double dStrikeDelta = (dSpot - dStrikeMin) / 2.0;
  double dStrikeMax = dSpot * 1.5;
  //OptionType optionType = Option_Put;
  OptionType optionType = Option_Call;

  size_t nNbYears = 2;
  Date maturityDate = valuationDate;
  for (size_t nIdx = 0; nIdx < nNbYears; nIdx++)
  {
    maturityDate.AddYears(1);

    double dStrike = dStrikeMin;
    while (dStrike <= dStrikeMax + 1.e-8)
    {

      shared_ptr<Option> pOption = InitOption(dStrike, maturityDate, optionType );
      pOption->SetSessionData(pSessionData);
      shared_ptr<finance::ModelOutput> output = pModel->Compute(*pOption);

      pOption->SetMarketPrice( output->GetPrice() );

      optionList.push_back(pOption);

      dStrike += dStrikeDelta;

    } // loop over strike

  } // loop over maturities

  return optionList;
}

void PriceUsingForward(std::list< shared_ptr<Option> > optionList, 
                       shared_ptr<hg::TheoreticalModel> pModel)
{
  // price using a forward
  ForwardOption forward( optionList );
  
  shared_ptr<MultiOutput> 
    forwardOutput( PriceForwardOption(*pModel, forward) );

  std::vector<double> pdPrices;
  pdPrices = forwardOutput->GetPrices();

  // Update the market prices
  size_t nIdxOption = 0;
  std::list< shared_ptr<Option> >::const_iterator iterOptions;
  for (iterOptions = optionList.begin();
       iterOptions != optionList.end();
       ++iterOptions)
  {

    double dOldMarketPrice = (*iterOptions)->GetMarketPrice();
    double dNewMarketPrice = pdPrices[nIdxOption];
    (*iterOptions)->SetMarketPrice( dNewMarketPrice );
    nIdxOption++;

    std::cout << "backward price = " << dOldMarketPrice
              << ", forward price = " << dNewMarketPrice
              << std::endl;
  }
}

std::list< shared_ptr<CDS> > 
ConstructCDSList(shared_ptr<hg::TheoreticalModel> pModel, 
                 shared_ptr<SessionData> pSessionData )
{
  // the return derivatives
  std::list< shared_ptr<CDS> > cdsList;

  // Get session data
  Date valuationDate = pSessionData->GetValuationDate();

  // Construct some CDSes  
  size_t nNbYears = 4;
  Date issueDate = valuationDate;
  Date firstDate = issueDate;
  firstDate.AddYears(1);
  Date lastDate = firstDate;
  lastDate.AddYears(3);
  double dSpread = 0.25;

  for (size_t nIdx = 0; nIdx < nNbYears; nIdx++)
  {    
    dSpread += nIdx * 0.02;
    shared_ptr<CDS> pCDS = InitCDS(issueDate, firstDate, lastDate, dSpread);
    pCDS->SetSessionData(pSessionData);
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCDS);
    double dPrice = output->GetPrice();
    pCDS->SetMarketPrice( dPrice );

    cdsList.push_back(pCDS);

    lastDate.AddYears(1);
  }

  return cdsList;
}

std::list< shared_ptr<EDS> > 
ConstructEDSList(shared_ptr<hg::TheoreticalModel> pModel, 
                 shared_ptr<SessionData> pSessionData )
{
  // the return derivatives
  std::list< shared_ptr<EDS> > edsList;

  // Get session data
  Date valuationDate = pSessionData->GetValuationDate();

  // Construct some EDSes
  size_t nNbYears = 4;
  Date issueDate = valuationDate;
  Date firstDate = issueDate;
  firstDate.AddYears(1);
  Date lastDate = firstDate;
  lastDate.AddYears(3);
  double dSpread = 0.25;
  double dBarrier = pSessionData->GetSpotSharePrice() * 0.30;
  
  for (size_t nIdx = 0; nIdx < nNbYears; nIdx++)
  {    
    dSpread += nIdx * 0.01;
    shared_ptr<EDS> pEDS = InitEDS(issueDate, firstDate, lastDate, dSpread, dBarrier);
    pEDS->SetSessionData(pSessionData);
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pEDS);
    double dPrice = output->GetPrice();
    pEDS->SetMarketPrice( dPrice );

    edsList.push_back(pEDS);

    lastDate.AddYears(1);
  }

  return edsList;
}

std::list< shared_ptr<VarianceSwap> >
ConstructVSList( shared_ptr<hg::TheoreticalModel> pModel, 
                 shared_ptr<SessionData> pSessionData )
{
  // the return derivatives
  std::list< shared_ptr<VarianceSwap> > list;

  Date startOfSamplingPeriod = pSessionData->GetValuationDate();
  size_t nNb = 3;
  for (size_t n = 0; n < nNb; n++)
  {
    Date maturityDate = startOfSamplingPeriod;
    maturityDate.AddMonths( 4 * ( n + 1) );
    double dVolatilityStrike = 0.1 + 0.02 * n;
    shared_ptr<VarianceSwapTerms>
      pTerms( new VarianceSwapTerms
                  ( maturityDate,
                    Swap_Variance, 
                    startOfSamplingPeriod,
                    Date::DaysDiff(startOfSamplingPeriod, maturityDate) ) );

    shared_ptr<VarianceSwap> 
      pVS( new VarianceSwap(pTerms, dVolatilityStrike) );

    pVS->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pVS);
    double dPrice = output->GetPrice();
    pVS->SetMarketPrice( dPrice );

    list.push_back(pVS);
  }

  return list;
}

double OutputSummary(shared_ptr<hg::TheoreticalModel> pModel,
                     const Derivatives::Elements& derivatives)
{
  size_t nIdx;
  shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess 
    = pModel->GetUnderlyingProcess();
  size_t nNbRegimes = pUnderlyingProcess->GetNbRegimes();

  // Output the calibrated volatilities
  std::vector<double> pdVols = pUnderlyingProcess->GetVolatilities();
  std::cout << "Calibrated volatilities" << std::endl;
  for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
    std::cout << "  " << pdVols[nIdx] << std::endl;

  std::cout << std::endl;
    
  // Output the calibrated default intensities
  std::vector<double> pdIntensities = pUnderlyingProcess->GetJumpsToDefault();
  std::cout << "Calibrated intensities" << std::endl;
  for (nIdx = 0; nIdx < nNbRegimes; nIdx++)  
    std::cout << "  " << pdIntensities[nIdx] << std::endl;

  std::cout << std::endl;

  // Output the calibrated jumps
  std::cout << "Jumps" << std::endl;
  size_t nIdx1, nIdx2;
  for (nIdx1 = 0; nIdx1 < nNbRegimes; nIdx1++)
  {
    for (nIdx2 = 0; nIdx2 < nNbRegimes; nIdx2++)
    {
      std::cout << "  " << nIdx1 << " --> " << nIdx2 << std::endl; 
      Jumps jumps = pUnderlyingProcess->GetJumps(nIdx1, nIdx2);

      Jumps::const_iterator iter;
      for (iter = jumps.begin(); iter != jumps.end(); ++iter)
      {
        std::cout << "    Intensity = " << iter->GetIntensity()
                  << ", amplitude = " << iter->GetAmplitude()
                  << std::endl;
      }

    } // loop over second regime
  } // loop over first regime
  std::cout << std::endl;

  // Price using the calibrated model
  std::cout << "Calibrated contracts" << std::endl;
  double dTotalError = 0.0;
  Derivatives::Elements::const_iterator iter;
  for (iter = derivatives.begin(); iter != derivatives.end(); ++iter)
  {
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*iter->first);

    double dComputedPrice = output->GetPrice();

    const double dMarketPrice = iter->first->GetMarketPrice();

    double dScale = fabs(dMarketPrice);
    if (dScale > 1.e-5)
      dScale = 1.0/dScale;
    else
      dScale = 1.0;

    double dRelError = fabs((dMarketPrice - dComputedPrice) * dScale);

    std::cout << "Market= " << dMarketPrice
              << ", calibrated= " << dComputedPrice
              << ", rel error= " << dRelError
              << std::endl;

    dTotalError += dRelError * dRelError;
  } // loop over the derivatives

  std::cout << "Sum of squared relative error = " << dTotalError << std::endl; 
  return dTotalError;
}


int TestGeneratedData()
{
  try
  {   
    shared_ptr<SessionData> pSessionData = InitSessionDataGenerated();

    size_t nNbRegimes = 2;

    std::vector<double> pdVols(nNbRegimes);
    std::vector<double> pdIntensities(nNbRegimes);

    size_t nIdx;
    for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
    {
      pdVols[nIdx] = 0.1 + nIdx * 0.2;
      pdIntensities[nIdx] = 0.01 + nIdx * 0.1;
    }

    shared_ptr<hg::UnderlyingProcess> 
      pUnderlyingProcess( new hg::UnderlyingProcess
                              (nNbRegimes, pdVols, pdIntensities) );

    std::vector<double> pdIntensity(1);
    std::vector<double> pdAmplitude(1);

    pdIntensity[0] = 0.1;
    pdAmplitude[0] = -0.2;
    pUnderlyingProcess->SetJumps(0, 0, pdIntensity, pdAmplitude);

    if (nNbRegimes > 1)
    {
      
      pdIntensity[0] = 0.15;
      pdAmplitude[0] = -0.25;
      pUnderlyingProcess->SetJumps(0, 1, pdIntensity, pdAmplitude);

      pdIntensity[0] = 0.11;
      pdAmplitude[0] = 0.35;
      pUnderlyingProcess->SetJumps(1, 0, pdIntensity, pdAmplitude);

      pdIntensity[0] = 0.2;
      pdAmplitude[0] = -0.5;
      pUnderlyingProcess->SetJumps(1, 1, pdIntensity, pdAmplitude);
 
    }

    shared_ptr<hg::TheoreticalModel> 
    pModel( new hg::TheoreticalModel(pUnderlyingProcess) );

    // Create some contracts to price
    std::list< shared_ptr<Option> > 
      optionList( ConstructOptionList(pModel, pSessionData) );
    
    std::list< shared_ptr<CDS> > 
      cdsList( ConstructCDSList(pModel, pSessionData) );
    
    std::list< shared_ptr<EDS> > 
      edsList( ConstructEDSList(pModel, pSessionData) );

    std::list< shared_ptr<VarianceSwap> > 
      vsList( ConstructVSList(pModel, pSessionData) );

    PriceUsingForward(optionList, pModel);

    //cdsList.clear();
    //edsList.clear();
    //optionList.clear();
    //vsList.clear();

    Derivatives derivatives;

    // for printing at the end...
    std::list< shared_ptr<CDS> >::const_iterator iterCDS;
    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives.AddWithWeight(*iterCDS, 1.0);

    std::list< shared_ptr<EDS> >::const_iterator iterEDS;
    for (iterEDS = edsList.begin(); iterEDS != edsList.end(); ++iterEDS)
      derivatives.AddWithWeight(*iterEDS, 1.0);

    std::list< shared_ptr<Option> >::const_iterator iterOp;
    for (iterOp = optionList.begin(); iterOp != optionList.end(); ++iterOp)
    {
      (*iterOp)->SetSessionData(pSessionData);
      derivatives.AddWithWeight(*iterOp, 0.05);
    }

    std::list< shared_ptr<VarianceSwap> >::const_iterator iterVS;
    for (iterVS = vsList.begin(); iterVS != vsList.end(); ++iterVS)
    {
      (*iterVS)->SetSessionData(pSessionData);
      derivatives.AddWithWeight(*iterVS, 10);
    }

    // Change the data so that calibration has something to do
    for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
    {
      pdVols[nIdx] = 0.2;
      pdIntensities[nIdx] = 0.05;
    }

    // Construct the underlying process
    shared_ptr<hg::UnderlyingProcess> 
      pUnderlyingProcessInitial( new hg::UnderlyingProcess
                                     (nNbRegimes, pdVols, pdIntensities) );

    size_t nIdx1, nIdx2;
    for (nIdx1 = 0; nIdx1 < nNbRegimes; nIdx1++)
    {
      for (nIdx2 = 0; nIdx2 < nNbRegimes; nIdx2++)
      {

        // Get the old jumps
        Jumps jumps = pModel->GetUnderlyingProcess()->GetJumps(nIdx1, nIdx2);
        Jumps::const_iterator iterJumps;

        // Make new jumps with default values
        Jumps newJumps;
        for (iterJumps = jumps.begin(); iterJumps != jumps.end(); ++iterJumps)
        {
          //double dAmplitude = (*iterJumps).GetAmplitude();
          //double dIntensity = (*iterJumps).GetIntensity();
          //newJumps.push_back( Jump(dIntensity, dAmplitude) );
          
          newJumps.push_back( Jump(0.05, -0.01) );
        }

        pUnderlyingProcessInitial->SetJumps(nIdx1, nIdx2, newJumps);
      }
    }

    // Now try calibration 
    hg::Parametrization parametrization(pUnderlyingProcessInitial);

    try
    {
      pUnderlyingProcess = parametrization.Calibrate(derivatives);
    }
    catch(ito33::Exception e)
    {
      std::cout << "Calibration failed!!!!!! Prices below are best guesses only" 
                << std::endl << std::endl;

      pUnderlyingProcess = parametrization.GetCalibratedUnderlyingProcess();
    }

    pModel = make_ptr( new hg::TheoreticalModel(pUnderlyingProcess) );
    
    OutputSummary(pModel, derivatives.GetAll());

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

int main()
{

  TestGeneratedData();

}

