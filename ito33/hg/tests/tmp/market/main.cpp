#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/finance/eds.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/option.h"
#include "ito33/finance/optiontype.h"
#include "ito33/finance/impliedcdsspreads.h"
#include "ito33/finance/impliededsspreads.h"
#include "ito33/finance/derivatives.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"
#include "ito33/hg/parametrization.h"
#include "ito33/hg/parametrization_cdsrecovery.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/link.h"

#include "marketdata.h"

ITO33_FORCE_LINK_MODULE(HGPriceEDS);
ITO33_FORCE_LINK_MODULE(HGPriceCDS);
ITO33_FORCE_LINK_MODULE(HGPriceOption);


using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

// options with values less than this are removed from calibration
double g_dOptionCut = 0.1;

double g_dOptionWeightCut = 1.0;
double g_dOptionWeightLarge = 0.01;
double g_dOptionWeightSmall = 0.001;

// the weighting of the cds and eds contracts
//double g_dCDSEDSWeight = 100.0; 
                         
double g_dCDSEDSWeight = 10.0; 

std::vector<double> g_pdOptionWeights;

void ComputeOptionWeights(std::list< shared_ptr<Option> > optionList)
{
  // allocate space
  g_pdOptionWeights.resize( optionList.size() );

  // loop through the options, computing weight for each one
  size_t nCounter = 0;
  std::list< shared_ptr<finance::Option> >::const_iterator iterOption;
  for (iterOption = optionList.begin(); iterOption != optionList.end(); ++iterOption)
  {
    // Equal weighting
    g_pdOptionWeights[nCounter] = 1.0;

    nCounter++;
    
  }
}


void OutputSurfaces(shared_ptr<hg::TheoreticalModel> pModel,
                     std::list< shared_ptr<Option> > optionList)
{
  std::string sImpliedVol("impliedvol.txt");
  std::ofstream sImpliedVolOut(sImpliedVol.c_str(), std::ios::out);
  ASSERT_MSG(sImpliedVolOut, "Problem opening implied vol file");

  std::string sImpliedVolComp("impliedvolcomputed.txt");
  std::ofstream sImpliedVolCompOut(sImpliedVolComp.c_str(), std::ios::out);
  ASSERT_MSG(sImpliedVolCompOut, "Problem opening implied vol computed file");

  std::string sMarket("market.txt");
  std::ofstream sMarketOut(sMarket.c_str(), std::ios::out);
  ASSERT_MSG(sMarketOut, "Problem opening market file");

  std::string sErrorRel("errorrel.txt");
  std::ofstream sErrorRelOut(sErrorRel.c_str(), std::ios::out);
  ASSERT_MSG(sErrorRelOut, "Problem opening relative error file");

  std::string sErrorAbs("errorabs.txt");
  std::ofstream sErrorAbsOut(sErrorAbs.c_str(), std::ios::out);
  ASSERT_MSG(sErrorAbsOut, "Problem opening absolute error file");

  std::list< shared_ptr<Option> >::const_iterator iter;

  Date valuationDate = optionList.front()->GetSessionData()->GetValuationDate();
  double dValuationTime = GetDoubleFrom(valuationDate);

  Date tmp = optionList.front()->GetMaturityDate();
  double dLastTime = GetDoubleFrom(tmp);

  for (iter = optionList.begin(); iter != optionList.end(); ++iter)
  {
    if ( (*iter)->GetMarketPrice() < g_dOptionCut )
      continue;

    // Assume the option list is ordered
    Date maturityDate = (*iter)->GetMaturityDate();
    double dMaturityTime = GetDoubleFrom(maturityDate);

    if ( fabs(dLastTime - dMaturityTime) > 1.e-6)
    {
      sImpliedVolOut << std::endl;
      sImpliedVolCompOut << std::endl;
      sMarketOut << std::endl;
      sErrorRelOut << std::endl;
      sErrorAbsOut << std::endl;
    }

    shared_ptr<finance::ModelOutput> output = pModel->Compute(**iter);

    double dComputedPrice = output->GetPrice();

    const double dMarketPrice = (*iter)->GetMarketPrice();

    sImpliedVolOut << dMaturityTime - dValuationTime
                   << " "
                   << (*iter)->GetStrike()
                   << " "
                   << (*iter)->GetImpliedVol()
                   << std::endl;

    sImpliedVolCompOut << dMaturityTime - dValuationTime
                       << " "
                       << (*iter)->GetStrike()
                       << " "
                       << (*iter)->GetImpliedVolFrom(dComputedPrice)
                       << std::endl;
    
    sMarketOut << dMaturityTime - dValuationTime
               << " "
               << (*iter)->GetStrike()
               << " "
               << dMarketPrice
               << std::endl;

    sErrorAbsOut << dMaturityTime - dValuationTime
                 << " "
                 << (*iter)->GetStrike()
                 << " "
                 << dMarketPrice - dComputedPrice 
                 << std::endl;

    double dScale = dMarketPrice;
    if (dMarketPrice < 1.e-6)
      dScale = 1.0;
    sErrorRelOut << dMaturityTime - dValuationTime
                 << " "
                 << (*iter)->GetStrike()
                 << " "
                 << (dMarketPrice - dComputedPrice) / dScale
                 << std::endl;

    dLastTime = dMaturityTime;
  }

}


double OutputPriceData(shared_ptr<hg::TheoreticalModel> pModel,
                       shared_ptr<finance::Derivative> pDerivative)
{

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pDerivative);

  double dComputedPrice = output->GetPrice();

  const double dMarketPrice = pDerivative->GetMarketPrice();

  double dScale = fabs(dMarketPrice);
  if (dScale > 1.e-6)
    dScale = 1.0/dScale;
  else
    dScale = 1.0;

  double dRelError = fabs((dMarketPrice - dComputedPrice) * dScale);

  std::cout << "Market= " << dMarketPrice
            << ", calib= " << dComputedPrice
            << ", rel err= " << dRelError
            << std::endl;

  return dRelError;
}

double OutputImpliedVolData(shared_ptr<hg::TheoreticalModel> pModel,
                            shared_ptr<finance::Option> pOption)
{

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pOption);

  double dComputedPrice = output->GetPrice();

  double dComputedImpliedVol = pOption->GetImpliedVolFrom(dComputedPrice);

  const double dMarketImpliedVol = pOption->GetImpliedVol();

  double dScale = fabs(dMarketImpliedVol);
  if (dScale > 1.e-6)
    dScale = 1.0/dScale;
  else
    dScale = 1.0;

  double dRelError = fabs((dMarketImpliedVol - dComputedImpliedVol) * dScale);

  std::cout << "Market= " << dMarketImpliedVol
            << ", calib= " << dComputedImpliedVol
            << ", rel err= " << dRelError
            << std::endl;

  return dRelError;
}

void ReportError(double dTotalError, double dTotalWeight, size_t nNbDerivs)
{

  std::cout << "Sum of weighted squared relative error = " << dTotalError << std::endl; 
  std::cout << "Square root of weighted squared relative error = " << sqrt(dTotalError) << std::endl; 
  std::cout << "Square root of (total error / sum of weights ) = " << sqrt(dTotalError/dTotalWeight) << std::endl; 
  std::cout << "Square root of (total error / N ) = " << sqrt(dTotalError/nNbDerivs) << std::endl; 
  std::cout << "1/N * square root of (total error) = " << sqrt(dTotalError)/nNbDerivs << std::endl; 
  std::cout << std::endl;
  std::cout << std::endl;

}

double OutputSummary(shared_ptr<hg::TheoreticalModel> pModel,
                     std::list< shared_ptr<Option> > optionList,
                     std::list< shared_ptr<CDS> > cdsList,
                     std::list< shared_ptr<EDS> > edsList)
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

  // Output the cds recovery, in case it was calibrated
  std::cout << "CDS recovery rate = " << cdsList.front()->GetRecoveryRate()
            << std::endl;

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
  double dTotalWeight = 0.0;
  size_t nNbDerivs = 0;

  size_t nCounter;
  std::list< shared_ptr<finance::Option> >::const_iterator iterOptions;
  for (iterOptions = optionList.begin(), nCounter = 0;
       iterOptions != optionList.end();
       ++iterOptions, nCounter++)
  {
    if ( (*iterOptions)->GetMarketPrice() < g_dOptionCut)
      continue;

    if ( (*iterOptions)->GetOptionType() == Option_Call )
      std::cout << "Call ";
    else
      std::cout << "Call ";

    if ( (*iterOptions)->GetExerciseType() == ExerciseType_European )
      std::cout << "Euro ";
    else
      std::cout << "Amer ";

    std::cout << (*iterOptions)->GetStrike() << " ";

    Date maturityDate = (*iterOptions)->GetMaturityDate();
    std::cout << maturityDate << " ";

    double dRelError = OutputPriceData(pModel, *iterOptions );
    //double dRelError = OutputImpliedVolData(pModel, *iterOptions );

    double dWeight = g_pdOptionWeights[nCounter];

    dTotalWeight += dWeight;
    dTotalError += dWeight * dRelError * dRelError; 
    nNbDerivs++;

  } // loop over the options

  std::cout << std::endl;
  ReportError(dTotalError, dTotalWeight, nNbDerivs);



  std::list< shared_ptr<finance::CDS> >::const_iterator iterCDS;
  for (iterCDS = cdsList.begin();
       iterCDS != cdsList.end();
       ++iterCDS)
  {
    std::cout << "CDS: ";

    std::cout << (*iterCDS)->GetMaturityDate() << " ";

    double dRelError = OutputPriceData(pModel, *iterCDS);

    // Convert into spreads
    Date valuationDate = (*iterCDS)->GetSessionData()->GetValuationDate();
    Date firstDate = (*iterCDS)->GetSpreadStream()->GetFirstPaymentDate();
    Date::DayCountConvention dayCount = 
      (*iterCDS)->GetSpreadStream()->GetDayCountConvention();
    finance::Frequency freq = 
      (*iterCDS)->GetSpreadStream()->GetPaymentFrequency();

    double dRecoveryRate = (*iterCDS)->GetRecoveryRate();
    
    ImpliedCDSSpreads impliedSpreads(valuationDate,
                                     firstDate, 
                                     dayCount,
                                     freq, 
                                     dRecoveryRate);

    std::vector<Date> pDates(1);
    pDates[0] = (*iterCDS)->GetMaturityDate();

    std::vector<double> pdSpreads = 
      impliedSpreads.Compute(pModel, (*iterCDS)->GetSessionData(), pDates);

    double dMarketSpread = (*iterCDS)->GetSpreadStream()->GetAnnualPaymentAmount();
    double dComputedSpread = pdSpreads[0];

    std::cout << "     market spread = " << dMarketSpread
              << " implied spread = " << dComputedSpread 
              << std::endl;

    dRelError = (dMarketSpread - dComputedSpread) / dMarketSpread;

    dTotalError += g_dCDSEDSWeight * dRelError * dRelError; 
    dTotalWeight += g_dCDSEDSWeight;
    nNbDerivs++;

  } // loop over the cds
  std::cout << std::endl;
  ReportError(dTotalError, dTotalWeight, nNbDerivs);


  std::list< shared_ptr<finance::EDS> >::const_iterator iterEDS;
  for (iterEDS = edsList.begin();
       iterEDS != edsList.end();
       ++iterEDS)
  {
    std::cout << "EDS: ";

    std::cout << (*iterEDS)->GetMaturityDate() << " ";

    double dRelError = OutputPriceData(pModel, *iterEDS);

    // Convert into spreads
    Date valuationDate = (*iterEDS)->GetSessionData()->GetValuationDate();
    Date firstDate = (*iterEDS)->GetSpreadStream()->GetFirstPaymentDate();
    Date::DayCountConvention dayCount = 
      (*iterEDS)->GetSpreadStream()->GetDayCountConvention();
    finance::Frequency freq = 
      (*iterEDS)->GetSpreadStream()->GetPaymentFrequency();

    double dRecoveryRate = (*iterEDS)->GetRecoveryRate();
    double dBarrier = (*iterEDS)->GetBarrier();
   
    ImpliedEDSSpreads impliedSpreads(valuationDate, firstDate, dayCount,
                                     freq, dBarrier, dRecoveryRate);

    std::vector<Date> pDates(1);
    pDates[0] = (*iterEDS)->GetMaturityDate();

    std::vector<double> pdSpreads = 
      impliedSpreads.Compute(pModel, (*iterEDS)->GetSessionData(), pDates);

    double dMarketSpread = (*iterEDS)->GetSpreadStream()->GetAnnualPaymentAmount();
    double dComputedSpread = pdSpreads[0] ;
    std::cout << "     market spread = " << dMarketSpread
              << " implied spread = " << dComputedSpread
              << std::endl;

    dRelError = (dMarketSpread - dComputedSpread) / dMarketSpread;

    dTotalError += g_dCDSEDSWeight * dRelError * dRelError; 
    dTotalWeight += g_dCDSEDSWeight;
    nNbDerivs++;

  } // loop over the eds
  std::cout << std::endl;
  ReportError(dTotalError, dTotalWeight, nNbDerivs);

  return dTotalError;
}


void ConvertToPrices(std::list< shared_ptr<Option> > options, 
                     shared_ptr<SessionData> pSessionData)
{

  std::list< shared_ptr<Option> >::iterator iterOptions;

  for (iterOptions = options.begin();
       iterOptions != options.end(); 
       ++iterOptions)
  {

    // Setup simple Black-Scholes like process
    size_t nNbRegimes = 1;
    double dImpliedVol = (*iterOptions)->GetMarketPrice();

    std::vector<double> pdVols(nNbRegimes, dImpliedVol);
    std::vector<double> pdIntensities(nNbRegimes, 0.0);

    shared_ptr<hg::UnderlyingProcess> pBSProcess
      ( new hg::UnderlyingProcess(nNbRegimes, pdVols, pdIntensities) );

    hg::TheoreticalModel model(pBSProcess); 

    // Now price and set as market price
    shared_ptr<finance::ModelOutput> output = model.Compute(**iterOptions);

    double dPrice = output->GetPrice();

    (*iterOptions)->SetMarketPrice(dPrice);
  } // end loop over options
}

int TestMarketData(Company company, bool bIsAdi)
{
  try
  {   
    shared_ptr<SessionData> pSessionData = InitSessionDataMarket(company);
    shared_ptr<hg::UnderlyingProcess>
      pInitialProcess = GetUnderlyingProcess(company);

    // Get the appropriate input file (date from Elie/Adi or polygon)
    std::string sInputFile = GetInputFile(company, bIsAdi);

    // Create some contracts to price        
    std::list< shared_ptr<Option> > optionList;
    if (bIsAdi)
      optionList = ReadOptionData2(sInputFile, pSessionData);
    else
      optionList = ReadOptionData(sInputFile, pSessionData);
      
    std::list< shared_ptr<CDS> > cdsList = GetCDSData(company, pSessionData);
    std::list< shared_ptr<EDS> > edsList = GetEDSData(company, pSessionData);

    // Copy to one list
    Derivatives derivatives;
/*
    std::list< shared_ptr<finance::CDS> >::const_iterator iterCDS;
    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives.AddWithWeight(*iterCDS, 1.0);

    std::list< shared_ptr<finance::EDS> >::const_iterator iterEDS;
    for (iterEDS = edsList.begin(); iterEDS != edsList.end(); ++iterEDS)
      derivatives.AddWithWeight(*iterEDS, 1.0);
*/
    std::vector<double> pdVegas(optionList.size());
    size_t nVegaCounter = 0;
    std::list< shared_ptr<finance::Option> >::const_iterator iterOption;
    for (iterOption = optionList.begin(); iterOption != optionList.end(); ++iterOption)
    {
      /*
      ihg::TheoreticalModel model;
      shared_ptr<ihg::ComputationalFlags> pFlags(new ihg::ComputationalFlags);
      pFlags->SetComputeVega(true);
      model.SetComputationalFlags(pFlags);

      double dVol = (*iterOption)->GetImpliedVol();
      shared_ptr<ihg::Volatility> pVolatility(new ihg::VolatilityFlat(dVol));

      shared_ptr<ihg::HazardRate> pHazardRate(new ihg::HazardRateFlat(0.0));
      
      model.SetVolatility(pVolatility);
      model.SetHazardRate(pHazardRate);

      shared_ptr<ihg::ModelOutput> pOutput = model.Compute(**iterOption);

      double dVega = pOutput->GetVega();
      std::cout << "vega = " << dVega << std::endl;

      dVega = std::min( 1.0/(dVega+1.e-12), 100.0);
      
      pdVegas[nVegaCounter++] = dVega;

      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives.AddWithWeight(*iterOption, dVega);
      */


      double dPrice = (*iterOption)->GetMarketPrice();
      if ( dPrice >= 1.0 )
        derivatives.AddWithWeight(*iterOption, 0.05);
      else if ( dPrice >= 0.1 )
        derivatives.AddWithWeight(*iterOption, 0.005);
      
      // Do not calibrate options with price less than 10 cents
    }

    // Now try calibration 
    hg::Parametrization parametrization(pInitialProcess);
    shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess;

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

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));

    OutputSummary(pModel, optionList, cdsList, edsList);
    OutputSurfaces(pModel, optionList);

    // Round two
    hg::Parametrization parametrization2(pUnderlyingProcess);
    Derivatives derivatives2;

    nVegaCounter = 0;
    for (iterOption = optionList.begin(); iterOption != optionList.end(); ++iterOption)
    {
      /*
      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives2.AddWithWeight(*iterOption, pdVegas[nVegaCounter]);

      nVegaCounter++;
      */

      double dPrice = (*iterOption)->GetMarketPrice();
      if ( dPrice >= 1.0 )
        derivatives2.AddWithWeight(*iterOption, 0.05);
      else if ( dPrice >= 0.1 )
        derivatives2.AddWithWeight(*iterOption, 0.005);
    }

    std::list< shared_ptr<finance::CDS> >::const_iterator iterCDS;
    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives2.AddWithWeight(*iterCDS, g_dCDSEDSWeight);

    try
    {
      pUnderlyingProcess = parametrization2.Calibrate(derivatives2);
    }
    catch(ito33::Exception e)
    {
      std::cout << "Calibration failed!!!!!! Prices below are best guesses only" 
                << std::endl << std::endl;
      pUnderlyingProcess = parametrization2.GetCalibratedUnderlyingProcess();
    }

    pModel = make_ptr( new hg::TheoreticalModel(pUnderlyingProcess) );

    OutputSummary(pModel, optionList, cdsList, edsList);



    // Round three
    hg::Parametrization parametrization3(pUnderlyingProcess);
    Derivatives derivatives3;

    nVegaCounter = 0;
    for (iterOption = optionList.begin(); iterOption != optionList.end(); ++iterOption)
    {
      /*
      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives3.AddWithWeight(*iterOption, pdVegas[nVegaCounter]);

      nVegaCounter++;
      */

      double dPrice = (*iterOption)->GetMarketPrice();
      if ( dPrice >= 1.0 )
        derivatives3.AddWithWeight(*iterOption, 0.05);
      else if ( dPrice >= 0.1 )
        derivatives3.AddWithWeight(*iterOption, 0.005);
    }


    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives3.AddWithWeight(*iterCDS, g_dCDSEDSWeight);

    std::list< shared_ptr<finance::EDS> >::const_iterator iterEDS;
    for (iterEDS = edsList.begin(); iterEDS != edsList.end(); ++iterEDS)
      derivatives3.AddWithWeight(*iterEDS, g_dCDSEDSWeight);

    try
    {
      pUnderlyingProcess = parametrization3.Calibrate(derivatives3);
    }
    catch(ito33::Exception e)
    {
      std::cout << "Calibration failed!!!!!! Prices below are best guesses only" 
                << std::endl << std::endl;
      pUnderlyingProcess = parametrization3.GetCalibratedUnderlyingProcess();      
    }


    pModel = make_ptr( new hg::TheoreticalModel(pUnderlyingProcess)  );

    OutputSummary(pModel, optionList, cdsList, edsList);

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


int TestMarketDataCDSRecovery(Company company, bool bIsAdi)
{
  try
  {   
    bool bUseSpreads = true;

    shared_ptr<SessionData> pSessionData = InitSessionDataMarket(company);
    shared_ptr<hg::UnderlyingProcess> 
      pInitialProcess = GetUnderlyingProcess(company);

    // Get the appropriate input file (date from Elie/Adi or polygon)
    std::string sInputFile = GetInputFile(company, bIsAdi);

    // Create some contracts to price        
    std::list< shared_ptr<Option> > optionList;
    if (bIsAdi)
      optionList = ReadOptionData2(sInputFile, pSessionData);
    else
      optionList = ReadOptionData(sInputFile, pSessionData);
      
    std::list< shared_ptr<CDS> > cdsList = GetCDSData(company, pSessionData);
    std::list< shared_ptr<EDS> > edsList = GetEDSData(company, pSessionData);

    // Copy to one list
    Derivatives derivatives;
    
    // compute option weights and copy to master list
    ComputeOptionWeights( optionList );

    size_t nCounter = 0;
    std::list< shared_ptr<finance::Option> >::const_iterator iterOption;
    for (iterOption = optionList.begin(); 
         iterOption != optionList.end(); 
         ++iterOption)
    {

      // Do not calibrate options with price less than some cutoff (eg 30 cents)
      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives.AddWithWeight(*iterOption, g_pdOptionWeights[nCounter]);

      nCounter++;
    }

    // Now try calibration 
    std::cout << "Round one" << std::endl;
    ParametrizationCDSRecovery parametrization(pInitialProcess);
    parametrization.CalibrateWithSpreads(bUseSpreads);

    double dCalibratedCDSRecovery = -1.0;
    shared_ptr<hg::UnderlyingProcess> pUnderlyingProcess;

    try
    {
      pUnderlyingProcess = parametrization.Calibrate(derivatives);
      dCalibratedCDSRecovery = parametrization.GetCalibratedCDSRecovery();
    }
    catch (const ito33::Exception& /* e */)
    {
      std::cout << "Calibration failed!!!! Prices below are best guesses only"
                << std::endl << std::endl;
      pUnderlyingProcess = parametrization.GetCalibratedUnderlyingProcess();
      dCalibratedCDSRecovery = parametrization.GetCalibratedCDSRecovery();
    }

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));

    // reconstruct cds list using the calibrated recovery value
    std::list< shared_ptr<CDS> > cdsListCalibrated;
    std::list< shared_ptr<finance::CDS> >::const_iterator iterCDS;
    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
    {
      shared_ptr<CDS> pCDS(new CDS(dCalibratedCDSRecovery,
                                  (*iterCDS)->GetSpreadStream() ));
      pCDS->SetSessionData( (*iterCDS)->GetSessionData() );
      pCDS->SetMarketPrice( (*iterCDS)->GetMarketPrice() );
      cdsListCalibrated.push_back(pCDS);
    }

    OutputSummary(pModel, optionList, cdsListCalibrated, edsList);
    OutputSurfaces(pModel, optionList);

    // Round two
    std::cout << "Round two" << std::endl;
    ParametrizationCDSRecovery parametrization2(pUnderlyingProcess);
    parametrization2.CalibrateWithSpreads(bUseSpreads);
    Derivatives derivatives2;

    nCounter = 0;
    for (iterOption = optionList.begin(); 
         iterOption != optionList.end(); 
         ++iterOption)
    {
      
      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives2.AddWithWeight(*iterOption, g_pdOptionWeights[nCounter]);

      nCounter++;
    }

    //std::list< shared_ptr<finance::CDS> >::const_iterator iterCDS;
    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives2.AddWithWeight(*iterCDS, g_dCDSEDSWeight);

    try
    {
      pUnderlyingProcess = parametrization2.Calibrate(derivatives2);
      dCalibratedCDSRecovery = parametrization2.GetCalibratedCDSRecovery();
    }
    catch (const ito33::Exception& /* e */)
    {
      std::cout << "Calibration failed!!!! Prices below are best guesses only"
                << std::endl << std::endl;
      pUnderlyingProcess = parametrization2.GetCalibratedUnderlyingProcess();
      dCalibratedCDSRecovery = parametrization2.GetCalibratedCDSRecovery();
    }

    cdsListCalibrated.clear();
    pModel = make_ptr( new hg::TheoreticalModel(pUnderlyingProcess) );

    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
    {
      shared_ptr<CDS> pCDS(new CDS(dCalibratedCDSRecovery,
                                  (*iterCDS)->GetSpreadStream() ));
      pCDS->SetSessionData( (*iterCDS)->GetSessionData() );
      pCDS->SetMarketPrice( (*iterCDS)->GetMarketPrice() );
      cdsListCalibrated.push_back(pCDS);
    }

    OutputSummary(pModel, optionList, cdsListCalibrated, edsList);



    // Round three
    std::cout << "Round three" << std::endl;
    ParametrizationCDSRecovery parametrization3(pUnderlyingProcess);
    parametrization3.CalibrateWithSpreads(bUseSpreads);
    Derivatives derivatives3;

    nCounter = 0;
    for (iterOption = optionList.begin(); 
         iterOption != optionList.end(); 
         ++iterOption)
    {
      
      double dPrice = (*iterOption)->GetMarketPrice();
      if (dPrice >= g_dOptionCut)
        derivatives3.AddWithWeight(*iterOption, g_pdOptionWeights[nCounter]);

      nCounter++;
    }


    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
      derivatives3.AddWithWeight(*iterCDS, g_dCDSEDSWeight);

    std::list< shared_ptr<finance::EDS> >::const_iterator iterEDS;
    for (iterEDS = edsList.begin(); iterEDS != edsList.end(); ++iterEDS)
      derivatives3.AddWithWeight(*iterEDS, g_dCDSEDSWeight);

    try
    {
      pUnderlyingProcess = parametrization3.Calibrate(derivatives3);
      dCalibratedCDSRecovery = parametrization3.GetCalibratedCDSRecovery();
    }
    catch (const ito33::Exception& /* e */)
    {
      std::cout << "Calibration failed!!!! Prices below are best guesses only"
                << std::endl << std::endl;
      pUnderlyingProcess = parametrization3.GetCalibratedUnderlyingProcess();      
      dCalibratedCDSRecovery = parametrization3.GetCalibratedCDSRecovery();
    }


    cdsListCalibrated.clear();
    pModel = make_ptr( new hg::TheoreticalModel(pUnderlyingProcess) );

    for (iterCDS = cdsList.begin(); iterCDS != cdsList.end(); ++iterCDS)
    {
      shared_ptr<CDS> pCDS(new CDS(dCalibratedCDSRecovery,
                                  (*iterCDS)->GetSpreadStream() ));
      pCDS->SetSessionData( (*iterCDS)->GetSessionData() );
      pCDS->SetMarketPrice( (*iterCDS)->GetMarketPrice() );
      cdsListCalibrated.push_back(pCDS);
    }

    OutputSummary(pModel, optionList, cdsListCalibrated, edsList);

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

  // pass true for adi data, false for bloomberg data
  
  //TestMarketData(Accor, true);
  //TestMarketData(Valeo, true);
  //TestMarketData(Renault, true);   // fits eds, close for cds
  //TestMarketData(FranceTelecom, true); // eds pretty good, cds worse
  //TestMarketData(Lafarge, true);

  //TestMarketDataCDSRecovery(Accor, true); // eds1 and cds1 both over
  //TestMarketDataCDSRecovery(Valeo, true); // terrible fit with options
  //TestMarketDataCDSRecovery(Carrefour, true); // terrible fit with option
  //TestMarketDataCDSRecovery(FranceTelecom, true); //terrible fit with options
  TestMarketDataCDSRecovery(Lafarge, true); // good fit to cds/eds
  //TestMarketDataCDSRecovery(BNBParibas, true); // not bad fit to cds/eds
  //TestMarketDataCDSRecovery(Endesea, true); // excellent fit to cds/eds
  //TestMarketDataCDSRecovery(Renault, true); // off for first cds

}

