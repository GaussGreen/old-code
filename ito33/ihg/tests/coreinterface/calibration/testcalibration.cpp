#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include <fstream>
#include "ito33/afterstd.h"

#include <cmath>

#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/calibrationprocess.h"

#include "ito33/finance/bondlike/convertiblelike.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/volatilitycallback.h"

#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratecallback.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/version.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_volflat_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volflat_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_volpower.h"
#include "ito33/ihg/parametrization_volpower_hrpower.h"
#include "ito33/ihg/parametrization_voltanh.h"
#include "ito33/ihg/parametrization_volwithtimecomponent.h"
#include "ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h"

#include "ihg/translator.h"
#include "ihg/calibratorgeneral.h"

#include "ito33/xml/write.h"
#include "ihg/xml/common.h"

using namespace std; 
using namespace ito33;
using namespace ito33::finance;

ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceEDS);
ITO33_FORCE_LINK_MODULE(IHGPriceOption);

void __stdcall parametricHR(double /* dTime */, const double *pdS, double *pdValues, 
                            size_t nNbS, int /* i */)
{
  // Implement the standard HR formula hr(S,t) = alpha * (S0/S)^beta

  // Set default parameters
  const double dAlpha = 0.2;
  const double dBeta = 0.5;
  const double dS0 = 100.0;

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    if (pdS[nIdx] < 1.e-16)
      pdValues[nIdx] = dAlpha * pow((dS0 / 1.e-16), dBeta);
    else
      pdValues[nIdx] = dAlpha * pow((dS0 / pdS[nIdx]), dBeta);
  }
}


void WriteToXML(std::string sFileName, 
                shared_ptr<SessionData> pSessionData,
                const finance::TermStructureCDS& tsCDS, 
                ihg::ParametrizationVolFlatHRWithSpotComponentPower& parametrization);

void MakeCDSList(ito33::finance::TermStructureCDS&);
void MakeEDSList(ito33::finance::TermStructureEDS&);
shared_ptr<SessionData> InitSessionData(Date issueDate);
shared_ptr<TermStructureCDS> MakeCDSListElie(shared_ptr<SessionData> pSessionData,
                                            Date issueDate,
                                            double dSpread1,
                                            double dSpread2, 
                                            double dPrice1,
                                            double dPrice2);

shared_ptr<SessionData> MakeSessionDataElie(Date valuationDate,
                                   double dSpot);

void HRWithTimeComponent();  
void VolFlatHRWithTimeComponent(); 
void TestImpliedCDSSpreads();   
void HRWithSpotComponentPower();
void VolFlatHRWithSpotComponentPower(); 
void VolFlatHRTimeComponentWithOptions();
void VolPower();    
void VolPowerHRPower();  
void HRWithTimeComponentEDS();   

void VolTanh(); // Calibrate vol tanh with 2 options

void VolWithTimeComponent();

void VolTimeonlyHRWithTimeComponent(); // general list, but newton

void GeneralCalibration(); // general list using calibrator general

//used to store information to
//run the tests suggested by Elie.
class Parameters
{

public:
 Date m_MaturityDate; 
 Date m_ValuationDate;
 Date m_IssueDate;
 double m_dStrike;           
 OptionType m_optionType; 
 ExerciseType m_exerciseType;
 double m_dOptionPrice;
 double m_dSpot;
 double m_dSpread1;
 double m_dSpread2; 
 double m_dCDSPrice1;
 double m_dCDSPrice2;
};

void RunTestElie(Parameters param);

void testElie1();
void testElie2();
void testElie3();
void testElie4();



int main()
{
  try {

    //HRWithTimeComponent();  

    //VolFlatHRWithTimeComponent();

    //TestImpliedCDSSpreads();

    //HRWithSpotComponentPower();  
    
    //VolFlatHRWithSpotComponentPower(); 
    
    //VolFlatHRTimeComponentWithOptions();

    //VolPower();   

    //VolPowerHRPower();

    //HRWithTimeComponentEDS();   

    //VolTanh();

    //VolWithTimeComponent();

    VolTimeonlyHRWithTimeComponent();

    //GeneralCalibration();

   //---------------------------------------------------------------------//
   // tests suggested by Elie
   //---------------------------------------------------------------------//

   //testElie1();
   //testElie2();
   //testElie3();
   //testElie4();


  }
  catch(const ito33::Exception& e)
  {
    cout << e.GetErrorMessage() << endl;
  }

  return 0;
}


void MakeCDSList(TermStructureCDS& tsCDS)
{
  // tsCDS.clear();

  size_t nNbCDS = 5, nIdx;
  
  for (nIdx = 0; nIdx < nNbCDS; nIdx++)
  {
    Date IssueDate(2003, ito33::Date::Dec, 15);
    Date FirstDate(2004, ito33::Date::Jan, 1);
    
    Date MaturityDate = FirstDate;
    MaturityDate.AddMonths( (nIdx + 2) * 6);

    shared_ptr<CashFlowStreamUniform> 
      pSpreadStream( new CashFlowStreamUniform
                         (
                           IssueDate,
                           FirstDate, 
                           MaturityDate,
                           0.2,
                           Date::DayCountConvention_Act365,
                           Frequency_BiMonthly
                         )
                   );

    shared_ptr<CDS> pCDS(new CDS(0.4, pSpreadStream) );

    tsCDS.Add(pCDS);
  } 
}

void MakeEDSList(TermStructureEDS& tsEDS)
{

  size_t nNbEDS = 5, nIdx;
  
  for (nIdx = 0; nIdx < nNbEDS; nIdx++)
  {
    Date IssueDate(2003, ito33::Date::Dec, 15);
    Date FirstDate(2004, ito33::Date::Jan, 1);
    
    Date MaturityDate = FirstDate;
    MaturityDate.AddMonths( (nIdx + 2) * 6);

    shared_ptr<CashFlowStreamUniform> 
      pSpreadStream( new CashFlowStreamUniform
                         (
                           IssueDate,
                           FirstDate, 
                           MaturityDate,
                           0.2,
                           Date::DayCountConvention_Act365,
                           Frequency_BiMonthly
                         )
                   );

    double dRecovery = 0.0;
    double dBarrier = 30.0;
    shared_ptr<EDS> pEDS(new EDS(dRecovery, pSpreadStream, dBarrier) );

    tsEDS.Add(pEDS);
  } 
}

void MakeOptionList(TermStructureOption& tsOption)
{
  size_t nNb = 5, nIdx;
  
  Date valuationDate("2004/2/1");
  double dStrike = 90;

  for (nIdx = 0; nIdx < nNb; nIdx++)
  {
   
    ito33::Date maturityDate = valuationDate;
    maturityDate.AddMonths( (nIdx + 2) * 6);

    shared_ptr<Option> 
      pOption(new Option(dStrike, maturityDate, 
                         Option_Put, ExerciseType_American) );

    tsOption.Add(pOption);
  } 
}

shared_ptr<SessionData> InitSessionData(Date valuationDate)
{
  // Setup the equity, and attach to session data
  //shared_ptr<Equity> pEquity( new Equity(33.5) );
  //shared_ptr<Equity> pEquity( new Equity(50.0) );
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity( new Equity(100.0, pCurrency) );
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.02) ) );
  
  shared_ptr<YieldCurve> pYC(new YieldCurveFlat(0.04));
 
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pYC);

  return 
    shared_ptr<SessionData>(new SessionData(pRateData, pEquity, valuationDate));
}

void WriteToXML(std::string sFileName, 
                shared_ptr<SessionData> pSessionData,
                const finance::TermStructureCDS& tsCDS, 
                const finance::Option& option,
                ihg::ParametrizationVolFlatHRWithSpotComponentPower& parametrization)
{

  std::ofstream ofs(sFileName.c_str());
  ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
  tagRoot.precision(10);
  tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

  // Need the braces so that the tags go out of scope.  This closes the
  // xml tags at the right spot
  {
    pSessionData->Dump(tagRoot);
  }

  //{
  //ito33::XML::Tag tagModel(XML_TAG_IHG_MODEL, tagRoot);
  //model.GetVolatility()->Dump(tagModel);
  //model.GetHazardRate()->Dump(tagModel);
  //}

  {
    ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
    parametrization.Dump(tagParam);
    tsCDS.Dump(tagParam);
    option.Dump(tagParam);
  }

}


void HRWithTimeComponent()
{ 

  std::cout << std::endl;
  std::cout << "ParametrizationHRWithTimeComponent" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.2)) );

  model.SetHazardRate( shared_ptr<ihg::HazardRate>
                       ( new ihg::HazardRateCallBack(parametricHR, 0) ) );
 
  // make CDS term structure
  TermStructureCDS tsCDS;
  MakeCDSList(tsCDS);

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
   
  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    CDSLike& cds = *(*ppCDS);

    cds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);

    // Take the computed prices as market price for testing
    cds.SetMarketPrice( output->GetPrice() );
  }

  ihg::ParametrizationHRWithTimeComponent parametrization;

  parametrization.SetSpotComponent
                  ( 
                    shared_ptr<ihg::SpotComponent>
                    ( new ihg::HRSpotComponentPower(0.4, 50) )
                  );

  parametrization.SetVolatility( model.GetVolatility() );

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.CalibrateWithCDSs(tsCDS);

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the hr in  the temporary model to the calibrated one
  model.SetHazardRate(pHR);

  // compared the prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }
}

void VolFlatHRWithTimeComponent()
{
  std::cout << std::endl;
  std::cout << "VolFlatHRWithTimeComponent" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.3)) );
  model.SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateCallBack(parametricHR, 0)) );
 
  // make a CDS term structure
  TermStructureCDS tsCDS;

  MakeCDSList(tsCDS);

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
   
  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    CDSLike& cds = *(*ppCDS);

    cds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);

    // Take the computed prices as market price for testing
    cds.SetMarketPrice( output->GetPrice() );
  }

  // Create an option 
  shared_ptr<Option> 
    pOption(new Option(45, Date("2005/02/01"), Option_Call, ExerciseType_European));

  // Price the option with the temporary model
  pOption->SetSessionData(pSessionData);

  shared_ptr<finance::ModelOutput> output = model.Compute(*pOption);

  // Take the computed prices as market price for testing
  pOption->SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolFlatHRWithTimeComponent parametrization;

  parametrization.SetSpotComponent
                  ( 
                    shared_ptr<ihg::SpotComponent>
                    ( new ihg::HRSpotComponentPower(0.4, 50) )
                  );

  parametrization.CalibrateWithOptionAndCDSs(*pOption, tsCDS);

  shared_ptr<ihg::VolatilityFlat>
    pVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.GetHazardRate();

  // Print out the calibrated flat vol
  std::cout << "the flat vol: " << pVol->GetValue() << std::endl;

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the vol and hr in  the temporary model to the calibrated one
  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);
  
  // Compare the cds prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }


  // Compare the option price
  output = model.Compute(*pOption);

  std::cout << "option" << std::endl;

  std::cout << "  Analytic (generated) price " << " = " 
            << pOption->GetMarketPrice() << std::endl;
  
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   
}


void VolFlatHRTimeComponentWithOptions()
{
  std::cout << std::endl;
  std::cout << "test vol flat hazard rate with time component with options" 
            << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.3)) );
  model.SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateFlat(0.02)) );
 
  // make an option term structure
  TermStructureOption tsOption;

  MakeOptionList(tsOption);

  const TermStructureOption::Elements& pElements = tsOption.GetAll();
  TermStructureOption::Elements::const_iterator ppOption;
   
  // Compute some test prices
  for ( ppOption = pElements.begin(); ppOption != pElements.end(); ++ppOption)
  {
    Option& option = *(*ppOption);

    option.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(option);

    // Take the computed prices as market price for testing
    option.SetMarketPrice( output->GetPrice() );
  }

  // Create an option 
  Option option(100, Date("2005/02/01"), Option_Call, ExerciseType_European);

  // Price the option with the temporary model
  option.SetSessionData(pSessionData);

  shared_ptr<finance::ModelOutput> output = model.Compute(option);

  // Take the computed prices as market price for testing
  option.SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolFlatHRWithTimeComponent parametrization;

  /*
  parametrization.SetSpotComponent
                  ( 
                    new ihg::HRSpotComponentPower(0.4, 50)
                  );
  */

  parametrization.CalibrateWithOptionAndOptions(option, tsOption);

  shared_ptr<ihg::VolatilityFlat>
    pVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.GetHazardRate();

  // Print out the calibrated flat vol
  std::cout << "the flat vol: " << pVol->GetValue() << std::endl;

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the vol and hr in  the temporary model to the calibrated one
  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);
  
  // Compare the cds prices
  for ( ppOption = pElements.begin(); ppOption != pElements.end(); ++ppOption)
  {
    const Option& option = *(*ppOption);

    shared_ptr<finance::ModelOutput> output = model.Compute(option);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << option.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }


  // Compare the option price
  output = model.Compute(option);

  std::cout << "option" << std::endl;

  std::cout << "  Analytic (generated) price " << " = " 
            << option.GetMarketPrice() << std::endl;
  
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   
}

void HRWithSpotComponentPower()
{
  std::cout << std::endl;
  std::cout << std::endl << "HRWithSpotComponentPower" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.4)) );
  model.SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateCallBack(parametricHR, 0)) );
  //model.SetHazardRate( new ihg::HazardRateFlat(0.4) );

  // make CDS term structure
  TermStructureCDS tsCDS;

  MakeCDSList(tsCDS);

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
   
  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    CDSLike& cds = *(*ppCDS);

    cds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);

    // Take the computed prices as market price for testing
    cds.SetMarketPrice( output->GetPrice() );
  }

  ihg::ParametrizationHRWithSpotComponentPower 
    parametrization(model.GetVolatility() );

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.CalibrateWithCDSs( tsCDS);

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues     = pHR->GetTimeComponentValues();

  std::cout << std::endl;

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the hr in  the temporary model to the calibrated one
  model.SetHazardRate(pHR);

  // compared the prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  } //end for


}//test4




void VolFlatHRWithSpotComponentPower()
{
  std::cout << std::endl;
  std::cout << "VolFlatHRWithSpotComponentPower" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.1)) );
  model.SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateCallBack(parametricHR, 0)) );
 
  // make a CDS term structure
  TermStructureCDS tsCDS;

  MakeCDSList(tsCDS);

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
   
  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    CDSLike& cds = *(*ppCDS);

    cds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);

    // Take the computed prices as market price for testing
    cds.SetMarketPrice( output->GetPrice() );
  }

  // Create an option 
  Option option(45, Date("2005/02/01"), Option_Call, ExerciseType_European);

  // Price the option with the temporary model
  option.SetSessionData(pSessionData);

  shared_ptr<finance::ModelOutput> output = model.Compute(option);

  // Take the computed prices as market price for testing
  option.SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolFlatHRWithSpotComponentPower parametrization;

  parametrization.CalibrateWithOptionAndCDSs(option, tsCDS);

  shared_ptr<ihg::VolatilityFlat>
    pVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.GetHazardRate();

  // Print out the calibrated flat vol
  std::cout << "the flat vol: " << pVol->GetValue() << std::endl;

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the vol and hr in  the temporary model to the calibrated one
  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);
  
  // Compare the cds prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }


  // Compare the option price
  output = model.Compute(option);

  std::cout << "option" << std::endl;

  std::cout << "  Analytic (generated) price " << " = " 
            << option.GetMarketPrice() << std::endl;
  
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   
} //test5


void VolPower()
{ 

  std::cout << std::endl;
  std::cout << "VolPower" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  double dS0 = pSessionData->GetSpotSharePrice();
  //model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityPower(0.4, 0.3, dS0)));
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityPower(0.35, 0.75, dS0)) );
  
  //model.SetHazardRate( new ihg::HazardRateCallBack(parametricHR, 0) );
  shared_ptr<ihg::HazardRate> pHR( new ihg::HazardRateFlat(0.2) );
  model.SetHazardRate( pHR );
 
  // make the two options for calibration
  Date maturityDate("2006/02/01");
  Option option1(dS0+5.0, maturityDate, Option_Call, ExerciseType_European);
  Option option2(dS0-5.0, maturityDate, Option_Call, ExerciseType_European);

  option1.SetSessionData(pSessionData);
  shared_ptr<finance::ModelOutput> output = model.Compute(option1);
  option1.SetMarketPrice( output->GetPrice() );

  option2.SetSessionData(pSessionData);
  output = model.Compute(option2);
  option2.SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolPower parametrization(pHR);

  parametrization.SetDebugOutputFile("testVolPower.xml");
  shared_ptr<ihg::VolatilityPower>
    pVol = parametrization.CalibrateWithOptions(option1, option2);

  // Output result
  std::cout << "alpha = " << pVol->GetAlpha()
            << ", beta = " << pVol->GetBeta()
            << ", dS0 = " << pVol->GetS0()
            << std::endl;

  // Set the vol in  the temporary model to the calibrated one
  model.SetVolatility(pVol);

  output = model.Compute(option1);
    
  std::cout << "  Analytic (generated) price " << " = " 
            << option1.GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   

  output = model.Compute(option2);
    
  std::cout << "  Analytic (generated) price " << " = " 
            << option2.GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   


}


void VolPowerHRPower()
{ 

  std::cout << std::endl;
  std::cout << "VolPowerHRPower" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  double dS0 = pSessionData->GetSpotSharePrice();
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityPower(0.25, -0.55, dS0)));

  shared_ptr<ihg::HazardRate> pHR( new ihg::HazardRatePower(0.1, 0.35, dS0) );
  model.SetHazardRate( pHR );
 
  // make the two options for calibration
  Date maturityDate("2006/02/01");
  shared_ptr<Option> 
    pOption1( new Option(dS0-2.0, maturityDate, Option_Call, ExerciseType_European) );
  
  shared_ptr<Option> 
    pOption2( new Option(dS0+10.0, maturityDate, Option_Call, ExerciseType_European) );

  // make the cds term structure for calibration
  TermStructureCDS tsCDS;
  MakeCDSList(tsCDS);

  // Compute some test prices
  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    CDSLike& cds = *(*ppCDS);

    cds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);

    // Take the computed prices as market price for testing
    cds.SetMarketPrice( output->GetPrice() );
  }


  pOption1->SetSessionData(pSessionData);
  shared_ptr<finance::ModelOutput> output = model.Compute( *pOption1 );
  pOption1->SetMarketPrice( output->GetPrice() );

  pOption2->SetSessionData(pSessionData);
  output = model.Compute( *pOption2);
  pOption2->SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolPowerHRPower parametrization;

  parametrization.SetDebugOutputFile("testVolPowerHRPower.xml");
  
  parametrization.CalibrateWithOptionsAndCDSs(*pOption1, *pOption2, tsCDS);
  
  shared_ptr<ihg::VolatilityPower> pCalVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateCombo> pCalHR = parametrization.GetHazardRate();

  shared_ptr<ihg::HRSpotComponentPower> 
    pCalHRSpot = parametrization.GetHRSpotComponentPower();

  // Output result
  std::cout << "vol alpha = " << pCalVol->GetAlpha()
            << ", vol beta = " << pCalVol->GetBeta()
            << ", hr beta = " << pCalHRSpot->GetBeta()
            << ", dS0 = " << pCalVol->GetS0()
            << std::endl;

  // Set the vol and hr in the temporary model to the calibrated ones
  model.SetVolatility(pCalVol);
  model.SetHazardRate(pCalHR);

  output = model.Compute( *pOption1);
  
  std::cout << "Option 1:"  << std::endl;
  std::cout << "  Analytic (generated) price " << " = " 
            << pOption1->GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   

  output = model.Compute( *pOption2);
    
  std::cout << "Option 2:"  << std::endl;
  std::cout << "  Analytic (generated) price " << " = " 
            << pOption2->GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   

  // Compare the cds prices
  std::cout << "CDS term structure:"  << std::endl;
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }

}


void HRWithTimeComponentEDS()
{ 

  std::cout << std::endl;
  std::cout << "HRWithTimeComponent using EDS" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);

  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.2)) );
  model.SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateCallBack(parametricHR, 0)) );
 
  // make EDS term structure
  TermStructureEDS tsEDS;
  MakeEDSList(tsEDS);

  const TermStructureEDS::Elements& pEDSElements = tsEDS.GetAll();
  TermStructureEDS::Elements::const_iterator ppEDS;
   
  // Compute some test prices
  for ( ppEDS = pEDSElements.begin(); ppEDS != pEDSElements.end(); ++ppEDS)
  {
    EDS& eds = *(*ppEDS);

    eds.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(eds);

    // Take the computed prices as market price for testing
    eds.SetMarketPrice( output->GetPrice() );
  }

  ihg::ParametrizationHRWithTimeComponent parametrization;

  parametrization.SetSpotComponent
     ( 
       shared_ptr<ihg::SpotComponent>
       (new ihg::HRSpotComponentPower(0.4, pSessionData->GetSpotSharePrice()))
     );

  parametrization.SetVolatility( model.GetVolatility() );

  parametrization.SetDebugOutputFile("testEDS.xml");

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.CalibrateWithEDSs(tsEDS);

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the hr in  the temporary model to the calibrated one
  model.SetHazardRate(pHR);

  // compared the prices
  for ( ppEDS = pEDSElements.begin(); ppEDS != pEDSElements.end(); ++ppEDS)
  {
    const EDS& eds = *(*ppEDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(eds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << eds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }
}


void VolTanh()
{ 
  std::cout << std::endl;
  std::cout << "test vol tanh" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  double dS0 = 100;

  pSessionData->GetEquity()->SetSpotSharePrice(dS0);

  // make the two options for calibration
  Date maturityDate("2005/02/01");

  Option option1(dS0 - 10.0, maturityDate, Option_Call, ExerciseType_European);
  option1.SetSessionData(pSessionData);

  Option option2(dS0 + 10.0, maturityDate, Option_Call, ExerciseType_European);
  option2.SetSessionData(pSessionData);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  
  shared_ptr<ihg::HazardRate> pHR( new ihg::HazardRateFlat(0.) );
  model.SetHazardRate( pHR );

  //model.SetVolatility( new ihg::VolatilityTanh(0.15, 0.35, dS0) ); 
  
  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.35)) );//
  
  shared_ptr<finance::ModelOutput> output = model.Compute(option1);
  option1.SetMarketPrice( output->GetPrice() );

  model.SetVolatility( shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(0.25)) );//

  output = model.Compute(option2);
  option2.SetMarketPrice( output->GetPrice() );

  ihg::ParametrizationVolTanh parametrization(pHR);

  // parametrization.SetDebugOutputFile("testVolTanh.xml");

  shared_ptr<ihg::VolatilityTanh>
    pVol = parametrization.CalibrateWithOptions(option1, option2);

  // Output result
  std::cout << "left = " << pVol->GetLeft()
            << ", right = " << pVol->GetRight()
            << ", dS0 = " << pVol->GetS0()
            << std::endl;

  // Set the vol in  the temporary model to the calibrated one
  model.SetVolatility(pVol);

  output = model.Compute(option1);
    
  std::cout << "  Analytic (generated) price " << " = " 
            << option1.GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   

  output = model.Compute(option2);
    
  std::cout << "  Analytic (generated) price " << " = " 
            << option2.GetMarketPrice() << std::endl;
    
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   
}

void VolWithTimeComponent()
{ 

  std::cout << std::endl;
  std::cout << "VolWithTimeComponent" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  double dS0 = pSessionData->GetSpotSharePrice();

  // make option term structure
  TermStructureOption tsOption;
  MakeOptionList(tsOption);

  // Generate prices using timeonly volatilities.  Should be able to
  // fit these prices exactly
  const TermStructureOption::Elements& pElements = tsOption.GetAll();
  size_t nNbElements = pElements.size();
  std::vector<Date> pDates(nNbElements);
  std::vector<double> pdValues(nNbElements);

  TermStructureOption::Elements::const_iterator ppElem;
  size_t nIdx;
  for ( ppElem = pElements.begin(), nIdx = 0; 
        ppElem != pElements.end(); 
        ++ppElem, nIdx++)
  {
    pDates[nIdx] = (*ppElem)->GetMaturityDate();
    pdValues[nIdx] = 0.1 + nIdx*0.04;
  }

  // Setup the model and generate prices
  ihg::TheoreticalModel model;

  shared_ptr<ihg::Volatility> 
    pVol( new ihg::VolatilityTimeOnly(pDates, pdValues));
  shared_ptr<ihg::HazardRate> pHR(new ihg::HazardRatePower(0.2, 0.3, dS0));

  model.SetVolatility( pVol );
  model.SetHazardRate( pHR );

  // Compute some test prices
  for ( ppElem = pElements.begin(); 
        ppElem != pElements.end(); 
        ++ppElem)
  {
    Option& option = *(*ppElem);

    option.SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> output = model.Compute(option);

    // Take the computed prices as market price for testing
    option.SetMarketPrice( output->GetPrice() );
  }


  // Setup the parametrization and calibrate
  ihg::ParametrizationVolWithTimeComponent parametrization(pHR);
  parametrization.SetDebugOutputFile("volwithtimecomponent.xml");

  shared_ptr<ihg::VolatilityWithTimeComponent>
    pCalVol = parametrization.CalibrateWithOptions(tsOption);

  // Print out the time component
  std::vector<Date> pMaturityDates = pCalVol->GetDates();
  std::vector<double> pdCalValues = pCalVol->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdCalValues[nIdx] << std::endl;

  // Set the hr in  the temporary model to the calibrated one
  model.SetVolatility(pCalVol);

  // compared the prices
  for ( ppElem = pElements.begin(); ppElem != pElements.end(); ++ppElem)
  {
    const Option& option = *(*ppElem);

    shared_ptr<finance::ModelOutput> output = model.Compute(option);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << option.GetMarketPrice() << std::endl;
    
    std::cout << "            Calibrated price " << " = " 
              << output->GetPrice() << std::endl;   
  }
}

//-----------------------------------------------------------------------------//

shared_ptr<SessionData> MakeSessionDataElie(Date valuationDate, double dSpot)
{
  // Setup the equity, and attach to Session data
  shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity( new Equity(dSpot, pCurrency) );

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>(new YieldCurveFlat(0.02)) );
 
  shared_ptr<YieldCurve> pYC(new YieldCurveFlat(0.02));
 
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pCurrency, pYC);

  return 
    shared_ptr<SessionData>(new SessionData(pRateData, pEquity, valuationDate));

} //MakeSessionDataElieconst Date valuationDate)

shared_ptr<TermStructureCDS> MakeCDSListElie(shared_ptr<SessionData> pSessionData,
                                            Date issueDate,
                                            double dSpread1,
                                            double dSpread2, 
                                            double dPrice1,
                                            double dPrice2)
{

  shared_ptr<TermStructureCDS> tsCDS(new TermStructureCDS);

  Date FirstDate(2005, ito33::Date::Jul, 3);

  Date MaturityDate1(2006,ito33::Date::Jul,3);
  Date MaturityDate2(2009,ito33::Date::Jul,3);

  double dRecoveryRate = 0;
 
  shared_ptr<CashFlowStreamUniform> 
    pSpreadStream1( new CashFlowStreamUniform
                        (
                          issueDate,
                          FirstDate, 
                          MaturityDate1,
                          dSpread1,
                          Date::DayCountConvention_Act365,
                          finance::Frequency_Annual
                        )
                  );

   shared_ptr<CashFlowStreamUniform> 
     pSpreadStream2( new CashFlowStreamUniform
                         (
                           issueDate,
                           FirstDate, 
                           MaturityDate2,
                           dSpread2,
                           Date::DayCountConvention_Act365,
                           finance::Frequency_Annual
                         )
                   );

    shared_ptr<CDS> pCDS1(new finance::CDS(dRecoveryRate, pSpreadStream1) );
    shared_ptr<CDS> pCDS2(new finance::CDS(dRecoveryRate, pSpreadStream2) );

    pCDS1->SetMarketPrice(dPrice1);
    pCDS2->SetMarketPrice(dPrice2);

    pCDS1->SetSessionData(pSessionData);
    pCDS2->SetSessionData(pSessionData);

    tsCDS->Add(pCDS1);
    tsCDS->Add(pCDS2);

    return tsCDS;

} // shared_ptr<TermStructureCDS> MakeCDSListElie


void RunTestElie(Parameters param)
{

  
  shared_ptr<Option> pOption(
    new finance::Option(
                        param.m_dStrike,
                        param.m_MaturityDate,
                        param.m_optionType,
                        param.m_exerciseType) );

  pOption->SetMarketPrice(param.m_dOptionPrice);

  shared_ptr<SessionData> pSessionData = 
                   MakeSessionDataElie
                        (
                        param.m_ValuationDate,
                        param.m_dSpot
                        );

  pOption->SetSessionData(pSessionData);

  
  shared_ptr<TermStructureCDS> tsCDS = 
                  MakeCDSListElie(
                                  pSessionData,
                                  param.m_IssueDate, 
                                  param.m_dSpread1,
                                  param.m_dSpread2, 
                                  param.m_dCDSPrice1,
                                  param.m_dCDSPrice2
                                  );

  ihg::ParametrizationVolFlatHRWithSpotComponentPower parametrization;

  parametrization.CalibrateWithOptionAndCDSs(*pOption, *tsCDS);

  shared_ptr<ihg::VolatilityFlat>
    pVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.GetHazardRate();

  // Print out the calibrated flat vol
  std::cout << "the flat vol: " << pVol->GetValue() << std::endl;

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // Set the vol and hr in  the temporary model to the calibrated one
  ihg::TheoreticalModel model;
  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);
  
  // Compare the cds prices
  const TermStructureCDS::Elements& pCDSElements = tsCDS->GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
  
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    shared_ptr<finance::ModelOutput> output = model.Compute(cds);
    
    std::cout << "  Analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "              Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }

   // Compare the option price
  shared_ptr<finance::ModelOutput> output = model.Compute(*pOption);

  std::cout << "option" << std::endl;

  std::cout << "  Analytic (generated) price " << " = " 
            << pOption->GetMarketPrice() << std::endl;
  
  std::cout << "              Computed price " << " = " 
            << output->GetPrice() << std::endl;   

  //std::string fileName("test.xml");
  //WriteToXML(fileName, pSessionData, *tsCDS, *pOption, parametrization);

}//RunTestElie(ElieParameter param)
/**
Option 	Strike 	Maturity	Price	Spot		
American Put	80	7/13/2005	3.29	100	

CDS
Maturity	First Coupon	Spread	Price
7/3/2006	7/3/2005	0.041	-0.0001
7/3/2009	7/3/2005	0.042	-0.0014

valuation date	12/1/2004
issue     date   12/1/2004

*/

void testElie1()
{
  std::cout << "test Elie 1 " << std::endl;
   
  Parameters param;

  param.m_MaturityDate    = Date(2005,Date::Jul,13);
  param.m_ValuationDate   = Date(2004,Date::Dec,1);
  param.m_IssueDate       = param.m_ValuationDate;
  param.m_dStrike         = 80.0;
  param.m_optionType      = finance::Option_Put;
  param.m_exerciseType    = ExerciseType_American;
  param.m_dOptionPrice    = 3.29;
  param.m_dSpot           = 100;
  param.m_dSpread1        = 0.041;
  param.m_dSpread2        = 0.042;
  param.m_dCDSPrice1      = -0.0001;
  param.m_dCDSPrice2      = -0.0014;

  RunTestElie(param);

} //void testElie1

/*
Option 	     Strike 	Maturity	Price	Spot		
American Call	65	   7/13/2005	1.31	50		
				
CDS
Maturity	First Coupon	Spread	 Price
7/3/2006	7/3/2005	    0.075	   -0.0007
7/3/2009	7/3/2005	   0.072	   0.0012

*/
void testElie2()
{

  std::cout << "test Elie 2 " << std::endl;

  Parameters param;

  param.m_MaturityDate    = Date(2005,Date::Jul,13);
  param.m_ValuationDate   = Date(2004,Date::Dec,1);
  param.m_IssueDate       = param.m_ValuationDate;
  param.m_dStrike         = 65.0;
  param.m_optionType      = finance::Option_Call;
  param.m_exerciseType    = ExerciseType_American;
  param.m_dOptionPrice    = 1.31;
  param.m_dSpot           = 50;
  param.m_dSpread1        = 0.075;
  param.m_dSpread2        = 0.072;
  param.m_dCDSPrice1      = -0.0007;
  param.m_dCDSPrice2      = 0.0012;

  RunTestElie(param);

} //testElie2()

/*
Option 	Strike 	Maturity	Price	Spot
American Call	100	7/13/2005	19.9750	100		

CDS
Maturity	First Coupon	Spread	Price
7/3/2006	7/3/2005	0.0485	0.0002
7/3/2009	7/3/2005	0.062	-0.0009

*/
void testElie3()
{
  std::cout << "test Elie 3 " << std::endl;

   Parameters param;

  param.m_MaturityDate    = Date(2005,Date::Jul,13);
  param.m_ValuationDate   = Date(2004,Date::Dec,1);
  param.m_IssueDate       = param.m_ValuationDate;
  param.m_dStrike         = 100.0;
  param.m_optionType      = finance::Option_Call;
  param.m_exerciseType    = ExerciseType_American;
  param.m_dOptionPrice    = 19.9750;
  param.m_dSpot           = 100;
  param.m_dSpread1        = 0.0485;
  param.m_dSpread2        = 0.062;
  param.m_dCDSPrice1      = 0.0002;
  param.m_dCDSPrice2      = -0.0009;

  RunTestElie(param);
}

/*
Option 	Strike 	Maturity	Price	Spot
American Put	30	7/13/2005	3.5779	40	

CDS
Maturity	First Coupon	Spread	Price
7/3/2006	7/3/2005	0.1	0.0072	
7/3/2009	7/3/2005	0.115	0.0029

*/
void testElie4()
{
  std::cout << "test Elie 4 " << std::endl;

  Parameters param;

  param.m_MaturityDate    = Date(2005,Date::Jul,13);
  param.m_ValuationDate   = Date(2004,Date::Dec,1);
  param.m_IssueDate       = param.m_ValuationDate;
  param.m_dStrike         = 30.0;
  param.m_optionType      = finance::Option_Put;
  param.m_exerciseType    = ExerciseType_American;
  param.m_dOptionPrice    = 3.5779;
  param.m_dSpot           = 40;
  param.m_dSpread1       = 0.1;
  param.m_dSpread2       = 0.115;
  param.m_dCDSPrice1      = 0.0072;
  param.m_dCDSPrice2      = 0.0029;

  RunTestElie(param);
}


void VolTimeonlyHRWithTimeComponent()
{

  std::cout << std::endl;
  std::cout << "VolTimeonlyHRWithTimeComponent calibration" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  
  shared_ptr<ihg::Volatility> pVol(new ihg::VolatilityFlat(0.3));
  shared_ptr<ihg::HazardRate> pHR(new ihg::HazardRateFlat(0.02));

  model.SetVolatility( pVol );
  model.SetHazardRate( pHR );

  // make CDS term structure and option term structure
  TermStructureCDS tsCDS;
  MakeCDSList(tsCDS);

  TermStructureOption tsOption;
  MakeOptionList(tsOption);

  // Setup
  Derivatives derivatives;

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
 
  shared_ptr<finance::ModelOutput> output;

  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    shared_ptr<CDSLike> pCDS = *ppCDS;

    pCDS->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> pOutput = model.Compute(*pCDS);

    pCDS->SetMarketPrice( pOutput->GetPrice() );

    derivatives.AddWithWeight(pCDS, 1.0);
  }

  const TermStructureOption::Elements& pOptionElements = tsOption.GetAll();
  TermStructureOption::Elements::const_iterator ppOption;

  for (ppOption = pOptionElements.begin();
       ppOption != pOptionElements.end();
       ++ppOption)
  {
    shared_ptr<Option> pOption = *ppOption;

    pOption->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> pOutput = model.Compute(*pOption);

    double dPrice = pOutput->GetPrice();

    pOption->SetMarketPrice( dPrice );

    derivatives.AddWithWeight(pOption, 0.05);
  }


  // Now do the calibration
  ihg::ParametrizationVolTimeOnlyHRWithTimeComponent calibrator;
  calibrator.SetDebugOutputFile("voltimeonly.xml");

  try
  {
    calibrator.CalibrateWithDerivatives(derivatives);
  }
  catch(ito33::Exception e)
  {
    std::cout << "Exception caught:" << std::endl  
              << e.GetFullMessage() << std::endl;

    std::cout << "Calibration failed!!!!!! Prices below are best guesses only" 
              << std::endl << std::endl;
  }

  //shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel() );

  model.SetVolatility( calibrator.GetVolatility());
  model.SetHazardRate( calibrator.GetHazardRate());

  // compare the prices
  std::cout << "cds" << std::endl;

  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    output = model.Compute(cds);
    
    std::cout << "  CDS analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "                  Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }

  for (ppOption = pOptionElements.begin();
       ppOption != pOptionElements.end();
       ++ppOption)
  {
    Option& option = *(*ppOption);

    output = model.Compute(option);    

    std::cout << "  Option analytic (generated) price " << " = " 
              << option.GetMarketPrice() << std::endl;
  
    std::cout << "                     Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }
}

void GeneralCalibration()
{

  std::cout << "General calibration" << std::endl;

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");
  shared_ptr<SessionData> pSessionData = InitSessionData(valuationDate);
  
  // Setup temporary model to compute some test prices
  ihg::TheoreticalModel model;
  
  shared_ptr<ihg::Volatility> pVol(new ihg::VolatilityFlat(0.3));
  shared_ptr<ihg::HazardRate> pHR(new ihg::HazardRateFlat(0.02));

  model.SetVolatility( pVol );
  model.SetHazardRate( pHR );

  // make CDS term structure and option term structure
  TermStructureCDS tsCDS;
  MakeCDSList(tsCDS);

  TermStructureOption tsOption;
  MakeOptionList(tsOption);

  // Setup
  Derivatives derivatives;

  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;
 
  shared_ptr<finance::ModelOutput> output;

  // Compute some test prices
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    shared_ptr<CDSLike> pCDS = *ppCDS;

    pCDS->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> pOutput = model.Compute(*pCDS);

    pCDS->SetMarketPrice( pOutput->GetPrice() );

    derivatives.AddWithWeight(pCDS, 1.0);
  }

  const TermStructureOption::Elements& pOptionElements = tsOption.GetAll();
  TermStructureOption::Elements::const_iterator ppOption;

  for (ppOption = pOptionElements.begin();
       ppOption != pOptionElements.end();
       ++ppOption)
  {
    shared_ptr<Option> pOption = *ppOption;

    pOption->SetSessionData(pSessionData);

    shared_ptr<finance::ModelOutput> pOutput = model.Compute(*pOption);

    double dPrice = pOutput->GetPrice();

    pOption->SetMarketPrice( dPrice );

    derivatives.AddWithWeight(pOption, 0.05);
  }

  // Now do the calibration
  shared_ptr<ihg::Translator> pTranslator(new 
    ihg::Translator(ihg::VolType_power, ihg::SpotType_power, 
                    derivatives, pSessionData->GetSpotSharePrice()));

  finance::CalibrationProcess calibrationProcess;
  ihg::CalibratorGeneral calibrator(calibrationProcess);

  shared_ptr<ihg::TheoreticalModel> pModel;
  try
  {
    pModel = calibrator.Calibrate(pTranslator.get(), derivatives);
  }
  catch(ito33::Exception e)
  {
    std::cout << "Calibration failed!!!!!! Prices below are best guesses only" 
              << std::endl << std::endl;
  }

  // compare the prices
  std::cout << "cds" << std::endl;

  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    output = pModel->Compute(cds);
    
    std::cout << "  CDS analytic (generated) price " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << "                  Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }

  for (ppOption = pOptionElements.begin();
       ppOption != pOptionElements.end();
       ++ppOption)
  {
    Option& option = *(*ppOption);

    output = pModel->Compute(option);    

    std::cout << "  Option analytic (generated) price " << " = " 
              << option.GetMarketPrice() << std::endl;
  
    std::cout << "                     Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }
}

