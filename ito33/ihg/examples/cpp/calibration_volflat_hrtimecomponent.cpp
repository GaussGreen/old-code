
// This example shows the calibration of price of an option and that of
// CDSs by flat volatility and hazard rate (time componenent)

#include <iostream>

#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_volflat_hrwithtimecomponent.h"

using namespace std; 
using namespace ito33;
using namespace ito33::finance;

int CalibrationVolFlatHRTimeComponent()
{
  try {

  // Setup the pricing machinery
  Date valuationDate("2004/02/01");

  // Setup the equity
  shared_ptr<Numeraire> pNumeraire( new Numeraire("EUR") );
  shared_ptr<Equity> pEquity( new Equity(100.0, pNumeraire) );
  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.02) ) );
  
  // Setup the rate data, and attach to the session data
  shared_ptr<RateData> pRateData( new RateData );
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.04) );
  pRateData->SetYieldCurve(pNumeraire, pyc);

  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );
 
  // make a CDS term structure
  TermStructureCDS tsCDS;

  size_t nNbCDS = 5, nIdx;

  double pdPrices[] = { 0.004890216031, 0.008064241158, 0.01144001564, 
                        0.01474396574, 0.01821804761 };

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
                           0.01,
                           Date::DayCountConvention_Act365,
                           Frequency_BiMonthly
                         )
                   );

    shared_ptr<CDS> pCDS(new CDS(0.2, pSpreadStream) );

    pCDS->SetSessionData(pSessionData);

    pCDS->SetMarketPrice(pdPrices[nIdx]);

    tsCDS.Add(pCDS);
  } 

  // Create an option 
  Option option(45, Date("2005/02/01"), Option_Call, ExerciseType_European);

  // Price the option with the temporary model
  option.SetSessionData(pSessionData);

  // Take the computed prices as market price for testing
  option.SetMarketPrice(55.65928949);

  ihg::ParametrizationVolFlatHRWithTimeComponent parametrization;

  // Uncomment the following code to user a hazard rate with a spot component
  /*
  parametrization.SetSpotComponent
                  ( new ihg::HRSpotComponentPower(0.4, 50) );
  */

  parametrization.CalibrateWithOptionAndCDSs(option, tsCDS);

  shared_ptr<ihg::VolatilityFlat>
    pVol = parametrization.GetVolatility();

  shared_ptr<ihg::HazardRateWithTimeComponent>
    pHR = parametrization.GetHazardRate();

  std::cout << std::endl;

  std::cout << "**************************************************\n";
  std::cout << "vol flat + Hr time component calibration\n";
  std::cout << "**************************************************\n";

  // Print out the calibrated flat vol
  std::cout << "the flat vol: " << std::endl;
  std::cout << pVol->GetValue() << std::endl;

  // Print out the time component
  std::vector<Date> pMaturityDates = pHR->GetDates();
  std::vector<double> pdValues = pHR->GetTimeComponentValues();

  std::cout << "Calibrated time component of the hazard rate" << std::endl;
  for (size_t nIdx = 0; nIdx < pdValues.size(); nIdx++)
    std::cout << pMaturityDates[nIdx] << '\t' << pdValues[nIdx] << std::endl;

  // construct a temporary model to check the calibrated prices
  ihg::TheoreticalModel model;

  model.SetVolatility(pVol);
  model.SetHazardRate(pHR);
  
  shared_ptr<finance::ModelOutput> output;

  // Compare the cds prices
  const TermStructureCDS::Elements& pCDSElements = tsCDS.GetAll();
  TermStructureCDS::Elements::const_iterator ppCDS;

  std::cout << "Comparison for CDSs\n";
  for ( ppCDS = pCDSElements.begin(); ppCDS != pCDSElements.end(); ++ppCDS)
  {
    const CDSLike& cds = *(*ppCDS);

    output = model.Compute(cds);
    
    std::cout << "  Market price  " << " = " 
              << cds.GetMarketPrice() << std::endl;
    
    std::cout << " Computed price " << " = " 
              << output->GetPrice() << std::endl;   
  }

  // Compare the option price
  output = model.Compute(option);

  std::cout << "Comparison for option\n";

  std::cout << "  Market price  " << " = " 
            << option.GetMarketPrice() << std::endl;
  
  std::cout << " Computed price " << " = " 
            << output->GetPrice() << std::endl;   
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n" << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught.\n";

    return 2;
  }

  return 0;
}
