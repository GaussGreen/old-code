#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "ihg/tests/testutils.h"

#include "ihg/numoutput.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

int TestPricing();

int TestMeshes();

shared_ptr<finance::Option> MakeOption(double dStrike, 
                                      Date maturityDate, 
                                      size_t nWhich)
{
  
  shared_ptr<Option> opt;
  if (nWhich == 0)
    opt = make_ptr( new Option(dStrike, maturityDate, 
                               Option_Call, ExerciseType_European) );
  else if (nWhich == 1)
    opt = make_ptr( new Option(dStrike, maturityDate,
                               Option_Put, ExerciseType_European) );
                                     
  return opt;
}

int TestMeshes()
{
  try
  {

    // Make a series of test options.  Do convergence testing on each.
    // The code can then be re-run with different mesh generation code
    // to see which one worked better on average

    std::vector<double> pdSpots(3);
    pdSpots[0] = 1.0;
    pdSpots[1] = 100.0;
    pdSpots[2] = 10000.0;

    std::vector<double> pdStrikeRatios(6);
    pdStrikeRatios[0] = 0.3;
    pdStrikeRatios[1] = 0.8;
    pdStrikeRatios[2] = 0.9;
    pdStrikeRatios[3] = 1.0;
    pdStrikeRatios[4] = 1.1;
    pdStrikeRatios[5] = 3.0;

    // number of months to add to valuation date for the maturity
    std::vector<size_t> pnMonths(3);
    pnMonths[0] = 6;
    pnMonths[1] = 12;
    pnMonths[2] = 60;

    for (size_t nSpotIdx = 0; nSpotIdx < pdSpots.size(); nSpotIdx++)
    {
      double dSpot = pdSpots[nSpotIdx];
      shared_ptr<finance::SessionData> pSessionData = MakeSessionData(dSpot);

      for (size_t nStrikeIdx=0; nStrikeIdx<pdStrikeRatios.size(); nStrikeIdx++)
      {
        double dStrike = pdStrikeRatios[nStrikeIdx] * pdSpots[nSpotIdx];

        for (size_t nMonthIdx = 0; nMonthIdx < pnMonths.size(); nMonthIdx++)
        {
          Date maturityDate = pSessionData->GetValuationDate();
          maturityDate.AddMonths( pnMonths[nMonthIdx] );

          // Create a pricing model
          shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
          shared_ptr<ihg::VolatilityFlat> pVol( new VolatilityFlat(0.2) ); 
          pModel->SetVolatility( pVol );
  
          shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(0.05) );
          pModel->SetHazardRate( pHR );

          // call 
          std::cout << "Spot = " << dSpot
                    << ", strike = " << dStrike
                    << ", maturity = " << maturityDate
                    << ", type = call " 
                    << std::endl;
          shared_ptr<finance::Option> pOption = MakeOption(dStrike, maturityDate, 0);
          pOption->SetSessionData(pSessionData);
          ShowConvergence(*pModel, *pOption, 5);
          std::cout << std::endl << std::endl;

          // put
          std::cout << "Spot = " << dSpot
                    << ", strike = " << dStrike
                    << ", maturity = " << maturityDate
                    << ", type = put " 
                    << std::endl;
          pOption = MakeOption(dStrike, maturityDate, 1);
          pOption->SetSessionData(pSessionData);
          ShowConvergence(*pModel, *pOption, 5);
          std::cout << std::endl << std::endl;

        }
      }
    }



    return 0;
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
}

int TestPricing()
{
  try
  {
    bool
      bComputeRho = true,
      bComputeVega = true,
      bComputeSurface = true;

    shared_ptr<SessionData> pSessionData = MakeSessionData(33.5);

    shared_ptr<Option> opt(new Option(45, 
                                      Date(2005, Date::Jan, 1),
                                      Option_Call,
                                      ExerciseType_European
                                     )
                          );
    opt->SetSessionData(pSessionData);

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);

    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
//    flags->SetAnalysisDate(Date(2004, Date::Jan, 1));

    opt->SetComputationalFlags(flags);    

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
    shared_ptr<ihg::VolatilityFlat> pVol( new VolatilityFlat(0.2) ); 
    pModel->SetVolatility( pVol );
  
    static const Date times[2] = {38000, 40000};
    static const double values[2] = { 0.0201, .0251};

    shared_ptr<ihg::HazardRateTimeOnly> 
      pHR( new ihg::HazardRateTimeOnly(times, values, 2) ); 
    pModel->SetHazardRate( pHR );

    pModel->SetDebugOutputFile("ihg_option.xml");

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*opt);
    std::cout.precision(10);
    std::cout << "Option price = " << output->GetPrice() << std::endl;
    std::cout << "       delta = " << output->GetDelta() << std::endl;
    std::cout << "       gamma = " << output->GetGamma() << std::endl;
    std::cout << "       theta = " << output->GetTheta() << std::endl;
    if ( output->HasRho() )
      std::cout << "         rho = " << output->GetRho() << std::endl;
    if ( output->HasVega() )
      std::cout << "        vega = " << output->GetVega() << std::endl;

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


int TestRho()
{
  try
  {
    bool
      bComputeRho = true,
      bComputeVega = false,
      bComputeSurface = false;

    shared_ptr<Option> opt(new Option(45, Date(2005, Date::Jan, 1),
                                     Option_Call, ExerciseType_European) );

    Date valuationDate(2003, Date::Feb, 1);


    double dLeft = 0.25;
    //double dLeft = 0.376;
    double dRight = 0.55;
    //double dRight = 0.378;
    double dStep = 0.001;

    double dOldPrice = 1.0;

    for (double dRate = dLeft; dRate < dRight; dRate += dStep)
    {
      shared_ptr<SessionData> pSessionData = MakeSessionData(33.5);

      opt->SetSessionData(pSessionData);
      
      shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
      flags->SetComputeRho(bComputeRho);
      flags->SetComputeVega(bComputeVega);
      flags->SetComputeSurface(bComputeSurface);

      opt->SetComputationalFlags(flags);

      shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
      static const Date times[2] = {38000, 40000};
      static const double values[2] = { 0.0201, .0251};

      shared_ptr<ihg::HazardRateTimeOnly> 
        pHR( new ihg::HazardRateTimeOnly(times, values, 2) ); 
      pModel->SetHazardRate( pHR );

      shared_ptr<ihg::VolatilityFlat> pVol( new VolatilityFlat(0.2) ); 
      
      pModel->SetVolatility( pVol );

      shared_ptr<finance::ModelOutput> output = pModel->Compute(*opt);

      double dPrice = output->GetPrice();

      double dRhoFD = (dPrice - dOldPrice) / dStep;

      double dRho = output->GetRho();

      // Should probably plot these results
      if (dRate > dLeft)
      {
      std::cout << dRate 
                << " " << dPrice 
                << " " << dRho
                << " " << dRhoFD
                << " " << dRho - dRhoFD
                << std::endl;
      }
      dOldPrice = dPrice;

    }
 
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
  int iReturnCode = 0;
  
  iReturnCode += TestPricing();

 // iReturnCode += TestMeshes();

  //iReturnCode += TestRho();

  return iReturnCode;
}
