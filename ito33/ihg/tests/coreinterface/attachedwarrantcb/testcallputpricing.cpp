#include <iostream>
#include <fstream>

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/surfacedouble.h"
#include "ito33/finance/domain.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratetimeonly.h"

#include "ito33/tests/showconvergence.h"

#include "utils.h"

using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{

  void OutputSurface(const finance::SharedSurface& surface,
    std::string sFileName, double dSpot);

void TestCallPutPricing()
{

  std::cout.precision(15);
  std::cout << "Test pricing with call and put on reset date" << std::endl << std::endl;

  try 
  {    
    bool bComputeRho = true;
    bool bComputeVega = true;
    bool bComputeSurface = true;
    
    shared_ptr<finance::AttachedWarrantConvertibleBond>
      //pAttachedWarrantCB_1D = InitAmericanExpressTest();
      pAttachedWarrantCB_1D = InitCarnivalTest(); 

    bool bShiftCallDate = true;
    shared_ptr<finance::AttachedWarrantConvertibleBond> 
      //  pAttachedWarrantCB_2D = InitAmericanExpressTest(bShiftCallDate);
      pAttachedWarrantCB_2D = InitCarnivalTest(bShiftCallDate);

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);

    pAttachedWarrantCB_1D->SetComputationalFlags(flags);
    
    pAttachedWarrantCB_2D->SetComputationalFlags(flags);

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

    //pModel->SetDebugOutputFile("ihgcbwithwarrants.xml");

    double dVol = 0.3;
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

    Date valuationDate = 
      pAttachedWarrantCB_1D->GetSessionData()->GetValuationDate();
    
    Date maturityDate = pAttachedWarrantCB_1D->GetMaturityDate();
    
    size_t nNbHRTimes = 2;
    static const Date pTimes[2] = { valuationDate , maturityDate  };
    static const double pdValues[2] = { 0.2, .2 };

    pModel->SetHazardRate
      ( shared_ptr<HazardRate>(new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes)) );

    /*
      double dLambda = .2; 
      pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dLambda)) );
    */

    shared_ptr<finance::ModelOutput> 
      output_1d = pModel->Compute(*pAttachedWarrantCB_1D);

    std::cout.precision(15);
    std::cout << "The 1d price is: " << output_1d->GetPrice() << std::endl;
    std::cout << "The 1d delta is: " << output_1d->GetDelta() << std::endl;
    std::cout << "The 1d gamma is: " << output_1d->GetGamma() << std::endl;

    //double dS = pAttachedWarrantCB_2D->GetSessionData()->GetSpotSharePrice();

    //OutputSurface(output_1d->GetPriceSurface(), "1d_price_surface", dS);

    if ( output_1d->HasRho() )
    {
      std::cout << "The 1d Rho  is: " << output_1d->GetRho() << std::endl;
      //OutputSurface(output_1d->GetRhoSurface(), "1d_rho_surface", dS);
    }

    if ( output_1d->HasVega() )
    {
      std::cout << "The 1d Vega is: " << output_1d->GetVega() << std::endl;
      //OutputSurface(output_1d->GetVegaSurface(), "1d_vega_surface", dS);
    }

    //ShowConvergence(*pModel, *pAttachedWarrantCB_1D, 4);

    shared_ptr<finance::ModelOutput> 
      output_2d = pModel->Compute(*pAttachedWarrantCB_2D);

    std::cout << std::endl;
    std::cout << "The 2d price is: " << output_2d->GetPrice() << std::endl;
    std::cout << "The 2d delta is: " << output_2d->GetDelta() << std::endl;
    std::cout << "The 2d gamma is: " << output_2d->GetGamma() << std::endl;

   //OutputSurface(output_2d->GetPriceSurface(), "2d_price_surface", dS);

    if ( output_2d->HasRho() )
    {
      std::cout << "The 2d Rho  is: " << output_2d->GetRho() << std::endl;
      //OutputSurface(output_2d->GetRhoSurface(), "2d_rho_surface", dS);
    }

    if ( output_2d->HasVega() )
    {
      std::cout << "The 2d Vega is: " << output_2d->GetVega() << std::endl;
      //OutputSurface(output_2d->GetVegaSurface(), "2d_vega_surface", dS);
    }

    //ShowConvergence(*pModel, *pAttachedWarrantCB_2D, 4);

  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
              << e.GetFullMessage() << std::endl;

  }
}


void OutputSurface(const finance::SharedSurface& surface,
    std::string sFileName, double dSpot)
{

  std::ofstream sOut(sFileName.c_str(), std::ios::out);
  ASSERT_MSG(sOut, "Problem opening output file");

  double dSMin = .1*dSpot;
  double dStep = .01*dSpot;
  double dSMax = 2.*dSpot;

  double dS    = dSMin;
  std::vector<double> pdS;
  
  while( dS < dSMax + dStep )
  {
    pdS.push_back(dS);
    dS   += dStep;
  }

  shared_ptr<finance::Domain> domain = surface->GetDomain();
  domain->SetUnderlyingSharePrices(pdS);
  finance::Domain::Spots spots = domain->GetUnderlyingSharePrices();
  finance::Domain::Dates dates = domain->GetDates();

  finance::Domain::Dates::const_iterator iterDates;

  size_t nIdxTime = 0;
  for (iterDates = dates.begin();
       iterDates != dates.end();
       ++iterDates)
  {
    // Get the values at this time
    finance::SurfaceDouble::Doubles values = surface->GetValuesAt(nIdxTime);
    nIdxTime++;

    for (size_t nIdx = 0; nIdx < spots.size(); nIdx++)
    {

        sOut << (*iterDates)
              << " "
              << spots[nIdx]
              << " "
              << values[nIdx]
              << std::endl;

    }
    sOut << std::endl;
  }

}

} //end namespace ito33
