#include <iostream>

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/domain.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ito33/tests/showconvergence.h"

#include "utils.h"

using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{

void TestBasicPricing()
{

  std::cout.precision(15);
  std::cout << "Test basic pricing and convergence" << std::endl << std::endl;

  try {
    bool bComputeRho = false;
    bool bComputeVega = false;
    bool bComputeSurface = false;
   
    shared_ptr<SessionData> pSessionData = InitSessionData();  

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

    //pModel->SetDebugOutputFile("ihgcbwithwarrants.xml");

    double dVol = 0.5;
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

    double dLambda = .2; 
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dLambda)) );

    shared_ptr<finance::AttachedWarrantConvertibleBond> 
      //pAttachedWarrantCB = InitAttachedWarrantCB(pSessionData);
      //pAttachedWarrantCB = InitAutobacsTest();
      pAttachedWarrantCB = InitGettyTest();

    pSessionData = pAttachedWarrantCB->GetSessionData();

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
    flags->SetComputeSurface(true);
    flags->SetAnalysisDate( pSessionData->GetValuationDate() );
    //Date maturityDate = pAttachedWarrantCB->GetMaturityDate();
    //flags->SetAnalysisDate( maturityDate.AddDays(-1) );

    pAttachedWarrantCB->SetComputationalFlags(flags);

    shared_ptr<finance::BondLikeOutput> 
      output( static_pointer_cast<finance::BondLikeOutput> 
              ( pModel->Compute(*pAttachedWarrantCB) ) );

    std::cout.precision(15);
    std::cout << "The price is: " << output->GetPrice() << std::endl;
    std::cout << "The delta is: " << output->GetDelta() << std::endl;
    std::cout << "The gamma is: " << output->GetGamma() << std::endl;

    //get the bond floor
    std::cout << "The bond floor is: " << output->GetBondFloor() << std::endl;

    if (output->HasRho())
      std::cout << "The rho is: " << output->GetRho() << std::endl;
    if (output->HasVega())
      std::cout << "The vega is: " << output->GetVega() << std::endl;
    
    if ( output->HasPriceSurface() )
    {
      std::cout << std::endl;
      std::cout << "Surface computed" << std::endl;

      finance::SharedSurface pSurface = output->GetPriceSurface();

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

    //ShowConvergence(*pModel, *pAttachedWarrantCB, 4);
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
              << e.GetFullMessage() << std::endl;

  }

  std::cout << std::endl << std::endl;

}


} //end namespace ito33
