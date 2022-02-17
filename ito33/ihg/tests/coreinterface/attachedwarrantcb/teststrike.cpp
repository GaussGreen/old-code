#include <iostream>

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/convertiblebond.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "utils.h"

using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{

void TestStrike()
{

  /*
    Basic idea: As the strike increases, the price should decrease. An infinite
                strike should give the same price as for a normal convertible
                bond.  A zero strike should give the maximum price (but what
                this price is unknown, and the price could still be increased
                by increasing the share factor).
  */

  std::cout.precision(15);
  std::cout << "Test the strike" << std::endl << std::endl;

  try {
    bool bComputeRho = false;
    bool bComputeVega = false;
    bool bComputeSurface = false;
   
    shared_ptr<SessionData> pSessionData = InitSessionData();  

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);    

    double dVol = 0.5;
    pModel->SetVolatility( shared_ptr<ihg::Volatility>(new VolatilityFlat(dVol)) );

    double dLambda = .2; 
    pModel->SetHazardRate( shared_ptr<ihg::HazardRate>(new ihg::HazardRateFlat(dLambda)) );

    shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
    //flags->SetAnalysisDate( Date("2003/02/01"));

    pCB->SetComputationalFlags(flags);

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);
    
    std::cout << "Normal convertible bond price = " << output->GetPrice() 
              << std::endl
              << std::endl;

    for (double dStrike = 0.0; dStrike < 200.0; dStrike += 20.0)
    {

      shared_ptr<finance::AttachedWarrantConvertibleBond> 
        pAttachedWarrantCB = InitAttachedWarrantCB(pSessionData, 1.0, dStrike);

      shared_ptr<finance::ModelOutput> output = pModel->Compute(*pAttachedWarrantCB);

      std::cout << "Strike = " << dStrike 
                << ", price = " << output->GetPrice() 
                << std::endl;
    
    } // loop increasing the share factor
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
              << e.GetFullMessage() << std::endl;

  }

  std::cout << std::endl << std::endl;
}


} //end namespace ito33
