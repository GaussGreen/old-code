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

void TestCap()
{

  /*
    Basic idea: As the cap increases, the price should increase. When the
                cap is the same as the initial conversion, the conversion
                ratio cannot change, and the price should be the same as
                for a normal CB.
  */

  std::cout.precision(15);
  std::cout << "Test the cap" << std::endl << std::endl;

  try {
    bool bComputeRho = false;
    bool bComputeVega = false;
    bool bComputeSurface = false;
   
    shared_ptr<SessionData> pSessionData = InitSessionData();  

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);    

    double dVol = 0.5;
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

    double dLambda = .2; 
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dLambda)) );

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

    for (double dCap = 1.0; dCap <= 1.21; dCap += 0.01)
    {

      shared_ptr<finance::AttachedWarrantConvertibleBond> 
        pAttachedWarrantCB = InitAttachedWarrantCB(pSessionData, 1.0, -1.0, 
                                                   dCap);

      shared_ptr<finance::ModelOutput> 
        output = pModel->Compute(*pAttachedWarrantCB);

      std::cout << "Cap = " << dCap 
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
