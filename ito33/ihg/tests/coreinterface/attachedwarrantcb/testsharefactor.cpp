#include <iostream>

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

void TestShareFactor()
{

  /*
    Basic idea: If the share factor is zero, it is just a normal convertible
                bond.  As the share factor increases, the price should
                increase.
  */

  std::cout.precision(15);
  std::cout << "Test the incremental share factor" << std::endl << std::endl;

  try {

    shared_ptr<SessionData> pSessionData = InitSessionData();  

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);    

    double dVol = 0.5;
    pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol)) );

    double dLambda = .2; 
    pModel->SetHazardRate( shared_ptr<HazardRate>(new HazardRateFlat(dLambda)) );

    shared_ptr<finance::ConvertibleBond> pCB = InitCB(pSessionData);

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);
    
    std::cout << "Normal convertible bond price = " << output->GetPrice() 
              << std::endl
              << std::endl;

    for (double dShareFactor = 0.0; dShareFactor <= 4.01; dShareFactor += 0.2)
    {
      // share factor cannot be zero, and it looks strange to have
      // factor of 0.1000001, 0.200001 etc
      double dRealShareFactor = dShareFactor;
      if (dRealShareFactor == 0.0)
        dRealShareFactor = 1.e-8;

      shared_ptr<finance::AttachedWarrantConvertibleBond> 
        pAttachedWarrantCB = InitAttachedWarrantCB(pSessionData, 
                                                   dRealShareFactor);

      shared_ptr<finance::ModelOutput> 
        output = pModel->Compute(*pAttachedWarrantCB);

      std::cout << "Share factor = " << dShareFactor 
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
