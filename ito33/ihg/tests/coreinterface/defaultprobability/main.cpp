#include <iostream>
#include <cmath>

#include "ito33/link.h"
#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/vector.h"
#include "ito33/dateutils.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/yieldcurve_flat.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "ihg/tests/testutils.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

using namespace ito33;
using namespace ito33::finance;

void DeveloperTest()
{

  // Create the session
  double dSpot = 50.0;
  shared_ptr<finance::SessionData> pSessionData = ihg::MakeSessionData(dSpot);

  // Create the model
  shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

  shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat(0.2) ); 
  pModel->SetVolatility( pVol );

  double dHR = 0.05;
  shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(dHR) );
  pModel->SetHazardRate( pHR );

  // Create the dates at which the probabilities are desired
  std::vector<Date> pDates( 3, pSessionData->GetValuationDate() );
  pDates[0].AddDays(6);
  pDates[1].AddYears(1);
  pDates[2].AddYears(2);
  //for (size_t nIdxDate = 0; nIdxDate < pDates.size(); nIdxDate++)
  //  pDates[nIdxDate].AddYears(nIdxDate+1);

  // Do the computation
  std::vector<double> pdProbs = 
    pModel->ComputeCumulativeDefaultProbability(*pSessionData, pDates);

  // Report
  std::cout << "Valuation date = " << pSessionData->GetValuationDate() 
            << std::endl;
  std::cout << std::endl;

  std::cout << "Computed values" << std::endl;
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    std::cout << pDates[nIdx] << ", prob = " << pdProbs[nIdx] << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Analytic values" << std::endl;
  for (size_t nIdx = 0; nIdx < pdProbs.size(); nIdx++)
  {
    double dValuationTime = GetDoubleFrom(pSessionData->GetValuationDate());
    double dDateTime = GetDoubleFrom(pDates[nIdx]);
    double dTime = dDateTime - dValuationTime;
    
    std::cout << pDates[nIdx] << ", prob = " << 1.0 - exp(-dHR*dTime) 
              << std::endl;
  }
}

int main()
{
  try
  {
    // For manual testing
    DeveloperTest();

    return 0;
  }
  catch (const ito33::Exception& e)
  {
    std::cerr << "Exception caught:" << std::endl;
    std::cerr << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught." << std::endl;

    return 2;
  }

}
