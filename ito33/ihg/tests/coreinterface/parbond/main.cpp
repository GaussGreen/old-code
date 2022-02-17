#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/parbond.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/dividends.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceParBond);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

int main()
{
 try
 {   
   bool
     bComputeRho = true,
     bComputeVega = true,
     bComputeSurface = true;

   Date valuationDate(2002, Date::Jul, 1);

   shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
   shared_ptr<Equity> pEquity(new Equity(45, pCurrency));

   shared_ptr<Dividends> pDividends = shared_ptr<Dividends> (new Dividends());
   pEquity->SetDividends(pDividends);

   shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.00));
     
   pEquity->SetBorrowCurve(pyf);
     
   double dRate = 0.05;
   shared_ptr<YieldCurve> pyc(new YieldCurveFlat(dRate));   
   shared_ptr<RateData> pRateData(new RateData);
   pRateData->SetYieldCurve(pCurrency, pyc);

   shared_ptr<SessionData> 
     pSessionData(new SessionData(pRateData, pEquity, valuationDate));

   shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

   pModel->SetDebugOutputFile("ihgparbond.xml");

   double dVol = 0.5;
   pModel->SetVolatility( shared_ptr<Volatility>(new VolatilityFlat(dVol) ));

//    size_t nNbHRTimes = 2;
//    static const Date pTimes[2] = { 38000, 40000};
//    static const double pdValues[2] = { 0.02, .02 };

   //pModel->SetHazardRate
   // ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );
   double dLambda = .2;
  // pModel->SetHazardRate( shared_ptr<HazardRate>(new ihg::HazardRateFlat(dLambda) ));

   pModel->SetHazardRate( shared_ptr<HazardRate>(new ihg::HazardRatePower(dLambda, 0.01, 100)));

   double dRecoveryRate = 0.5;

   Date issueDate(2001, Date::Feb, 1);

   finance::ParBond parbond(pSessionData->GetValuationDate(),
                            24,
                            0.03,
                            0.005,
                            Frequency_Quarterly,
                            Date::DayCountConvention_30360,
                            dRecoveryRate);

   parbond.SetSessionData(pSessionData);

   shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
   flags->SetComputeRho(bComputeRho);
   flags->SetComputeVega(bComputeVega);
   flags->SetComputeSurface(bComputeSurface);

   parbond.SetComputationalFlags(flags);

   shared_ptr<finance::ModelOutput> output = pModel->Compute(parbond);

   std::cout.precision(10);
    
   std::cout << parbond.GetCouponRate() << std::endl;
   std::cout << "The parbond price is: " << output->GetPrice() << std::endl;
   std::cout << "The parbond delta is: " << output->GetDelta() << std::endl;
   if ( output->HasRho() )
     std::cout << "The parbond Rho is: " << output->GetRho() << std::endl;
   if ( output->HasVega() )
     std::cout << "The parbond Vega is: " << output->GetVega() << std::endl;
    
 //  ShowConvergence(*pModel, parbond);

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

