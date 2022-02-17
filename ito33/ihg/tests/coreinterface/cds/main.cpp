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
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/dividends.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCDS);

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

   shared_ptr<Dividends> pDividends( new Dividends() );
   pEquity->SetDividends(pDividends);

   shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
     
   pEquity->SetBorrowCurve(pyf);
     
   double dRate = 0.05;
   shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );
   
   shared_ptr<RateData> pRateData(new RateData);
   pRateData->SetYieldCurve(pCurrency, pyc);

   shared_ptr<SessionData> 
     pSessionData(new SessionData(pRateData, pEquity, valuationDate));

   shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

   pModel->SetDebugOutputFile("ihgcds.xml");

   double dVol = 0.5;
   shared_ptr<Volatility> pVol(new VolatilityFlat(dVol));
   pModel->SetVolatility(pVol);

   double dLambda = .2;
   shared_ptr<HazardRate> pHazardRate(new HazardRateFlat(dLambda) );

//    size_t nNbHRTimes = 2;
//    static const Date pTimes[2] = { 38000, 40000};
//    static const double pdValues[2] = { 0.02, .02 };

   //pModel->SetHazardRate
   // ( new HazardRateTimeOnly(&pTimes[0], &pdValues[0], nNbHRTimes) );

   
   pModel->SetHazardRate(pHazardRate);

  // pModel->SetHazardRate( new ihg::HazardRatePower(dLambda, 0.000001, 100));

   double dRecoveryRate = 0.5;
   double dSpread       = 0.25;

   Date issueDate(2001, Date::Feb, 1);
   Date firstSpreadDate(2003, Date::May, 1);
   Date lastSpreadDate(2003, Date::May, 1);

   shared_ptr<CashFlowStreamUniform>
     pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstSpreadDate,
                          lastSpreadDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                        )
                   );

   shared_ptr<finance::CDS>
     pCDS( new ito33::finance::CDS(dRecoveryRate, pSpreadStream) );

   pCDS->SetSessionData(pSessionData);

   shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
   flags->SetComputeRho(bComputeRho);
   flags->SetComputeVega(bComputeVega);
   flags->SetComputeSurface(bComputeSurface);

   pCDS->SetComputationalFlags(flags);

   shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCDS);

    std::cout.precision(10);
    std::cout << "The cds price is: " << output->GetPrice() << std::endl;
    if ( output->HasRho() )
      std::cout << "The cds Rho is: " << output->GetRho() << std::endl;
    if ( output->HasVega() )
      std::cout << "The cds Vega is: " << output->GetVega() << std::endl;

   ShowConvergence(*pModel, *pCDS);

   // this is just to test computing reference CDS works
   shared_ptr<finance::ReferenceCDS>
     pReferenceCDS(new finance::ReferenceCDS(48,
                                             Frequency_Quarterly,
                                             Date::DayCountConvention_30360,
                                             0.4));

   pReferenceCDS->SetSpread(0.02);
   pReferenceCDS->SetSessionData(pSessionData);
   output = pModel->Compute(*pReferenceCDS);

   std::cout << "The reference cds price is: " << output->GetPrice() << std::endl;

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
