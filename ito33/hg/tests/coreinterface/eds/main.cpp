#include "ito33/beforestd.h"
#include <iostream>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/eds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/domain.h"

#include "ito33/numeric/predicatedouble.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/jumps.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#ifdef ITO33_TEST_MODENV
#include "ito33/numeric/schemetype.h"
#include "ito33/numeric/numparams_reference.h"
#include "ito33/numeric/numparams_modifyreference.h"
#endif

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceEDS);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

int main()
{
 try
 {   
   Date valuationDate(2002, Date::Jul, 1);

   shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

   shared_ptr<Equity> pEquity(new Equity(31., pCurrency));

   shared_ptr<Dividends> pDividends( new Dividends() );
   pEquity->SetDividends(pDividends);

   shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
     
   pEquity->SetBorrowCurve(pyf);
     
   double dRate = 0.05;
   shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

   shared_ptr<RateData> pRateData(new RateData);
   pRateData->SetYieldCurve( pCurrency, pyc);

   shared_ptr<SessionData> 
     pSessionData(new SessionData(pRateData, pEquity, valuationDate));

   size_t nNbRegimes = 2;
   std::vector<double> pdVols(nNbRegimes);
   std::vector<double> pdDefaultIntensities(nNbRegimes);
  
   size_t nIdx;
   for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
   {
     pdVols[nIdx] = 0.2 + nIdx * 0.1;
     pdDefaultIntensities[nIdx] = 0.05;
   }

   shared_ptr<hg::UnderlyingProcess>
     pUnderlyingProcess(new hg::UnderlyingProcess
                            ( nNbRegimes, pdVols, pdDefaultIntensities ));

   if (nNbRegimes > 1)
   {
     std::vector<double> pdIntensity(1);
     std::vector<double> pdAmplitude(1);

     size_t nIdx1, nIdx2;
     for (nIdx1 = 0; nIdx1 < nNbRegimes; nIdx1++)
     {
       for (nIdx2 = 0; nIdx2 < nNbRegimes; nIdx2++)
       {
         std::vector<double> pdIntensity(1);
         std::vector<double> pdAmplitude(1);

         pdIntensity[0] = 0.1;
         pdAmplitude[0] = 0.0;

         if (nIdx1 != nIdx2)
           pUnderlyingProcess->SetJumps(nIdx1, nIdx2, pdIntensity, pdAmplitude);
       }
     }
   } // if more than one regime

   shared_ptr<hg::TheoreticalModel> 
     pModel( new hg::TheoreticalModel(pUnderlyingProcess) );

   pModel->SetDebugOutputFile("hg_eds.xml");    

   double dRecoveryRate = 0.5;
   double dSpread       = 0.25;

   Date issueDate(2001, Date::Feb, 1);
   Date firstDate(2003, Date::May, 1);
   Date lastDate(2003, Date::May, 1);

   shared_ptr<CashFlowStreamUniform>
     pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                        )
                  );

   double dBarrier = 30.0;
   shared_ptr<finance::EDS>
     pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

   pEDS->SetSessionData(pSessionData);
  
   shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);

   pFlags->ActivateAllSensitivities(true);

   //pFlags->SetSensitivityMethod(1);

   Date analysisDate = valuationDate;
   analysisDate.AddMonths(2);
   pFlags->SetAnalysisDate(analysisDate);

   //Date tmp(2003,Date::Apr, 28);
   //pFlags->SetAnalysisDate(tmp);

   pFlags->SetComputeSurface(true);

   pEDS->SetComputationalFlags(pFlags);

   shared_ptr<finance::ModelOutput> output = pModel->Compute(*pEDS);

   std::cout.precision(10);
   std::cout << "The EDS price is: " << output->GetPrice() << std::endl;
   
   shared_ptr<hg::NumOutput>
     pNumOutput( static_pointer_cast<hg::NumOutput>(output->GetNumOutput()) );

   // Output sensitivities, if computed
   if ( pNumOutput->HasSensitivities() )
   {
      std::cout << std::endl;
      std::cout << "Sensitivities:" << std::endl;
      std::vector<double> pdSensitivities = pNumOutput->GetSensitivities();
      for (size_t nIdx = 0; nIdx < pdSensitivities.size(); nIdx++)
        std::cout << "index = " << nIdx << ", value = " << pdSensitivities[nIdx] << std::endl;
   }
   std::cout << std::endl;

   // Compute derivatives by finite differences
   if ( pNumOutput->HasSensitivities() )
   {
     double dShift = 1.e-8;
      
     std::vector<double> 
       sensitivityByFD( hg::ComputeSensitivity(*pModel, *pEDS, dShift) );

     std::cout << "sensitivity by finite difference" << std::endl;

     for (size_t nIdx = 0; nIdx < sensitivityByFD.size(); nIdx++)
       std::cout << "index = " << nIdx << ", value = "
                 << sensitivityByFD[nIdx] << std::endl;
   }
   
   bool bTestExtraPoint = true;
   if ( bTestExtraPoint && output->HasPriceAtAnalysisDate() )
   {
     std::cout.precision(12);
     std::cout << "Analysis date data" << std::endl;
     finance::Values pdSpots = output->GetSpotsAtAnalysisDate();
     finance::Values pdPrices = output->GetPricesAtAnalysisDate();
     finance::Values pdDeltas = output->GetDeltasAtAnalysisDate();
     finance::Values pdGammas = output->GetGammasAtAnalysisDate();

     for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
       std::cout << pdSpots[nIdx] 
                 << " " << pdPrices[nIdx] 
                 << " " << pdDeltas[nIdx] 
                 << " " << pdGammas[nIdx] 
                 << std::endl;
     std::cout << std::endl;
   }
   

   if ( bTestExtraPoint && output->HasPriceSurface() )
   {
     std::cout << "Surface at valuation date" << std::endl;

     finance::Domain::Spots pdSpots(20);
     pdSpots[0] = dBarrier - 1.0;
     for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
       pdSpots[nIdxS] = pdSpots[nIdxS - 1] + 0.5;

     finance::SharedSurface priceSurface = output->GetPriceSurface();
     finance::SharedSurface deltaSurface = output->GetDeltaSurface();
     finance::SharedSurface gammaSurface = output->GetGammaSurface();

     priceSurface->GetDomain()->SetUnderlyingSharePrices(pdSpots);

     finance::Domain::Dates pdDates = priceSurface->GetDomain()->GetDates();

     for (size_t nIdxDate = 0; nIdxDate < pdDates.size(); nIdxDate++)
     {
       std::cout << "date = " << pdDates[nIdxDate] << std::endl;
       finance::SurfaceDouble::Doubles pdPrices 
         = priceSurface->GetValuesAt(nIdxDate);

       finance::SurfaceDouble::Doubles pdDeltas 
         = deltaSurface->GetValuesAt(nIdxDate);

       finance::SurfaceDouble::Doubles pdGammas 
         = gammaSurface->GetValuesAt(nIdxDate);

       for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
         std::cout << nIdxDate << " " << pdSpots[nIdx] 
                   << " " << pdPrices[nIdx] 
                //   << " " << pdDeltas[nIdx] 
                //   << " " << pdGammas[nIdx] 
                   << std::endl;
       
       std::cout << std::endl;
     }
   }

   // ShowConvergence(*pModel, *pEDS, 3);

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
