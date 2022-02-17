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
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/dividends.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/jumps.h"

#include "ito33/tests/showconvergence.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceCDS);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

int main()
{
 try
 {   
   Date valuationDate(2003, Date::Apr, 20);
   shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

   shared_ptr<Equity> pEquity(new Equity(45, pCurrency) );

   shared_ptr<Dividends> pDividends( new Dividends() );
   pEquity->SetDividends(pDividends);

   shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
     
   pEquity->SetBorrowCurve(pyf);
     
   double dRate = 0.04;

   shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

   shared_ptr<RateData> pRateData(new RateData);
   pRateData->SetYieldCurve( pCurrency, pyc );

   shared_ptr<SessionData> 
     pSessionData(new SessionData(pRateData, pEquity, valuationDate));

   size_t nNbRegimes = 2;

   double dVol = 0.3;
   std::vector<double> pdVols;
   pdVols.resize(nNbRegimes, dVol);

   double dIntensity = 0.02;
   std::vector<double> pdDefaultIntensities;
   pdDefaultIntensities.resize(nNbRegimes, dIntensity);
   pdDefaultIntensities[1] = 0.3;
  
   shared_ptr<hg::UnderlyingProcess> 
     pUnderlyingProcess(new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities));

   hg::Jumps jumps;
   jumps.push_back(hg::Jump(0.2, -0.8));

   pUnderlyingProcess->SetJumps(0, 1, jumps);

   jumps.clear();
   jumps.push_back(hg::Jump(0.3, -0.9));

   pUnderlyingProcess->SetJumps(1, 0, jumps); 

   shared_ptr<hg::TheoreticalModel> 
     pModel(new hg::TheoreticalModel(pUnderlyingProcess));
  


   double dRecoveryRate = 0.5;
   double dSpread       = 0.25;

   Date issueDate(2001, Date::Feb, 1);
   Date firstDate(2003, Date::May, 1);
   Date lastDate(2008, Date::May, 1);

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

   shared_ptr<finance::CDS>
     pCDS( new ito33::finance::CDS(dRecoveryRate, pSpreadStream) );

   pCDS->SetSessionData(pSessionData);
   
   shared_ptr<ComputationalFlags>
     pFlags(new ComputationalFlags);

   pFlags->SetComputeSurface(true);
   pFlags->SetComputeRho(true);

   pFlags->ActivateAllSensitivities(true);

   pCDS->SetComputationalFlags(pFlags);

   shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCDS);

   std::cout.precision(10);
   std::cout << "The cds price is: " << output->GetPrice() << std::endl;
   std::cout << "        delta is: " << output->GetDelta() << std::endl;
   std::cout << "        gamma is: " << output->GetGamma() << std::endl;
   std::cout << "        theta is: " << output->GetTheta() << std::endl;
      
   if ( output->HasRho() )
     std::cout << "         rho is: " << output->GetRho() << std::endl;

   shared_ptr<hg::NumOutput>
     pNumOutput( dynamic_pointer_cast<hg::NumOutput>( output->GetNumOutput() ) );
   
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
        sensitivityByFD( hg::ComputeSensitivity(*pModel, *pCDS, dShift) );

      std::cout << "sensitivity by finite difference" << std::endl;

      for (size_t nIdx = 0; nIdx < sensitivityByFD.size(); nIdx++)
        std::cout << "index = " << nIdx << ", value = "
                  << sensitivityByFD[nIdx] << std::endl;
    }

   //ShowConvergence(*pModel, *pCDS);

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
