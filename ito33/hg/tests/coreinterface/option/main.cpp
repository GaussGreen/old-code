#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/domain.h"

#include "ito33/hg/theoreticalmodel.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_MODULE(HGPriceOption);

#include "ito33/timer.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;


int main()
{
  try
  {
    Date valuationDate(2003, Date::Feb, 1);

    shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

    double dSpot = 60.0;
    shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
    pEquity->SetBorrowCurve(pyf);

    shared_ptr<Dividends> pDividends( new Dividends );
    pDividends->AddYield(Date(2003, Date::Jul, 1), 0.02);
    pEquity->SetDividends(pDividends);
    
    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.0801) );

    shared_ptr<RateData> pRateData(new RateData);
    pRateData->SetYieldCurve( pCurrency, pyc);
    
    shared_ptr<SessionData> 
      pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

    shared_ptr<Option> opt(new Option(67.5, 
                                     Date(2004, Date::Jul, 1),
                                     Option_Put,
                                     ExerciseType_European)
                         );

    opt->SetSessionData(pSessionData);
    shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);
    
    pFlags->SetComputeSurface(true);

    pFlags->SetAnalysisDate(valuationDate);

    pFlags->ActivateAllSensitivities(true);

    pFlags->SetComputeRho(true);

    opt->SetComputationalFlags(pFlags);

    size_t nNbRegimes = 2;

    std::vector<double> pdVols;
    pdVols.resize(nNbRegimes, 0.1);

    std::vector<double> pdDefaultIntensities;
    pdDefaultIntensities.resize(nNbRegimes, 0.05);

    shared_ptr<hg::UnderlyingProcess>
      pUnderlyingProcess( new hg::UnderlyingProcess
                              (nNbRegimes, pdVols, pdDefaultIntensities) );

    hg::Jumps jumps;
 
    jumps.push_back(hg::Jump(0.1, -0.2));
    pUnderlyingProcess->SetJumps(0, 0, jumps);
    
    if (nNbRegimes > 1)
    {
      
      jumps.clear();
      jumps.push_back(hg::Jump(0.15, -0.25));
      pUnderlyingProcess->SetJumps(0, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.2, -0.5));
      pUnderlyingProcess->SetJumps(1, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.11, 0.35));
      pUnderlyingProcess->SetJumps(1, 0, jumps); 
      
    }

    if (nNbRegimes > 2)
    {
      jumps.clear();
      jumps.push_back(hg::Jump(0.1, -0.2));
      pUnderlyingProcess->SetJumps(0, 2, jumps); 
 
      jumps.clear();
      jumps.push_back(hg::Jump(0.1, -0.3));
      pUnderlyingProcess->SetJumps(1, 2, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.3, -0.2));
      pUnderlyingProcess->SetJumps(2, 0, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.2, -0.4));
      pUnderlyingProcess->SetJumps(2, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.05, -0.25));
      pUnderlyingProcess->SetJumps(2, 2, jumps); 
    }

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));

    // pModel->SetDebugOutputFile("hg.xml");

    StopWatch sw;

    Time::GetCurrentTicks();
    
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*opt);
    
    std::cout << "Running time: " << sw() << std::endl;

    std::cout.precision(10);
    std::cout << "Option price = " << output->GetPrice() << std::endl;    
    std::cout << "       delta = " << output->GetDelta() << std::endl;    
    std::cout << "       gamma = " << output->GetGamma() << std::endl;    
    std::cout << "       theta = " << output->GetTheta() << std::endl;    

    if (output->HasRho())
      std::cout << "         rho = " << output->GetRho() << std::endl;    

    std::cout << std::endl;
    std::cout << "Total time (ms): " << sw() << std::endl;

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
        sensitivityByFD( hg::ComputeSensitivity(*pModel, *opt, dShift) );

      std::cout << "sensitivity by finite difference" << std::endl;

      for (size_t nIdx = 0; nIdx < sensitivityByFD.size(); nIdx++)
        std::cout << "index = " << nIdx << ", value = "
                  << sensitivityByFD[nIdx] << std::endl;
      std::cout << std::endl;
    }

    bool bVerbose = false;
    if ( bVerbose && output->HasRhoAtAnalysisDate() )
    {
      std::cout << "Rhos at analysis date" << std::endl;

      std::vector<double> pdSpots = output->GetSpotsAtAnalysisDate();
      std::vector<double> pdRhos = output->GetRhosAtAnalysisDate();

      size_t nNbS = pdSpots.size();
      for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
      {
        std::cout << pdSpots[nIdx] << " " << pdRhos[nIdx] << std::endl;
      }
      std::cout << std::endl;
    }

   if ( bVerbose && output->HasPriceSurface() )
   {
     std::cout << "surfaces" << std::endl;

     finance::Domain::Spots pdSpots(20);
     pdSpots[0] = dSpot/2.0;
     double dDiff = (dSpot - pdSpots[0]) / 10.0;
     for (size_t nIdxS = 1; nIdxS < 20; nIdxS++)
       pdSpots[nIdxS] = pdSpots[nIdxS - 1] + dDiff;

     finance::SharedSurface priceSurface = output->GetPriceSurface();
     finance::SharedSurface deltaSurface = output->GetDeltaSurface();
     finance::SharedSurface gammaSurface = output->GetGammaSurface();
     
     finance::SharedSurface rhoSurface;
     if ( output->HasRhoSurface() )
       rhoSurface = output->GetRhoSurface();

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

       finance::SurfaceDouble::Doubles pdRhos;
       if ( rhoSurface )
         pdRhos = rhoSurface->GetValuesAt(nIdxDate);

       for (size_t nIdx = 0; nIdx < pdSpots.size(); nIdx++)
       {
         std::cout << nIdxDate 
                   << " " << pdSpots[nIdx] 
                   << " " << pdPrices[nIdx] 
                   << " " << pdDeltas[nIdx] 
                   << " " << pdGammas[nIdx];

         if ( rhoSurface )
           std::cout << " " << pdRhos[nIdx];

         std::cout << std::endl;
       }
       std::cout << std::endl;

     } // loop over dates
   }

   
    // ShowConvergence(*pModel, *opt, 3);

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
