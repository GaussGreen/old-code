#include "ito33/beforestd.h"
#include <iostream>
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"
ITO33_FORCE_LINK_MODULE(HGPriceOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

int main()
{
  try
  {
    Date valuationDate(2003, Date::Feb, 1);

    shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

    shared_ptr<Equity> pEquity(new Equity(33.5, pCurrency));

    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.05));
    
    pEquity->SetBorrowCurve(pyf);
    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.01) );

    shared_ptr<RateData> pRateData( new RateData() );
    pRateData->SetYieldCurve( pCurrency, pyc );

    shared_ptr<SessionData> 
      pSessionData = new SessionData(pRateData, pEquity, valuationDate);

    shared_ptr<Option> opt(new Option(45, 
                                     Date(2005, Date::Jan, 1),
                                     Option_Put,
                                     ExerciseType_European)
                         );

    opt->SetSessionData(pSessionData);

    size_t nNbRegimes = 2; 

    std::vector<double> pdVols;
    pdVols.resize(nNbRegimes, 0.2);

    std::vector<double> pdDefaultIntensities;
    pdDefaultIntensities.resize(nNbRegimes, 0.2);
   
    shared_ptr<hg::UnderlyingProcess>
      pUnderlyingProcess(new hg::UnderlyingProcess
                            (nNbRegimes, pdVols, pdDefaultIntensities) );

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));
    
    shared_ptr<hg::ComputationalFlags> flags(new hg::ComputationalFlags);
    flags->SetComputeSensitivities(true);
    pModel->SetComputationalFlags(flags);

    std::ofstream sAmplitude("c:\\ito33\\output\\amplitude.txt");
    
    double dAmplitude = 1.7;
    std::vector<double> pdAmplitude;
    std::vector<double> pdPrices;
    std::vector<double> pdSens;

    for (size_t i = 0; dAmplitude < 1.8; i++)
    {
      dAmplitude += 1.e-3;

      hg::Jumps jumps;
      jumps.push_back(hg::Jump(0.2, dAmplitude));

      pUnderlyingProcess->SetJumps(0, 1, jumps);//

      jumps.clear();
      jumps.push_back(hg::Jump(0.2, -0.1));

      pUnderlyingProcess->SetJumps(1, 0, jumps);//

      shared_ptr<hg::ModelOutput> output = pModel->Compute(*opt);

      std::cout.precision(10);
      
      pdAmplitude.push_back(dAmplitude);
      pdPrices.push_back( output->GetPrice() );

      pdSens.push_back( output->GetSensitivities()[5] );

      std::cout << dAmplitude << std::endl;
    }

    for (size_t nIdx = 0; nIdx < pdAmplitude.size() - 1; nIdx++)
      sAmplitude << pdAmplitude[nIdx] << '\t' << (pdPrices[nIdx + 1] - pdPrices[nIdx]) * 1.e3 
                 << " " << pdSens[nIdx] << " " << pdPrices[nIdx] << '\n';
      //sAmplitude << pdAmplitude[nIdx] << '\t' << (pdPrices[nIdx + 1] - pdPrices[nIdx]) * 1.e3 << '\n';

    sAmplitude.close();

    //  ShowConvergence(*pModel, *opt, 3); 

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
