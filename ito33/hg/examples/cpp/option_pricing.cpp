#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/underlyingprocess.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

int OptionPricing()
{
  try
  {
    bool
      bComputeRho = true,
      bComputeVega = true,
      bComputeSurface = true;

    ///////////////////////////////////////////////////////////////////////////
    /// financial part ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    Date valuationDate(2003, Date::Feb, 1);

    // initialization of equity
    shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));
    shared_ptr<Equity> pEquity(new Equity(33.5, pNumeraire));
    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.005));    
    pEquity->SetBorrowCurve(pyf);

    // initialization of rate data
    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.03) );
    shared_ptr<RateData> pRateData( new RateData );
    pRateData->SetYieldCurve(pNumeraire, pyc);
    
    // create the session data
    shared_ptr<SessionData> 
      pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

    // create an european call
    shared_ptr<Option> opt(new Option(45, 
                                      Date(2005, Date::Jan, 1),
                                      Option_Call,
                                      ExerciseType_European
                                     )
                          );

    // associate session to the option
    opt->SetSessionData(pSessionData);

    // computational flags
    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
    flags->SetAnalysisDate(Date(2004, Date::Jan, 1));

    opt->SetComputationalFlags(flags);

    ///////////////////////////////////////////////////////////////////////////
    /// model part ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    size_t nNbRegimes = 2;

    std::vector<double> pdVols;
    pdVols.resize(nNbRegimes, 0.1);

    std::vector<double> pdDefaultIntensities;
    pdDefaultIntensities.resize(nNbRegimes, 0.05);

    shared_ptr<hg::UnderlyingProcess>
      pUnderlyingProcess( new hg::UnderlyingProcess
                              (nNbRegimes, pdVols, pdDefaultIntensities) );

    Jumps jumps;
    
    jumps.push_back(Jump(0.1, -0.2));
    pUnderlyingProcess->SetJumps(0, 0, jumps);
    
    if (nNbRegimes > 1)
    {
      jumps.clear();
      jumps.push_back(Jump(0.15, -0.25));
      pUnderlyingProcess->SetJumps(0, 1, jumps); 

      jumps.clear();
      jumps.push_back(Jump(0.2, -0.5));
      pUnderlyingProcess->SetJumps(1, 1, jumps); 

      jumps.clear();
      jumps.push_back(Jump(0.11, 0.35));
      pUnderlyingProcess->SetJumps(1, 0, jumps); 
    }

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));

    // output file
    pModel->SetDebugOutputFile("hg.xml");

    ///////////////////////////////////////////////////////////////////////////
    //// Compute //////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*opt);

    std::cout.precision(10);
    std::cout << "Option price = " << output->GetPrice() << std::endl;
    std::cout << "       delta = " << output->GetDelta() << std::endl;
    std::cout << "       theta = " << output->GetTheta() << std::endl;
    if(output->HasRho())
      std::cout << "         rho = " << output->GetRho() << std::endl;
    if(output->HasVega())
      std::cout << "        vega = " << output->GetVega() << std::endl;

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

