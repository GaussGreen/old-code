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
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

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

    ///////////////////////////////////////////////////////////////////////////
    /// model part ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

    // volatility
    pModel->SetVolatility( shared_ptr<ihg::Volatility>
                           ( new VolatilityFlat(0.2) ) );
  
    // hazard rate
    static const Date dates[2] = {Date(2003, Date::Jul, 1),
                                  Date(2003, Date::Dec, 1)};
    static const double values[2] = { 0.02, .025};

    pModel->SetHazardRate( shared_ptr<ihg::HazardRate>
                            ( new HazardRateTimeOnly(dates, values, 2) ) );

    // computational flags
    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
    flags->SetAnalysisDate(Date(2004, Date::Jan, 1));

    opt->SetComputationalFlags(flags);

    // output file
    pModel->SetDebugOutputFile("ihg_option.xml");

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
