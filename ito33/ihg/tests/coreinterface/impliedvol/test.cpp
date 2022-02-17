#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

ITO33_FORCE_LINK_MODULE(IHGPriceOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(ihg::TheoreticalModel);
}

int main()
{
  try
  {

    Date valuationDate(2003, Date::Feb, 1);

    shared_ptr<Numeraire> pCurrency(new Numeraire("EUR"));
    shared_ptr<Equity> pEquity(new Equity(33.5, pCurrency));

    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.01));
    
    pEquity->SetBorrowCurve(pyf);
    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.01) );
    shared_ptr<RateData> pRateData(new RateData);
    pRateData->SetYieldCurve(pCurrency, pyc);
    
    shared_ptr<SessionData> 
      pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

    shared_ptr<Option> opt(new Option(45, 
                                      Date(2005, Date::Jan, 1),
                                      Option_Call,
                                      ExerciseType_European
                                     )
                          );
    opt->SetSessionData(pSessionData);

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
    
    double dVol = .2;
    std::cout << "Volatility: " << dVol << std::endl;

    shared_ptr<ihg::VolatilityFlat> pVol( new ihg::VolatilityFlat( dVol ) ); 
    pModel->SetVolatility( pVol );
  
    shared_ptr<ihg::HazardRateFlat> pHR( new ihg::HazardRateFlat(.0) ); 
    pModel->SetHazardRate( pHR );

    pModel->SetDebugOutputFile("ihg.xml");

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*opt);
    std::cout.precision(10);
    std::cout << "Option price = " << output->GetPrice() << std::endl;
    std::cout << "       delta = " << output->GetDelta() << std::endl;
    std::cout << "       theta = " << output->GetTheta() << std::endl;
  
    // implied price
    opt->SetImpliedVol( dVol );

    std::cout << "Implied Price = " << opt->GetMarketPrice() << std::endl;

    //implied vol
    std::cout << "Implied Volatility = " 
              << opt->GetImpliedVolFrom( output->GetPrice() ) << std::endl;

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
