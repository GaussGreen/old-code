//
//                       UBS/DAIMLER 4.35% 2010
//

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/link.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/finance/bondlike/bondliketerms.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratepower.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;


int GeneralizedPEPSLikePricing()
{
  try
  {
    bool
      bComputeRho = true,
      bComputeVega = true,
      bComputeSurface = false;

    ///////////////////////////////////////////////////////////////////////////
    /// financial part ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    Date valuationDate(2005, Date::Mar, 1);

    // initialization of equity
    shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));
    shared_ptr<Equity> pEquity(new Equity(50, pNumeraire));
    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.005));    
    pEquity->SetBorrowCurve(pyf);

    // initialization of rate data
    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.03) );
    shared_ptr<RateData> pRateData( new RateData );
    pRateData->SetYieldCurve(pNumeraire, pyc);
    
    // create the session data
    shared_ptr<SessionData> 
      pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

    // generalized peps-like
    Date
      issueDate = Date(2005, Date::Feb, 18),
      maturityDate = Date(2010, Date::Feb, 18);

    double
      dIssuePrice = 1.,
      dNominal = 50,
      dRecoveryRate = 0.;

    shared_ptr<BondLikeTerms> 
      pblt( new BondLikeTerms(issueDate, dIssuePrice,
                              maturityDate, dNominal, dRecoveryRate) );

    Date firstDate(2006, ito33::Date::Feb, 18);
    shared_ptr<CashFlowStream>
      pInterests( new CashFlowStreamUniform
                      (
                        issueDate,
                        firstDate, 
                        maturityDate,
                        0.0435,
                        Date::DayCountConvention_Act365,
                        Frequency_Annual
                      )
                );

    pblt->SetCashDistribution(pInterests);

    double
      dConversionRatio = 1;
    double
      dLowerStrike = 45,
      dHigherStrike = 72.5;
    shared_ptr<finance::GeneralizedPEPSLike> 
      pPEPS( new finance::GeneralizedPEPSLike
                              ( pblt,
                                dConversionRatio,
                                dLowerStrike,
                                dConversionRatio,
                                dHigherStrike)
	   );

    pPEPS->EnableOptionalConversion();

    pPEPS->SetSessionData(pSessionData);
    ///////////////////////////////////////////////////////////////////////////
    /// model part ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);

    // volatility
    pModel->SetVolatility( shared_ptr<Volatility>( new VolatilityFlat(0.2) ) );
  
    // hazard rate
    {
      static const Date dates[2] = {Date(2003, Date::Jul, 1),
                                    Date(2003, Date::Dec, 1)};
      static const double values[2] = { 0.02, .025};

      pModel->SetHazardRate( shared_ptr<HazardRate>
                             ( new HazardRateTimeOnly(dates, values, 2) ) );
    }

    // computational flags
    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetComputeRho(bComputeRho);
    flags->SetComputeVega(bComputeVega);
    flags->SetComputeSurface(bComputeSurface);
    flags->SetAnalysisDate(Date(2004, Date::Jan, 1));

    pPEPS->SetComputationalFlags(flags);

    // output file
    pModel->SetDebugOutputFile("ihg_peps.xml");

    ///////////////////////////////////////////////////////////////////////////
    //// Compute //////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pPEPS);

    std::cout.precision(10);
    std::cout << "pepslike price = " << output->GetPrice() << std::endl;
    std::cout << "         delta = " << output->GetDelta() << std::endl;
    std::cout << "         theta = " << output->GetTheta() << std::endl;
    if(output->HasRho())
      std::cout << "           rho = " << output->GetRho() << std::endl;
    if(output->HasVega())
      std::cout << "          vega = " << output->GetVega() << std::endl;

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
