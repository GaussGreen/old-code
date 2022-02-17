#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/frequency.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"

#include "ito33/hg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"

#include "hg/computesensitivity.h"
#include "hg/numoutput.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceCB);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

shared_ptr<TermStructureCDS> 
MakeCDSList(const shared_ptr<SessionData>& pSessionData)
{
  shared_ptr<TermStructureCDS> tsCDS( new TermStructureCDS() );

  size_t nNbCDS = 5, nIdx;
  
  for (nIdx = 0; nIdx < nNbCDS; nIdx++)
  {
    Date IssueDate(2003, ito33::Date::Dec, 15);
    Date FirstDate(2004, ito33::Date::Jan, 1);
    
    Date MaturityDate = FirstDate;
    MaturityDate.AddMonths( (nIdx + 2) * 6);

    shared_ptr<CashFlowStreamUniform> 
      pSpreadStream( new CashFlowStreamUniform
                         (
                           IssueDate,
                           FirstDate, 
                           MaturityDate,
                           0.02,
                           Date::DayCountConvention_Act365,
                           Frequency_BiMonthly
                         )
                   );

    shared_ptr<CDS> pCDS(new CDS(0.4, pSpreadStream) );
    pCDS->SetSessionData(pSessionData);
    pCDS->SetMarketPrice(0);
    tsCDS->Add(pCDS);
  } 

  return tsCDS;
}
shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<SessionData>& pSessionData)
{
  Date
    issueDate = Date(2003, Date::Feb, 1),
    maturityDate = Date(2004, Date::May, 1);

  Date
    conversionStart = Date(2003, Date::Apr, 1),
    conversionEnd = Date(2004, Date::May, 1);

  double
    dIssuePrice = 1,
    dParValue = 110,
    dRedemptionRate = 1,
    dRecoveryRate = 0.;

  shared_ptr<BondTerms> 
    bc( new BondTerms(issueDate, dIssuePrice,
                      maturityDate, dParValue, dRedemptionRate,
                      dRecoveryRate) );

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  paymentDates.push_back(Date(2003, Date::May, 1));
  paymentDates.push_back(maturityDate);
  paymentRates.push_back(.02);
  paymentRates.push_back(.02);

  shared_ptr<CashFlowStream> 
    pInterests( new CashFlowStreamGeneral(issueDate, paymentDates, paymentRates,
    Date::DayCountConvention_Act365,Frequency_SemiAnnual) );
  
  bc->SetCashDistribution(pInterests);

  shared_ptr<ConversionSchedule> conv( new ConversionSchedule() );
  conv->SetKeepAccrued(true);

  shared_ptr<ConversionPeriod> 
    pConvPeriod( new ConversionPeriod(conversionStart, conversionEnd, 1) );
  conv->AddConversionPeriod(pConvPeriod);

  shared_ptr<finance::ConvertibleBond> 
    pConvertibleBond( new finance::ConvertibleBond(bc, conv) );

  Date startCallDate(issueDate);

  Date endCallDate(maturityDate);

  shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
  
  size_t nNoticePeriod = 20;

//    pCallSchedule->SetNoticePeriod(nNoticePeriod);
//    startCallDate.AddMonths(4);
//    endCallDate.AddMonths(-4);
        
  shared_ptr<CallPeriod> 
    pCallPeriod( CallPeriod::CreateWithStrike
                             (startCallDate, endCallDate, 1.) );

  pCallSchedule->AddCallPeriod( pCallPeriod);

        
  pConvertibleBond->SetCallSchedule(pCallSchedule);// 

  /*
  shared_ptr<finance::TermStructureCDS> tsCDS;

  pConvertibleBond->SetExchangeable(MakeCDSList(pSessionData), false);
  */

  pConvertibleBond->SetSessionData(pSessionData);
  
  return pConvertibleBond;
}


shared_ptr<SessionData> InitSessionData()
{
  // Setup the pricing machinery
  Date valuationDate("2003/02/01");

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

  // Setup the equity, and attach to session data
  shared_ptr<Equity> pEquity(new Equity(100, pCurrency));

  pEquity->SetBorrowCurve( shared_ptr<YieldCurve>( new YieldCurveFlat(0.0) ) );
 
  // Setup the issuer, and attach to the session data
  shared_ptr<RateData> pRateData( new RateData() );
  pRateData->SetYieldCurve(pCurrency, 
    shared_ptr<YieldCurve>( new YieldCurveFlat(0.04) ) );
  
  shared_ptr<SessionData> 
    pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

  return pSessionData;
}

//*************************************************************************

int main()
{
  try {
   
  bool
    bComputeRho = true,
    bComputeSurface = false;
   
  shared_ptr<SessionData> pSessionData = InitSessionData();  

  size_t nNbRegimes = 2; 

  std::vector<double> pdVols;
  pdVols.resize(nNbRegimes, 0.2);

  std::vector<double> pdDefaultIntensities;
  pdDefaultIntensities.resize(nNbRegimes, 0.02);
 
  shared_ptr<hg::UnderlyingProcess>
    pUnderlyingProcess(new hg::UnderlyingProcess
                           (nNbRegimes, pdVols, pdDefaultIntensities) );

  hg::Jumps jumps;
  jumps.push_back(hg::Jump(0.2, -0.8));

  pUnderlyingProcess->SetJumps(0, 1, jumps);

  jumps.clear();
  jumps.push_back(hg::Jump(0.2, -0.1));

  pUnderlyingProcess->SetJumps(1, 0, jumps);
  /**/

  shared_ptr<hg::TheoreticalModel> 
    pModel(new hg::TheoreticalModel(pUnderlyingProcess));
  
  pModel->SetDebugOutputFile("hgcb.xml");

  shared_ptr<ConvertibleBond> pCB = InitCB(pSessionData);
  
  shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);
  pFlags->SetComputeRho(bComputeRho);
  pFlags->ActivateAllSensitivities(true);
  pFlags->SetComputeSurface(bComputeSurface);

  // pFlags->SetAnalysisDate( Date("2003/02/01"));

  pCB->SetComputationalFlags(pFlags);

  shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);

  std::cout.precision(15);
  std::cout << "The cb price is: " << output->GetPrice() << std::endl;
  std::cout << "       delta is: " << output->GetDelta() << std::endl;
  std::cout << "       gamma is: " << output->GetGamma() << std::endl;
  std::cout << "       theta is: " << output->GetTheta() << std::endl;
  if (output->HasRho())
    std::cout << "         rho is: " << output->GetRho() << std::endl;

  // Output sensitivities, if computed
  shared_ptr<hg::NumOutput>
    pNumOutput( static_pointer_cast<hg::NumOutput>(output->GetNumOutput()) );
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
    double dShift = 1.e-6;
      
    std::vector<double> 
      sensitivityByFD( hg::ComputeSensitivity(*pModel, *pCB, dShift) );

    std::cout << "sensitivity by finite difference" << std::endl;

    for (size_t nIdx = 0; nIdx < sensitivityByFD.size(); nIdx++)
      std::cout << "index = " << nIdx << ", value = "
                << sensitivityByFD[nIdx] << std::endl;
    std::cout << std::endl;
  }

  //ShowConvergence(*pModel, *pCB);// 

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
