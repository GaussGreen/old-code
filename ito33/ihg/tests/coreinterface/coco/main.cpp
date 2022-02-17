
#include "ito33/beforestd.h"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"


#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/cocotype.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"

#include "utils.h"
#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(IHGPriceCB);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;
using namespace std;

//--------------------------------------------------------------------------
// Create the different windows of constraints
//--------------------------------------------------------------------------
Date issueDate     = Date(2003,Date::Feb, 1);
Date maturityDate  = Date(2004,Date::May, 1);

Date startCallDate = Date(2003,Date::Mar, 1);
Date endCallDate   = Date(2003,Date::Dec, 1);

//Date startConvDate = issueDate;
//Date endConvDate   = maturityDate;

Date startConvDate = Date(2003,Date::Mar, 1);
Date endConvDate   = Date(2003,Date::Dec, 1);



//--------------------------------------------------------------------------
// Risk parameters
//--------------------------------------------------------------------------
double dVol = 0.5;
shared_ptr<ihg::Volatility> pVolatility(new ihg::VolatilityFlat(dVol) );
double dLambda = .02;
shared_ptr<ihg::HazardRate> pHazardRate(new ihg::HazardRateFlat(dLambda) );



int main()
{

  try
  {
    shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);
    shared_ptr<finance::ConvertibleBond> 
      pCB = InitCB(pSessionData,issueDate,maturityDate);

    //shared_ptr<CallSchedule>       pCall = new CallSchedule();
    //shared_ptr<PutSchedule>        pPut  = new PutSchedule();
    shared_ptr<ConversionSchedule> pConv( new ConversionSchedule() );

    // Add coco period
    double dTriggerRate = 1.0;

    // normal window. No trigger. Issue date to maturity date
    //shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(issueDate, maturityDate, 1.0) );

    // normal window. No trigger
//    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );

/*
    // anytime/checkdate
    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );
    conversionPeriod->SetCoCo(dTriggerRate, CoCoType_CheckAnyTimeAndConvertOnCheckDate, 0.0, dTriggerRate);
*/
/*
    // anytime/as of check date
    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );
    conversionPeriod->SetCoCo(dTriggerRate, CoCoType_CheckAnyTimeAndConvertAsOfCheckDate, 0.0, dTriggerRate);
*/
/*
    // quarterly/as of check date
    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );
    conversionPeriod->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertAsOfCheckDate, 0.0, dTriggerRate);
*/

    // quarterly/next quarter
    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );
    conversionPeriod->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertDuringNextQuarter, 0.0, dTriggerRate);

/*
    // mixed window
    shared_ptr<ConversionPeriod> conversionPeriod(new ConversionPeriod(startConvDate, endConvDate, 1.0) );
    conversionPeriod->SetCoCo(dTriggerRate, CoCoType_CheckQuarterlyAndConvertAsOfCheckDate, 0.0, dTriggerRate);

    Date startConvDate2 = endConvDate;
    Date endConvDate2 = endConvDate;
    startConvDate2.AddMonths(2);
    endConvDate2.AddMonths(4);
    shared_ptr<ConversionPeriod> conversionPeriod2(new ConversionPeriod(startConvDate2, endConvDate2, 1.0) );
    conversionPeriod2->SetCoCo(dTriggerRate, CoCoType_CheckAnyTimeAndConvertOnCheckDate, 0.0, dTriggerRate);
*/
    pConv->SetKeepAccrued(true);
    pConv->AddConversionPeriod(conversionPeriod);
    //pConv->AddConversionPeriod(conversionPeriod2);

    //pCB->SetCallSchedule(pCall);
    //pCB->SetPutSchedule(pPut);
    pCB->SetConversionSchedule(pConv);



    // Contingent call section
    size_t nDays = 5;
    shared_ptr<CallSchedule> pCallSchedule( new CallSchedule() );
    pCallSchedule->SetTriggerCheckPeriod(nDays, 0);
    pCallSchedule->SetKeepAccrued(true);
    shared_ptr<CallPeriod> 
      pCallPeriod( CallPeriod::CreateWithStrike
                    (startCallDate, endCallDate, 1.0) ) ;
    pCallPeriod->SetTrigger(1.0);

    pCallSchedule->AddCallPeriod( pCallPeriod);
  
    pCB->SetCallSchedule(pCallSchedule);

    shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
    flags->SetAnalysisDate( issueDate );

    pCB->SetComputationalFlags(flags);

    shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
    pModel->SetVolatility( pVolatility );
    pModel->SetHazardRate( pHazardRate );

    shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);

    std::cout.precision(15);
    std::cout << "testcocofinancialinterface" << std::endl;
    std::cout << "  Conversion trigger window: " 
              << startConvDate << " to " << endConvDate << std::endl;
    std::cout << "  Trigger rate: " << dTriggerRate << std::endl;
    std::cout << "  CB price = " << output->GetPrice() << std::endl;
    std::cout << std::endl;

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
