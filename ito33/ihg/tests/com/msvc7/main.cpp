#include <iostream>

#include "ito33/common.h"
#include "ito33/string.h"

#include "ito33/com/safearray.h"

#pragma hdrstop

#import "../../../../common/msvc7/itocom/itocom.tlb"
#import "../../msvc7/ihgcom/ihg.tlb"


class COMInit
{
  public:
    COMInit() { ::CoInitialize(NULL); }
    ~COMInit() { ::CoUninitialize(); }
};

int main()
{
  using namespace std;

  COMInit initCOM;

  try
  {
    IHG::IGlobalIHGPtr ihg("IHG.GlobalIHG");
    ITO33::IITO33Ptr common("ITO33.ITO33");

    ito33::COM::BasicString currencyCode("EUR");

    ITO33::INumerairePtr 
      currency(common->NewNumeraire(currencyCode.Detach()));
  
    ITO33::IRateDataPtr rateData("ITO33.RateData");
    rateData->SetYieldCurve(currency, common->NewYieldCurveFlat(0.01));

    ITO33::IEquityPtr equity(common->NewEquity(33.5, currency));

    equity->BorrowCurve = common->NewYieldCurveFlat(0.01);
    
    ITO33::ISessionDataPtr 
      sessionData(common->NewSessionData(rateData, equity, 37653));

    IHG::ITheoreticalModelPtr tm("IHG.TheoreticalModel");
    tm->Volatility = ihg->NewVolatilityFlat(0.2);
    ito33::COM::SafeArray<DATE> aTimes(2);
    ito33::COM::SafeArray<double> aValues(2);
    {
      ito33::COM::SafeArrayAccessor<DATE> times(aTimes);
      ito33::COM::SafeArrayAccessor<double> values(aValues);

      times[0] = 38000;
      times[1] = 40000;
      values[0] = 0.02;
      values[1] = 0.025;
    }
    tm->HazardRate = ihg->NewHazardRateTimeOnly(aTimes, aValues);

    tm->DebugOutputFile = "./ihg.xml";

    ITO33::IComputationalFlagsPtr flags("ITO33.ComputationalFlags");
    flags->ComputeRho = VARIANT_TRUE;
    flags->ComputeVega = VARIANT_TRUE;

    bool option_price = true;

    if ( option_price == true )
    {
      ITO33::IOptionPtr option(common->NewOption(
            50,     // Strike
            38353,  // MaturityDate
            ITO33::Option_Call, 
            ITO33::ExerciseType_American));
      option->MarketPrice = 4.35;
      option->SessionData = sessionData;

      option->ComputationalFlags = flags;

      ITO33::IModelOutputPtr output = tm->Compute(
          reinterpret_cast<ITO33::IDerivative*>(option.GetInterfacePtr()));

      cout << "IHG Output Option:" << endl
           << "========================================" << endl
           << "Price:\t\t" << output->Price << endl  // ->get_Price(&x) for
           << "Gamma:\t\t" << output->Gamma << endl  // more compatibility.
           << "Delta:\t\t" << output->Delta << endl
           << "Vega: \t\t" << output->Vega  << endl;
    }
    else 
    {
      ITO33::ICashFlowStreamUniformPtr 
        pCashFlow(common->NewCashFlowStreamUniform
            (
             37653,
             37742,
             40664,
             0.025,
             ITO33::DayCountConvention_Act365,
             ITO33::Frequency_BiMonthly
            )
            );

      ITO33::ICDSPtr cds( common->NewCDS(.01, pCashFlow) );

      cds->SessionData = sessionData;

      cds->ComputationalFlags = flags;

      ITO33::IModelOutputPtr output = tm->Compute(
          reinterpret_cast<ITO33::IDerivative*>(cds.GetInterfacePtr()));

      cout << "IHG Output CDS:" << endl
           << "========================================" << endl
           << "Price:\t\t" << output->Price << endl  // ->get_Price(&x) for
           << "Gamma:\t\t" << output->Gamma << endl  // more compatibility.
           << "Delta:\t\t" << output->Delta << endl // ;
           << "Vega: \t\t" << output->Vega  << endl;
    }

  } //end try
  catch ( const _com_error& ce )
  {
    cout << "COM exception: " << ce.Error() << endl;

    IErrorInfo *pErrorInfo;
   
    pErrorInfo = ce.ErrorInfo();
    
    if ( pErrorInfo )
    {
      BSTR bstr;
      pErrorInfo->GetDescription(&bstr);
      if ( *bstr )
      {
        cout << ": " << (char *)_bstr_t(bstr);
      }

      pErrorInfo->Release();
    }

    if ( SUCCEEDED(::GetErrorInfo(0, &pErrorInfo)) && pErrorInfo )
    {
      BSTR bstr;
      pErrorInfo->GetDescription(&bstr);
      if ( *bstr )
      {
        cout << ": " << (char *)_bstr_t(bstr);
      }

      pErrorInfo->Release();
    }

    cout << endl;
  }

  return 0;
}
