/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xll/ihgxll.cpp
// Purpose:     main file of IHG XLL implementation
// Author:      Vadim Zeitlin
// Created:     2006-04-11
// RCS-ID:      $Id: ihgxll.cpp,v 1.78 2006/08/23 16:31:52 willy Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/XL/addin.h"
#include "ito33/XL/sheetobject.h"

#include "ito33/link.h"
#include "ito33/optionalvalue.h"
#include "ito33/numeric/interpolation.h"

#include "ito33/finance/computationalflags.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/referencecds.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/floatingrates.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/option.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_swap.h"

#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/convertiblebond.h"
#include "ito33/finance/bondlike/cboption.h"
#include "ito33/finance/bondlike/reset.h"
#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/bondterms.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/pepsaveragingperiod.h"
#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/finance/bondlike/utils.h"

#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/exercisetype.h"
#include "ito33/xml/finance/optiontype.h"

#include "ito33/ihg/version.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/perfect_hedge_ratios.h"

#include "ito33/ihg/bondlikeoutput.h"

#include <algorithm>

ITO33_FORCE_LINK_MODULE(IHGPriceReset);
ITO33_FORCE_LINK_MODULE(IHGPriceCB);
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceOption);
ITO33_FORCE_LINK_MODULE(IHGPriceBond);
ITO33_FORCE_LINK_MODULE(IHGPriceCBOption);
ITO33_FORCE_LINK_MODULE(IHGPriceAttachedWarrantConvertibleBond);
ITO33_FORCE_LINK_MODULE(IHGPriceGeneralizedPEPSLike);

using namespace ito33;

// ----------------------------------------------------------------------------
// helper classes
// ----------------------------------------------------------------------------

typedef shared_ptr<finance::ModelOutput> ModelOutputPtr;

struct OutputData
{
  shared_ptr<finance::Derivative> Derivative;
  shared_ptr<finance::SessionData> Sessiondata;
  ModelOutputPtr Output;
  shared_ptr<ihg::TheoreticalModel> Model;
};

typedef XL::SheetObjectManager<OutputData> OutputDataManager;

typedef shared_ptr<ihg::PerfectHedgeRatios> PerfectHedgeRatiosPtr;

typedef XL::SheetObjectManager<PerfectHedgeRatiosPtr> HedgeRatiosManager;

typedef XL::SheetObjectManager< shared_ptr<finance::YieldCurve> >
        YieldCurvesManager;

// objects created by ITO33CurrencyData() formulae
struct CurrencyData
{
  shared_ptr<finance::Numeraire> currency;
  shared_ptr<finance::YieldCurve> yieldcurve;
};

typedef XL::SheetObjectManager<CurrencyData> CurrenciesManager;

// objects created by ITO33Equity() formulae
struct EquityData
{
  shared_ptr<finance::Numeraire> currency;
  shared_ptr<finance::YieldCurve> yieldcurve;
  shared_ptr<finance::Equity> equity;
};

typedef XL::SheetObjectManager<EquityData> EquitiesManager;

// objects created by ITO33Cashflow(), ITO33UniformCashflow() or ITO33FloatCashflow()
struct CashflowData
{
  // used for non floating cashflows
  shared_ptr<finance::CashFlowStream> cashflow;

  // use for the floaters
  shared_ptr<finance::FloatingRates> floaters;


  // default ctor only because we have some other ones -- it does nothing
  CashflowData()
  {
  }

  // ctor for non floaters
  CashflowData(const shared_ptr<finance::CashFlowStream> cashflow_)
    : cashflow(cashflow_)
  {
  }

  // another ctor for non floaters: takes ownership of the pointer
  CashflowData(finance::CashFlowStream *cashflow_)
    : cashflow(cashflow_)
  {
  }

  // use this method instead of modifying floaters directly because it also
  // resets the cashflow a only one of the pointers above can be non NULL
  void SetFloater(const shared_ptr<finance::FloatingRates>& floaters_)
  {
    cashflow.reset();
    floaters = floaters_;
  }
};

typedef XL::SheetObjectManager<CashflowData> CashflowsManager;

// objects created by ITO33OID() and ITO33CashPayToZero()
struct OIDData
{
  OIDData() { hasAccrRate = false; grossY2M = 0.; }

  // use these methods instead of changing the corresponding fields directly
  void SetOID(double grossY2M_)
  {
    hasAccrRate = false;
    grossY2M = grossY2M_;
  }

  void SetCashPayToZero(double accrRate_)
  {
    hasAccrRate = true;
    accrRate = accrRate_;
  }


  shared_ptr<finance::CashFlowStream> cashflow;

  bool hasAccrRate;
  union
  {
    double grossY2M;
    double accrRate;
  };

  finance::Frequency frequency;
  Date::DayCountConvention dcc;
};

typedef XL::SheetObjectManager<OIDData> OIDsManager;

struct CrossCurrencyData
{
  double spotFXRate;
  std::string baseCurrency;
  std::string foreignCurrency;

  // Needed parameter for fixed quanto cb
  bool isFixedQuanto;
  double volatility;
  double correlation;
  EquityData equityData;

  CrossCurrencyData():spotFXRate(0.),isFixedQuanto(false)
  {
  }

  void SetCrossCurrency(std::string base, std::string foreign, double fxRate)
  {
    spotFXRate =  fxRate;
    baseCurrency = base;
    foreignCurrency = foreign;
  }

  void SetFixedQuanto(double dVol, double dCorrel,const EquityData& eq)
  {
    volatility = dVol;
    correlation = dCorrel;
    equityData = eq;

    isFixedQuanto = true;
  }

  bool CompareEquityDataWith(const EquityData& eq)
  {
    return 
      equityData.currency->GetCode() == eq.currency->GetCode() &&
      equityData.equity == eq.equity &&
      equityData.yieldcurve == eq.yieldcurve;
  }
};

typedef XL::SheetObjectManager< shared_ptr<CrossCurrencyData> >
        CrossCurrenciesManager;

// objects created by ITO33Conversion() and ITO33Conversions()
typedef XL::SheetObjectManager< shared_ptr<finance::ConversionSchedule> >
        ConversionsManager;

typedef XL::SheetObjectManager< shared_ptr<finance::ResetConversionSchedule> >
        ResetsManager;

typedef XL::SheetObjectManager< shared_ptr<finance::ShareDependentConversion> >
        ShareDepConversionsManager;

// objects created by ITO33BondTerms()
typedef XL::SheetObjectManager< shared_ptr<finance::BondTerms> > BondTermsManager;

struct CBTerms
{
  bool newShare;
  finance::TriggerAsPercentageOf triggerAsPercentageOf;
  finance::TriggerInCurrencyOf triggerInCurrencyOf;
  shared_ptr<finance::Issuer> CBIssuer;
  bool IsExchangeUponDefault;

  CBTerms():newShare(false),
            triggerAsPercentageOf(finance::TriggerAsPercentageOf_Principal),
            triggerInCurrencyOf(finance::TriggerInCurrencyOf_Derivative),
            IsExchangeUponDefault(false)
  {
  }
};

//objects created by ITO33CBTerms()
typedef XL::SheetObjectManager< shared_ptr<CBTerms> > CBTermsManager;

// objects created by ITO33PEPSTerms()

struct PEPSTerms
{
  double downConvRatio;
  double lowerStrike;
  double upConvRatio;
  double higherStrike;
  bool hasOptionalConv;
  shared_ptr<finance::PEPSAveragingPeriod> PEPSAveraging;
};

typedef XL::SheetObjectManager< shared_ptr<PEPSTerms> > PEPSTermsManager;

typedef XL::SheetObjectManager<shared_ptr<finance::GeneralizedPEPSLikeCall> > PEPSLikeCallManager;

typedef XL::SheetObjectManager<shared_ptr<finance::PEPSAveragingPeriod> > PEPSAveragingManager;

// objects created by ITO33CB(), ITO33Option(), ITO33CDS() &c formulae
struct DerivativeData
{
  shared_ptr<finance::Numeraire> eqCurrency;
  shared_ptr<finance::YieldCurve> eqYieldcurve;
  shared_ptr<finance::Equity> equity;
  
  shared_ptr<finance::Derivative> derivative;
  shared_ptr<finance::Numeraire> derCurrency;
  shared_ptr<finance::YieldCurve> derYieldcurve;
  
  //Fx Rate for cross currency derivative
  shared_ptr<finance::SpotFXRates> FXRate;

  bool IsCBOption;
  
  void InitFrom(const EquityData& equityData)
  {
    ResetCrossCurrency();
    eqCurrency = equityData.currency;
    eqYieldcurve = equityData.yieldcurve;
    equity = equityData.equity;
  }
  //use it to set a cross currency derivative
  void SetCrossCurrency(const shared_ptr<finance::Numeraire> & derivCurr,
                        const shared_ptr<finance::YieldCurve>& derivYC,
                        double dFXRate)
  {
    derCurrency = derivCurr;
    derYieldcurve = derivYC;
    FXRate = make_ptr( new finance::SpotFXRates );
    FXRate->SetFXRate(eqCurrency,derCurrency,dFXRate);
  }
  
  void ResetAll()
  {
    eqCurrency.reset(); 
    eqYieldcurve.reset(); 
    equity.reset();
    ResetCrossCurrency();
  }

  //Unset cross currency object to NULL
  void ResetCrossCurrency()
  {
    derivative.reset();
    derCurrency.reset();
    derYieldcurve.reset();
    FXRate.reset();
  }
  //Init the derivative data
  DerivativeData():IsCBOption( false ){}
};

typedef XL::SheetObjectManager<DerivativeData> DerivativesManager;

// this is a somewhat degenerate case: we don't create any objects in
// ITO33FloatingRates() but just remember its parameters in this struct to
// allow creating FloatingRates in ITO33FloatCashflow() later -- the only
// benefit of having it is that it allows us to have less parameters in
// ITO33FloatCashflow() by moving some of them to ITO33FloatingRates()
struct FloaterData
{
  OptionalValue<double> 
    cap,
    floor,
    multiplier;

  OptionalValue<int> fixingDelay;
  
  OptionalValue<Date::DayCountConvention> dcc;

  finance::Frequency frequency;
  
  double margin;
  
  Date startAccr,
       firstUnknown,
       lastUnknown;

  void ResetAll()
  {
    cap.Reset();
    floor.Reset();
    multiplier.Reset();
    fixingDelay.Reset();
    dcc.Reset();
  }
};

typedef XL::SheetObjectManager<FloaterData> FloatersManager;

struct AdvCallFeatures
{
  int noticePeriod;
  int triggerPeriod;
  int triggerHistory;
};

typedef XL::SheetObjectManager<AdvCallFeatures> AdvCallFeaturesManager;

// objects created by ITO33Calls[ByYield]() and ITO33Puts()
typedef XL::SheetObjectManager< shared_ptr<finance::CallSchedule> > CallsManager;
typedef XL::SheetObjectManager< shared_ptr<finance::PutSchedule> > PutsManager;

// ----------------------------------------------------------------------------
// main add-in class
// ----------------------------------------------------------------------------

class IHGAddIn : public XL::AddIn
{
public:
  IHGAddIn() : XL::AddIn("ITO33 IHG Add-in", "ITO33 IHG") { }

  // called by the library to register all the exported functions
  virtual void RegisterAllFunctions();

  // accessors for the collections of various objects we use
  AdvCallFeaturesManager& AdvCallFeats() { return m_advcall; }
  BondTermsManager& BondTerms() { return m_bondterms; }
  CBTermsManager& CBTerms() { return m_cbterms; }
  CallsManager& Calls() { return m_calls; }
  CashflowsManager& Cashflows() { return m_cashflows; }
  ConversionsManager& Conversions() { return m_conversions; }
  CrossCurrenciesManager& CrossCurrencies() { return m_crosscurrency; }
  CurrenciesManager& Currencies() { return m_currencies; }
  DerivativesManager& Derivatives() { return m_derivatives; }
  EquitiesManager& Equities() { return m_equities; }
  FloatersManager& Floaters() { return m_floaters; }
  HedgeRatiosManager& HedgeRatios() { return m_hedge; }
  OIDsManager& OIDs() { return m_oids; }
  OutputDataManager& Outputs() { return m_outputs; }
  PEPSAveragingManager& PEPSAveragingPeriods() { return m_pepsaveraging; }
  PEPSLikeCallManager& PEPSLikeCalls() { return m_pepscalls; }
  PEPSTermsManager& PEPSTerms() { return m_pepsterms; }
  PutsManager& Puts() { return m_puts; }
  ResetsManager& Resets() {return m_reset;}
  ShareDepConversionsManager& ShareDepConversions() {return m_sharedepconvs;}
  YieldCurvesManager& YieldCurves() { return m_yieldcurves; }

private:
  AdvCallFeaturesManager m_advcall;
  BondTermsManager m_bondterms;
  CallsManager m_calls;
  CashflowsManager m_cashflows;
  CBTermsManager m_cbterms;
  CrossCurrenciesManager m_crosscurrency;
  ConversionsManager m_conversions;
  CurrenciesManager m_currencies;
  DerivativesManager m_derivatives;
  EquitiesManager m_equities;
  FloatersManager m_floaters;
  HedgeRatiosManager m_hedge;
  OIDsManager m_oids;
  OutputDataManager m_outputs;
  PEPSAveragingManager m_pepsaveraging;
  PEPSLikeCallManager m_pepscalls;
  PEPSTermsManager m_pepsterms;
  PutsManager m_puts;
  ResetsManager m_reset;
  ShareDepConversionsManager m_sharedepconvs;
  YieldCurvesManager m_yieldcurves;


  NO_COPY_CLASS(IHGAddIn);
};

XL_IMPLEMENT_ADDIN(IHGAddIn)


// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

// this is like RestoreEnumValueFromName() except that it is case-insensitive
// and understands non-ambiguous abbreviations
//
// note that N template parameter has to be explicitly specified, the compiler
// is unable to deduce it for us unfortunately, this is why it comes first
template <int N, typename E>
static bool
GetEnumFromAbbrevName(E *e,
                      const char *str,
                      const EnumValuesNames<E> names[N])
{
  // keep track of the names which match str so far
  bool matches[N];

  // initially everything matches
  std::fill(matches, matches + N, true);
  int numMatches = N;

  // check for each character in the string if the names still match it
  for ( int c = 0; *str != '\0'; c++, str++ )
  {
    char ch = *str;
    if ( isupper(ch) )
      ch = static_cast<char>(tolower(ch));

    for ( int n = 0; n < N; n++ )
    {
      if ( matches[n] )
      {
        if ( names[n].name[c] != ch )
        {
          matches[n] = false;
          if ( !--numMatches )
          {
            // nothing matches str
            return false;
          }
        }
      }
      //else: this one already doesn't match, don't care about it
    }
  }

  // find the unique value corresponding to it
  bool hasValue = false;
  for ( int n = 0; n < N; n++ )
  {
    if ( matches[n] )
    {
      if ( hasValue )
      {
        // check that it is the same
        if ( names[n].value != *e )
        {
          // ambiguous string
          return false;
        }
      }
      else // first match
      {
        *e = names[n].value;
        hasValue = true;
      }
    }
  }

  return true;
}

// throws an XL::Exception with the message of the form "prefix: format"
static void ThrowValueError(const char *prefix, const char *format, va_list ap)
{
  throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              String::Printf
              (
                "%s: %s",
                prefix,
                String::PrintfV(format, ap).c_str()
              )
            );
}

// convert a string of the form "30/360" to a DayCountConvention value, throws
// if the string is invalid
static Date::DayCountConvention DccFromExcel(const char *dcc)
{
  // first get the DCC value, for this first remove any slashes (i.e. accept
  // "act/act" as "actact") from the provided string
  std::string dccNormalized;
  while ( *dcc )
  {
    const char ch = *dcc++;
    if ( ch == '/' )
      continue;

    dccNormalized += isupper(ch) ? static_cast<char>(tolower(ch)) : ch;
  }

  Date::DayCountConvention dccVal;
  if ( !RestoreEnumValueFromName
        (
          dccNormalized.c_str(),
          dccVal,
          SIZEOF(g_dayCountConventions),
          g_dayCountConventions)
     )
  {
    throw XL::EXCEPTION_MSG
              (
                xlerrValue,
                String::Printf("Not a valid day count string \"%s\"", dcc)
              );
  }

  return dccVal;
}

static inline Date::DayCountConvention DccFromExcel(const std::string& dcc)
{
  return DccFromExcel(dcc.c_str());
}

// convert to Frequency enum value, throw if the parameter is invalid
static finance::Frequency FreqFromExcel(const char *freqStr)
{
  finance::Frequency freq;
  unsigned long freqVal;
  if ( String::ToULong(freqStr, &freqVal) )
  {
    // if frequency is a number, it should be given as a value of our enum so
    // no conversion is needed but check for invalid ones
    freq = static_cast<finance::Frequency>(freqVal);
    if ( !finance::IsValid(freq) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Not a valid payment frequency %d", freqVal)
                );
    }
  }
  else // not a number
  {
    static const EnumValuesNames<finance::Frequency> frequencies[] =
    {
      "annual",       finance::Frequency_Annual,
      "semiannual",   finance::Frequency_SemiAnnual,
      "semi-annual",  finance::Frequency_SemiAnnual,
      "quarterly",    finance::Frequency_Quarterly,
      "bimonthly",    finance::Frequency_BiMonthly,
      "bi-monthly",   finance::Frequency_BiMonthly,
      "monthly",      finance::Frequency_Monthly,
    };

    if ( !GetEnumFromAbbrevName<SIZEOF(frequencies)>(&freq, freqStr, frequencies) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown payment frequency \"%s\"", freqStr)
                );
    }
  }

  return freq;
}

static inline finance::Frequency FreqFromExcel(const std::string& freq)
{
  return FreqFromExcel(freq.c_str());
}

// convert to Dividend::Type enum value, throw if the parameter is invalid
static finance::Dividend::Type DividendTypeFromExcel(const char* divtype)
{
  finance::Dividend::Type type;

  static const EnumValuesNames<finance::Dividend::Type> dividendTypes[] =
    {
      { "cash",   finance::Dividend::Cash   },
      { "yield",  finance::Dividend::Yield  },
    };

    if ( !GetEnumFromAbbrevName<SIZEOF(dividendTypes)>
          (
            &type,
            divtype,
            dividendTypes
          ) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown dividend type \"%s\"", divtype)
                );
    }
    

  return type;
}

static inline finance::Dividend::Type DividendTypeFromExcel(const std::string& divtypeStr)
{
  return DividendTypeFromExcel( divtypeStr.c_str() );
}
// convert to TimeUnit enum value, throw if the parameter is invalid
static finance::TimeUnit TimeUnitFromExcel(const char *unit)
{
  static const EnumValuesNames<finance::TimeUnit> timeUnits[] =
  {
    "days",   finance::TimeUnit_Day,
    "months", finance::TimeUnit_Month,
    "years",  finance::TimeUnit_Year,
  };

  finance::TimeUnit unitVal;
  if ( !GetEnumFromAbbrevName<SIZEOF(timeUnits)>(&unitVal, unit, timeUnits) )
  {
    throw XL::EXCEPTION_MSG
              (
                xlerrValue,
                String::Printf("Unknown time unit \"%s\"", unit)
              );
  }

  return unitVal;
}

static inline finance::TimeUnit TimeUnitFromExcel(const std::string& unit)
{
  return TimeUnitFromExcel(unit.c_str());
}

static finance::TriggerAsPercentageOf TriggerPercentFromExcel(const char* triggerPercentStr )
{
  static const EnumValuesNames<finance::TriggerAsPercentageOf> TriggerPercentValues[] =
  {
    { "principal", finance::TriggerAsPercentageOf_Principal },
    { "issue price", finance::TriggerAsPercentageOf_IssuePrice },
    { "claim", finance::TriggerAsPercentageOf_Claim }
  };
  finance::TriggerAsPercentageOf triggerAsPercentage;   
  if ( 
       !GetEnumFromAbbrevName<SIZEOF(TriggerPercentValues)>
       (
         &triggerAsPercentage,
         triggerPercentStr,
         TriggerPercentValues
       ) 
     )
  {
    throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("Unknown trigger percent of \"%s\"", triggerPercentStr)
          );
  }
  return triggerAsPercentage;
}

static inline finance::TriggerAsPercentageOf TriggerPercentFromExcel(const std::string& triggerPercentStr)
{
  return TriggerPercentFromExcel( triggerPercentStr.c_str() );
}
// converts to MakeWholeType, returns MakeWholeType_Max if not specified or
// none (because there is no MakeWholeType_None)
//
// if MakeWholeType_Coupon is returned, isCouponPresentValue is filled with the
// appropriate argument for Call::SetCouponMakeWhole()
static finance::MakeWholeType
MakeWholeTypeFromExcel(const char *str, bool *isCouponPresentValue)
{
  // all valid strings ("CA", "CP", "P" and "N" or empty) are handled inside
  // this switch, the invalid ones break out of it and throw an exception below
  switch ( str ? toupper(*str) : '\0' )
  {
    case 'C':
      {
        const int ch = toupper(str[1]);
        if ( ch == 'P' )
          *isCouponPresentValue = true;
        else if ( ch == 'A' )
          *isCouponPresentValue = false;
        else
          break;

        if ( str[2] != '\0' )
          break;
      }
      return finance::MakeWholeType_Coupon;

    case 'P':
      if ( str[1] != '\0' )
        break;
      return finance::MakeWholeType_Premium;

    case 'N':
      if ( str[1] != '\0' )
        break;
      // fall through

    case '\0':
      // empty or not specified
      return finance::MakeWholeType_Max;
  }

  throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              String::Printf("Unknown make whole type \"%s\"", str)
            );
}

// return the hazard rate object corresponding to the given parameter which can
// be either a scalar (flat hazard rate) or an array/range (time-dependent HR)
//
// if the parameter can't be coerced to any hazard rate form, throw an exception
static shared_ptr<ihg::HazardRate> HRFromExcel(const XL::Oper *xloHR)
{
    shared_ptr<ihg::HazardRate> hr;

    std::vector<Date> dates;
    std::vector<double> values;
    double hrFlat;
    if ( xloHR->AsTable(dates, values) )
    {
      hr = make_ptr( new ihg::HazardRateTimeOnly(dates, values) );
    }
    else if ( xloHR->As(&hrFlat) )
    {
      hr = make_ptr( new ihg::HazardRateFlat(hrFlat) );
    }

    if ( !hr )
    {
      throw XL::EXCEPTION_MSG(xlerrValue, "Unrecognized hazard rate format");
    }

    return hr;
}

// return the valuation date object from the given parameter or today if the
// parameter is missing
//
// throws if the parameter is not missing but can't be converted to date
static Date ValDateFromExcel(const XL::Oper *xloDate)
{
  Date date;
  if ( !XL::Oper::IsMissing(xloDate) )
  {
    if ( !xloDate->As(&date) )
    {
      throw XL::EXCEPTION_MSG(xlerrValue, "Invalid valuation date format");
    }
  }
  else // reference date is today by default
  {
    date = Date::Today();
  }

  return date;
}

// ============================================================================
// exported functions
// ============================================================================

// ----------------------------------------------------------------------------
// miscellaneous helpers
// ----------------------------------------------------------------------------

/*
   Template for a new function:

XLLEXPORT XLOPER *
ITO33()
{
  try
  {
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

*/

XLLEXPORT const char *
ITO33IHGVersion()
{
  return ITO33_IHG_VERSION_DOT_STRING;
}

XLLEXPORT XLOPER *ITO33T(const XL::Oper *xlo)
{
  try
  {
    // retrieve the array data
    XL::Oper xloArray;
    if ( !xlo->AsArray(&xloArray) )
      return XL::Oper::ErrValue();

    unsigned w, h;
    const XL::Oper *src = xloArray.GetAsArray(w, h);

    const unsigned size = w*h;
    XL::Oper *data = new XL::Oper[size];

    // transpose it into the new array
    XL::Oper *dst = data;
    for ( const XL::Oper * const srcEnd = src + w; src < srcEnd; ++src )
    {
      const XL::Oper *col = src;
      for ( unsigned y = 0; y < h; ++y, ++dst, col += w )
      {
        *dst = *col;
      }
    }

    return XL::Oper::ReturnArray(h, w, data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *ITO33Split(const XL::Oper* xlo, int inbxCol)
{
  try
  {
    // retrieve the array data
    XL::Oper xloArray;
    if ( !xlo->AsArray(&xloArray) )
      return XL::Oper::ErrValue();

    if ( !inbxCol )
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("second parameter can not be null")
                );
    
    unsigned w, h;
    const XL::Oper *src = xloArray.GetAsArray(w, h);

    if ( h != 1 || (w % inbxCol) )
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("invalid range")
                );

    unsigned 
      newh(w / inbxCol),
      neww( inbxCol );
    const unsigned size = w;
    XL::Oper *data = new XL::Oper[size];

    XL::Oper *cpdata = data;
    for(const XL::Oper * const srcEnd = src + w ; src != srcEnd  ; ++src,++cpdata)
    {
      *cpdata = *src; 
    }

    return XL::Oper::ReturnArray(neww, newh, data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// yield cuve objects creation
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33ZeroCouponYC(XL::Oper *xloLegs, XL::Oper *xloDate)
{
  try
  {
    // check/extract parameters
    std::vector<int> daysOrig;
    std::vector<double> rates;
    if ( !xloLegs || !xloLegs->AsTable(daysOrig, rates) )
    {
      return XL::Oper::ErrValue();
    }

    const size_t ndays = daysOrig.size();
    std::vector<size_t> days(ndays);
    for ( size_t n = 0; n < ndays; n++ )
    {
      if ( daysOrig[n] < 0 )
        return XL::Oper::ErrValue();

      days[n] = static_cast<size_t>(daysOrig[n]);
    }

    // create the yield curve object and initialize it
    shared_ptr<finance::YieldCurveAnnuallyCompounded>
      yc(new finance::YieldCurveAnnuallyCompounded(ValDateFromExcel(xloDate)));
    yc->SetLegs(days, rates);

    // and remember it
    return GetTheIHGAddIn().YieldCurves().SetForThisCell(yc);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// used by ITO33YC() to signal errors with swaps data
static void ThrowSwapFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid swaps data format", format, ap);

  va_end(ap);
}

static void ThrowCashFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid cash rate data format", format, ap);

  va_end(ap);
}

XLLEXPORT XLOPER *
ITO33YC(XL::Oper *xloSwapRates,
        XL::Oper *xloSwapDCC,
        XL::Oper *xloCashRates,
        XL::Oper *xloCashDCC,
        XL::Oper *xloDate)
{
  try
  {
    // parse swaps rates data
    enum
    {
      Col_SwapDuration,
      Col_SwapTimeUnit,
      Col_SwapRate,
      Col_SwapFreq,
      Col_SwapMax
    };

    finance::SwapRates swapRates;
    if ( !XL::Oper::IsMissing(xloSwapRates) )
    {
      unsigned w, h;
      const XL::Oper *data = xloSwapRates->GetAsArray(w, h);
      if ( w != Col_SwapMax )
        ThrowSwapFormatError("incorrect number of columns");

      for ( unsigned n = 0; n < h; ++n, ++data )
      {
        int duration;
        if ( !data->As(&duration) )
          ThrowSwapFormatError("bad swap duration format in row %u", n + 1);

        std::string unit;
        if ( !(++data)->As(&unit) )
          ThrowSwapFormatError("bad swap time unit format in row %u", n + 1);

        double rate;
        if ( !(++data)->As(&rate) )
          ThrowSwapFormatError("bad swap rate format in row %u", n + 1);

        std::string freq;
        if ( !(++data)->As(&freq) )
          ThrowSwapFormatError("bad swap frequency format in row %u", n + 1);

        swapRates.AddLeg
                  (
                    rate,
                    duration,
                    TimeUnitFromExcel(unit),
                    FreqFromExcel(freq)
                  );
      }

      Date::DayCountConvention dcc;
      if ( !XL::Oper::IsMissing(xloSwapDCC) )
      {
        std::string dccStr;
        if ( !xloSwapDCC->As(&dccStr) )
          ThrowSwapFormatError("bad swaps dcc format in row %u", n + 1);

        dcc = DccFromExcel(dccStr);
      }
      else // default dcc
      {
        dcc = Date::DayCountConvention_ActAct;
      }

      swapRates.SetBasis(dcc);
    }
    else // it doesn't make sense to use swaps dcc if we don't have them
    {
      if ( !XL::Oper::IsMissing(xloSwapDCC) )
        ThrowSwapFormatError("swaps dcc must be omitted if no swaps");
    }


    // parse cash rates data
    enum
    {
      Col_CashDuration,
      Col_CashTimeUnit,
      Col_CashRate,
      Col_CashMax
    };

    finance::CashRates cashRates;
    if ( !XL::Oper::IsMissing(xloCashRates) )
    {
      unsigned w, h;
      const XL::Oper *data = xloCashRates->GetAsArray(w, h);
      if ( w != Col_CashMax )
        ThrowCashFormatError("incorrect number of columns");

      for ( unsigned n = 0; n < h; ++n, ++data )
      {
        int duration;
        if ( !data->As(&duration) )
          ThrowCashFormatError("bad cash rate duration format in row %u", n + 1);

        std::string unit;
        if ( !(++data)->As(&unit) )
          ThrowCashFormatError("bad cahs rate time unit format in row %u", n + 1);

        double rate;
        if ( !(++data)->As(&rate) )
          ThrowCashFormatError("bad cash rate format in row %u", n + 1);

        cashRates.AddLeg
                  (
                    rate,
                    duration,
                    TimeUnitFromExcel(unit)
                  );
      }

      Date::DayCountConvention dcc;
      if ( !XL::Oper::IsMissing(xloCashDCC) )
      {
        std::string dccStr;
        if ( !xloCashDCC->As(&dccStr) )
          ThrowCashFormatError("bad swaps dcc format in row %u", n + 1);

        dcc = DccFromExcel(dccStr);
      }
      else // default dcc
      {
        dcc = Date::DayCountConvention_Act360;
      }

      cashRates.SetBasis(dcc);
    }
    else // it doesn't make sense to use cash rates dcc if we don't have them
    {
      if ( !XL::Oper::IsMissing(xloCashDCC) )
        ThrowCashFormatError("cash rates dcc must be omitted if no cash rates");
    }


    // create the yield curve object
    shared_ptr<finance::YieldCurveSwap>
      yc(new finance::YieldCurveSwap(ValDateFromExcel(xloDate)));
    if ( !swapRates.GetAll().empty() )
      yc->SetSwapRates(swapRates);
    if ( !cashRates.GetAll().empty() )
      yc->SetCashRates(cashRates);

    // remember the yield curve
    return GetTheIHGAddIn().YieldCurves().SetForThisCell(yc);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// miscellaneous utility objects creation
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33CurrencyData(const char *name, XL::Oper *xloYC)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    shared_ptr<finance::YieldCurve> yieldcurve;
    YieldCurvesManager::Handle hyc;
    if ( xloYC->As(&hyc) )
    {
      const shared_ptr<finance::YieldCurve> * const
        yc = addin.YieldCurves().FromHandle(hyc);
      if ( yc )
        yieldcurve = *yc;
    }

    if ( !yieldcurve )
    {
      double ycFlat;
      if ( !xloYC->As(&ycFlat) )
        return XL::Oper::ErrValue();

      yieldcurve = make_ptr( new finance::YieldCurveFlat(ycFlat) );
    }


    XLL_TRACE("Creating new currency data \"%s\"", name);

    CurrenciesManager& manager = addin.Currencies();
    CurrencyData& data = manager.GetForThisCell();
    data.currency = make_ptr( new finance::Numeraire(name) );
    data.yieldcurve = yieldcurve;

    return manager.HandleOf(data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// cashflow objects creation
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33UniformCashflow(double issue,
                     double first,
                     double last,
                     double premium,
                     const char *dcc,
                     const char *freq,
                     const XL::Oper* xloLastPaymentType)
{
  try
  {
    // create the cashflow stream (this can throw)
    const Date dateFirst(Date::FromExcel(first));
    const Date dateLast(Date::FromExcel(last));
    shared_ptr<finance::CashFlowStream> cashflow;

    if (!XL::Oper::IsEmpty(xloLastPaymentType))
    {
      static const EnumValuesNames<finance::LastPaymentType> lastpaymentTypes[] =
      {
        "short",    finance::LastPaymentType_Short,
        "long",    finance::LastPaymentType_Long
      };
      
      std::string lastpaymentStr;
      if ( !xloLastPaymentType->As(&lastpaymentStr))
        return XL::Oper::ErrValue();
      
      
      finance::LastPaymentType lastpaymentType;

      if ( !GetEnumFromAbbrevName<SIZEOF(lastpaymentTypes)>
          (
            &lastpaymentType,
            lastpaymentStr.c_str(),
            lastpaymentTypes
          ) )
      {
        throw XL::EXCEPTION_MSG
                  (
                    xlerrValue,
                    String::Printf("Unknown last payment type\"%s\"", lastpaymentStr.c_str())
                  );
      }
      cashflow = make_ptr(new finance::CashFlowStreamUniform
                            (
                              Date::FromExcel(issue),
                              dateFirst,
                              dateLast,
                              premium,
                              DccFromExcel(dcc),
                              FreqFromExcel(freq),
                              lastpaymentType
                            ));
    }
    else //regular uniform coupon
    {
      cashflow = make_ptr(new finance::CashFlowStreamUniform
                            (
                              Date::FromExcel(issue),
                              dateFirst,
                              dateLast,
                              premium,
                              DccFromExcel(dcc),
                              FreqFromExcel(freq)
                            ));
    }
    // and store it
    XLL_TRACE("Creating new cashflow with premium=%g, freq=%d, from %s to %s",
              premium,
              freq,
              dateFirst.Format().c_str(),
              dateLast.Format().c_str());

    return GetTheIHGAddIn().Cashflows().SetForThisCell(cashflow);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Cashflow(double contract,
              XL::Oper *xloPayments,
              const char *dcc,
              const char *freq)
{
  try
  {
    std::vector<Date> dates;
    std::vector<double> amounts;
    if ( !xloPayments || !xloPayments->AsTable(dates, amounts) )
    {
      return XL::Oper::ErrValue();
    }

    return GetTheIHGAddIn().Cashflows().SetForThisCell
           (
             new finance::CashFlowStreamGeneral
             (
                Date::FromExcel(contract),
                dates,
                amounts,
                DccFromExcel(dcc),
                FreqFromExcel(freq)
             )
           );
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33FloatingRates(double margin,
                   double startAccr,
                   double firstUnknown,
                   double lastUnknown,
                   const char *freq,
                   XL::Oper *xloDCC,
                   XL::Oper *xloFloor,
                   XL::Oper *xloCap,
                   XL::Oper *xloMultiplier,
                   XL::Oper *xloFixingDelay
                   )
{
  try
  {
    // just remember the values passed to us in a sheet object, they're going
    // to be used later from ITO33FloatCashflow()
    FloatersManager& manager = GetTheIHGAddIn().Floaters();
    FloaterData& data = manager.GetForThisCell();
    data.ResetAll();

    data.margin = margin;

    data.frequency = FreqFromExcel(freq);

    data.startAccr.SetExcel(static_cast<unsigned long>(startAccr));
    data.firstUnknown.SetExcel(static_cast<unsigned long>(firstUnknown));
    data.lastUnknown.SetExcel(static_cast<unsigned long>(lastUnknown));
    
    double 
      cap,
      floor,
      multiplier;
    
    int fixingdelay;

    std::string dcc;

    if ( xloCap && !xloCap->IsMissing() )
    {
      if (xloCap->As(&cap))
        data.cap.Set(cap);
      else
        return XL::Oper::ErrValue();
    }

    if ( xloFloor && !xloFloor->IsMissing() )
    {
      if (xloFloor->As(&floor))
        data.floor.Set(floor);
      else
        return XL::Oper::ErrValue();
    }

    if ( xloMultiplier && !xloMultiplier->IsMissing() )
    {
      if (xloMultiplier->As(&multiplier))
        data.multiplier.Set(multiplier);
      else
        return XL::Oper::ErrValue();
    }
    
    if ( xloFixingDelay && !xloFixingDelay->IsMissing() )
    {
      if (xloFixingDelay->As(&fixingdelay))
        data.fixingDelay.Set(fixingdelay);
      else
        return XL::Oper::ErrValue();
    }

    if ( xloDCC && !xloDCC->IsMissing() )
    {
      if (xloDCC->As(&dcc))
        data.dcc.Set( DccFromExcel(dcc) );
      else
        return XL::Oper::ErrValue();
    }
        
    return manager.HandleOf(data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33FloatCashflow(XL::SheetObjectHandle hfloater,
                   XL::Oper *xloPayments
                  )
{
  try
  {
    // check the parameters
    IHGAddIn& addin = GetTheIHGAddIn();

    const FloaterData& floaterData = addin.Floaters()[hfloater];

    std::vector<Date> dates;
    std::vector<double> amounts;
    if ( xloPayments &&
          !xloPayments->IsMissing() &&
            !xloPayments->AsTable(dates, amounts) )
    {
      return XL::Oper::ErrValue();
    }
         
    // construct and initialize the FloatingRates object
    shared_ptr<finance::FloatingRates> rates(new finance::FloatingRates
                                                (
                                                  floaterData.margin,
                                                  floaterData.startAccr,
                                                  floaterData.firstUnknown,
                                                  floaterData.lastUnknown,
                                                  floaterData.frequency
                                                ));
    
    if ( floaterData.dcc )
      rates->SetDayCountConvention(*floaterData.dcc);
    if ( floaterData.multiplier )
      rates->SetMultiplier(*floaterData.multiplier);
    if ( floaterData.floor )
      rates->SetFloor(*floaterData.floor);
    if ( floaterData.cap )
      rates->SetCap(*floaterData.cap);
    if ( floaterData.fixingDelay )
      rates->SetFixingDelay(*floaterData.fixingDelay);

    if ( !dates.empty() )
      rates->SetKnownPaymentStream(dates, amounts);


    // now remember it
    CashflowsManager& manager = addin.Cashflows();
    CashflowData& data = manager.GetForThisCell();
    data.SetFloater(rates);

    return manager.HandleOf(data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33OID(double grossY2M, const char *frequency, const XL::Oper *xloDCC , const XL::Oper *xloCashflow)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();
    OIDData oidData;

    oidData.SetOID(grossY2M);
    oidData.frequency = FreqFromExcel(frequency);
    if ( !xloDCC || xloDCC->IsMissing() )
      oidData.dcc = Date::DayCountConvention_ActAct;
    else
    {
      std::string dccStr;
      if (!xloDCC->As(&dccStr))
        return XL::Oper::ErrValue();
      oidData.dcc = DccFromExcel(dccStr);
    }

    if ( xloCashflow && !xloCashflow->IsMissing() )
    {
      CashflowsManager::Handle hcashflow;
      if ( !xloCashflow->As(&hcashflow) )
      {
        throw XL::EXCEPTION_MSG(xlerrValue, "Invalid cashflow handle");
      }

      CashflowData& cashflowData = addin.Cashflows()[hcashflow];
      if ( cashflowData.floaters )
      {
        throw XL::EXCEPTION_MSG
                  (
                    xlerrValue,
                    "Floating rates cashflow not supported with OIDs"
                  );
      }

      oidData.cashflow = cashflowData.cashflow;
    }
    //else: non OID zero-coupon bond

    return addin.OIDs().SetForThisCell(oidData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33CashPayToZero(double accrRate,
                   const char *frequency,
                   const XL::Oper *xloCashflow)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();
    OIDData oidData;

    oidData.SetCashPayToZero(accrRate);
    oidData.frequency = FreqFromExcel(frequency);

    if ( xloCashflow && !xloCashflow->IsMissing() )
    {
      CashflowsManager::Handle hcashflow;
      if ( !xloCashflow->As(&hcashflow) )
      {
        throw XL::EXCEPTION_MSG(xlerrValue, "Invalid cashflow handle");
      }

      CashflowData& cashflowData = addin.Cashflows()[hcashflow];
      if ( cashflowData.floaters )
      {
        throw XL::EXCEPTION_MSG
                  (
                    xlerrValue,
                    "Floating rates cashflow not supported with OIDs"
                  );
      }

      oidData.cashflow = cashflowData.cashflow;
    }
    //else: zero-coupon bond

    return addin.OIDs().SetForThisCell(oidData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// put and call provisions of convertible bonds
// ----------------------------------------------------------------------------

// used by CallsFromExcelRange() to signal an error
static void ThrowCallsFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid calls data format", format, ap);

  va_end(ap);
}

// Construct an object that stores advanced call feature 
// like notice period and trigger period 
XLLEXPORT XLOPER * 
ITO33AdvCallFeatures(int noticePeriod,
                     int triggerPeriod,
                     int triggerHistory)
{
  try
  {
    AdvCallFeatures AdvFeat;
    AdvFeat.noticePeriod = noticePeriod;
    AdvFeat.triggerPeriod = triggerPeriod;
    AdvFeat.triggerHistory = triggerHistory;

    return GetTheIHGAddIn().AdvCallFeats().SetForThisCell(AdvFeat);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// set notice preiod and trigger period for all type af call
static void SetAdvCallFeatures(shared_ptr<finance::Call> call,
                               const XL::Oper *xloAdvCallFeat)
{
  if ( !XL::Oper::IsEmpty( xloAdvCallFeat ) )
  {
    AdvCallFeaturesManager::Handle hadvcall;
    if (!xloAdvCallFeat->As(&hadvcall))
      throw XL::Oper::ErrValue();

    AdvCallFeatures AdvFeat( GetTheIHGAddIn().AdvCallFeats()[hadvcall] );
    if ( AdvFeat.noticePeriod )
      call->SetNoticePeriod(static_cast<size_t>(AdvFeat.noticePeriod));

    if ( AdvFeat.triggerPeriod )
      call->SetTriggerCheckPeriod(AdvFeat.triggerPeriod, AdvFeat.triggerHistory); 
  }
}

// construct a CallSchedule from the data in the given range, may throw
static shared_ptr<finance::CallSchedule>
CallsFromExcelRange(XL::Oper *xloCalls, double yield = 0.)
{
  // the columns in the calls range
  enum
  {
    Col_Start,
    Col_End,
    Col_Value,    // not used when using yield
    Col_Trigger,
    Col_Max
  };

  const bool hasYield = yield != 0.;
  const unsigned colMin = hasYield ? Col_End : Col_Value,
                 colMax = hasYield ? Col_Trigger : Col_Max;

  unsigned w, h;
  const XL::Oper *data = xloCalls->GetAsArray(w, h);
  if ( !data || w <= colMin || w > colMax )
  {
    ThrowCallsFormatError("bad number of columns");
  }

  // if we have the trigger column, this is a soft call, otherwise no triggers
  const bool isSoftCall = w == colMax;
  shared_ptr<finance::CallSchedule> calls(new finance::CallSchedule());
  for ( unsigned n = 0; n < h; ++n, ++data )
  {
    Date start,
         end;

    if ( !data->As(&start) || !(++data)->As(&end) )
    {
      ThrowCallsFormatError("bad date in row %u", n + 1);
    }

    shared_ptr<finance::CallPeriod> call;
    if ( hasYield )
    {
      call = finance::CallPeriod::CreateWithYield(start, end, yield);
    }
    else // must have strike
    {
      double strike;
      if ( !(++data)->As(&strike) )
        ThrowCallsFormatError("bad strike in row %u", n + 1);

      call = finance::CallPeriod::CreateWithStrike(start, end, strike);
    }

    // even when the range has Col_Max columns, we could still have some
    // hard calls as well, check for this
    if ( isSoftCall && !(++data)->IsMissing() )
    {
      double trigger;
      if ( !data->As(&trigger) )
        ThrowCallsFormatError("bad trigger in row %u", n + 1);

      call->SetTrigger(trigger);
    }

    calls->AddCallPeriod(call);
  }

  return calls;
}

// common implementation of ITO33Calls() and ITO33CallsByYield()
static XLOPER *
CreateCallSchedule(XL::Oper *xloCalls,
                   double yield,
                   XL::Oper *xloKeepAccr,
                   bool forfeitCoupon,
                   const char *makeWholeType,
                   double makeWholePremium,
                   const XL::Oper *xloTriggerPercentOf,
                   const XL::Oper *xloAdvCallFeat)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    if ( !xloCalls )
      return XL::Oper::ErrValue();

    // construct the call schedule from the calls data
    shared_ptr<finance::CallSchedule>
      calls(CallsFromExcelRange(xloCalls, yield));

    // deal with the rest of the parameters
    if ( xloKeepAccr && !xloKeepAccr->IsMissing() )
    {
      bool keepAccrued;
      if ( !xloKeepAccr->As(&keepAccrued) )
        return XL::Oper::ErrValue();

      calls->SetKeepAccrued(keepAccrued);
    }
    //else: default is fine, it's true by default

    calls->SetForfeitCoupon(forfeitCoupon);

    if ( !XL::Oper::IsEmpty(xloTriggerPercentOf) )
    {
      std::string triggerPercentStr;
      if ( !xloTriggerPercentOf->As(&triggerPercentStr ) )
        return XL::Oper::ErrValue();
      calls->SetTriggerAsPercentageOf (TriggerPercentFromExcel(triggerPercentStr) );
    }
    
    SetAdvCallFeatures(calls, xloAdvCallFeat);
    
    bool isCouponPresentValue;
    switch ( MakeWholeTypeFromExcel(makeWholeType, &isCouponPresentValue) )
    {
      case finance::MakeWholeType_Coupon:
        calls->SetCouponMakeWhole(isCouponPresentValue);
        break;

      case finance::MakeWholeType_Premium:
        calls->SetPremiumMakeWhole(makeWholePremium);

        // make the check below pass
        makeWholePremium = 0.;
        break;

      case finance::MakeWholeType_Max:
        // no make whole protection
        break;

      default:
        FAIL("unknown make whole type value");
    }

    if ( makeWholePremium != 0. )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  "Premium is ignored if the make whole type is not \"P\""
                );
    }

    // now store the call schedule
    return addin.Calls().SetForThisCell(calls);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Calls(XL::Oper *xloCalls,
           XL::Oper *xloKeepAccr,
           bool forfeitCoupon,
           const char *makeWholeType,
           double makeWholePremium,
           XL::Oper *xloTriggerPercentOf,
           XL::Oper *xloAdvCallFeat)
{
  return CreateCallSchedule(xloCalls, 0., xloKeepAccr, forfeitCoupon,
                            makeWholeType, makeWholePremium, xloTriggerPercentOf, xloAdvCallFeat);
}

XLLEXPORT XLOPER *
ITO33CallsByYield(XL::Oper *xloCalls,
                  double yield,
                  XL::Oper *xloKeepAccr,
                  bool forfeitCoupon,
                  const char *makeWholeType,
                  double makeWholePremium,
                  XL::Oper *xloTriggerPercentOf,
                  XL::Oper *xloAdvCallFeat)
{
  if ( yield == 0. )
  {
    // must have the yield here
    return XL::Oper::ErrValue();
  }

  return CreateCallSchedule(xloCalls, yield, xloKeepAccr, forfeitCoupon,
                            makeWholeType, makeWholePremium, xloTriggerPercentOf, xloAdvCallFeat);
}

// construct a PutSchedule from the data in the given range, may throw
static shared_ptr<finance::PutSchedule>
PutsFromExcelRange(XL::Oper *xloPuts, double yield = 0.)
{
  // the columns in the puts range
  enum
  {
    Col_Date,
    Col_Value,
    Col_Max
  };

  bool hasYield = yield != 0.;
  
  // check the parameters and create the put schedule object
  unsigned w, h;
  const XL::Oper *data = xloPuts->GetAsArray(w, h);
  if ( !data)
  {
    data = xloPuts;
    w = h = 1;
  }
   
  if ( XL::Oper::IsMissing( data ) || 
      ( w != Col_Max && !hasYield ) || 
      ( w != Col_Max - 1 && hasYield  ) )
  {
    throw XL::EXCEPTION_MSG(xlerrRef, "Invalid put data format");
  }

  shared_ptr<finance::PutSchedule> puts(new finance::PutSchedule());
  for ( unsigned n = 0; n < h; ++n, ++data )
  {
    Date date;
    double strike;

    if ( !data->As(&date) )
    {
      throw XL::EXCEPTION_MSG(xlerrValue, "Invalid put data format");
    }

    if (!hasYield )
    {
      if ( !(++data)->As(&strike) )
        throw XL::EXCEPTION_MSG(xlerrValue, "Invalid put data format");
      puts->AddPutWithStrike(date, strike);
    }
    else
      puts->AddPutWithYield(date, yield);
  }

  return puts;
}

// common implementation of ITO33Puts() and ITO33PutsByYield()
static XLOPER *
CreatePutSchedule(XL::Oper *xloPuts,double yield, XL::Oper *xloKeepAccr, bool forfeitCoupon)
{
  try
  {
    if ( !xloPuts )
      return XL::Oper::ErrValue();

    // construct the put schedule from the puts data
    shared_ptr<finance::PutSchedule>
      puts(PutsFromExcelRange(xloPuts, yield));
  
    if ( xloKeepAccr && !xloKeepAccr->IsMissing() )
    {
      bool keepAccrued;
      if ( !xloKeepAccr->As(&keepAccrued) )
        return XL::Oper::ErrValue();

      puts->SetKeepAccrued(keepAccrued);
    }
    //else: default is fine, it's true by default

    puts->SetForfeitCoupon(forfeitCoupon);

    // now store the put schedule
    return GetTheIHGAddIn().Puts().SetForThisCell(puts);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Puts(XL::Oper *xloPuts, XL::Oper *xloKeepAccr, bool forfeitCoupon)
{
  return CreatePutSchedule(xloPuts,0.,xloKeepAccr, forfeitCoupon);
}

XLLEXPORT XLOPER *
ITO33PutsByYield(XL::Oper *xloPuts, double yield, XL::Oper *xloKeepAccr, bool forfeitCoupon)
{
  if ( yield == 0. )
  {
    // must have the yield here
    return XL::Oper::ErrValue();
  }

  return CreatePutSchedule(xloPuts,yield,xloKeepAccr, forfeitCoupon);
}

// ----------------------------------------------------------------------------
// convertible bonds conversions
// ----------------------------------------------------------------------------

// used by ConvsFromExcelRange() and MakeConversion() to signal an error
static void ThrowConvsFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid conversions data format", format, ap);

  va_end(ap);
}

static finance::CoCoType CoCoTypeFromExcel(const char *coco)
{
  finance::CoCoType cocoType = finance::CoCoType_Max;
  
  static const EnumValuesNames<finance::CoCoType> cocoTypes[] =
  {
    {"AP",    finance::CoCoType_CheckAnyTimeAndConvertAsOfCheckDate},
    {"AF",    finance::CoCoType_CheckAnyTimeAndConvertOnCheckDate},
    {"QC",    finance::CoCoType_CheckQuarterlyAndConvertAsOfCheckDate},
    {"QQ",    finance::CoCoType_CheckQuarterlyAndConvertDuringNextQuarter}
  };
  
  if (!RestoreEnumValueFromName(coco,cocoType,SIZEOF(cocoTypes),cocoTypes))
    ThrowConvsFormatError("unknown coco type \"%s\"", coco);

  return cocoType;
}

static inline finance::CoCoType CoCoTypeFromExcel(std::string coco)
{
  return CoCoTypeFromExcel(coco.c_str() );
}

// create a conversion period from the provided data, used by
// ConvsFromExcelRange() and ITO33Conversion()
static shared_ptr<finance::ConversionPeriod>
MakeConversion(const Date& start, const Date& end, double ratio,
               const XL::Oper *xloCoCoType, const XL::Oper *xloTrigger,
               const XL::Oper *xloTrigRate, const XL::Oper *xloTrigEnd,
               const XL::Oper *xloCash)
{
  shared_ptr<finance::ConversionPeriod> conv(
      new finance::ConversionPeriod(start, end, ratio));

  // is this a coco?
  if ( !XL::Oper::IsMissing(xloCoCoType) )
  {
    std::string cocoTypeStr;
    if ( !xloCoCoType->As(&cocoTypeStr) )
      ThrowConvsFormatError("invalid coco type");

    finance::CoCoType cocoType = CoCoTypeFromExcel(cocoTypeStr);
    
    double trigger = 0;
    if ( !xloTrigger || !xloTrigger->As(&trigger) )
      ThrowConvsFormatError("invalid trigger value");

    // does the trigger vary?
    const bool constTrig = XL::Oper::IsMissing(xloTrigRate);
    if ( constTrig != XL::Oper::IsMissing(xloTrigEnd) )
    {
      ThrowConvsFormatError("trigger rate and ceiling must be both present or not");
    }

    if ( constTrig )
    {
      double trigRate,
             trigEnd;
      if ( !xloTrigRate->As(&trigRate) )
        ThrowConvsFormatError("invalid trigger rate");

      if ( !xloTrigEnd->As(&trigEnd) )
        ThrowConvsFormatError("invalid trigger ceiling");

      conv->SetCoCo(trigger, cocoType, trigRate, trigEnd);
    }
    else // trigger doesn't change, use 0 change rate
    {
      conv->SetCoCo(trigger, cocoType, 0, trigger);
    }
  }
  else // trigger-related stuff doesn't make sense if not a CoCo
  {
    if ( !XL::Oper::IsMissing(xloTrigger) ||
            !XL::Oper::IsMissing(xloTrigRate) ||
              !XL::Oper::IsMissing(xloTrigEnd) )
    {
      ThrowConvsFormatError("trigger parameters only valid for cocos");
    }
  }

  // do we have cash?
  if ( xloCash && !xloCash->IsMissing() )
  {
    double cash;
    if ( !xloCash->As(&cash) )
      ThrowConvsFormatError("invalid cash value");

    conv->SetCash(cash);
  }

  return conv;
}
// used by ConvsFromExcelRange() and MakeConversion() to signal an error
static void ThrowResetFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid reset data format", format, ap);

  va_end(ap);
}
// create a conversion price reset from the provided data, used by
// ITO33Reset()
static shared_ptr<finance::ConversionPriceReset>
MakeReset(const Date& resetDate, double floor,
          const XL::Oper *xloMultiplier, const XL::Oper *xloCap
          )
{
  shared_ptr<finance::ConversionPriceReset> reset (
    new finance::ConversionPriceReset( resetDate, floor ) );
  
  if ( xloMultiplier && !xloMultiplier->IsMissing() )
  {
    double dMultiplier;
    if (!xloMultiplier->As(&dMultiplier))
      ThrowResetFormatError("error on the multiplier");
    reset->SetMultiplier(dMultiplier);
  }
  if ( xloCap && !xloCap->IsMissing() )
  {
    double dCap;
    if (!xloCap->As(&dCap) )
      ThrowResetFormatError("error on the cap");
    reset->SetCap(dCap);
  }
  return reset; 
}

static shared_ptr<finance::ConversionSchedule>
ConvsFromExcelRange(XL::Oper *xloConvs)
{
  // the columns in the conversions range
  enum
  {
    Col_Start,
    Col_End,
    Col_Ratio,
    Col_CoCoType,
    Col_Trigger,
    Col_TrigRate,
    Col_TrigEnd,
    Col_Cash,
    Col_Max
  };

  // check the parameters and create the conversions schedule object
  unsigned w, h;
  const XL::Oper *data = xloConvs->GetAsArray(w, h);
  if ( !data || (w != Col_CoCoType && w != Col_TrigRate &&
                  w != Col_Cash && w != Col_Max) )
  {
    ThrowConvsFormatError("invalid number of columns in conversions data range");
  }

  const bool hasCoCo = w > Col_CoCoType;
  const bool hasTrigRate = w > Col_TrigRate;
  const bool hasCash = w > Col_Cash;
  shared_ptr<finance::ConversionSchedule> convs(new finance::ConversionSchedule());
  for ( unsigned n = 0; n < h; ++n, ++data )
  {
    Date start, end;
    double ratio;

    if ( !data->As(&start) || !(++data)->As(&end) )
      ThrowConvsFormatError("invalid conversion dates");

    if ( !(++data)->As(&ratio) )
      ThrowConvsFormatError("invalid conversion ratio");

    const XL::Oper *xloCoCoType,
                   *xloTrigger;
    if ( hasCoCo )
    {
      xloCoCoType = ++data;
      xloTrigger = ++data;
    }
    else // no CoCo fields
    {
      xloCoCoType =
      xloTrigger = NULL;
    }

    const XL::Oper *xloTrigRate,
                   *xloTrigEnd;
    if ( hasTrigRate )
    {
      xloTrigRate = ++data;
      xloTrigEnd = ++data;
    }
    else // trigger doesn't vary
    {
      xloTrigRate =
      xloTrigEnd = NULL;
    }

    const XL::Oper *xloCash = hasCash ? ++data : NULL;

    convs->AddConversionPeriod(MakeConversion(start, end, ratio,
                                              xloCoCoType, xloTrigger,
                                              xloTrigRate, xloTrigEnd,
                                              xloCash));
  }

  return convs;
}

// helper of ITO33Conversion() and ITO33Conversions(): sets keep accrued flag for
// the conversions object which is not completely trivial because it is true by
// default in IHG but false if omitted in XL so we need to only set it if it
// was really specified as false
static void
ConvSetKeepAccrued(finance::Conversion *convs,
                   const XL::Oper *xloKeepAccr)
{
  if ( xloKeepAccr && !xloKeepAccr->IsMissing() )
  {
    bool keepAccrued;
    if ( !xloKeepAccr->As(&keepAccrued) )
    {
      ThrowConvsFormatError("invalid keep accrued flag value");
    }

    convs->SetKeepAccrued(keepAccrued);
  }
  //else: default is fine, it's true by default
}

XLLEXPORT XLOPER *
ITO33Conversion(double start, double end, double ratio,
                const XL::Oper *xloCoCoType, const XL::Oper *xloTrigger,
                const XL::Oper *xloTrigRate, const XL::Oper *xloTrigEnd,
                const XL::Oper *xloCash,
                const XL::Oper *xloKeepAccr,
                bool forfeitCoupon)
{
  try
  {
    shared_ptr<finance::ConversionPeriod>
      conv(MakeConversion(Date::FromExcel(start), Date::FromExcel(end), ratio,
                          xloCoCoType, xloTrigger,
                          xloTrigRate, xloTrigEnd,
                          xloCash));

    // create the conversion schedule consisting of this single conversion
    shared_ptr<finance::ConversionSchedule> convs(new finance::ConversionSchedule());
    convs->AddConversionPeriod(conv);

    ConvSetKeepAccrued(convs.get(), xloKeepAccr);

    convs->SetForfeitCoupon(forfeitCoupon);

    // now store the conversion schedule
    return GetTheIHGAddIn().Conversions().SetForThisCell(convs);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Conversions(XL::Oper *xloConvs, XL::Oper *xloKeepAccr,
                 bool forfeitCoupon,bool triggerCondMet)
{
  try
  {
    if ( !xloConvs )
      return XL::Oper::ErrValue();

    // construct the conversions schedule from the conversion data
    shared_ptr<finance::ConversionSchedule> convs(ConvsFromExcelRange(xloConvs));

    ConvSetKeepAccrued(convs.get(), xloKeepAccr);

    convs->SetForfeitCoupon(forfeitCoupon);

    convs->SetIsLastTriggerConditionMet(triggerCondMet);
    // now store the conversion schedule
    return GetTheIHGAddIn().Conversions().SetForThisCell(convs);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Reset(double start, double end,double initialConvPrice,
           double currentConvPrice,const char* flooredBy,
           XL::Oper *xloReset,XL::Oper *xloCash,
           XL::Oper *xloKeepAccr,bool forfeitCoupon)
{
  try
  {
    if (!xloReset )
      return XL::Oper::ErrValue();
    
    finance::ResetFlooredBy flooredByType;

    static const EnumValuesNames<finance::ResetFlooredBy> flooredTypes[] =
    {
      { "prevailing",   finance::ResetFlooredBy_PrevailingConversionPrice   },
      { "initial",  finance::ResetFlooredBy_InitialConversionPrice  },
    };

    if ( !GetEnumFromAbbrevName<SIZEOF(flooredTypes)>
          (
            &flooredByType,
            flooredBy,
            flooredTypes
          ) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown reset floored by type \"%s\"", flooredBy)
                );
    }

     // construct a reset schedule
    shared_ptr<finance::ResetConversionSchedule> reset (
      new finance::ResetConversionSchedule(Date::FromExcel(start),
                                           Date::FromExcel(end),
                                           initialConvPrice,
                                           currentConvPrice,
                                           flooredByType) );
    // the columns in the conversions range
    enum
    {
      Col_ResetDate,
      Col_Floor,
      Col_Multiplier,
      Col_Cap,
      Col_Max
    };

    unsigned w,h;
    const XL::Oper *data = xloReset->GetAsArray(w, h);
    if (w < Col_Multiplier || w > Col_Max )
      ThrowResetFormatError("invalid number of columns in conversion price reset data range");
    
    bool hasMultiplier = w > Col_Multiplier;
    bool hasCap = w > Col_Cap;
    for ( unsigned n = 0; n < h; ++n, ++data )
    {
      Date resetDate;
      double floor;
      if ( data->IsMissing() || !data->As(&resetDate))
        ThrowResetFormatError("invalid reset date");
      if (!(++data)->As(&floor) || data->IsMissing() )
        ThrowResetFormatError("invalid floor rate format");

      const XL::Oper *xloMultiplier,
               *xloCap;
      if (hasMultiplier)
      {
        xloMultiplier = ++data;
      }
      else
      {
        xloMultiplier = NULL;
      }
      if (hasCap)
      {
        xloCap = ++data;
      }
      else
        xloCap =NULL;
      reset->AddConversionPriceReset( MakeReset(resetDate,floor,xloMultiplier,xloCap) );
    }
    ConvSetKeepAccrued(reset.get(), xloKeepAccr);

    reset->SetForfeitCoupon(forfeitCoupon);
    
    //cash upon conversion
    if ( xloCash && !xloCash->IsMissing() )
    {
      double cash;
      if ( !xloCash->As(&cash) )
        ThrowResetFormatError("invalid cash value");
      reset->SetCash(cash);
    }
    // now store the reset conversion schedule
    return GetTheIHGAddIn().Resets().SetForThisCell(reset);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33AttachedWarrant(double start, double end,double baseRatio,
                     double incrementalFactor,
                     XL::Oper *xloCapRatio,
                     XL::Oper *xloResetDate,
                     XL::Oper *xloFixedStrike,
                     XL::Oper *xloCurrentRatio,
                     XL::Oper *xloKeepAccrued,
                     bool forfeitcoupon)
{
  try
  {
    shared_ptr<finance::ShareDependentConversion>
      dependentConv (new finance::ShareDependentConversion(Date::FromExcel(start),
                                                           Date::FromExcel(end),
                                                           baseRatio,
                                                           incrementalFactor) );

    if (xloCapRatio && !xloCapRatio->IsMissing() )
    {
      double cap;
      if ( !xloCapRatio->As(&cap))
        return XL::Oper::ErrValue();
      dependentConv->SetCapRatio(cap);
    }
    
    if (xloResetDate && !xloResetDate->IsMissing() )
    {
      Date resetDate;
      if ( !xloResetDate->As(&resetDate))
        return XL::Oper::ErrValue();
      dependentConv->SetResetDate(resetDate);
    }
    
    if (xloFixedStrike && !xloFixedStrike->IsMissing() )
    {
      double fixedStrike;
      if ( !xloFixedStrike->As(&fixedStrike))
        return XL::Oper::ErrValue();
      dependentConv->SetFixedStrike(fixedStrike);
    }
    
    if (xloCurrentRatio && !xloCurrentRatio->IsMissing() )
    {
      double ratio;
      if ( !xloCurrentRatio->As(&ratio))
        return XL::Oper::ErrValue();
      dependentConv->SetCurrentRatio(ratio);
    }

    if (xloCurrentRatio && !xloCurrentRatio->IsMissing() )
    {
      double ratio;
      if ( !xloCurrentRatio->As(&ratio))
        return XL::Oper::ErrValue();
      dependentConv->SetCurrentRatio(ratio);
    }

    ConvSetKeepAccrued(dependentConv.get(), xloKeepAccrued);

    dependentConv->SetForfeitCoupon(forfeitcoupon);
    
    //now store the sharedependent conversion
    return GetTheIHGAddIn().ShareDepConversions().SetForThisCell(dependentConv);

  }
  catch( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33AttachedWarrantWithCoCo(XL::SheetObjectHandle hatWarrant,
                             char * cocoType,
                             double trigger,
                             XL::Oper *xloTrigRate,
                             XL::Oper *xloTrigEnd,
                             bool isTrigMet)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    shared_ptr<finance::ShareDependentConversion> 
      dependentConv = addin.ShareDepConversions()[hatWarrant];

    shared_ptr<finance::ShareDependentConversion> 
      cocoDepConv( new finance::ShareDependentConversion(*dependentConv) );
    
    double trigRate;

    if ( xloTrigRate && !xloTrigRate->IsMissing() )
    {
      if ( !xloTrigRate->As(&trigRate))
        ThrowConvsFormatError("invalid trigger rate");
    }
    else
      trigRate = 0;
    
    double trigEnd;
    if ( xloTrigEnd && !xloTrigEnd->IsMissing() )
    {
      if ( !xloTrigEnd->As(&trigEnd))
        ThrowConvsFormatError("invalid trigger ceil");
    }
    else
      trigEnd = trigger;
    
    if ( XL::Oper::IsMissing(xloTrigRate) != XL::Oper::IsMissing(xloTrigEnd) )
      ThrowConvsFormatError("trigger rate and Trigger ceil must be set both or not");

    cocoDepConv->SetCoCo(trigger,CoCoTypeFromExcel(cocoType),trigRate,trigEnd,isTrigMet);

    return addin.ShareDepConversions().SetForThisCell(cocoDepConv);
  }
  catch( ... )
  {
    return XL::ReturnException();
  }

}
// ----------------------------------------------------------------------------
// equities and derivatives creation
// ----------------------------------------------------------------------------

// used by ITO33Equity()
static void ThrowDividendsFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid dividends data format", format, ap);

  va_end(ap);
}

XLLEXPORT XLOPER *
ITO33Equity(XL::SheetObjectHandle hcurrency,
            double borrow,
            XL::Oper *xloDividends,
            int fiscalYearStartMonth,
            int fiscalYearStartDay)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    // check the currency data handle
    const CurrencyData& currencyData = addin.Currencies()[hcurrency];

    // deal with the dividends
    shared_ptr<finance::Dividends> dividends;
    if ( xloDividends && !xloDividends->IsMissing() )
    {
      enum
      {
        Col_Date,
        Col_Amount,
        Col_Type,
        Col_Max
      };

      unsigned w, h;
      const XL::Oper *data = xloDividends->GetAsArray(w, h);
      if ( !data || (w != Col_Type && w != Col_Max) )
        ThrowDividendsFormatError("invalid number of columns in the range");

      const bool hasType = w > Col_Type;
      dividends = make_ptr( new finance::Dividends );
      for ( unsigned n = 0; n < h; ++n, ++data )
      {
        Date date;
        double value;
        std::string s;
        if ( !data->As(&date) )
          ThrowDividendsFormatError("invalid date in row %u", n + 1);

        if ( !(++data)->As(&value) )
          ThrowDividendsFormatError("invalid amount in row %u", n + 1);

        finance::Dividend::Type type;
        if ( hasType )
        {
          std::string typeStr;
          if ( !(++data)->As(&typeStr) )
            ThrowDividendsFormatError("invalid type format in row %u", n + 1);
          
          type =  DividendTypeFromExcel( typeStr );
          
        }
        else // by default dividends are in %
        {
          type = finance::Dividend::Yield;
        }

        dividends->Add(type, date, value);
      }
    }

    // by default fiscal year starts at January 1
    if ( !fiscalYearStartMonth )
      fiscalYearStartMonth = 1;
    if ( !fiscalYearStartDay )
      fiscalYearStartDay = 1;

    // create the temporary issuer (we only use it internally)
    const Date fiscalYearStart(2006,
                               static_cast<Date::Month>(fiscalYearStartMonth),
                               fiscalYearStartDay);
    shared_ptr<finance::Issuer> issuer(new finance::Issuer);
    issuer->SetFiscalYearStartDate(fiscalYearStart);


    // finally create the equity object itself
    XLL_TRACE("Creating new equity for currency \"%s\"",
               currencyData.currency->GetCode().c_str());

    EquityData& equityData = addin.Equities().GetForThisCell();
    equityData.currency = currencyData.currency;
    equityData.yieldcurve = currencyData.yieldcurve;

    // NB: we must use a valid spot here but we don't have it, so just invent
    //     something ridiculously small to keep IHG happy
    finance::Equity *equity = new finance::Equity(0.0001, currencyData.currency);
    //Setting dividends
    if ( dividends )
      equity->SetDividends( dividends );
    equityData.equity = make_ptr(equity);
    equity->SetIssuer(issuer);
    equity->SetBorrowCurve(
        shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(borrow)));

    return addin.Equities().HandleOf(equityData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

//Test if the last "Freq number" of dividends hold in a year 
static bool IsUnregularDividend(finance::Frequency Freq,const XL::Oper *xloDividends, unsigned weight)
{
  ito33::Date Previousdate;
  const XL::Oper * data = xloDividends;
  ito33::Date limitdate;
  for ( unsigned n = 0; n < static_cast<unsigned>( Freq ) ; ++n ,data += weight)
  {

    ito33::Date date;
    if ( !data->As( &date ) )
      //TO BE CHANGED    
      ThrowDividendsFormatError(
        "invalid date %u rows before the end of range",
        static_cast<unsigned>( Freq ) - ( n + 1) );
    
    if ( Previousdate.IsValid() && Previousdate >= date )
      ThrowDividendsFormatError(
        "non increasing date in the %u last input date",
        static_cast<unsigned>( Freq ) );
        
    if (!n) 
    {
      limitdate = date;
      limitdate.AddYears(1);      
    }
    
    if ( date >= limitdate ) return false;
    Previousdate = date;    
  }


  return true;
}

XLLEXPORT XLOPER *
ITO33GenerateDividends(XL::Oper *xloDividends,double end,
                        double accretion, const char *frequency,
                        XL::Oper * xloIsAdjusted)
{
  try
  {
    if ( !xloDividends || xloDividends->IsMissing() )
      ThrowDividendsFormatError("null range is not handdle");

    enum
    {
      Col_Date,
      Col_Amount,
      Col_Type,
      Col_Max
    };

    unsigned w, h;
    const XL::Oper 
      *data = xloDividends->GetAsArray(w, h),
      *dataInit = data;
    if ( !data || (w != Col_Max && w!= Col_Type) )
      ThrowDividendsFormatError("invalid number of columns in the range");

    const bool hasType = w > Col_Type;       
    ito33::Date 
      StartDate ,
      EndDate ( Date::FromExcel( end ) );
    
    finance::Frequency Freq = FreqFromExcel( frequency );
    bool bIsUnregularDiv(false);
    int MonthInc = 12 / static_cast<int> (Freq);

    double amount;
    std::string typeStr;

    //Make it point on the last dividend data
    data = dataInit + w * (h - 1) ;

    std::vector<ito33::Date> dates;
    std::vector<double> values;
    std::vector<std::string> typesStr;
    
    bool bIsAdjusted(false);
    if ( xloIsAdjusted && !xloIsAdjusted->IsMissing() )
      xloIsAdjusted->As(&bIsAdjusted);

    if ( bIsAdjusted && h >= static_cast<unsigned> (Freq) ) 
      bIsUnregularDiv = IsUnregularDividend(Freq, data - (Freq - 1)*w , w);
    
    if ( bIsAdjusted != bIsUnregularDiv )
      ThrowDividendsFormatError("input data incompatible with the last parameter");

    if (!bIsAdjusted)
    {
     
      //Get the last dividend input data
      if ( !data->As( &StartDate ) )
        ThrowDividendsFormatError("invalid date in last row");
      if ( !(++data)->As( &amount ) )
        ThrowDividendsFormatError("invalid cash amount in last row");

      if ( hasType )
      {

        if ( !(++data)->As( &typeStr ) )
          ThrowDividendsFormatError("invalid string in last row");
    
        //Restore the full name of typeStr from the abreviation
        //Throw if typeStr is not a valid string  
        typeStr = (DividendTypeFromExcel(typeStr) != finance::Dividend::Cash)? "yield":"cash" ;
      }
      else
      {
        typeStr = "yield";
      }
              
      int i =0;
      ito33::Date::Year_t PreviousYear = StartDate.GetYear();
      StartDate.AddMonths(MonthInc);
      //dividend accretes only after the first payment of a year
      for(ito33::Date PaymentDate = StartDate; PaymentDate <= EndDate; PaymentDate.AddMonths(MonthInc),i++)
      {
        if (PreviousYear != PaymentDate.GetYear() ) amount *= 1 + accretion ;
        dates.push_back( PaymentDate );
        typesStr.push_back( typeStr );
        values.push_back( amount );
        PreviousYear = PaymentDate.GetYear();
      }
    }
    else
    {

      data = data - (Freq -1)*w;
      ito33::Date DivDate;
      //Loop for the first years of generated dividends
      for( unsigned n = 0; n < static_cast<unsigned>( Freq ) ; ++n, ++data )
      {
        if ( !data->As( &DivDate ) )
          ThrowDividendsFormatError(
            "invalid date %u rows before the end of range",
            static_cast<unsigned>( Freq ) - ( n + 1) );
        if ( !(++data)->As( &amount ) )
          ThrowDividendsFormatError(
            "invalid cash amount %u rows before the end of range",
            static_cast<unsigned>( Freq ) - ( n + 1) );
        if ( hasType)
        {
          if ( !(++data)->As( &typeStr ) )
            ThrowDividendsFormatError(
              "invalid string %u rows before the end of range",
              static_cast<unsigned>( Freq ) - ( n + 1) );
          //Restore the full name of typeStr from the abreviation
          //Throw if typeStr is not a valid string  
          typeStr = (DividendTypeFromExcel(typeStr) != finance::Dividend::Cash)? "yield":"cash" ;
        }
        else
        {
          typeStr = "yield";
        }
        amount *= 1 + accretion ;
        DivDate.AddYears(1);
        dates.push_back( DivDate );
        typesStr.push_back( typeStr );
        values.push_back( amount );
      }
             
      std::vector<ito33::Date>::const_iterator iterDate  = dates.begin();        
      std::vector<double>::const_iterator iterAmount = values.begin();
      std::vector<std::string>::const_iterator iterTypeStr = typesStr.begin();
      
      DivDate = *iterDate;
      //Loop for all the remaining dividends we have to generated      
      for(DivDate.AddYears(1); DivDate <= EndDate; DivDate.AddYears(1))
      {    
        typeStr = *iterTypeStr;
        amount = *iterAmount * ( 1 + accretion );
        
        dates.push_back( DivDate );
        typesStr.push_back( typeStr );
        values.push_back( amount );
        
        iterDate = dates.end() - Freq;
        iterAmount = values.end() - Freq;
        iterTypeStr = typesStr.end() - Freq;
        DivDate = *iterDate;
      }     
    }
    // and return the result as an array
    const unsigned size = static_cast<unsigned>(dates.size());
    XL::Oper *array = new XL::Oper[3*(size + h)];
    
    XL::Oper *p = array;

    // 1st loop: p hold the input dividends 
    for ( unsigned n = 0; n < h ; ++n,++dataInit )
    {
      *p++ = *dataInit;
      *p++ = *(++dataInit);
      if (hasType)
      {
        std::string str;
        if (!(++dataInit)->As(&str)) 
          ThrowDividendsFormatError("invalid dividend type in the input dividend range");
        *p++ = (DividendTypeFromExcel(str) != finance::Dividend::Cash)? "yield":"cash" ;
      }
      else
      {
        *p++ = "yield";
      }
    }
    //2nd loop: p hold the generated dividends
    for ( unsigned n = 0; n < size; ++n )
    {
      *p++ = dates[n];
      *p++ = values[n];
      *p++ = typesStr[n];
    }

    return XL::Oper::ReturnArray(3, size + h , array);    
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Option(XL::SheetObjectHandle hequity,
            double strike,
            double maturity,
            const char *type,
            const char *exercise)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    const EquityData& equityData = addin.Equities()[hequity];

    finance::OptionType optType;
    if ( !GetEnumFromAbbrevName<SIZEOF(g_optionTypes)>
          (
            &optType,
            type,
            g_optionTypes
          ) )
      return XL::Oper::ErrValue();

    finance::ExerciseType optExe;
    if ( !GetEnumFromAbbrevName<SIZEOF(g_exerciseTypes)>
          (
            &optExe,
            exercise,
            g_exerciseTypes
          ) )
      return XL::Oper::ErrValue();


    const Date dateMat(Date::FromExcel(maturity));

    XLL_TRACE("Creating new option with strike=%g, maturity %s, %s, %s",
              strike, dateMat.Format().c_str(), type, exercise);

    DerivativeData& derivData = addin.Derivatives().GetForThisCell();
    derivData.InitFrom(equityData);
    derivData.IsCBOption = false;
    derivData.derivative = make_ptr( new finance::Option
                                         (
                                           strike,
                                           dateMat,
                                           optType,
                                           optExe
                                         ) );

    return addin.Derivatives().HandleOf(derivData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33CDS(XL::SheetObjectHandle hequity,
         XL::SheetObjectHandle hcashflow,
         XL::Oper *xloRecovery)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    // check/dereference the handles
    const EquityData& equityData = addin.Equities()[hequity];
    const CashflowData& cashflowData = addin.Cashflows()[hcashflow];

    const shared_ptr<finance::CashFlowStreamUniform> cashflowUniform =
      dynamic_pointer_cast<finance::CashFlowStreamUniform>(cashflowData.cashflow);
    if ( !cashflowUniform )
      return XL::Oper::ErrValue();


    // check the recovery parameter if any
    double recovery;
    if ( !xloRecovery || xloRecovery->IsMissing() )
    {
      // default value for most of the CDS
      recovery = 0.40;
    }
    else if ( !xloRecovery->As(&recovery) )
    {
      return XL::Oper::ErrValue();
    }


    // do create the new CDS
    XLL_TRACE("Creating new CDS with recovery=%g%%", recovery*100.);

    DerivativesManager& manager = addin.Derivatives();
    DerivativeData& derivData = manager.GetForThisCell();
    derivData.InitFrom(equityData);
    derivData.IsCBOption = false;
    derivData.derivative = make_ptr( new finance::CDS
                                         (recovery, cashflowUniform) );

    return manager.HandleOf(derivData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

/*XLLEXPORT XLOPER *
ITO33ReferenceCDS(int duration,
                  char *units,
                  double spread,
                  XL::Oper *xloDate,
                  XL::Oper *xloRecovery,
                  XL::Oper *xloFreq,
                  XL::Oper *xloDCC)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    // check the recovery parameter if any
    double recovery;
    if ( !xloRecovery || xloRecovery->IsMissing() )
    {
      // default value for most of the CDS
      recovery = 0.40;
    }
    else if ( !xloRecovery->As(&recovery) )
    {
      return XL::Oper::ErrValue();
    }

    // check the frequency parameter if any
    finance::Frequency freq;
    std::string freqStr;
    
    if ( !xloFreq || xloFreq->IsMissing() )
    {
      // default frequency for most of the CDS
      freq = finance::Frequency_Quarterly;
    }
    else 
    {
      if ( !xloFreq->As(&freqStr) )
        return XL::Oper::ErrValue();
      else      
        freq =FreqFromExcel(freqStr);
    }
    
    

    // check the frequency parameter if any
    Date::DayCountConvention dcc;
    std::string dccStr;
    
    if ( !xloDCC || xloDCC->IsMissing() )
    {
      // default frequency for most of the CDS
      dcc = Date::DayCountConvention_30360;
    }
    else 
    {
      if ( !xloDCC->As(&dccStr) )
        return XL::Oper::ErrValue();
      else      
        dcc =DccFromExcel(dccStr);
    }

    switch ( TimeUnitFromExcel(units) )
    {
      case finance::TimeUnit_Day:
        throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("maturity in days not allowed")
          );
        // break; -- unnecessary after throw
      
      case finance::TimeUnit_Month:
        // nothing done, duration is good
        break;
      case finance::TimeUnit_Year:
        duration *= 12;
        break;

      default:
        FAIL("unknown maturity unit value");
    }

    DerivativesManager& manager = addin.Derivatives();
    DerivativeData& derivData = manager.GetForThisCell();
    
    //derivData.InitFrom(equityData);
    Date dateVal = ValDateFromExcel(xloDate);

    derivData.derivative = new finance::ReferenceCDS( dateVal,
                                                     duration,
                                                     spread,
                                                     dcc,
                                                     freq,
                                                     recovery);

    return manager.HandleOf(derivData);
    
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}*/

XLLEXPORT XLOPER *
ITO33BondTerms(double issue,
               double maturity,
               double iprice,
               double nominal,
               double redemption,
               double recovery,
               XL::Oper *xloCashflow)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    const Date dateIssue(Date::FromExcel(issue));
    const Date dateMat(Date::FromExcel(maturity));

    XLL_TRACE("Creating new bond terms %s..%s with nominal %g",
              dateIssue.Format().c_str(),
              dateMat.Format().c_str(),
              nominal);

    shared_ptr<finance::BondTerms> terms
                                  (
                                    new finance::BondTerms
                                        (
                                          dateIssue,
                                          iprice,
                                          dateMat,
                                          nominal,
                                          redemption,
                                          recovery
                                        )
                                  );

    if ( xloCashflow && !xloCashflow->IsMissing() )
    {
      CashflowsManager::Handle hcashflow;
      if ( !xloCashflow->As(&hcashflow) )
        return XL::Oper::ErrValue();

      const CashflowData * const
        cashflowData = addin.Cashflows().FromHandle(hcashflow);
      if ( cashflowData )
      {
        if ( cashflowData->floaters )
        {
          terms->SetFloatingRates(cashflowData->floaters);
        }
        else // not a floater
        {
          terms->SetCashDistribution(cashflowData->cashflow);
        }
      }
      else // must be an OID
      {
        const OIDData& oidData = addin.OIDs()[hcashflow];
        if ( oidData.hasAccrRate )
        {
          if ( oidData.cashflow  )
            terms->SetCashDistribution(oidData.cashflow );
          terms->SetCashPayToZero(oidData.accrRate);          
        }
        else
        {
          if ( oidData.cashflow  )
            terms->SetCashDistribution( oidData.cashflow );
          terms->SetAccretingBond(oidData.grossY2M);
          terms->SetYieldDayCountConvention(oidData.dcc);
          terms->SetYieldCompoundingFrequency(oidData.frequency);
        }       
      }
    }
    //else: zero-coupon bond

    return GetTheIHGAddIn().BondTerms().SetForThisCell(terms);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// helper function used by ITO33HazardRate() to signal errors in input
static void ThrowCDSTSFormatError(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  ThrowValueError("Invalid CDS term structure format", format, ap);

  va_end(ap);
}

//throw if errors
shared_ptr<finance::TermStructureCDS>
CreateTermStructureCDS( const XL::Oper *xloCDS,
                        const  shared_ptr<finance::SessionData>& sessionData)
{
    
    // create the CDSs from the range data
    enum
    {
      Col_MaturityDuration,
      Col_MaturityUnits,
      Col_Spread,
      Col_Recovery,
      Col_Frequency,
      Col_DCC,
      Col_Max
    };
    
    unsigned w, h;
    const XL::Oper *data = xloCDS->GetAsArray(w, h);
    if ( w <= Col_Spread || w >= Col_Max )
      ThrowCDSTSFormatError("incorrect number of columns %u", w);

    shared_ptr<finance::TermStructureCDS> ts ( new finance::TermStructureCDS );
    for ( unsigned n = 0; n < h; ++n, ++data )
    {
      // first the mandatory columns
      int duration;
      if ( !data->As(&duration) )
        ThrowCDSTSFormatError("bad maturity duration in row %u", n + 1);

      std::string units;
      if ( !(++data)->As(&units) )
        ThrowCDSTSFormatError("bad maturity units in row %u", n + 1);

      double spread;
      if ( !(++data)->As(&spread) )
        ThrowCDSTSFormatError("bad spread format in row %u", n + 1);

      spread /= 10000.; // convert from bp to absolute value


      // next the optional ones
      double recovery;
      if ( w > Col_Recovery && !(++data)->IsMissing() )
      {
        if ( !data->As(&recovery) )
          ThrowCDSTSFormatError("bad recovery format in row %u", n + 1);
      }
      else // recovery not specified
      {
        recovery = 0.40;
      }

      finance::Frequency freq;
      if ( w > Col_Frequency && !(++data)->IsMissing() )
      {
        std::string freqStr;
        if ( !data->As(&freqStr) )
          ThrowCDSTSFormatError("bad frequency format in row %u", n + 1);

        freq = FreqFromExcel(freqStr);
      }
      else // frequency not specified
      {
        freq = finance::Frequency_Quarterly;
      }

      Date::DayCountConvention dcc;
      if ( w > Col_DCC && !(++data)->IsMissing() )
      {
        std::string dccStr;
        if ( !data->As(&dccStr) )
          ThrowCDSTSFormatError("bad day counting convention in row %u", n + 1);

        dcc = DccFromExcel(dccStr);
      }
      else // dcc not specified
      {
        dcc = Date::DayCountConvention_30360;
      }

      switch ( TimeUnitFromExcel(units) )
      {
        case finance::TimeUnit_Day:
          ThrowCDSTSFormatError("maturity in days not allowed in row %u", n + 1);
          // break; -- unnecessary after throw
        
        case finance::TimeUnit_Month:
          // nothing done, duration is good
          break;

        case finance::TimeUnit_Year:
          duration *= 12;
          break;

        default:
          FAIL("unknown maturity unit value");
      }

      shared_ptr<finance::ReferenceCDS> cds(new finance::ReferenceCDS( duration,
                                                                      freq,
                                                                      dcc,                                                                      
                                                                      recovery));

      cds->SetSpread( spread );
      
      cds->SetSessionData(sessionData);

      cds->SetMarketPrice(0.); // market price of ref CDS is always 0

      ts->Add(cds);
    }
    return ts;
  
}

XLLEXPORT XLOPER *
ITO33CrossCurrencyData(double fxRate, 
                       const char * baseCurr,
                       const char * foreignCurr, 
                       XL::Oper *xloVolFixedQuanto,
                       XL::Oper *xloCorrelFixedQuanto,
                       XL::Oper *xloEquity)
{
  try
  {

    IHGAddIn& addin = GetTheIHGAddIn();

    shared_ptr<CrossCurrencyData> crossCurrency (new CrossCurrencyData);
    crossCurrency->SetCrossCurrency(baseCurr, foreignCurr, fxRate);
    
    bool
      hasVolatility,
      hasCorrelation,
      hasEquityHandle;

    hasVolatility = 
    hasCorrelation = 
    hasEquityHandle = false;

    double dVol(0.),dCorrel(0.);
    EquitiesManager::Handle hequity;

    EquityData equityData;

    if ( xloVolFixedQuanto && !xloVolFixedQuanto->IsMissing())
    {
      if ( xloVolFixedQuanto->As(&dVol) )
        hasVolatility = true;
      else
        return XL::Oper::ErrValue();
    }

    if ( xloCorrelFixedQuanto && !xloCorrelFixedQuanto->IsMissing())
    {
      if ( xloCorrelFixedQuanto->As(&dCorrel) )
        hasCorrelation= true;
      else
        return XL::Oper::ErrValue();
    }

    if ( xloEquity && !xloEquity->IsMissing())
    {
      if ( xloEquity->As(&hequity) )
      {
        equityData  = addin.Equities()[hequity];
        hasEquityHandle = true;
      }
      else
        return XL::Oper::ErrValue();
    }

    if (hasCorrelation == hasVolatility &&
        hasVolatility == hasEquityHandle )
    {
      if (hasVolatility)
      {
        crossCurrency->SetFixedQuanto(dVol, dCorrel, equityData);
      }
    }
    else 
      throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("The three parameters \"VolFixedQuanto\", \"CorrelFixedQuanto\","
                           " \"Equity\" must been set if the bond is fixed quanto,"
                           " if it not the case, they must not be set")
          );

    return addin.CrossCurrencies().SetForThisCell(crossCurrency);

  }catch( ... )
  {
    return XL::ReturnException();
  }

}

XLLEXPORT XLOPER *
ITO33CBTerms(XL::Oper *xloNewShare,
             XL::Oper *xloTriggerPercentOf,
             XL::Oper *xloTriggerCurrencyOf,
             XL::Oper *xloHr,
             XL::Oper *xloExchUponDef)
{
  try
  {
    shared_ptr<CBTerms> cbTerms (new CBTerms);
    if ( xloNewShare && !xloNewShare->IsMissing() )
      if ( !xloNewShare->As(&cbTerms->newShare) )
         return XL::Oper::ErrValue();
  
    if ( xloTriggerPercentOf && !xloTriggerPercentOf->IsMissing() )
    {
      std::string triggerPercentStr;
      if ( !xloTriggerPercentOf->As(&triggerPercentStr ) )
        return XL::Oper::ErrValue();
      cbTerms->triggerAsPercentageOf = TriggerPercentFromExcel(triggerPercentStr);
    }
    
    if ( xloTriggerCurrencyOf && !xloTriggerCurrencyOf->IsMissing() )
    {
      std::string triggerCurrencyStr;
      if ( !xloTriggerCurrencyOf->As(&triggerCurrencyStr ) )
        return XL::Oper::ErrValue();

      static const EnumValuesNames<finance::TriggerInCurrencyOf> TriggerCurrencyValues[] =
          {
           { "derivative", finance::TriggerInCurrencyOf_Derivative },
           { "underlying", finance::TriggerInCurrencyOf_Underlying }          
         };
    
      if ( 
            !GetEnumFromAbbrevName<SIZEOF(TriggerCurrencyValues)>
              (
                &cbTerms->triggerInCurrencyOf,
                triggerCurrencyStr.c_str(),
                TriggerCurrencyValues
              ) 
         )
      {
        throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown trigger currency of \"%s\"", triggerCurrencyStr)
                );
      }

    }
      
    if ( xloHr && !xloHr->IsMissing() )
    {
      shared_ptr<finance::Issuer> issuer ( new finance::Issuer() );
      shared_ptr<ihg::HazardRate> Hr ( HRFromExcel( xloHr ) );
      
      shared_ptr<ihg::HazardRateTimeOnly>
        HrTC( dynamic_pointer_cast<ihg::HazardRateTimeOnly>(Hr) );
      if (!HrTC)
        throw XL::EXCEPTION_MSG
              (
                xlerrValue,
                String::Printf("Hazard Rate must be time dependent")
              );
      std::vector<double> pvalues ( HrTC->GetValues() );
      std::vector<Date> pdates ( HrTC->GetDates() );
      
      if ( pdates.size() )
      {
        issuer->SetDefaultIntensity( pdates ,pvalues);
      
        cbTerms->CBIssuer = issuer;
      }
      else 
        throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              String::Printf("Hazard Rate must not be flat, It must be Time dependent")
            );
    }

    if ( xloExchUponDef && !xloExchUponDef->IsMissing() )
      if ( !xloExchUponDef->As(&cbTerms->IsExchangeUponDefault) )
        return XL::Oper::ErrValue();

    return GetTheIHGAddIn().CBTerms().SetForThisCell(cbTerms);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
  
}

//Pass to any ConvertibleLike value of CBTerms
//return true if the operation succeded
static bool CBTermsToConvertibleLike(const shared_ptr<finance::ConvertibleLike> & cbLike,
                                     const shared_ptr<CBTerms> & cbTerms )
{
  if (!cbLike)
    return false;
  
  cbLike->SetConvertIntoNewShare(cbTerms->newShare);
  if ( cbTerms->CBIssuer )
  {
    cbLike->SetIssuer( cbTerms->CBIssuer );
    cbLike->SetExchangeable( cbTerms->IsExchangeUponDefault );
  }
  cbLike->SetConversionTriggerAsPercentageOf( cbTerms->triggerAsPercentageOf );
  cbLike->SetTriggerInCurrencyOf( cbTerms->triggerInCurrencyOf );
  
  return true;
}

XLLEXPORT XLOPER *
ITO33CB(XL::SheetObjectHandle hbondterms,
        XL::SheetObjectHandle hequity,
        XL::Oper *xloCBTermsHandler,
        XL::Oper *xloConvs,
        XL::Oper *xloCalls,
        XL::Oper *xloPuts,
        XL::Oper *xloCurrency,
        XL::Oper *xloCrossCurrency,
        XL::Oper *xloFixedExchangeRate)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    // check/dereference the parameters
    const shared_ptr<finance::BondTerms>& terms = addin.BondTerms()[hbondterms];

    const EquityData& equityData = addin.Equities()[hequity];
    
    shared_ptr<CBTerms> cbterms;
    if ( xloCBTermsHandler && !xloCBTermsHandler->IsMissing() )
    {
      CBTermsManager::Handle hcbterms;
      if ( !xloCBTermsHandler->As(&hcbterms) )
        return XL::Oper::ErrValue();

      cbterms = addin.CBTerms()[hcbterms];
    }
    
    // conversions can be either a scalar value indicating a single conversion
    // during the entire bond lifetime or a handle of conversion object
    shared_ptr<finance::ConversionSchedule> conversions;
    shared_ptr<finance::ResetConversionSchedule> reset;
    shared_ptr<finance::ShareDependentConversion> sharedepconv;
    
    if ( xloConvs && !xloConvs->IsMissing() )
    {
      ConversionsManager::Handle hconv;
      if ( xloConvs->As(&hconv) )
      {
        const shared_ptr<finance::ConversionSchedule> * const
          convs = addin.Conversions().FromHandle(hconv);
        if ( convs )
          conversions = *convs;
      }
      //check for reset handle
      if ( !conversions )
      {
        ResetsManager::Handle hres;
        if ( xloConvs->As(&hres) )
        {
          const shared_ptr<finance::ResetConversionSchedule> * const
            res = addin.Resets().FromHandle(hres);
          if ( res )
            reset = *res;
        }
        
        if ( !reset )
        {
          ShareDepConversionsManager::Handle hsharedepconv;
          if ( xloConvs->As(&hsharedepconv) )
          {
            const shared_ptr<finance::ShareDependentConversion> * const
              sharedep = addin.ShareDepConversions().FromHandle(hsharedepconv);
            if ( sharedep )
              sharedepconv = *sharedep;
          }
        }

      }

      if ( !conversions && !reset && !sharedepconv )
      {
        if ( xloConvs->GetType() == XL::Oper::Type_Num )
        {
          double ratio;
          if ( xloConvs->As(&ratio) )
          {
            conversions = make_ptr( new finance::ConversionSchedule );
            conversions->AddConversionPeriod
                         (
                            terms->GetIssueDate(),
                            terms->GetMaturityDate(),
                            ratio
                         );
          }
        }
        else // must be a range passed directly
        {
          conversions = ConvsFromExcelRange(xloConvs);
        } 
      }
    }
    
    if ( !conversions && !reset && !sharedepconv )
      return XL::Oper::ErrValue();
    
    // calls can be either missing, be a handle to an existing call schedule
    // object or just a reference to the call data
    shared_ptr<finance::CallSchedule> calls;
    if ( xloCalls && !xloCalls->IsMissing() )
    {
      CallsManager::Handle hcalls;
      if ( xloCalls->As(&hcalls) )
      {
        const shared_ptr<finance::CallSchedule> * const
          callData = addin.Calls().FromHandle(hcalls);
        if ( callData )
          calls = *callData;
      }

      if ( !calls )
        calls = CallsFromExcelRange(xloCalls);
    }

    // puts are handled in the same way as calls
    shared_ptr<finance::PutSchedule> puts;
    if ( xloPuts && !xloPuts->IsMissing() )
    {
      PutsManager::Handle hputs;
      if ( xloPuts->As(&hputs) )
      {
        const shared_ptr<finance::PutSchedule> * const
          putData = addin.Puts().FromHandle(hputs);
        if ( putData )
          puts = *putData;
      }

      if ( !puts )
        puts = PutsFromExcelRange(xloPuts);
    }

    // support for cross-currency CBs
    shared_ptr<finance::Numeraire> currency;
    shared_ptr<finance::YieldCurve> yieldcurve;
    if ( xloCurrency && !xloCurrency->IsMissing() )
    {
      CurrenciesManager::Handle hcurrency;
      if ( !xloCurrency->As(&hcurrency) )
        return XL::Oper::ErrValue();

      // FIXME: should create RateData here!
      CurrencyData BondCurrency = addin.Currencies()[hcurrency];
      currency = BondCurrency.currency;
      yieldcurve = BondCurrency.yieldcurve;
    }
    
    bool hasFXRate( false );
    double dFXrate(0.);
    const shared_ptr<CrossCurrencyData> *crossCurrency = NULL;

    if ( !XL::Oper::IsEmpty(xloCrossCurrency) )
    {
      CrossCurrenciesManager::Handle hcrosscurrency;
      if ( !xloCrossCurrency->As(&hcrosscurrency) )
        return XL::Oper::ErrValue();
      
      crossCurrency = addin.CrossCurrencies().FromHandle(hcrosscurrency);
      
      if ( !crossCurrency )
        if ( !xloCrossCurrency->As(&dFXrate))
          return XL::Oper::ErrValue();
      
      hasFXRate = true;
    }
    
    bool hasFixedExchangeRate( false );
    double dFixedExchangeRate(0.);
    if ( xloFixedExchangeRate && !xloFixedExchangeRate->IsMissing() )
    {
      if ( !xloFixedExchangeRate->As(&dFixedExchangeRate) )
        return XL::Oper::ErrValue();
      hasFixedExchangeRate = true;
    }

    shared_ptr<finance::CBBase> cbBase;
    if ( conversions )
    {
      XLL_TRACE("Creating new CB with nominal %g and %lu conversions",
              terms->GetNominal(),
              static_cast<unsigned long>(conversions->GetAll().size()));
    
      cbBase = make_ptr( new finance::ConvertibleBond
                             ( terms, conversions ) );
    }
    else if ( reset )
      cbBase = make_ptr( new finance::Reset(terms,reset) );
    else
      cbBase = make_ptr( new finance::AttachedWarrantConvertibleBond
                             (terms,sharedepconv) );

    if ( calls )
      cbBase->SetCallSchedule(calls);

    if ( puts )
      cbBase->SetPutSchedule(puts);

    if ( currency )
      cbBase->SetNumeraire(currency);

    if ( cbterms )
    {
      CBTermsToConvertibleLike(cbBase,cbterms);
    }
    else
    {
      CBTermsToConvertibleLike(cbBase, shared_ptr<CBTerms> (new CBTerms()) );
    }
    
    DerivativesManager& manager = addin.Derivatives();
    DerivativeData& derivData = manager.GetForThisCell();
    derivData.InitFrom(equityData);
    derivData.IsCBOption = false;
    if ( currency && *currency != *equityData.currency )
    {
      if ( hasFXRate)
      {
        if ( crossCurrency )
        {
          // test the validity of the object
          if ( 
               ((*crossCurrency)->baseCurrency == *currency ||
                (*crossCurrency)->baseCurrency == *equityData.currency) &&
               ((*crossCurrency)->foreignCurrency == *currency ||
                (*crossCurrency)->foreignCurrency == *equityData.currency) &&
               ((*crossCurrency)->baseCurrency != (*crossCurrency)->foreignCurrency)
             )
          {
            // Set the FXrate, as the 1FX = FXrate *Base and the ito33 standard choose 
            // as base currency the currency of the bond,
            // we set it as 1 / FXrate if the base is the equity currency 
            if ( (*crossCurrency)->baseCurrency != *currency ) 
              dFXrate = 1 / (*crossCurrency)->spotFXRate;
            else
              dFXrate = (*crossCurrency)->spotFXRate;

            derivData.SetCrossCurrency( currency, yieldcurve, dFXrate);
          }
          else // throw if cross currency data is incompatible with the bond and the equity
            throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              String::Printf("Cross currency data are incompatible with the bond and equity currencies!")
            );
          if ( (*crossCurrency)->isFixedQuanto )
          {
            if ( (*crossCurrency)->CompareEquityDataWith( equityData ) )
              cbBase->SetFixedQuanto( (*crossCurrency)->volatility, (*crossCurrency)->correlation );
            else
              throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              String::Printf("Equity Handle must be the same between the cross currency data and the bond!")
            ); 
          }
        }
        else // we have enter a double directly.
          derivData.SetCrossCurrency( currency, yieldcurve, dFXrate );
      }
      else
        throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("For a cross currency bond, FX Rate parameter can not be null")
          );
    }
    else
    {
      if ( hasFXRate)
        throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("FX Rate must be null for a bond which is not a cross currency one!")
          );
    }
    if ( hasFixedExchangeRate )
    cbBase->SetFixedFXRate( dFixedExchangeRate );

    derivData.derivative = cbBase;

    return manager.HandleOf(derivData);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33CBOption(XL::SheetObjectHandle hcb,
              XL::SheetObjectHandle hfloatingrates,
              double maturity,
              const char * notionalType)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();

    const DerivativeData& cbderivData = addin.Derivatives()[hcb];
    shared_ptr<finance::ConvertibleBond> 
      cb( dynamic_pointer_cast<finance::ConvertibleBond>
          (cbderivData.derivative) );

    const CashflowData& cashflowData = addin.Cashflows()[hfloatingrates];
    
    if ( !cb )
      throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("The  derivative must be a convertible bond")
          );

    if ( !cashflowData.floaters )  
      throw XL::EXCEPTION_MSG
          (
            xlerrValue,
            String::Printf("the cash flow must be floating rates")
          );

    const Date maturityDate(Date::FromExcel(maturity));

    static const EnumValuesNames<finance::ASWNotionalIs>
    g_ASWNotionalIs[] =
      {
        { "issue price",
          finance::ASWNotionalIs_IssuePrice },
        { "put price",
          finance::ASWNotionalIs_PutPrice }
      };
    
    finance::ASWNotionalIs aswNotional;

    if ( !GetEnumFromAbbrevName<SIZEOF(g_ASWNotionalIs)>(&aswNotional, notionalType , g_ASWNotionalIs) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown asset swap notional type \"%s\"", notionalType)
                );
    }

    DerivativeData& derivData = addin.Derivatives().GetForThisCell();
    derivData.eqCurrency = cbderivData.eqCurrency;
    derivData.eqYieldcurve = cbderivData.eqYieldcurve;
    derivData.equity = cbderivData.equity;
    derivData.derCurrency = cbderivData.derCurrency;
    derivData.derYieldcurve = cbderivData.derYieldcurve;
    derivData.FXRate = cbderivData.FXRate;
    derivData.IsCBOption = true;
    
    shared_ptr<finance::CBOption> 
      cbOpt( new finance::CBOption( cb,
                                    cashflowData.floaters,
                                    maturityDate )
            ); 
    cbOpt->SetASWNotionalIs(aswNotional);
    
    derivData.derivative = cbOpt;

    return addin.Derivatives().HandleOf(derivData);
    
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33PEPSCall(double start,
              double end,
              double trigger,
              const char * calltype,
              const XL::Oper *xloKeepAccr,
              bool forfeitCoupon,
              const char *makeWholeType,
              double makeWholePremium,
              const XL::Oper *xloTriggerPercentOf,
              const XL::Oper *xloAdvCallFeat)
{
  try
  {
    static const EnumValuesNames<finance::GeneralizedPEPSLikeCallType>
      g_GeneralizedPEPSLikeCallType[] =
      {
        { "fixed shares",
          finance::GeneralizedPEPSLikeCallType_FixedShare},
        { "variable shares",
          finance::GeneralizedPEPSLikeCallType_VariableShare }
      };
    
    finance::GeneralizedPEPSLikeCallType pepscallType;

    if ( !GetEnumFromAbbrevName<SIZEOF(g_GeneralizedPEPSLikeCallType)>(&pepscallType, calltype, g_GeneralizedPEPSLikeCallType) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown generalized PEPS call type \"%s\"", calltype)
                );
    }

    shared_ptr<finance::GeneralizedPEPSLikeCall> 
      PEPScall ( new finance::GeneralizedPEPSLikeCall(Date::FromExcel(start),
                                                      Date::FromExcel(end),
                                                      trigger,
                                                      pepscallType) );

       // deal with the rest of the parameters
    if ( xloKeepAccr && !xloKeepAccr->IsMissing() )
    {
      bool keepAccrued;
      if ( !xloKeepAccr->As(&keepAccrued) )
        return XL::Oper::ErrValue();

      PEPScall->SetKeepAccrued(keepAccrued);
    }
    //else: default is fine, it's true by default

    PEPScall->SetForfeitCoupon(forfeitCoupon);

    SetAdvCallFeatures(PEPScall, xloAdvCallFeat);

    if ( !XL::Oper::IsEmpty(xloTriggerPercentOf) )
    {
      std::string triggerPercentStr;
      if ( !xloTriggerPercentOf->As(&triggerPercentStr ) )
        return XL::Oper::ErrValue();
      PEPScall->SetTriggerAsPercentageOf (TriggerPercentFromExcel(triggerPercentStr) );
    }
    
 
    bool isCouponPresentValue;
    switch ( MakeWholeTypeFromExcel(makeWholeType, &isCouponPresentValue) )
    {
      case finance::MakeWholeType_Coupon:
        PEPScall->SetCouponMakeWhole(isCouponPresentValue);
        break;

      case finance::MakeWholeType_Premium:
        PEPScall->SetPremiumMakeWhole(makeWholePremium);

        // make the check below pass
        makeWholePremium = 0.;
        break;

      case finance::MakeWholeType_Max:
        // no make whole protection
        break;

      default:
        FAIL("unknown make whole type value");
    }

    if ( makeWholePremium != 0. )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  "Premium is ignored if the make whole type is not \"P\""
                );
    }

    // now store the call schedule
    return GetTheIHGAddIn().PEPSLikeCalls().SetForThisCell(PEPScall);;
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// creates a enum to reuse GetEnumFromAbbrevName function
static enum PEPSAverageType
{
  PEPSAverage_Stock,
  PEPSAverage_Ratio,
  PEPSAverage_Max
};

XLLEXPORT XLOPER *
ITO33PEPSAveraging(double start,
                   double end,
                   int nbsampling,
                   const char * averageType,
                   double currentAveraging,
                   int nbsampleused)
{
  try
  {
    
    PEPSAverageType pepsaveragetype; 
      
    static const EnumValuesNames<PEPSAverageType>
    g_PEPSAverageType[] =
    {
      { "stock",
        PEPSAverage_Stock},
      { "ratio",
        PEPSAverage_Ratio}
    };
    
    if ( !GetEnumFromAbbrevName<SIZEOF(g_PEPSAverageType)>(&pepsaveragetype, averageType, g_PEPSAverageType) )
    {
      throw XL::EXCEPTION_MSG
                (
                  xlerrValue,
                  String::Printf("Unknown PEPS averaging type \"%s\"", averageType)
                );
    }
    
    shared_ptr<finance::PEPSAveragingPeriod> PEPSAveraging;

    if ( pepsaveragetype == PEPSAverage_Stock )
    {
      PEPSAveraging = finance::PEPSAveragingPeriod::CreateWithStock(Date::FromExcel(start),
                                                                    Date::FromExcel(end),
                                                                    nbsampling);
      
      if ( currentAveraging )
      {
        PEPSAveraging->SetCurrentStockAverage(currentAveraging,nbsampleused);
      }
    }
    else
    {
      PEPSAveraging = finance::PEPSAveragingPeriod::CreateWithConversionRatio(Date::FromExcel(start),
                                                                              Date::FromExcel(end),
                                                                              nbsampling);
      
      if ( currentAveraging )
      {
        PEPSAveraging->SetCurrentConversionRatioAverage (currentAveraging,nbsampleused);
      }
    }

    return GetTheIHGAddIn().PEPSAveragingPeriods().SetForThisCell(PEPSAveraging);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33PEPSTerms(double downCR,
               double lowerStrike,
               double upCR,
               double higherStrike,
               bool hasOptionalConv,
               const XL::Oper* xloAveraging)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();
    
    shared_ptr<PEPSTerms> pepsterms (new PEPSTerms);
    
    pepsterms->downConvRatio = downCR;
    pepsterms->lowerStrike = lowerStrike;
    pepsterms->upConvRatio = upCR;
    pepsterms->higherStrike = higherStrike;
    pepsterms->hasOptionalConv = hasOptionalConv;
    
    if ( !XL::Oper::IsEmpty(xloAveraging) )
    {
      PEPSAveragingManager::Handle haveraging;
      if ( xloAveraging->As(&haveraging) )
      {
        shared_ptr<finance::PEPSAveragingPeriod> 
          pepsAveraging( GetTheIHGAddIn().PEPSAveragingPeriods()[haveraging] );
        pepsterms->PEPSAveraging = pepsAveraging;
      }
      else
        return XL::Oper::ErrValue();
    }
            
    return addin.PEPSTerms().SetForThisCell(pepsterms);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33PEPS(XL::SheetObjectHandle hbondterms,
          XL::SheetObjectHandle hpepsterms,
          XL::SheetObjectHandle hequity,
          XL::Oper * xloCalls)
{
  try
  {
    IHGAddIn& addin = GetTheIHGAddIn();
    
    // check/dereference the parameters
    const shared_ptr<finance::BondTerms>& bondterms = addin.BondTerms()[hbondterms];
    const shared_ptr<PEPSTerms> pepsterms = addin.PEPSTerms()[hpepsterms];
    const EquityData& equityData = addin.Equities()[hequity];
    
    shared_ptr<finance::GeneralizedPEPSLike> peps (new finance::GeneralizedPEPSLike( bondterms,
                                                                  pepsterms->downConvRatio,
                                                                  pepsterms->lowerStrike,
                                                                  pepsterms->upConvRatio,
                                                                  pepsterms->higherStrike) );

    if ( pepsterms->hasOptionalConv )
      peps->EnableOptionalConversion();

    DerivativeData& derivData = addin.Derivatives().GetForThisCell();
    derivData.InitFrom( equityData );
    derivData.IsCBOption = false;
    derivData.derivative = peps;
    
    shared_ptr<finance::CallSchedule> fixedcalls;
    shared_ptr<finance::GeneralizedPEPSLikeCall> pepslikecall;

    if ( !XL::Oper::IsEmpty(xloCalls) ) 
    {
      CallsManager::Handle hcalls;
      if ( xloCalls->As(&hcalls) )
      {
        const shared_ptr<finance::CallSchedule> * const
          calls = addin.Calls().FromHandle(hcalls);
        if ( calls )
          peps->SetCallFixedCash(*calls);
        else
        {
          const shared_ptr<finance::GeneralizedPEPSLikeCall> * const
            pepscall = addin.PEPSLikeCalls().FromHandle(hcalls);
          if ( pepscall )
            peps->SetGeneralizedPEPSLikeCall(*pepscall);
          else 
            throw XL::EXCEPTION_MSG
            (
              xlerrValue,
              "invalid handle for calls"
            );
        }
      }
      else
        return XL::Oper::ErrValue();
    }

    if ( pepsterms->PEPSAveraging )
      peps->SetAveragingPeriod( pepsterms->PEPSAveraging );
    

    return addin.Derivatives().HandleOf(derivData);    
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// calibration
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33HazardRate(const XL::Oper *xloCDS,
                XL::SheetObjectHandle hcurrency,
                const XL::Oper *xloDate)
{
  try
  {
    // setup the common pricing objects
    const Date dateVal = ValDateFromExcel(xloDate);

    const CurrencyData& currencyData = GetTheIHGAddIn().Currencies()[hcurrency];
    const shared_ptr<finance::Numeraire>& currency = currencyData.currency;

    // the spot used here doesn't matter as we're calibrating a time-only
    // hazard rate, i.e. nothing depend on spot anyhow
    shared_ptr<finance::Equity> equity(new finance::Equity(1.0, currency));

    shared_ptr<finance::RateData> rateData(new finance::RateData);
    rateData->SetYieldCurve(currency, currencyData.yieldcurve);

    shared_ptr<finance::SessionData>
      sessionData(new finance::SessionData(rateData, equity, dateVal));
       
    shared_ptr<finance::TermStructureCDS> ts ( CreateTermStructureCDS( xloCDS, sessionData )) ;

    // do calibrate
    ihg::ParametrizationHRWithTimeComponent param;
    shared_ptr<ihg::HazardRateWithTimeComponent> hr(param.CalibrateWithCDSs(*ts));

    const std::vector<Date>& dates = hr->GetDates();
    const std::vector<double>& values = hr->GetTimeComponentValues();

    // and return the result as an array
    const unsigned size = static_cast<unsigned>(dates.size());
    XL::Oper *array = new XL::Oper[2*size];

    XL::Oper *p = array;
    for ( unsigned n = 0; n < size; ++n )
    {
      *p++ = dates[n];
      *p++ = values[n];
    }

    return XL::Oper::ReturnArray(2, size, array);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

static shared_ptr<finance::SessionData> CreateSessionData(const DerivativeData& derivData,
                                                         double spotPrice,
                                                         XL::Oper *xloDate)
{
  const shared_ptr<finance::Equity> equity = derivData.equity;

  equity->SetSpotSharePrice(spotPrice);

  // create everything we need for the computation
  shared_ptr<finance::RateData> rateData(new finance::RateData);
  rateData->SetYieldCurve(derivData.eqCurrency, derivData.eqYieldcurve);
    
  if (derivData.derCurrency &&
      derivData.derYieldcurve && 
      *derivData.eqCurrency != *derivData.derCurrency)
  {
    rateData->SetYieldCurve(derivData.derCurrency, derivData.derYieldcurve);
    rateData->SetSpotFXRates(derivData.FXRate);
    //throw XL::EXCEPTION_MSG(xlerrValue, "Cross Currency is not supported yet");
  }

  Date valDate = ValDateFromExcel(xloDate);
  //Set the SessionData
  shared_ptr<finance::SessionData> sessionData(new finance::SessionData
                                                    (
                                                      rateData,
                                                      equity,
                                                      valDate
                                                    ));
  return sessionData;
}

// Because of the the validation test that ensure the to session of the cb and
// the cboption be the same, we have to call this function as the CB can have 
//a session previously set

static void SetSessionDataToTheUnderlyingCBOfCBOption( finance::Derivative& deriv,
                                                            const shared_ptr<finance::SessionData>& sessionData )
{
  finance::CBOption* cbOpt = dynamic_cast<finance::CBOption *> ( &deriv );
  if ( cbOpt )
    cbOpt->GetConvertibleBond()->SetSessionData(sessionData);
}

XLLEXPORT XLOPER *
ITO33BrownianVolFromPrice(const XL::Oper *xloHR,
                          XL::SheetObjectHandle hderiv,
                          double derivPrice,
                          double spotPrice,
                          XL::Oper *xloDate,
                          XL::Oper *xloDebugStr)                          
{
  try
  {
    // set up the scene for computation
    const DerivativeData& derivData = GetTheIHGAddIn().Derivatives()[hderiv];
    finance::Derivative& deriv = *derivData.derivative;
    
    deriv.SetMarketPrice(derivPrice);

    shared_ptr<finance::SessionData> 
      sessionData( CreateSessionData(derivData, spotPrice,xloDate) );

    // Do a special session setting for the CBOption
    // as the underlying CB must have the same session
    if ( derivData.IsCBOption )
      SetSessionDataToTheUnderlyingCBOfCBOption(deriv, sessionData );
        
    deriv.SetSessionData(sessionData);

    shared_ptr<ihg::TheoreticalModel> model(new ihg::TheoreticalModel);
    model->SetHazardRate(HRFromExcel(xloHR));

    //Set DebugOutput to true if a debug output file is specified
    if ( xloDebugStr && !xloDebugStr->IsMissing() )
    {
      std::string debugfileStr;
      if ( !xloDebugStr->As(&debugfileStr) )
        return XL::Oper::ErrValue();
      
      model->EnableDebugOutput();
      model->SetDebugOutputFile( debugfileStr );
    }

    // do compute and return
    return XL::Oper::Return(model->ComputeImpliedBrownianVol(deriv));
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33BrownianVolFromBSVol(const XL::Oper *xloHR,
                          XL::SheetObjectHandle hderiv,
                          double implVol,
                          double spotPrice,
                          XL::Oper *xloDate,
                          XL::Oper *xloDebugStr)
{
  try
  {
    // compute the implied volatility
    shared_ptr<ihg::TheoreticalModel> model(new ihg::TheoreticalModel);
    model->SetHazardRate(HRFromExcel(xloHR));

    const DerivativeData& derivData = GetTheIHGAddIn().Derivatives()[hderiv];
    finance::Option * const
      option = dynamic_cast<finance::Option *>(derivData.derivative.get());
    if ( !option ) 
    {
      throw XL::EXCEPTION_MSG(xlerrValue,
                              "Implied volatility is for options only");
    }

    option->SetImpliedVol(implVol);
    
    option->SetSessionData(CreateSessionData(derivData, spotPrice,xloDate));

    //Set DebugOutput to true if a debug output file is specified
    if ( xloDebugStr && !xloDebugStr->IsMissing() )
    {
      std::string debugfileStr;
      if ( !xloDebugStr->As(&debugfileStr) )
        return XL::Oper::ErrValue();
      
      model->EnableDebugOutput();
      model->SetDebugOutputFile( debugfileStr );
    }

    return XL::Oper::Return(model->ComputeImpliedBrownianVol(*option));
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// computation
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33Output(XL::SheetObjectHandle hderiv,
            double spot,
            double vol,
            XL::Oper *xloHR,
            XL::Oper *xloDate,
            XL::Oper *xloDebugStr)
{
  try
  {
    // get the derivative to compute and the associated equity
    const DerivativeData& derivData = GetTheIHGAddIn().Derivatives()[hderiv];

    finance::Derivative& deriv = *derivData.derivative;

    Date valDate = ValDateFromExcel(xloDate);
    
    shared_ptr<finance::SessionData> 
      sessionData( CreateSessionData(derivData, spot,xloDate) );

    // Do a special session setting for the CBOption
    // as the underlying CB must have the same session
    if ( derivData.IsCBOption )
      SetSessionDataToTheUnderlyingCBOfCBOption(deriv,sessionData);
      
    deriv.SetSessionData(sessionData);
        
    shared_ptr<ihg::TheoreticalModel> model(new ihg::TheoreticalModel);
    model->SetVolatility(shared_ptr<ihg::Volatility>(new ihg::VolatilityFlat(vol)));
    model->SetHazardRate(HRFromExcel(xloHR));

    // enable computation of the greeks which are turned off by default as we
    // have no way to know if the user is, or not, going to ask for them using
    // =ITO33Vega() &c later
    shared_ptr<finance::ComputationalFlags> flags(new finance::ComputationalFlags);
    flags->SetComputeRho(true);
    flags->SetComputeVega(true);
    flags->SetComputeFugit(true);
    flags->SetAnalysisDate(valDate);
    deriv.SetComputationalFlags(flags);

    //Set DebugOutput to true if a debug output file is specified
    if ( xloDebugStr && !xloDebugStr->IsMissing() )
    {
      std::string debugfileStr;
      if ( !xloDebugStr->As(&debugfileStr) )
        return XL::Oper::ErrValue();
      
      model->EnableDebugOutput();
      model->SetDebugOutputFile( debugfileStr );
    }

    // finally do compute output (this can throw)
    XLL_TRACE("Computing the output with spot=%g, vol=%g", spot, vol);
    shared_ptr<finance::ModelOutput> mo = model->Compute(deriv);

    // create an OutputData object that store inputs and computation result
    //we do it to reuse it later in the hedging functions

    OutputData outputs;
    outputs.Output = mo;
    outputs.Model = model;
    outputs.Sessiondata = sessionData;
    outputs.Derivative = derivData.derivative;

    // remember the computation result and the inputs to allow other cells to access it
    return GetTheIHGAddIn().Outputs().SetForThisCell(outputs);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33HedgeOutput(XL::SheetObjectHandle houtput,
                 XL::SheetObjectHandle hdefhedgingderiv)           
{
  try
  {
    // get the default hedging derivative to compute.
    const DerivativeData& derivData = GetTheIHGAddIn().Derivatives()[hdefhedgingderiv];

    finance::Derivative& defhedgingderiv = *derivData.derivative;

    // get the model and the target derivative to compute the  associated hedge ratio.  
    const OutputData& outputData = GetTheIHGAddIn().Outputs()[houtput];

    finance::Derivative& targetderiv = *outputData.Derivative;

    // Do a special session setting for the CBOption
    // as the underlying CB must have the same session
    finance::CBOption* cbOpt = dynamic_cast<finance::CBOption *> ( &targetderiv );
    if ( cbOpt )
      cbOpt->GetConvertibleBond()->SetSessionData(outputData.Sessiondata);
    
    targetderiv.SetSessionData(outputData.Sessiondata);
    
    // Do a special session setting for the CBOption
    // as the underlying CB must have the same session
    if ( derivData.IsCBOption )
      SetSessionDataToTheUnderlyingCBOfCBOption(defhedgingderiv,outputData.Sessiondata);
   
    defhedgingderiv.SetSessionData(outputData.Sessiondata);

    ihg::TheoreticalModel& model = *outputData.Model;

    ihg::PerfectHedgeRatios ratios = model.ComputePerfectHedgeRatios( targetderiv , defhedgingderiv );

    return GetTheIHGAddIn().HedgeRatios().SetForThisCell( PerfectHedgeRatiosPtr( new ihg::PerfectHedgeRatios(ratios) ) );
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}


XLLEXPORT XLOPER *
ITO33SpotScenario(XL::SheetObjectHandle houtput,
                  double perturbation)
{
  try
  {
    std::vector<double> 
      pdSpots,
      pdValues;
        
    // get the model and the target derivative to compute the  associated hedge ratio.  
    const OutputData& outputData = GetTheIHGAddIn().Outputs()[houtput];

    ModelOutputPtr output = outputData.Output;

    double dSpot = outputData.Sessiondata->GetSpotSharePrice();
    
    pdSpots = output->GetSpotsAtAnalysisDate();
    pdValues = output->GetPricesAtAnalysisDate();
        
    double 
      pdSpotsInterpolated[2],
      pdValuesInterpolated[2];
    
    //LinearInterpolate accept only sorted values so we do it
    if (perturbation > 0)
    {
      pdSpotsInterpolated[0]= dSpot;
      pdSpotsInterpolated[1]= dSpot * (1 + perturbation);
    }
    else
    {
      pdSpotsInterpolated[0]= dSpot * (1 + perturbation);
      pdSpotsInterpolated[1]= dSpot;
    }
    
    numeric::LinearInterpolate(&pdSpots[0],&pdValues[0],pdSpots.size(),pdSpotsInterpolated,pdValuesInterpolated,2);
    
    //return result with the right sign
    if (perturbation > 0)
      return XL::Oper::Return( (pdValuesInterpolated[1] - pdValuesInterpolated[0])/pdValuesInterpolated[0]);
    else 
      return XL::Oper::Return( (pdValuesInterpolated[0] - pdValuesInterpolated[1])/pdValuesInterpolated[1]);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

// ----------------------------------------------------------------------------
// ModelOutput accessors
// ----------------------------------------------------------------------------

static XLOPER *
GetOutputField(XL::SheetObjectHandle houtput,
               double (finance::ModelOutput::*func)() const)
{
  try
  {
    const OutputData& outputdata = GetTheIHGAddIn().Outputs()[houtput];
    const ModelOutputPtr& output = outputdata.Output; 

    return XL::Oper::Return((output.get()->*func)());
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

static XLOPER *
GetHedgeOutputField(XL::SheetObjectHandle hhedgeoutput,
               double (ihg::PerfectHedgeRatios::*func)() const )
{
  try
  {
    const PerfectHedgeRatiosPtr& ratio = GetTheIHGAddIn().HedgeRatios()[hhedgeoutput];

    return XL::Oper::Return((ratio.get()->*func)());
  }
  catch( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *
ITO33Price(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetPrice);
}

XLLEXPORT XLOPER *
ITO33Delta(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetDelta);
}

XLLEXPORT XLOPER *
ITO33Gamma(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetGamma);
}

XLLEXPORT XLOPER *
ITO33Theta(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetTheta);
}

XLLEXPORT XLOPER *
ITO33Rho(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetRho);
}

XLLEXPORT XLOPER *
ITO33UnderlyingRho(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetUnderlyingRho );
}

XLLEXPORT XLOPER *
ITO33Vega(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetVega);
}

XLLEXPORT XLOPER *
ITO33Fugit(XL::SheetObjectHandle houtput)
{
  return GetOutputField(houtput, &finance::ModelOutput::GetFugit);
}

XLLEXPORT XLOPER *
ITO33BondFloor(XL::SheetObjectHandle houtput)
{
  try
  {
    const OutputData& outputdata = GetTheIHGAddIn().Outputs()[houtput];
    const ModelOutputPtr& output = outputdata.Output; 
        
    const shared_ptr<finance::BondLikeOutput>
      bondoutput( dynamic_pointer_cast<finance::BondLikeOutput>(output) );

    if ( bondoutput )
      return XL::Oper::Return(bondoutput->GetBondFloor());
    else
      throw XL::EXCEPTION_MSG(xlerrValue,
                              "The bond floor is not defined for this kind of derivative");
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
  
}

//Template function to avoid repetition of code if function
//ITO33YieldToMaturity and ITO33YieldToPut. U parameter must 
//be a pointer to function: 
//[double (T::*func)(double,finance::frequency,Date::dcc)]
template<class T,class U> static  XLOPER *
GetFixedIncomeOutput(T arg,
                     U func,
                     double price)
{
  try
  {
    finance::Frequency freq;
    Date::DayCountConvention dcc;

    const shared_ptr<finance::BondTerms> 
      bondterms = arg->GetBondTerms();                            //point of customization for typename T
      
    freq = 
      GetYieldCompoundingFrequencyEvenIfUndefined( *bondterms );

    dcc = 
      GetYieldDayCountConventionEvenIfUndefined( *bondterms );

    if ( freq == finance::Frequency_Undefined )
      freq = finance::Frequency_Annual;
        
    if ( dcc == Date::DayCountConvention_Max )
      dcc = Date::DayCountConvention_ActAct;
    // no coupon and no yield
            
    return XL::Oper::Return( (arg.get()->*func)(price,
                                                freq,
                                                dcc
                                               )                  //point of customization for typename U
                           );
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

//Compute the yield from the pricing date to the maturity
XLLEXPORT XLOPER *
ITO33YieldToMaturity(XL::SheetObjectHandle houtput, double price)
{
  try 
  {
    const shared_ptr<finance::CBBase> cbbase = 
      dynamic_pointer_cast<finance::CBBase>
      ( GetTheIHGAddIn().Outputs()[houtput].Derivative );
    
    if ( cbbase )
      return GetFixedIncomeOutput(cbbase, &finance::CBBase::ComputeYieldToMaturity, price);
    
    const shared_ptr<finance::Bond> bond =
      dynamic_pointer_cast<finance::Bond>
      ( GetTheIHGAddIn().Outputs()[houtput].Derivative );
    
    if ( bond )
      return GetFixedIncomeOutput(bond, &finance::Bond::ComputeYieldToMaturity, price);
    
    throw XL::EXCEPTION_MSG(xlerrValue,
                              "The yield to maturity is not defined for this kind of derivative");
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

//Compute the yield from the pricing date to the first put date
XLLEXPORT XLOPER *
ITO33YieldToPut(XL::SheetObjectHandle houtput, double price)
{
  try 
  {
    const shared_ptr<finance::CBBase> cbbase = 
      dynamic_pointer_cast<finance::CBBase>
      ( GetTheIHGAddIn().Outputs()[houtput].Derivative );
    
    if ( cbbase )
      return GetFixedIncomeOutput(cbbase, &finance::CBBase::ComputeYieldToPut, price);
    
    const shared_ptr<finance::Bond> bond =
      dynamic_pointer_cast<finance::Bond>
      ( GetTheIHGAddIn().Outputs()[houtput].Derivative );
    
    if ( bond )
      return GetFixedIncomeOutput(bond, &finance::Bond::ComputeYieldToPut, price);
    
    throw XL::EXCEPTION_MSG(xlerrValue,
                              "The yield to put is not defined for this kind of derivative");
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

//Compute the accrued interest for the ConvertibleLike
//TODO: Make it Works for simple bond
XLLEXPORT XLOPER *
ITO33AccruedInterest(XL::SheetObjectHandle houtput)
{
  try
  {
    const shared_ptr<finance::ConvertibleLike> cblike = 
      dynamic_pointer_cast<finance::ConvertibleLike>
      ( GetTheIHGAddIn().Outputs()[houtput].Derivative );
    if ( cblike )
      return XL::Oper::Return( cblike->GetAccruedInterestValue() );
    else
      throw XL::EXCEPTION_MSG(xlerrValue,
                              "The accrued is not defined for this kind of derivative");
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

//
XLLEXPORT XLOPER *
ITO33StraightBond(XL::SheetObjectHandle houtput)
{
  try
  {    
    // get the model and the cb like to compute the  associated straight bond.  
    const OutputData& outputData = GetTheIHGAddIn().Outputs()[houtput];

    shared_ptr<finance::CBBase> 
      cb( dynamic_pointer_cast<finance::CBBase>(outputData.Derivative) );

    if ( cb )
    {
      ihg::TheoreticalModel& model = *outputData.Model;

      cb->SetSessionData(outputData.Sessiondata);

      shared_ptr<finance::ModelOutput> mo = model.Compute( *cb->GetStraightBond() );

      return XL::Oper::Return( mo->GetPrice() );

    }
    else
      throw XL::EXCEPTION_MSG(xlerrValue,
                              "Function can not be called for this derivative");
     
  }
  catch( ... )
  {
    return XL::ReturnException();
  }

}

// ----------------------------------------------------------------------------
// HedgeOutput accessors
// ----------------------------------------------------------------------------

XLLEXPORT XLOPER *
ITO33UnderlyingRatio(XL::SheetObjectHandle hhedgeoutput)
{
  return GetHedgeOutputField(hhedgeoutput, &ihg::PerfectHedgeRatios::GetUnderlyingHedgeRatio);
}

XLLEXPORT XLOPER *
ITO33CreditRatio(XL::SheetObjectHandle hhedgeoutput)
{
  return GetHedgeOutputField(hhedgeoutput, &ihg::PerfectHedgeRatios::GetDefaultHedgeRatio);
}

// ============================================================================
// IHGAddIn implementation
// ============================================================================

void IHGAddIn::RegisterAllFunctions()
{
  XL_REGISTER_0
  (
    ITO33IHGVersion,
    "Returns the version of the IHG library"
  );

  XL_REGISTER_1
  (
    ITO33T,
    "Transposes the given range to be passed to other functions",
    "Range",
    "A (usually horizontal) range to be transposed"
  );

  XL_REGISTER_2
  (
    ITO33Split,
    "Break a line and transforme it to a matrix",
    "Range",
    "A horizontal range with one line",
    "Column number",
    "Number of colum of the output range"
  );

  XL_REGISTER_2
  (
    ITO33ZeroCouponYC,
    "Creates a zero coupon yield curve and returns its handle",
    "Legs", "Range containing dates and corresponding rates",
    "Date", "Reference date from which the (today by default)"
  );

  XL_REGISTER_5
  (
    ITO33YC,
    "Creates an interbank yield curve and returns its handle",
    "SwapRates",  "The swap rates",
    "SwapDCC",    "Day counting convention for swap rates (ACT/ACT by default)",
    "CashRates",  "The cash rates (may be omitted)",
    "CashDCC",    "Day counting convention for cash rates (ACT/360 by default)",
    "Date",       "Reference date (today by default)"
  );

  XL_REGISTER_2
  (
    ITO33CurrencyData,
    "Creates a currency data object and returns its handle",
    "Name",       "Name of the currency, e.g. \"EUR\"",
    "YieldCurve", "Yield curve for this currency or a flat rate"
  );

  XL_REGISTER_5
  (
    ITO33Equity,
    "Creates a equity object and returns its handle",
    "currencyData",   "Currency of this equity",
    "BorrowRate",     "Borrow rate (0 by default)",
    "Dividends",      "Range containing dividends dates, amounts and types",
    "FiscalMonth",    "Fiscal year start month (January by default)",
    "FiscalDay",      "Fiscal year start day (1st by default)"
  );

  XL_REGISTER_4
  (
    ITO33Cashflow,
    "Creates a cashflow object with the specified payments",
    "ContractDate", "The contracting date",
    "Payments",     "The coupon dates and values (as % of the nominal)",
    "DCC",          "The day counting convention (ACT/ACT by default)",
    "Frequency",    "The payment frequency (quaterly by default)"
  );

  XL_REGISTER_7
  (
    ITO33UniformCashflow,
    "Creates a cashflow with the regular payments and returns its handle",
    "SettlementDate",   "The settlement date",
    "FirstDate",        "The first payment date",
    "LastDate",         "The last payment date",
    "AnnualAmount",     "The annual payment amount in % of the nominal",
    "DCC",              "The day counting convention (ACT/ACT by default)",
    "Frequency",        "The payment frequence (quaterly by default)",
    "LastPaymentType",  "The last payment type. It must be use if the last payment"
                        " is short or long (Optional, value in \"short\" or \"long\""  
  );

  XL_REGISTER_10
  (
    ITO33FloatingRates,
    "Creates an object containing floating rates information and returns its handle",
    "Margin",       "Margin between the floating and reference rate",
    "StartAccrued", "The start of accrued date",
    "FirstUnknown", "The first unknown payment date",
    "LastUnknown",  "The last unknown payment date",
    "Frequency",    "The payment frequency",
    "DCC",          "The day counting convention (ACT/360 by default)",
    "Floor",        "The minimal amount of the annual coupon (can be ommitted)",
    "Cap",          "The maximal amount of the coupon (can be omitted)",
    "Multiplier",   "Multiplier for the reference rate (default is 1)",
    "FixingDelay",  "The fixing delay (can be omitted)"
    /* "lastButOne", "The last but one unknown payment date", */
  );

  XL_REGISTER_2
  (
    ITO33FloatCashflow,
    "Creates a cashflow with fixed and floating rates and returns its handle",
    "Floaters",   "The floating rates object",
    "Payments",   "The known coupon dates and values"
  );

  XL_REGISTER_4
  (
    ITO33OID,
    "Creates a cashflow for an OID bond and returns its handle",
    "GrossY2M",             "The gross yield to maturity of the bond",
    "CompoundingFrequency", "The compounding frequency used for the "
                            "computation of the accreted principal",
    "DCC",                  "The day counting convention (ACT/ACT by default)",
    "Cashflow",             "The handler of the cashflow object representing "
                            "the coupons paid by the bond."
  );

  XL_REGISTER_3
  (
    ITO33CashPayToZero,
    "Creates a cashflow for a bond paying coupons until a certain date only",
    "AccretionRate",        "Accreetion rate from the date of the last "
                            "coupon up to maturity",
    "CompoundingFrequency", "The compounding frequency used for the "
                            "computation of the accreted principal",
    "Cashflow",             "The handler of the cashflow object representing "
                            "the coupons paid by the bond."
  );

  XL_REGISTER_3
  (
    ITO33AdvCallFeatures,
    "Creates an object to define the special call feature (the computation will be slowed!)",
    "NoticePeriod",   "The number of days before the call during which the "  
                      "holder of the bond may convert into shares. If zero, " 
                      "there is no notice period (0 by default)",             
    "TriggerPeriod",  "The number of succesive days that the spot price "     
                      "must cross the trigger in order to activate the call " 
                      "provision (0, meaning none, by default)",              
    "TriggerHistory", "The number of days that the spot share price has "     
                      "already crossed the trigger. Only meaningfull if "     
                      "TriggerPeriod is used."    
  );
  
  // ITO33Calls() and ITO33CallsByYield() share all tail parameters so define
  // this macro to avoid duplicating them
  #define IHGCALLS_PARAMS                                                               \
    "KeepAccrued",              "Keep accrued flag (true by default)",                  \
    "ForfeitCoupon",            "Forfeit coupon flag (false by default)",               \
    "MakeWholeType",            "Make whole protection type: \"CA\" for absolute "      \
                                "coupon, \"CP\" for present value coupon, \"P\" for "   \
                                "premium or \"N\" or omitted for none",                 \
    "Premium",                  "Make whole premium used if make whole type is \"P\"",  \
    "TriggerPercentageOf",      "can be either \"Principal\",\"Issue Price\",\"Clam\""  \
                                "(default \"Principal\")",                              \
    "AdvancedCallFeatures",   "The handler of the advanced call feature"

  IHGAddIn& addin = GetTheIHGAddIn();
  addin.Register
  (
    ITO33Calls, "ITO33Calls",
    "Creates a call schedule object and returns its handle",
    "Calls",          "Call provisions range with start, end, strike and, "
                      "optionally if there are any soft calls, trigger columns",
    IHGCALLS_PARAMS
  );

  addin.Register
  (
    ITO33CallsByYield, "ITO33CallsByYield",
    "Creates a call schedule with strikes defined by yield and returns its handle",
    "Calls",          "Call provisions range with start, end, and, "
                      "optionally if there are any soft calls, trigger columns",
    "Yield",          "The yield applied to compute the call strike.",
    IHGCALLS_PARAMS
  );

  #undef IHGCALLS_PARAMS

  XL_REGISTER_3
  (
    ITO33Puts,
    "Creates a put schedule object and returns its handle",
    "Puts",           "Put provisions range with date",
    "KeepAccrued",    "Keep accrued flag (true by default)",
    "ForfeitCoupon",  "Forfeit coupon flag (false by default)"
  );

  XL_REGISTER_4
  (
    ITO33PutsByYield,
    "Creates a put schedule with yield and returns its handle",
    "Puts",           "Put provisions range with date and strike columns",
    "Yield",          "The yield applied to compute the call strike",
    "KeepAccrued",    "Keep accrued flag (true by default)",
    "ForfeitCoupon",  "Forfeit coupon flag (false by default)"
  );

  XL_REGISTER_10
  (
    ITO33Conversion,
    "Creates a conversion schedule consisting of a single conversion",
    "StartDate",      "Conversion start date",
    "EndDate",        "Conversion end date",
    "Ratio",          "Conversion ratio",
    "CoCoType",       "Contingent conversion type (may be omitted)",
    "Trigger",        "Contingent conversion trigger (may be omitted)",
    "TriggerRate",    "Trigger change rate (may be omitted)",
    "TriggerCeil",    "Trigger ceiling value (may be omitted)",
    "Cash",           "Cash paid on conversion (may be omitted)",
    "KeepAccrued",    "Keep accrued flag (true by default)",
    "ForfeitCoupon",  "Forfeit coupon flag (false by default)"
  );

  XL_REGISTER_4
  (
    ITO33Conversions,
    "Creates a conversion schedule object and returns its handle",
    "Conversions",    "Range containing conversions data with the columns "
                      "corresponding to ITO33Conversion() arguments",
    "KeepAccrued",    "Keep accrued flag (true by default)",
    "ForfeitCoupon",  "Forfeit coupon flag (false by default)",
    "TriggerCondMet", "Trigger condition met flag (false by default)"
  );
  
  XL_REGISTER_9
  (
    ITO33Reset,
    "Creates a reset schedule object and returns its handle",
    "StartDate",          "Reset conversion start date",
    "EndDate",            "Reset conversion end date",
    "initial CP",         "Initial conversion price",
    "current CP",         "Current conversion price",
    "Reset Floored By",   "Must be \"prevaling\" or \"initial\" conversion price",
    "Reset range",        "range containing the reset date",
    "Cash",               "Cash paid on conversion (may be ommited)",
    "KeepAccrued",        "Keep accrued flag (true by default)",
    "ForfeitCoupon",      "Forfeit coupon flag (false by default)"
  );
  
  XL_REGISTER_10
  (
    ITO33AttachedWarrant,
    "Creates a share dependent conversion object and returns its handle",
    "StartDate",          "The conversion start date",
    "EndDate",            "the conversion end date",
    "BaseRatio",          "The base ratio",
    "IncrementalFactor",  "The incremental factor",
    "CapRatio",           "The cap ratio of the share dependent ratio",
    "ResetDate",          "The reset date (optional)",
    "FixedStrike",        "The fixed strike, it make sense only if a reset date is specified (optional)",
    "CurrentRatio",       "The current ratio, (optional, must be set if the reset date is before the valuation date)",
    "KeepAccrued",        "Keep accrued flag (true by default)",
    "ForfeitCoupon",      "Forfeit coupon flag (false by default)"
  );

  XL_REGISTER_6
  (
    ITO33AttachedWarrantWithCoCo,
    "Creates a share dependent conversion with CoCo object and returns its handle",
    "AttachedWarrant",      "The handler of an attached warrant object",
    "CoCoType",             "Contingent conversion type",
    "Trigger",              "Contingent conversion trigger",
    "TriggerRate",          "Trigger change rate (may be omitted)",
    "TriggerCeil",          "Trigger ceiling value (may be omitted)",
    "IsLastTriggerCondMet", "flag specifying if the last trigger condition have been met (default false)"
  );

  XL_REGISTER_5
  ( 
    ITO33GenerateDividends,
    "Generates Dividends and returns a range",
    "Dividends",      "A dividend type range",
    "MaxDate",        "Maximun Date of dividend payment generation",
    "GrowthRate",     "Annual growth rate",
    "Frequency",      "Dividend payments frequency",
    "IsAdjusted",     "Unregular dividend flag (false by default)"
  );

  XL_REGISTER_5
  (
    ITO33Option,
    "Creates an option object and returns its handle",
    "Equity",     "The underlying equity",
    "Strike",     "The strike of the option",
    "Maturity",   "The expiration date",
    "Type",       "The contract type: \"call\" or \"put\"",
    "Exercise",   "The exercise type: \"american\" or \"european\""
  );

  XL_REGISTER_3
  (
    ITO33CDS,
    "Creates a CDS object and returns its handle",
    "Equity",     "The associated equity",
    "Cashflow",   "The uniform cashflow object",
    "Recovery",   "The recovery value (40% by default)"
  );

  XL_REGISTER_7
  (
    ITO33BondTerms,
    "Creates an object containing the bond terms and returns its handle",
    "IssueDate",  "Issue date",
    "Maturity",   "Maturity date",
    "IssuePrice", "Issue price",
    "Nominal",    "Nominal value",
    "Redemption", "Redemption rate",
    "Recovery",   "The recovery value",
    "Cashflow",   "The cashflow or OID object (may be omitted for "
                  "zero coupon bonds)"
  );
  
  XL_REGISTER_5
  (
    ITO33CBTerms,
    "Create an object containing the common terms with all CBlike Instrument and return its handle",
    "NewShare",               "New share flag (default false)",
    "TriggerPercentOf",       "can be either \"Principal\",\"Issue Price\",\"Clam\" (default \"Principal\")",
    "TriggerCurrencyOf",      "can be either \"Derivative\",\"Underlying\" (default \"Derivative\")",
    "IssuerHazardRate",       "range of hazard rate by \"ITO33HazardRate\" function (optional)",
    "ExchangeUponDefault",    "Exchange upon default flag (default false)"
  );

  XL_REGISTER_6
  (
    ITO33CrossCurrencyData,
    "Create an object containing all data to price a cross currency bond",
    "SpotFXRate",             "The spot exchange rate (i.e the number of Base currency"
                              "needed to purchase one unit of Foreign currency )",
    "Base",                   "The name of the base currency",
    "Foreign",                "The name of the foreign currency",
    "VolFixedQuanto",         "The volatility of the FX Rate (optional, must" 
                              " be set only if the derivative is fixed quanto)",
    "CorrelFixedQuanto",      "The correlation of the FX Rate with the equity (optional, must" 
                              " be set only if the derivative is fixed quanto)",
    "Equity",                 "The handle of the equity (optional, must" 
                              " be set only if the derivative is fixed quanto)"
  );

  XL_REGISTER_9
  (
    ITO33CB,
    "Creates a convertible bond object and returns its handle",
    "BondTerms",              "The bond terms object",
    "Equity",                 "The underlying equity",
    "CBTerms(optional)",      "The CBterms object handler",
    "Conversions",            "Conversions object, range containing the conversion "
                              "periods or constant conversion rate if there is a single "
                              "conversion running from issue date to the maturity",
    "Calls",                  "Calls object or range with calls data",
    "Puts",                   "Puts object or range with puts data",
    "CurrencyData",           "Bond currency data for cross-currency bonds",
    "CrossCurrencyData",      "The Cross Currency Data Handler or the"
                              " number of units of the derivative currency needed to purchase "
                              "one unit of the equity currency",
    "FixedExchangeRate",      "The fixed exchange rate (optional)"       
  );

  XL_REGISTER_4
  (
    ITO33CBOption,
    "Creates a CB option object and returns its handle",
    "CB",               "The CB handler",
    "FloatingRates",    "The Floating Rates handler",
    "Maturity",         "The maturity of the asset swap",
    "NotionalType",     "The notional type (\"issue price\" or \"put price\")"
  );

  XL_REGISTER_10
  (
    ITO33PEPSCall,
    "Creates a PEPS like call and returns its handler",
    "StartDate",                "The start date of the call period",
    "EndDate",                  "The end date of the call period",
    "Trigger",                  "The trigger that activates the call feature",
    "CallType",                 "the call type, can be either \"fixed shares\" or \"variable shares\"",
    "KeepAccrued",              "Keep accrued flag (true by default)",
    "ForfeitCoupon",            "Forfeit coupon flag (false by default)",
    "MakeWholeType",            "the make whole",
    "Premium",                  "The premium",
    "TriggerPercentOf",         "can be either \"Principal\",\"Issue Price\""
                                ",\"Clam\" (default \"Principal\")",
    "AdvancedCallFeatures",     "The handler of the advanced call feature"
  );

  XL_REGISTER_6
  (
    ITO33PEPSAveraging,
    "Creates a PEPS Averaging object an return its handler",
    "AveragingStartDate",     "The average period start date",
    "AveragingEndDate",       "The average period end date",
    "NbOfSampling",           "The number of sampling",
    "AverageType",            "The type of the averaging period. It can be either \"stock\" or \"ratio\"",
    "CurrentAverage",         "The current average",
    "NbOfSamplesUsed",        "The number of samples used to compute the current average"
  );
  
  XL_REGISTER_6
  (
    ITO33PEPSTerms,
    "Creates an object containing description for a peps like instrument",
    "DownConversionRatio",      "The conversion ratio that will be aplied if the "
                                "underlying share price is less then the LowerStrike",
    "LowerStrike",              "The lower level of the underlying stock price "
                                "that will be used in the formula for conversion ratio",
    "UpConversionRatio",        "The conversion ratio that will be applied if "
                                "the underlying stock price is greater then HigherStrike",
    "HighStrike",               "The higher level of the underlying stock "
                                "price that will be used in the formula for conversion ratio",
    "HasOptionalConversion",    "If this argument is tru the holder of the instrument may " 
                                "convert into ordinary shares before maturity. "
                                "The number of shares will be computed using the "
                                "same formula as for the mandatory conversion at maturity",
    "AveragingPeriod",          "The number of days used to compute the average stock price "
                                "at maturity that will be used as input for the payoff formula. "
                                "If zero no averaging will be used"
  );

  XL_REGISTER_4
  (
    ITO33PEPS,
    "Creates a peps like object and returns its handle",
    "BondTerms",    "The handler for the bond of this structure as returned by the ITO33BondTerms function",
    "PEPSTerms",    "The handler of the PEPSTerms",
    "Equity",       "The underlying equity handle",
    "Call",         "The call handler"
  );

  XL_REGISTER_3
  (
    ITO33HazardRate,
    "Computes the hazard rate from the given CDS term structure. "
    "This is an array function and must be entered in a 2 by number of CDS "
    "range using Ctrl-Shift-Enter",
    "TermStructure",  "A range containing CDS duration and units (y/m), "
                      "spreads and optionally recovery (40% by default), "
                      "frequency (quaterly) and day count (30/360)",
    "CurrencyData",   "Handle of the currency data object",
    "ValuationDate",  "The valuation date (today by default)"
  );

  XL_REGISTER_6
  (
    ITO33BrownianVolFromPrice,
    "Computes the brownian volatility corresponding to the given hazard rate "
    "and market price.",
    "HazardRate",           "Hazard rate: constant value or a range with date and rate "
                            "columns",
    "Derivative",           "The derivative object to be used for volatility "
                            "computation",
    "MarketPrice",          "The market price of the derivative",
    "SpotPrice",            "The spot price of the underlying",
    "Date",                 "The valuation date (today by default)",
    "DebugFile(optional)",  "The xml debug file"
  );

  XL_REGISTER_6
  (
    ITO33BrownianVolFromBSVol,
    "Computes the brownian volatility corresponding to the given hazard rate "
    "and implied Black-Scholes volatility.",
    "HazardRate",           "Hazard rate: constant value or a range with date and rate "
                            "columns",
    "Option",               "The option object to be used for volatility computation",
    "BSVol",                "The implied Black-Scholes volatility of the given option",
    "SpotPrice",            "The spot price of the underlying",
    "Date",                 "The valuation date (today by default)",
    "DebugFile(optional)",  "The xml debug file"
  );


  XL_REGISTER_6
  (
    ITO33Output,
    "Compute all model outputs for the given derivative",
    "Derivative",           "The derivative object to be computed",
    "Spot",                 "The spot price to use for the computation",
    "Vol",                  "Volatility to use for the computation",
    "HazardRate",           "Hazard rate: constant value or a range with date and rate "
                            "columns",
    "Date",                 "The valuation date (today by default)",
    "DebugFile(optional)",  "The xml debug file"
  );

  XL_REGISTER_2
  (
    ITO33HedgeOutput,
    "Compute the perfect hedge ratios given the default hedging derivative, return a handler",
    "Output",               "The output handler",
    "Hedging derivative",   "The hedging derivative"           
  );

  XL_REGISTER_1(ITO33Price,
                "Return the computed price",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33Delta,
                "Return the computed delta",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33Gamma,
                "Return the computed gamma",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33Theta,
                "Return the computed theta",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33Rho,
                "Return the computed rho",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33UnderlyingRho,
                "Return the computed underlying rho (must be used for cross currency instrument)",
                "Output",   "The handle of the computation output object");
  
  XL_REGISTER_1(ITO33Vega,
                "Return the computed vega",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33Fugit,
                "Return the computed fugit",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33BondFloor,
                "Return the computed bond Floor",
                "Output",   "The handle of the computation output object");
  
  XL_REGISTER_2(ITO33YieldToMaturity,
                "Return the yield to maturity of a bond via his output object",
                "Output",   "The handle of the computation output object",
                "Price",    "The bond target price (dirty)"
                );

  XL_REGISTER_2(ITO33YieldToPut,
                "Return the yield to put of a bond via his output object",
                "Output",   "The handle of the computation output object",
                "Price",    "The bond target price (dirty)"
                );

  XL_REGISTER_1(ITO33AccruedInterest,
                "Return the accrued interest of a bond",
                "Output",   "The handle of the computation output object");

  XL_REGISTER_1(ITO33StraightBond,
                "Return the straight bond for a cb like instrument",
                "Output",   "The handle of the computation output object");
  
  XL_REGISTER_1(ITO33UnderlyingRatio,
                "Return the underlying hedge ratio",
                "HedgeOutput",   "The handle of the hedge output object"
                );
  XL_REGISTER_1(ITO33CreditRatio,
                "Return the credit hedge ratio",
                "HedgeOutput",   "The handle of the hedge output object"
                );

  XL_REGISTER_2(ITO33SpotScenario,
                "Return the convexity",
                "Output",        "The handle of the output",
                "Perturbation",  "The perturbation"
                );
}

