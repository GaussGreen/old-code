///////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/error.cpp
// Purpose:     implementation of ito33::finance::Error
// Author:      ITO 33 Canada
// Created:     March 30, 2005
// RCS-ID:      $Id: error.cpp,v 1.19 2006/08/22 12:28:42 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
///////////////////////////////////////////////////////////////////////////////


#include "ito33/finance/error.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{
namespace finance
{
  ITO33_IMPLEMENT_ERROR_CLASS;
}
}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------

using ito33::finance::Error;

extern const Error ITO33_OUTPUT_NOT_AVAILABLE
("Required output is not available.");

extern const Error ITO33_DOMAIN_SPOTS_NOT_SET
("Underlying share prices must be set in Domain before getting values at.");

extern const Error ITO33_INVALID_DERIVATIVE
("Null derivative pointer."); 

extern const Error ITO33_INVALID_DERIVATIVEWEIGHT
("Derivative weights must be strictly positive.");

extern const Error ITO33_MATURITYBEFOREVALUATION
("Maturity date before valuation date.");

extern const Error ITO33_ISSUEDATE_AFTER_MATURITYDATE
("Maturity date before issue date.");

extern const Error ITO33_INVALID_FISCAL_YEAR
("Issuer definition: Invalid fiscal year start date.");

extern const Error ITO33_VALUATION_BEFORE_ISSUE
("Valuation date is before issue date.");

extern const Error ITO33_INVALID_RECOVERYRATE_1
("Setting invalid recovery rate value: %f.");

extern const Error ITO33_INVALID_BARRIER
("Barrier must be positive.");

//------------------------------------------------------------------------------
// ENUM
//------------------------------------------------------------------------------
extern const Error ITO33_INVALID_FREQUENCY
("Invalid payment frequency specified.");

extern const Error ITO33_INVALID_TIMEUNIT
("Invalid TimeUnit value.");

//------------------------------------------------------------------------------
// CASH FLOW
//------------------------------------------------------------------------------
extern const Error ITO33_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT
("Annual payment of uniform cash flow stream must be positive.");

extern const Error ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE
("First payment date before contraction date.");

extern const Error ITO33_CASHFLOWSTREAM_FIRST_PAYMENT_DATE_AFTER_LAST
("First payment date after last payment date.");

extern const Error ITO33_CASHFLOWSTREAM_INCONSISTENT_DATES
("First payment date is not consistent with "
  "last or last but one payment date.");

extern const Error ITO33_INVALID_LASTPAYMENTTYPE
("Invalid last payment type (for payment date generation) specified.");

extern const Error ITO33_NULL_COMPUTATIONALFLAGS
("Invalid computational flags in TheoreticalModel.");

//------------------------------------------------------------------------------
// MODEL
//------------------------------------------------------------------------------

extern const Error ITO33_POST_DEFAULT_VOLATILITY_NEGATIVE
("The post default volatility must be non negative.");

extern const Error ITO33_POST_DEFAULT_VOLATILITY_TOO_LARGE
("The post default volatility must be less than 500%.");

extern const Error ITO33_INVALID_UNDERLYINGPROCESS
("Invalid underlying process.");

//------------------------------------------------------------------------------
// DERIVATIVE
//------------------------------------------------------------------------------

//___________________________ GENERAL      _____________________________________
extern const Error ITO33_DERIVATIVE_INVALID_SESSIONDATA
("Derivative definition: setting invalid SessionData.");

extern const Error ITO33_DERIVATIVE_INVALID_NUMERAIRE
("Derivative definition: setting invalid Numeraire.");

extern const Error ITO33_DERIVATIVE_UNDEFINED_SESSIONDATA
("Derivative definition: SessionData is undefined.");

extern const Error ITO33_DERIVATIVE_VISITOR_NOT_IMPLEMENTED
("Modifying visitor is not supported for this derivative.  "
"Please report this error to the ITO33 support team.");

//___________________________ MARKET PRICE _____________________________________
extern const Error ITO33_INVALID_MARKETPRICE_NEGATIVE
("Market price must be strictly positive.");

extern const Error ITO33_INVALID_MARKETPRICE_TOOBIGTOBESET_1
("Derivative definition: setting enormous market price value %f.");

extern const Error ITO33_NO_MARKETPRICE
("Market price or data is required, but is not set in the derivative.");

extern const Error ITO33_PREVIOUS_SHARE_PRICE_UNDEFINED
("The previous share price must be set if the "
"valuation date is past or equal to the start of the sampling period.");

extern const Error ITO33_START_SHARE_PRICE_UNDEFINED
("The share price at the start of a period is required for payoff if the "
"valuation date is past or equal to the start of the period.");


//------------------------------------------------------------------------------
// CALIBRATION
//------------------------------------------------------------------------------
extern const Error ITO33_CALIBRATION_FAIL
("Calibration failed, the given market data can't be matched.");

extern const Error ITO33_CALIBRATION_CANCELLED
("Calibration has been cancelled by user.");

//------------------------------------------------------------------------------
// DERIVATIVES, TERMSTRUCTURE
//------------------------------------------------------------------------------
extern const Error ITO33_INCONSISTENT_SESSION_DATA
("The term structure elements do not have the same session data.");

extern const Error ITO33_TERMSTRUCTURE_INCOMPLETE
("The term structure is incomplete.");

extern const Error ITO33_TERMSTRUCTURE_ADD_NULL
("Adding null to term structure.");

extern const Error ITO33_TERMSTRUCTURE_ADD_DUPLICATED
("The term structure has already an instrument having the same maturity "
"as the one to be added.");

extern const Error ITO33_EMPTY_CDSCURVE 
("Empty CDS curve.");

extern const Error ITO33_EMPTY_PARBONDCURVE
("Empty ParBond curve.");
  
extern const Error ITO33_EMPTY_EDSCURVE 
("Empty EDS curve.");

extern const Error ITO33_EMPTY_OPTIONCURVE 
("Empty option curve.");

extern const Error ITO33_EMPTY_DERIVATIVELIST
("Empty derivative list.");

extern const Error ITO33_EMPTY_DERIVATIVECURVE
("Empty derivative curve.");

//------------------------------------------------------------------------------
// SWAP/YIELD CURVE
//------------------------------------------------------------------------------
extern const Error ITO33_YIELDCURVE_REFERENCEDATE_UNDEF
("Yield curve definition: reference date must be specified.");

extern const Error ITO33_SWAP_RATE_MATURITY_UNIT
("Time unit part of maturity term of a swap rate must be either 'Month' or 'Year'.");

extern const Error ITO33_SWAP_RATE_WRONG_FREQUENCY
("Payment frequency and the maturity term of the swap rate are not consistent.");

extern const Error ITO33_YIELDCURVE_SWAP_NODATA
("Swap curve definition: neither cash rates nor swap rates is defined.");

extern const Error ITO33_YIELDCURVE_SWAP_CASHRATES_MATURITES
("Swap curve definition: a maturity term for cash rates defined more than twice.");

extern const Error ITO33_YIELDCURVE_SWAP_SWAPRATES_MATURITES
("Swap curve definition: a maturity term for swap rates defined more than twice.");

extern const Error ITO33_YIELDCURVE_SWAP_INCONSISTENCY_SWAPCASH
("Swap curve definition: any maturity term for cash rate should be less than "
" any maturity term for swap rates.");

extern const Error ITO33_YIELDCURVELEG_RATE 
("Rate must be greater than or equal to zero.");

extern const Error ITO33_YIELDCURVELEG_MATURITY
("The maturity duration of a yield curve must be strictly positive.");

extern const Error ITO33_NEG_FLATRATE
("Negative flat rate.");

extern const Error ITO33_BOOTSTRAPPING_FAILED
("Bootstrap method failed. It is normally due to the existence of irregular "
"values in the given swap curve.");

extern const Error ITO33_YIELDCUVE_ANNUALLYCCOMPOUNDED
("Invalid or empty YieldCurveAnnuallyCompounded.");

extern const Error ITO33_YIELDCURVE_INVALID_RATE 
("YieldCurveAnnuallyCompounded definition: rate value must be positive.");

extern const Error ITO33_YIELDCURVE_SETLEGS_BADSIZE 
("YieldCurveAnnuallyCompounded definition: "
"Arrays of dates and values must have the same number of elements.");

//------------------------------------------------------------------------------
// SESSION DATA CONCERNED (ISSUER, RATEDATA etc.)
//------------------------------------------------------------------------------
extern const Error ITO33_INVALID_RATEDATA
("Invalid RateData specified.");

extern const Error ITO33_RATEDATA_INVALID_YIELD_CURVE
("RateData definition: invalid YieldCurve specified.");

extern const Error ITO33_RATEDATA_INVALID_SPOTFXRATES 
("RateData definition: invalid SpotFXRates specified.");

extern const Error ITO33_NUMERAIRE_INVALID
("Invalid Numeraire specified.");

extern const Error ITO33_RATEDATA_MISSING_YIELDCURVE
("The yield curve for the specified currency (numeraire) could not "
"be found.");

extern const Error ITO33_SPOTFXRATES_NEGATIVERATE 
("SpotFXRates definition: exchange rates must be strictly positive.");

extern const Error ITO33_SPOTFXRATES_MISSING_RATE
("The exchange rate for the specified currencies (numeraires) could "
"not be found.");

extern const Error ITO33_SESSIONDATA_INVALID_RATEDATA 
("SessionData definition: invalid RateData specified.");

extern const Error ITO33_SESSIONDATA_INVALID_EQUITY
("SessionData definition: invalid Equity specified.");

extern const Error ITO33_ISSUER_INVALID_DEFAULT_INTENSITY
("Issuer definition: setting invalid default intensity.\n"
" Dates must be "
"increasing, values must be non negative, and both must have same non "
"zero size");

extern const Error ITO33_ISSUER_NO_DEFAULT_INTENSITY
("The default intensity of the issuer is required (by exchangeable bond) "
"but is not available.");

//------------------------------------------------------------------------------
// TIME FUNCTION
//------------------------------------------------------------------------------
extern const Error ITO33_TIMEFUNCTION_EMPTY
("TimeFunction defintion: given value arrays is empty.");

extern const Error ITO33_TIMEFUNCTION_INVALID_SIZE
("TimeFunction defintion: date and value arrays must have same size.");

extern const Error ITO33_TIMEFUNCTION_INVALID_DATES
("TimeFunction defintion: dates must be increasing.");

//------------------------------------------------------------------------------
// EQUITY
//------------------------------------------------------------------------------
extern const Error ITO33_NEG_SPOT
("The spot share price must be strictly positive.");

extern const Error ITO33_UNDEF_SPOT
("The spot share price is not defined.");

extern const Error ITO33_NEG_PREVIOUS_SHARE_PRICE
("The previous share price must be strictly positive.");

extern const Error ITO33_EQUITY_INVALID_ISSUER
("Equity definition: invalid Issuer specified.");

extern const Error ITO33_EQUITY_INVALID_NUMERAIRE
("Equity definition: invalid Numeraire specified.");

extern const Error ITO33_EQUITY_INVALID_DIVIDENDS
("Equity definition: invalid Dividends specified.");

extern const Error ITO33_EQUITY_INVALID_BORROWCURVE
("Equity definition: invalid borrow curve specified.");

extern const Error ITO33_MONEYMARKET_INVALID_YIELD_CURVE
("MoneyMarket definition: invalid YieldCurve specified");

extern const Error ITO33_MONEYMARKET_INVALID_NUMERAIRE
("MoneyMarket definition: invalid Numeraire specified");

//------------- PAR BOND --------------------------------------------------
extern const Error ITO33_PARBOND_INVALID_CONTRACTING_DATE
("Par bond definition: invalid contracting date.");

extern const Error ITO33_PARBOND_NULL_MATURITY
("Par bond definition: maturity should not be 0.");

extern const Error ITO33_PARBOND_NEGATIVE_YTM
("Par bond definition: setting negative risk-free yield to maturity.");

extern const Error ITO33_PARBOND_NEGATIVE_SPREAD
("Par bond definition: setting negative spread.");

extern const Error ITO33_PARBOND_INCONSISTENT_FREQUENCY_MATURITY
("Par bond definition: maturity and payment frequency are inconsistent.");

//__________________________ Averaging Period _____________________________
extern const Error ITO33_AVERAGING_PERIOD_WRONG_START_END
("Averaging period: end date before start date.");

extern const Error ITO33_AVERAGING_PERIOD_INVALID_FREQUENCY
("Averaging period: invalid frequency selected.");

extern const Error ITO33_AVERAGING_NEGATIVE_CURRENT_AVERAGE
("Averaging period: current average negative.");

extern const Error ITO33_AVERAGINGPERIOD_ALREADY_SET
("Averaging period: period already set.");

extern const Error ITO33_AVERAGINGPERIOD_NULL
("Averaging period: null pointer.");

//__________________________Floating Rates________________________________
extern const Error ITO33_FLOATINGRATES_INCOHERENT_PAYMENTSTREAM
("Floating payment stream definition: paymentDates and paymentRates "
"must have the same size.");

extern const Error ITO33_FLOATINGRATES_INVALID_PAYMENTSTREAM
("Floating payment stream definition: Known payment rates must be "
"before unknown ones");

extern const Error ITO33_FLOATINGRATES_NO_UNKNOWN_PAYMENTS
("Floating rates data: Required data for unknown payments not available");

extern const Error ITO33_FLOATINGRATES_FIRSTUNKNOWN_AFTER_LASTUNKNOWN
("FloatingRates definition: first unknown payment date must not be after "
 "last unknown payment date.");

extern const Error
ITO33_FLOATINGRATES_FIRSTKNOWNPAYMENT_BEFORE_STARTOFACCRUED
("FloatingRates definition: first known payment date must be greater "
"than the start of accrued date.");

extern const Error ITO33_NULL_FLOATING_RATES
("Invalid floating rates.");


//_________________________CDS________________________________________
extern const Error ITO33_CDS_INVALID_SPREADSTREAM
("CDS definition: invalid spread stream.");

extern const Error ITO33_REFERENCECDS_MATURITY
("The maturity of the reference CDS must be positive.");

extern const Error ITO33_REFERENCECDS_NONZERO_PRICE
("The price of a reference CDS must be zero.");

extern const Error ITO33_REFERENCECDS_SPREAD_NOT_AVAILABLE
("The spread of the reference CDS is not set.");

extern const Error ITO33_REFERENCECDS_SPREADSTREAM_NOT_AVAILABLE
("The spread stream is not available. The spread must be set first.");

extern const Error ITO33_REFERENCECDS_SPREADSTREAM_UNDEFSESSIONDATA
("The spread stream is not available. The session data of the reference"
 " CDS must be defined first.");

extern const Error ITO33_REFERENCECDS_INVALID_SPREAD
("The spread of a reference CDS must be strictly positive.");

extern const Error ITO33_REFERENCECDS_MATURITYDATE_NOT_AVAILABLE
("The maturity date of the reference CDS is not available. "
 "The session data of the reference CDS must be defined first.");

extern const Error ITO33_REFERENCECDS_FIRSTPAYMENTDATE_NOT_AVAILABLE
("The first payment date of the reference CDS is not available. "
 "The session data of the reference CDS must be defined first.");
//_________________________Variance swaps_____________________________
extern const Error ITO33_VARIANCESWAP_INVALID_MATURITY
("Variance swap definition: invalid maturity date.");

extern const Error ITO33_VARIANCESWAP_INVALID_STARTOFSAMPLINGPERIOD
("Variance swap definition: invalid start of sampling period.");

extern const Error ITO33_VARIANCESWAP_START_AFTER_MATURITY
("Variance swap definition: the start of the sampling period must be "
"before the maturity date.");

extern const Error ITO33_VARIANCESWAP_NEGATIVE_VOLATILITYSTRIKE
("Variance swap definition: volatility strike cannot be negative.");

extern const Error ITO33_VARIANCESWAP_NULL_TERMS
("Variance swap definition: null terms.");

extern const Error ITO33_VARIANCESWAP_NEGATIVE_CURRENTVOLATILITY
("Variance swap definition: current volatility cannot be negative.");

extern const Error ITO33_VARIANCESWAP_NOCURRENTVOLATILITY
("Variance swap definition: current volatility must be set if the "
"valuation date is greater than the start of the sampling "
"period.");

extern const Error ITO33_VARIANCESWAP_INCONSISTENT_CURRENT
("Variance swap definition: if the current volatility is set, "
"the start of the sampling period must be before the valuation date.");

extern const Error ITO33_VARIANCESWAP_INVALID_NBSAMPLINGRETURNS
("Variance swap definition: must have at least one return to sample."); 

extern const Error ITO33_VARIANCESWAP_INVALID_NB_SAMPLES_USED
("Variance swap definition: the number of samples used to compute the "
"current volatility must be less than the total number of sampling "
"returns.");

extern const Error ITO33_VARIANCESWAP_INVALID_CURRENT_COUNT
("Variance swap definition: the number of days within the corridor "
"between the sampling start date and the valuation date must be greater"
"than or equal to zero and less than the total number of sampling returns.");

extern const Error ITO33_VARIANCESWAP_ZERO_NB_SAMPLES_USED
("Variance swap definition: the number of samples used to compute the "
"current volatility must be at least one.");

extern const Error ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_SMALL
("Variance swap definition: the cap multiplier must be greater or equal "
"to one.");

extern const Error ITO33_VARIANCESWAP_CAP_MULTIPLIER_TOO_LARGE
("Variance swap definition: the cap multiplier must be less or equal "
"to ten.");

extern const Error ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_VANILLA
("Variance swap definition: can only compute implied volatility strike "
"for variance swap without cap multiplier.");

extern const Error ITO33_VARIANCESWAP_IMPLIEDSTRIKE_NOT_IMPLEMENTED
("Variance swap definition: implied volatility strike is not supported or "
 "not yet implemented.");

extern const Error ITO33_VARIANCESWAP_IMPLIEDSTRIKE_FOR_STARTED
("Variance swap definition: can only compute implied volatility strike "
"for variance swap terms when the sampling is not yet started.");

extern const Error ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_ANNUAL_RETURN_FREQUENCY
("Variance swap definition: the number of business days in a year must "
"be strictly positive.");

extern const Error ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_CORRIDOR
("Variance swap definition: the corridor barrier must be strictly positive.");

extern const Error ITO33_VARIANCESWAP_UP_CORRIDOR_INVALID
("Variance swap definition: the up corridor barrier must be less than "
"the down corridor barrier.");

extern const Error ITO33_VARIANCESWAP_DOWN_CORRIDOR_INVALID
("Variance swap definition: the down corridor barrier must be greater than "
"the up corridor barrier.");

extern const Error ITO33_VARIANCESWAP_NON_STRICTLY_POSITIVE_START_PRICE
("Variance swap definition: the share price at the start of the sampling "
 "period must be strictly positive.");

extern const Error ITO33_VARIANCESWAP_START_SHARE_PRICE_UNDEFINED
("Variance swap definition: the share price at the start of the sampling "
"period must be set if the valuation date is equal to or greater than "
"the start of the sampling period.");

extern const Error ITO33_VARIANCESWAP_CURRENT_COUNT_UNDEFINED
("Variance swap definition: the current conditional count must be "
"set if the valuation date is equal to or greater than the start of "
"the sampling period.");

extern const Error ITO33_VARIANCESWAP_FORWARD_STARTING_GAMMA
("Variance swap definition: forward starting gamma variance swaps "
"are not currently supported.");

extern const Error ITO33_VARIANCESWAP_IHG_CONDITIONAL
("Variance swap definition: conditional variance swaps are not currently "
"supported by the IHG model.");

extern const Error ITO33_VARIANCESWAPTION_NULL_VSTERMS
("Variance swaption definition: null varaince swap terms.");

extern const Error ITO33_VARIANCESWAPTION_MATURITY_AFTER_STARTOFVSSAMPLING
("Variance swaption definition: The start of sampling return period of the "
 "underlying variance swap of this swaption must not be before the maturity "
 "date of the swaption.");

extern const Error ITO33_VARIANCESWAP_INVALID_OPTIONTYPE
("Option variance swap definition: invalid option type.");

extern const Error ITO33_VARIANCESWAP_LOGFORMULA_NOT_APPLICABLE
("Variance swap definition: only vanilla variance swap with log return of "
 "variance can be computed using the log formula.");
