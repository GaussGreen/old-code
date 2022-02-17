///////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/bonderror.cpp
// Purpose:     implementation of ito33::Error for bondlike
// Author:      ZHANG Yunzhi
// Created:     06.14.04
// RCS-ID:      $Id: bonderror.cpp,v 1.96 2006/07/28 17:46:36 dave Exp $
// Copyright:   (c) 2004- Trilemma LLP
///////////////////////////////////////////////////////////////////////////////


#include "ito33/finance/bondlike/bonderror.h"

// ----------------------------------------------------------------------------
// Implement Error class
// ----------------------------------------------------------------------------
namespace ito33
{
namespace finance
{
  ITO33_IMPLEMENT_GIVEN_ERROR_CLASS(BondError);
}
}

// ----------------------------------------------------------------------------
// definitions of the Error objects
// ----------------------------------------------------------------------------
using ito33::finance::BondError;

//________________________puts____________________________________

extern const BondError ITO33_BONDLIKE_PUTPERIOD_NOT_STRIKE
("Strike value of PutPeriod is not defined.\n"
"Suggestion: PutPeriod::HasYield() should have been checked before.");

extern const BondError ITO33_BONDLIKE_PUTPERIOD_NOT_YIELD
("Guaranteed yield of PutPeriod is not defined.\n"
"Suggestion: PutPeriod::HasYield() should have been checked before.");

extern const BondError ITO33_BONDLIKE_PUT_SCHEDULE_DUPLICATEDDATE
("PutSchedule definition: this put date has already been "
"added to the put schedule.");

extern const BondError ITO33_BONDLIKE_PUT_SCHEDULE_BADSTRIKE
("PutSchedule definition: the strike must be positive.");

extern const BondError ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_SMALL
("PutSchedule definition: guaranteed yield must be -10% or greater.");

extern const BondError ITO33_BONDLIKE_PUT_SCHEDULE_GUARANTEED_YIELD_TOO_LARGE
("PutSchedule definition: guaranteed yield must be 20% or less.");

extern const BondError ITO33_BONDLIKE_PUT_SCHEDULE_INCONSISTENT_YIELD
("PutSchedule definition: if a yield has been specified, KeepAccrued "
"must be true and ForfeitCoupon must be false.");

extern const BondError ITO33_BONDLIKE_NO_PUT
("No put at given date.");

//_______________Calls_______________________________________________

extern const BondError ITO33_BONDLIKE_CALLS_SOFTAFTERHARD
("Calls definition: The hard call period should "
"follow the soft call period.");

extern const BondError ITO33_BONDLIKE_CALLS_WRONGSCHEDULE
("Calls definition: Wrong call schedule, two call periods "
"overlap.");

extern const BondError ITO33_BONDLIKE_SOFTCALL_BADTRIGGER
("Soft call definition : The call trigger must "
"be a strictly positive number.");

extern const BondError ITO33_BONDLIKE_CALL_BADSTRIKE
("Call definition: The call strike must be positive.");

extern const BondError ITO33_BONDLIKE_CALL_WRONG_START_END
("Call definition: Start date must be for the end date.");

extern const BondError ITO33_BONDLIKE_CALL_NOTICEPERIOD
("Call definition: Call notice period has to be at least one day and "
"less than or equal to 120.");

extern const BondError ITO33_BONDLIKE_CALL_TRIGGERPERIOD
("Call definition: Trigger period must not be longer than 30 days.");

extern const BondError ITO33_BONDLIKE_CALL_INCONSISTENT_YIELD
("Call schedule definition: If any call period has a yield, KeepAccrued "
"must be true and ForfeitCoupon must be false.");

extern const BondError ITO33_BONDLIKE_CALL_STRIKE_AND_YIELD
("Call definition: The call strike must be 1.0 if a yield is defined.");

extern const BondError ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_SMALL
("Call definition: Guaranteed yield to call must be -10% or greater.");

extern const BondError ITO33_BONDLIKE_CALL_GUARANTEED_YIELD_TOO_LARGE
("Call definition: Guaranteed yield must be 20% or less.");

extern const BondError ITO33_BONDLIKE_CALLPERIOD_NOT_STRIKE
("Strike value of CallPeriod is not defined.\n"
"Suggestion: CallPeriod::HasYield() should have been checked before.");

extern const BondError ITO33_BONDLIKE_CALLPERIOD_NOT_YIELD
("Guaranteed yield of CallPeriod is not defined.\n"
"Suggestion: CallPeriod::HasYield() should have been checked before.");

extern const BondError ITO33_BONDLIKE_NO_CALL
("No call at given date.");

//______________make-whole_________________________________________________

extern const BondError ITO33_BONDLIKE_MAKEWHOLETYPE_UNDEFINED
("The bond has no make-whole provision.");

extern const BondError ITO33_BONDLIKE_MAKEWHOLETYPE_NOTPREMIUM
("Only valid for coupon make-whole.");

extern const BondError ITO33_BONDLIKE_MAKEWHOLETYPE_NOTCOUPON
("Only valid for premium make-whole.");

extern const BondError ITO33_BONDLIKE_MAKEWHOLETYPE_INITIAL_PREMIUM
("Initial make-whole premium must be positive.");

//_______________Conversions_____________________________________________

extern const BondError ITO33_BONDLIKE_NO_CONVERSION
("Conversion must be defined for a convertible bond.");

extern const BondError ITO33_BONDLIKE_CONVERSIONS_WRONGSCHEDULE
("Conversions definition: Invalid conversion schedule, "
"two conversion periods overlap.");

extern const BondError ITO33_BONDLIKE_CONVERSIONS_NULL_PERIOD
("Conversions definition: Invalid conversion period.");

extern const BondError ITO33_BONDLIKE_CONVERSION_WRONG_START_END 
("Conversion definition: Start date must be before end date.");

extern const BondError ITO33_BONDLIKE_CONVERSION_BADRATIO 
("Conversion definition: The conversion ratio must be a positive number.");

extern const BondError ITO33_BONDLIKE_COCO_BADTRIGGER
("Contingent Conversion: The conversion trigger must be a positive number.");

extern const BondError ITO33_BONDLIKE_COCO_TYPE 
("Contingent conversion type has not been specified.");

extern const BondError ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE
("Invalid extreme trigger level.");

extern const BondError ITO33_BONDLIKE_COCO_CURRENTLYACTIVE
("The last trigger condition met flag is inconsistent with the CoCo type of "
"the conversion period containing the valuation date.");

extern const BondError ITO33_BONDLIKE_NO_COCO
("Contingent conversion is not specified for this conversion period.");

//__________________Terms______________________________________________

extern const BondError ITO33_BONDLIKE_TERMS_COUPON_ISSUEDATE
("Invalid cash distribution: The contracting date must be equal to "
"the issue date of the instrument");

extern const BondError ITO33_BONDLIKE_TERMS_ISSUEPRICE
("Issue price is negative or too high(Please note that issue price is "
"expressed as percentage of nominal).");

extern const BondError ITO33_BONDLIKE_TERMS_BAD_NOMINAL
("Bond terms definition: Invalid nominal.");

extern const BondError ITO33_BONDLIKE_TERMS_INCONSISTENT_FREQUENCIES
("Yield compounding frequency is not the same as the payement frequency "
"of coupons.");

extern const BondError ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_COMPOUNDING_FREQUENCY 
("Yield compounding frequency is not defined (might be needed for "
"accreting zero coupon or guaranteed yield for call or put).");

extern const BondError ITO33_BONDLIKE_TERMS_INCONSISTENT_DCCS
("Yield day count convention is not the same as the day count convention "
"of coupons.");

extern const BondError ITO33_BONDLIKE_TERMS_UNKNOWN_YIELD_DCC
("Yield day count convention is not defined (might be needed for "
"accreting zero coupon or guaranteed yield for call or put).");

extern const BondError ITO33_BONDLIKE_TERMS_MUSTBE_ACCRETING
("Failed to get accreting bond yield for non-accreting bond."
"Suggestion: The function IsAccretingBond() should have been "
"checked before.");

extern const BondError ITO33_BONDLIKE_TERMS_MUSTBE_CASHPAYTOZERO
("Failed to get AccretionRateOfCashPayToZero for non "              
"cash-pay-to-zero bond.\n"                                         
"Suggestion: The function IsCashPayToZeroBond() should have been " 
"checked before.");

extern const BondError ITO33_BONDLIKE_TERMS_EITHER_OID_OR_CASHPAYTOZERO
("Bond terms definition: Can't call both SetAccretingBond and "
"SetCashPayToZero.");

//__________________Bond___________________________________________________
extern const BondError ITO33_BONDLIKE_BOND_SOFTCALL
("Soft call is not permitted for a non convertible bond.");

//_________________Cross Currency__________________________________________
extern const BondError ITO33_BONDLIKE_CROSSCURRENCY_FIXEDFXRATE
("The fixed FX rate must be positive.");

extern const BondError ITO33_BONDLIKE_FIXEDFXRATE_NOT_SET
("The fixed FX rate is required but not set before!");

extern const BondError ITO33_BONDLIKE_NOT_CROSSCURRENCY
("Information specific to a cross currency is required but the security "
"is not a cross currency.");

//__________________ConvertibleLike________________________________________

extern const BondError ITO33_BONDLIKE_NULL_BONDLIKETERMS
("ConvertibleLike definition: invalid BondLikeTerms.");

extern const BondError ITO33_BONDLIKE_NOT_EXCHANGEABLE
("Querying properties relevant to exchangeable bond for a non exchangeable bond.");

extern const BondError ITO33_EXCHANGEABLE_NO_ISSUER_DEFAULT
("Exchangeable bond requires the default of its issuer, but either the "
"issuer is not set or the default intensity is not set for the issuer.");

extern const BondError ITO33_BONDLIKE_NEWSHARE_AND_EXCHANGEABLE
("Exchangeable bond with new share feature is not supported!");

extern const BondError ITO33_BONDLIKE_TRIGGERASPERCENTAGEOF 
("The trigger must be represented as a percentage "
"of accreted principal or fixed price.");

extern const BondError ITO33_BONDLIKE_TRIGGERINCURRENCYOF 
("The trigger must be in the derivative currency "
"or the underlying share currency");

extern const BondError ITO33_BONDLIKE_FIXEDQUANTO_FXVOLATILITY 
("The FX volatility can't be negative or greater than 500%!");

extern const BondError ITO33_BONDLIKE_FIXEDQUANTO_CORRELATION
("The correlation must be between -100% and 100%.");

extern const BondError ITO33_BONDLIKE_NOT_FIXEDQUANTO 
("Information specific to a fixed quanto is required but the security "
"is not a fixed quanto.");

extern const BondError ITO33_BONDLIKE_FIXED_QUANTO_NOT_CROSSCURRENCY
("Fixed quanto must be a cross-currency instrument.");

//___________________bonds_________________________________________________
extern const BondError ITO33_BONDLIKE_NULL_BONDTERMS
("Invalid BondTerms.");

extern const BondError ITO33_BONDLIKE_NULL_CALLSCHEDULE
("Invalid call schedule.");

extern const BondError ITO33_BONDLIKE_NULL_PUTSCHEDULE
("Invalid put schedule.");

extern const BondError ITO33_BONDLIKE_MISSING_CASHFLOWS
("At least one payment (cashflow) is required if accretion rate is set.");

extern const BondError ITO33_BONDLIKE_ISSUEPRICE_TOO_SMALL
("The issue price must be at least 94% for non-OID (including partial "
"OID) bonds.");

extern const BondError ITO33_BONDLIKE_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION
("First unknown coupon date of the instrument must be greater than the "
"valuation date.");

extern const BondError ITO33_BONDLIKE_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE
("Last known coupon date or start of accrued date of the floating rates "
"of the instrument must be greater than the valuation date + fixing delay.");

//___________________CBBase________________________________________________
extern const BondError ITO33_BONDLIKE_CANNOT_READ_BONDLIKETERMS
("Accessing BondLikeTerms of a bond security. Please use BondTerms instead.");

extern const BondError ITO33_BONDLIKE_STRAIGHTBOND_FOR_EXCHANGEABLE
("Straight bond for exchangeable CB does not make sense.");

//___________________Convertible bond______________________________________
extern const BondError ITO33_BONDLIKE_NULL_CONVERSION
("ConvertibleBond definition: Invalid conversion schedule.");

//_________________Reset___________________________________________________
extern const BondError ITO33_BONDLIKE_RESET_AND_COCO
("Contingent conversion is not allowed for a reset convertible bond.");

extern const BondError ITO33_BONDLIKE_RESET_AND_COCALL
("Contingent calls are not allowed for a reset.");

extern const BondError ITO33_BONDLIKE_RESET_BEFORE_ISSUE
("Reset definition: Conversion start date cannot be before the issue date.");

extern const BondError ITO33_BONDLIKE_RESET_AFTER_MATURITY
("Reset definition: Conversion end date cannot be after the maturity date.");

extern const BondError ITO33_BONDLIKE_RESET_DATE_AFTER_MATURITY
("Reset definition: Cannot have a reset date on or after the maturity date.");

extern const BondError ITO33_BONDLIKE_EMPTY_RESET_DATES
("Reset definition: Reset dates not specified.");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_CAP_RATE
("Reset definition: Invalid cap rate. 100% <= cap <= 1000%.");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_FLOOR_RATE
("Reset definition: Invalid floor rate. 0 <= floor <= 100%.");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_MULTIPLIER
("Reset definition: Invalid multiplier. 0 <= multiplier <= 1000%.");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_PERIOD
("Reset schedule definition: Invalid conversion period. The start date must "
"be before the end date");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_DATE
("Reset schedule definition: Invalid reset date. The reset dates must "
"be unique and within the conversion start and end dates");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_INITIAL_CONV_PRICE
("Reset schedule definition: Invalid initial conversion price." 
" Price must be positive");

extern const BondError ITO33_BONDLIKE_INVALID_RESET_CURRENT_CONV_PRICE
("Reset schedule definition: Invalid current price. Price must be positive");

extern const BondError ITO33_BONDLIKE_RESET_INCONSISTENCY_INITIAL_CURRENT_PRICE
("Reset definition: Current conversion price should be equal to "
"initial conversion price "
"as valuation date is before the first reset date.");

extern const BondError ITO33_BONDLIKE_WITH_FLOATING_RATES_AND_CASHFLOWS
("Bond-like security terms definition: Floating rate bond must not have " 
"a cash flow stream.");

extern const BondError ITO33_BONDLIKE_WITH_FLOATING_RATES_AFTER_MATURITY
("Bond-like definition: last payment date of the floating rate must be "
"lower than (or equal to) the  maturity date.");

extern const BondError ITO33_BONDLIKE_NULL_RESET_CONVERSION_SCHEDULE
("Reset definition: Invalid reset schedule.");

extern const BondError ITO33_BONDLIKE_RESET_FIXEDRATE_NOT_SET
("Reset definition: Reset is cross currency, as such a fixed rate "
"must be specified.");

//___________________________mandatories___________________________________
extern const BondError ITO33_PEPSLIKE_CALL_ALREADY_SET 
("PEPS-like definition: Call provision has already been defined. \n" 
"Note that either fixed cash or fixed share call provision is allowed.");

extern const BondError ITO33_PEPSLIKE_CALLFIXEDSHARE_RATIO_NEGATIVE
("PEPS-like definition: ratio of call FixedShare must be positive.");

extern const BondError ITO33_PEPSLIKE_CALL_FIXED_SHARE_RATIO
("PEPS-like definition: the ratio of optional conversion at issuer's "
"option must be the minimum conversion ratio of the security.");

extern const BondError ITO33_PEPSLIKE_CALL_FIXED_SHARE_TRIGGER 
("PEPS-like definition: the trigger of optional conversion at issuer's "
"option must be greater than or equal to 100%.");

extern const BondError ITO33_MANDATORY_GUARANTEED_YIELD_TO_CALL_INVALID
("Guaranteed yield to call is not supported.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_CALL_TRIGGER
("Generalized PEPS-like call definition: trigger must be greater than or "
"equal to 100%.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_CALL_ALREADY_SET
("Generalized PEPS-like definition: "
"Call provision has already been defined. \n"
"Note that either fixed cash call or generalized PEPS-like "
"call is allowed.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_CALL_TYPE
("Generalized PEPS-like definition: Invalid call type.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_DOWNSIDE_CONVERSION_RATIO
("PEPS definition: The downside conversion ratio at maturity must be" 
"a positive number.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_UPSIDE_CONVERSION_RATIO
("PEPS definition: The upside base conversion ratio must be a postive "
"number.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_LOWER_HIGHER_STRIKE
("PEPS definition: lower strike should be positive and be less than "
"higher strike.");

extern const BondError ITO33_NULL_GENERALIZEDPEPSLIKE_CALL
("Invalid generalized PEPS-like call.");

extern const BondError ITO33_PERCSLIKE_CAP_PRICE
("PERCS definition: The cap price must be a positive number.");

extern const BondError ITO33_PERCSLIKE_MAX_CONVERSION_RATIO
("PERCS definition: The maximum conversion ratio must be a positive number.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_TOO_MANY_DAYS
("PEPS definition: The number of days for the monitoring period "
"is too large.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_NOT_STOCK_AVERAGING
("PEPS definition: Requesting stock average while monitoring "
"the conversion ratio.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_NOT_CONVERSION_RATIO_AVERAGING
("PEPS definition: Requestion conversion ratio average while "
"monitoring the stock average.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_NEGATIVE_AVERAGE
("PEPS definition: Setting a negative average.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_NO_CURRENT_AVERAGE_SET
("PEPS definition: Current Average should be set.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_INVALID_NBSAMPLINAVERAGES
("PEPS definition: There must be at least one sampling average.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_AVG_PERIOD_START_DATE_AFTER_END_DATE
("PEPS definition: The averaging period start date is after "
"the averaging period end date.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_INVALID_NB_SAMPLES_USED
("PEPS definition: the number of samples used to compute the " 
"current average must be less than the total number of sampling " 
"averages.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_ZERO_NB_SAMPLES_USED
("PEPS definition: the number of samples used to compute the " 
"current average must be at least one.");

extern const BondError ITO33_GENERALIZEDPEPSLIKE_AVERAGINGPERIOD_AVG_END_DATE_INVALID
("PEPS definition: the number of days between the end of the average "
"period and the maturity date can not be greater than five.");

//________________________attached warrant convertible bond________________
extern const BondError ITO33_AWCB_RESET_BEFORE_ISSUE
("Attached warrant convertible bond definition: "
"Conversion start date cannot be before the issue date.");

extern const BondError ITO33_AWCB_RESET_AFTER_MATURITY
("Attached warrant convertible bond definition: "
"Conversion end date cannot be after the maturity date.");

extern const BondError ITO33_AWCB_CROSSCURRENCY_UNDEFINED_STRIKE
("Attached warrant convertible bond definition: "
"Fixed strike must be set in case of cross currency.");

extern const BondError ITO33_AWCB_NULL_CONVERSION
("Attached warrant convertible bond definition: Invalid share dependent " 
"conversion specification.");

extern const BondError ITO33_AWCB_INVALID_FUNCTION_CALL_TO_GET_RESET_DATE
("Attached warrant convertible bond definition: Reset date has not been "
"defined, invalid function call.");

extern const BondError ITO33_AWCB_INVALID_CALL_NOTICE
("Attached warrant convertible bond definition: Call notice is not "
"supported when there is a reset date.");

extern const BondError ITO33_AWCB_NO_CURRENT_CONVERSION_RATIO_SET
("Attached warrant convertible bond definition: current conversion ratio "
"has not been set.");

//_____________________share dependent conversion__________________________
extern const BondError ITO33_BONDLIKE_SHAREDEPENDENT_BADCAP
("Conversion definition: Cap ratio must be greater than or equal to the " 
"initial ratio.");

extern const BondError ITO33_BONDLIKE_SHAREDEPENDENT_BADRESET 
("Conversion definition: Reset date must be after the start date and " 
"before the end of the conversion date.");

extern const BondError ITO33_BONDLIKE_SHAREDEPENDENT_BADSTRIKE
("Conversion definition: The share strike must be a positive number.");

extern const BondError ITO33_BONDLIKE_SHAREDEPENDENT_BADINCREMENTALSHAREFACTOR
("Conversion definition: The share factor must be strictly positive.");

//___________________________yield computation_____________________________
extern const BondError ITO33_BONDLIKE_COMPUTE_YTP_INEXIST
("Can't compute Yield-to-put: price too high or too low");

extern const BondError ITO33_BONDLIKE_COMPUTE_YTP_NOPUT
("Can't compute yield-to-put: the contract doesn't have put clause.");

extern const BondError ITO33_BONDLIKE_COMPUTE_YTP_NOPUT_AFTER_VALUATIONDATE
("Can't compute yield-to-put: no put after valuation date.");

extern const BondError ITO33_BONDLIKE_COMPUTE_YTM_INEXIST
("Can't compute yield-to-maturity: price too high or too low");

//___________________________CB Option_____________________________________
extern const BondError ITO33_CBOPTION_INVALID_FLOATINGRATES
("CBOption definition: Invalid floating rates.");

extern const BondError ITO33_CBOPTION_WITH_FLOATING_RATES_AFTER_MATURITY
("CBOption definition: maturity date must be greater "
"than (or equal to) the last payment date of the floating leg.");

extern const BondError ITO33_CBOPTION_INCOHERENT_ASWNOTIONAL
("CBOption definition: ASWNotional can't be the put price since there "
"is no put at the maturity date of the cb option.");

extern const BondError ITO33_CBOPTION_FIRST_UNKNOWN_PAYMENT_BEFORE_VALUATION
("First unknown payment date of the swap of the cb option must be "
"greater than the valuation date.");

extern const BondError ITO33_CBOPTION_PREVIOUS_REF_DATE_BEFORE_VALUATION_DATE
("Last known payment date or start of accrued date of the floating rates "
"of the swap of the cb option must be greater than the valuation date + " 
"fixing delay.");

extern const BondError ITO33_CBOPTION_AND_CB_NOT_SAME_CURRENCY
("CBOption should have the same currency as the underlying "
"convertible bond.");

extern const BondError ITO33_CBOPTION_AND_CB_NOT_SAME_SESSIONDATA
("CBOption should have the same session data as the underlying "
"convertible bond.");

extern const BondError ITO33_CBOPTION_WITH_BAD_MATURITY
("CBOption definition: The maturity date of the option must be either "
"the maturity date of the cb or a put date.");

extern const BondError ITO33_CBOPTION_WITH_INVALID_CB
("CBOption definition: Invalid convertible bond.");

extern const BondError ITO33_CBOPTION_WITH_FLOATING_CB
("CBOption definition: The convertible bond of the option "
"must be a non floating one.");

