#ifndef IRX_RATE_H
#define IRX_RATE_H

#include "irxflow.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Special values of rate type */
extern IrxTRateType IRX_CONTINUOUS_RATE;
extern IrxTRateType IRX_ANNUAL_RATE;
extern IrxTRateType IRX_SIMPLE_RATE;
extern IrxTRateType IRX_DISCOUNT_RATE;
extern IrxTRateType IRX_COMPOUND_RATE(int frequency);
    
/**
 * Converts a rate type from a string.
 *
 * S=simple, C=continuous, D=discount rate, Fn=frequency n
 */
IrxTRateType* irxRateTypeFromString (char* str);

/**
 * Converts a rate type from an integer.
 *
 * 0=simple, >365=continuous, 1-365 = given frequency
 */
IrxTRateType* irxRateTypeFromLong (long freq);

/**
 * Validates a rate type.
 *
 * For the frequency rate type, the frequency must be 1,2,4,12 or 365.
 * For other rate types, the frequency is ignored and is required to be 0.
 */
int irxRateTypeValidate (IrxTRateType *rateType);

/**
 * Converts from a rate to a discount factor given dates, day count
 * convention and rate type.
 */
extern int irxRateToDiscount (
    double           rate,
    IrxTDate         startDate,
    IrxTDate         endDate,
    IrxTDayCountConv dcc,
    IrxTRateType     rateType,
    double          *discount);
    
/**
 * Converts from a discount factor to a rate given dates, day count
 * convention and rate type.
 *
 * Fails when the discount is not strictly positive or the startDate and
 * endDate are the same.
 */
extern int irxDiscountToRate (
    double           discount,
    IrxTDate         startDate,
    IrxTDate         endDate,
    IrxTDayCountConv dcc,
    IrxTRateType     rateType,
    double          *rate);
    
    
/**
 * Converts from a discount factor to a rate given year fraction and 
 * rate type.
 *
 * Fails when the discount is not strictly positive or the yearFraction
 * is zero.
 */
extern int irxDiscountToRateYearFrac (
    double       discount,
    double       yearFraction,
    IrxTRateType rateType,
    double      *rate);

/**
 * Converts from a rate to a discount facto given year fraction and rate type.
 */
extern int irxRateToDiscountYearFrac(
    double       rate,
    double       yearFraction,
    IrxTRateType rateType,
    double      *discount);

/**
 * Validates a rate. Rates can sometimes be negative, but this will depend
 * on the rate type. For example, an annually compounded rate cannot be
 * less than -1.
 */
int irxRateValid(
    const char      *routine,
    double           rate,  
    IrxTDate         startDate,
    IrxTDate         endDate,
    IrxTDayCountConv rateDayCountConv,
    IrxTRateType     rateType);

/**
 * Validates a rate. Rates can sometimes be negative, but this will depend
 * on the rate type. For example, an annually compounded rate cannot be
 * less than -1.
 */
int irxRateValidYearFrac(
    const char  *routine,
    double       rate,
    double       yearFraction,
    IrxTRateType rateType);


#ifdef __cplusplus
}
#endif

#endif
