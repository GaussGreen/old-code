#ifndef USD_STRUCT_SWAP_H
#define USD_STRUCT_SWAP_H

#include "utDates.h"
#include "utgreeks.h"

Err BMA_Option(
    double         dToday,      /* Today's date (in number format) */
    double         dStartDate1, /* The start date of the first Libor fixing (in number format) */
    double         dEndDate1,   /* The end date of the first Libor fixing (in number format) */
    double         dStartDate2, /* The start date of the last Libor fixing (in number format) */
    double         dEndDate2,   /* The end date of the last Libor fixing (in number format) */
    double         dBMAfwd,     /* The BMA forward for the option period */
    double         dBMAratio,   /* The ratio BMA / Libor for the option period */
    double         dLevel,      /* The level for the option (Disc*Coverage for a caplet) */
    double         dStrike,     /* The unadjusted strike of the option */
    SrtCallPutType srtCallPut,  /* Whether it is pay/rec */
    int            bIsACap,     /* 1 for a cap, 0 for a swaption (not implemented) */
    int            bUseNormalModel, /* 1 for combining normal vols, 0 for using lognormal vols */
    char*          szVolCurveName,  /* VolCurve name passed to GetVol function */
    char*          szRefRate,       /* RefRate passed to GetVol function */
    Err (*GetVol)(
        char*, double, double, double, int, char*, double*, double*), /* Function to get vol */
    double* dPrice,
    double* out_dVega); /* Return value of the option price */

Err getExerciseDate(long lStart, char* szRefRate, SrtBusDayConv enumBusDayConv, long* lExercise);

Err getEndDate(long lStart, char* szRefRate, SrtBusDayConv enumBusDayConv, long* lEnd);

// ---------------------------------------------------------------------------------- //
// A MidCurve option pricer
//
// return value:  MidCurve Option Price
//
//	mkt:				for getting fwds and vols
//	lOptionStart		delivery date of the option.  exercise date calculated using lag
//	dStrike				strike of the option
//	lFutureStart		start date for the Libor underlying the midcurve option
//	dFuturePrice		market price of the future
//	iFutureTenorInMonths	length of underlying Libor in months (generally 1 or 3)
//	iConvexityModel		0:	midcurve strike adjusted by forward-future difference
//	iVolatilityModel	0:  stationary normal volatility
//	enumCallPut			Call = long
//	enumGreek			PREMIUM, DELTA, etc
//	out_dVol			vol used to price.  Will depend on the volatility model
//
char* MidCurveCaller(  // For fwds and vols
    long  lToday,
    char* szYieldCurveName,
    char* szVolCurve,
    char* szRefRate,
    char* (*getCashVol)(
        char*   szVolCurve,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),
    char* (*getDF)(char* szYieldCurve, double dStart, double dEnd, double* dDF),
    char* (*getSpread)(long start_date, long end_date, const char* szRefRate, double* dSpread),
    long           lOptionExpiry,
    double         dStrike,
    long           lFutureExpiry,
    long           lFutureStart,
    double         dFuturePrice,
    int            iConvexityModel,
    int            iVolatilityModel,
    SrtCallPutType enumCallPut,  // Call = long
    SrtGreekType   enumGreek,
    double         dLongCashFwd,
    double         dShortCashFwd,
    double         dLongCashVol,
    double         dShortCashVol,
    double*        out_dPrice,
    double*        out_dVol);

char* MidCurveOption(
    long           lToday,
    long           lOptionExpiry,
    double         dStrike,
    long           lFutureExpiry,
    double         dFuturePrice,
    int            iConvexityModel,
    int            iVolatilityModel,
    SrtCallPutType enumCallPut,  // Call = long
    SrtGreekType   enumGreek,
    double         dLongCashFwd,
    double         dLongSpread,
    double         dLongCashVol,
    double         dShortCashVol,
    double*        out_dPrice,
    double*        out_dVol);

#endif