#ifndef USD_MID_CURVE_H
#define USD_MID_CURVE_H
#include "utDates.h"

// ---------------------------------------------------------------------------------- //
// A Eurodollar option pricer
//
// return value:  Eurodollar Option Price and implied vol
//
//	mkt:				for getting fwds and vols
//	lOptionExpiry		expiry date of the option.  
//	dStrike				strike of the option
//	lFutureExpiry		start date for the Libor underlying the option
//	dFuturePrice		market price of the future
//	enumCallPut			Call = long
//	enumGreek			PREMIUM, DELTA, etc
//	out_dVol			implied vol 
//
char* EurodollarOptionCaller(  // For fwds and vols
						long lToday,
						char* szYieldCurveName,
						char* szVolCurve,
						char* szRefRate,
						char* (*getCashVol)(	char *szVolCurve,
													double	start_date, 
													double	end_date,
													double	cash_strike,
													int		zero,
													char	*ref_rate_name,
													double	*vol,
													double	*power),
						char* (*getDF)(char *szYieldCurve, double dStart, double dEnd, double *dDF),
						char* (*getSpread)(long start_date,long end_date,const char *szRefRate, double *dSpread),
						long lOptionExpiry,
						double dStrike,
						long lFutureExpiry,
						double dFuturePrice,
						int iConvexityModel,
						SrtDiffusionType enumDiffType,
						SrtCallPutType enumCallPut, // Call = long
						SrtGreekType enumGreek,
						double dDifVol,
						double dJumpAvg,
						double dJumpVol,
						double dJumpInt,
						long lIterNo,
						double* out_dPrice,
						double* out_dVol
						);

char* EurodollarOption( 
						long lToday,
						double dOptionExpiry,
						double dCashStrike,
						double dCashRate,
						SrtDiffusionType enumDiffType,
						SrtCallPutType enumCallPut, // Call = long
						SrtGreekType enumGreek,
						double dDifVol,
						double dJumpAvg,
						double dJumpVol,
						double dJumpInt,
						long lIterNo,
						double dOptionExpiryDf,
						double* out_dPrice,
						double* out_dVol
						);

#endif