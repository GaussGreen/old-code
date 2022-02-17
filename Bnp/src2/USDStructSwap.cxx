/* -------------------------------------------------------------------------------------------------------------------------------------------------
 */
/* USDStructSwap.cxx */
/* Contains functions: */
/*		BMA_Option() */

/* Include files */
#include "USDStructSwap.h"
#include "math.h"
#include "opFnctns.h"
#include "srt_h_all.h"
#include "swp_h_all.h"
#include "swp_h_ccy_param.h"
#include "swp_h_spread.h"
#include "utTypes.h"
#include "utallhdr.h"

Err BMA_Option(double dToday, double dStartDate1, double dEndDate1,
               double dStartDate2, double dEndDate2, double dBMAfwd,
               double dBMAratio, double dLevel, double dStrike,
               SrtCallPutType srtCallPut, int bIsACap, int bUseNormalModel,
               char *szVolCurveName, char *szRefRate,
               Err (*GetVol)(char *, double, double, double, int, char *,
                             double *, double *),
               double *out_dPrice, double *out_dVega) {
  /* Variable declaration */
  double dVol1, dVol2, dEffStrike, dLibor, dExpiry1, dExpiry2, dPower, dEffVol;
  long lExercise;
  SrtDiffusionType srtPricing, srtVol;
  Err err;
  SrtBusDayConv enumBusDayConv;
  SrtCcyParam *ptrCcyParam;
  char szCcy[10] = "USD";

  /* Data checking */
  if (bIsACap != 1)
    return "BMA_Option is not implemented for swaptions.";

  if (dToday > dEndDate2) {
    *out_dPrice = 0.0;
    return 0;
  }

  /* get the parameters */
  ptrCcyParam = new_CcyParam();
  if (err = swp_f_get_CcyParam_from_CcyStr(szCcy, &ptrCcyParam)) {
    if (ptrCcyParam)
      free_CcyParam(ptrCcyParam);
    return err;
  }
  enumBusDayConv = ptrCcyParam->swap_bus_day_conv;
  free_CcyParam(ptrCcyParam);

  /* Calculate the expiry times */
  if (err = getExerciseDate((long)dStartDate1, szRefRate, enumBusDayConv,
                            &lExercise))
    return err;
  dExpiry1 = (lExercise - dToday) / 365.0;
  dExpiry1 = dExpiry1 < 0.0 ? 0.0 : dExpiry1;

  if (err = getExerciseDate((long)dStartDate2, szRefRate, enumBusDayConv,
                            &lExercise))
    return err;
  dExpiry2 = (lExercise - dToday) / 365.0;
  dExpiry2 = dExpiry2 < 0.0 ? 0.0 : dExpiry2;

  /* Calculate the effective strike */
  dEffStrike = dStrike / dBMAratio;
  dLibor = dBMAfwd / dBMAratio;

  /* Get the volatilities */
  if (dStartDate1 > dToday) {
    if (err = (*GetVol)(szVolCurveName, dStartDate1, dEndDate1, dEffStrike,
                        bIsACap, szRefRate, &dVol1, &dPower))
      return err;
  } else
    dVol1 = 0.0;

  if (dStartDate2 > dToday) {
    if (err = (*GetVol)(szVolCurveName, dStartDate2, dEndDate2, dEffStrike,
                        bIsACap, szRefRate, &dVol2, &dPower))
      return err;
  } else
    dVol2 = 0.0;

  /* Get the vol types */
  if (bUseNormalModel == 1)
    srtPricing = SRT_NORMAL;
  else
    srtPricing = SRT_LOGNORMAL;

  if (dPower == 1)
    srtVol = SRT_LOGNORMAL;
  else
    srtVol = SRT_NORMAL;

  /* Convert the vols to the requested model type */
  srt_f_optsarbvol(dLibor, dEffStrike, dExpiry1, dVol1, 0.0, dPower, 0.0,
                   srtVol, srtPricing, &dVol1);
  srt_f_optsarbvol(dLibor, dEffStrike, dExpiry2, dVol2, 0.0, dPower, 0.0,
                   srtVol, srtPricing, &dVol2);

  /* Calculate the effective vol*/
  dEffVol = sqrt(2.0 * dVol1 * dVol1 * dExpiry1 / 3.0 / dExpiry2 +
                 dVol2 * dVol2 / 3.0);

  /* calculate the price */
  if (bUseNormalModel == 1)
    *out_dPrice = srt_f_optblknrm(dLibor, dEffStrike, dEffVol, dExpiry2, dLevel,
                                  srtCallPut, PREMIUM);
  else
    *out_dPrice = srt_f_optblksch(dLibor, dEffStrike, dEffVol, dExpiry2, dLevel,
                                  srtCallPut, PREMIUM);
  *out_dPrice *= dBMAratio;

  /* Bump the vols        , recalculate the price and hence obtain the vega */
  if (bUseNormalModel == 1) {
    dVol1 = dVol1 + 0.01 * dLibor;
    dVol2 = dVol2 + 0.01 * dLibor;
    dEffVol = sqrt(2.0 * dVol1 * dVol1 * dExpiry1 / 3.0 / dExpiry2 +
                   dVol2 * dVol2 / 3.0);
    *out_dVega = srt_f_optblknrm(dLibor, dEffStrike, dEffVol, dExpiry2, dLevel,
                                 srtCallPut, PREMIUM);
  } else {
    dVol1 = dVol1;
    dVol2 = dVol2;
    dEffVol = sqrt(2.0 * dVol1 * dVol1 * dExpiry1 / 3.0 / dExpiry2 +
                   dVol2 * dVol2 / 3.0);
    *out_dVega = srt_f_optblksch(dLibor, dEffStrike, dEffVol, dExpiry2, dLevel,
                                 srtCallPut, PREMIUM);
  }
  *out_dVega = dBMAratio * (*out_dVega) - *out_dPrice;

  /* end */
  return 0;
}

Err getExerciseDate(long lStart, char *szRefRate, SrtBusDayConv enumBusDayConv,
                    long *lExercise) {
  // Get the lag from the refrate
  int iLag;
  Err err;
  if (err = srt_f_get_spot_lag_from_refrate(szRefRate, &iLag))
    return err;

  // Calculate the exercise date
  *lExercise = add_unit(lStart, -iLag, SRT_BDAY, enumBusDayConv);
  return 0;
}

Err getEndDate(long lStart, char *szRefRate, SrtBusDayConv enumBusDayConv,
               long *lEnd) {
  // Get the compounding from the refrate
  Err err;
  SrtBasisCode basis;
  SrtCompounding compounding;
  if (err = swp_f_get_ref_rate_details(szRefRate, &basis, &compounding))
    return err;

  *lEnd = add_unit(lStart, 12 / compounding, SRT_MONTH, enumBusDayConv);
  return 0;
}

Err getCoverage(long lStart, long lEnd, char *szRefRate, double *dCoverage) {
  // Get the compounding from the refrate
  Err err;
  SrtBasisCode basis;
  SrtCompounding compounding;
  if (err = swp_f_get_ref_rate_details(szRefRate, &basis, &compounding))
    return err;

  *dCoverage = coverage(lStart, lEnd, basis);
  return 0;
}

char *getForward(long lStart, char *szYieldCurveName, char *szRefRate,
                 char *(*getDF)(char *szYieldCurve, double dStart, double dEnd,
                                double *dDF),
                 char *(*getSpread)(long start_date, long end_date,
                                    const char *szRefRate, double *dSpread),
                 double *out_dCashFwd, double *out_dSpread) {
  // Variable declarations
  SrtCurvePtr srtCurvePtr;
  SrtCcyParam *pCcyParam;
  char *err;
  long lEnd;
  double dRatioDF, dCvg;

  // Get the yield curve pointer
  if (!(srtCurvePtr = swp_f_lookup_curve(szYieldCurveName)))
    return "Could not find yield curve";

  // Get the market conventions
  if (err = swp_f_get_CcyParam_from_CcyStr(srtCurvePtr->curve_ccy, &pCcyParam))
    return err;

  // Get the end date
  if (err = getEndDate(lStart, szRefRate, pCcyParam->swap_bus_day_conv, &lEnd))
    return err;

  // Calculate the spread
  if (err = (*getSpread)(lStart, lEnd, szRefRate, out_dSpread))
    return err;

  // Calculate the discount factors
  if (err = (*getDF)(szYieldCurveName, (double)lStart, (double)lEnd, &dRatioDF))
    return err;

  // Get the coverage
  if (err = getCoverage(lStart, lEnd, szRefRate, &dCvg))
    return err;

  // Calculate the cash forward
  *out_dCashFwd = (dRatioDF - 1.0) / dCvg;

  return 0;
}

char *getCashVolwithType(
    char *szVolCurve, long lToday, long lExpiry, long lStart, double lEnd,
    double dCashFwd, double dCashStrike, char *szRefRate, int iOutVol,
    char *(*getCashVol)(char *szVolCurve, double start_date, double end_date,
                        double cash_strike, int zero, char *ref_rate_name,
                        double *vol, double *power),
    double *out_dCashVol) {
  double dPower;
  SrtDiffusionType enumInVol, enumOutVol;
  char *err;

  // Get the vol
  if (err = getCashVol(szVolCurve, (double)lStart, (double)lEnd, dCashStrike, 0,
                       szRefRate, out_dCashVol, &dPower))
    return err;

  // Check that the volatility is consistent with the model        , otherwise
  // convert it as required
  if (dPower == 0.0)
    enumInVol = SRT_NORMAL;
  else
    enumInVol = SRT_LOGNORMAL;

  if (iOutVol == 0.0)
    enumOutVol = SRT_NORMAL;
  else
    enumOutVol = SRT_LOGNORMAL;

  if (enumInVol != enumOutVol) {
    double dExpiry;
    dExpiry = (lExpiry - lToday) / 365.0;
    if (err =
            srt_f_optsarbvol(dCashFwd, dCashStrike, dExpiry, *out_dCashVol, 0.0,
                             dPower, 0.0, enumInVol, enumOutVol, out_dCashVol))
      return err;
  }

  return 0;
}

// ----------------------------------------------------------------------------------
// // A MidCurve wrapper
//
// return value:  MidCurve Option Price
//
//	mkt:				for getting fwds and vols
//	lOptionStart		delivery date of the option.  exercise date
// calculated using lag 	dStrike				strike of the
// option 	lFutureStart		start date for the Libor underlying the midcurve
// option 	dFuturePrice		market price of the future
// iFutureTenorInMonths length of underlying Libor in months (generally 1 or 3)
// iConvexityModel 0:	midcurve strike adjusted by forward-future difference
// iVolatilityModel 0:  stationary normal volatility 	enumCallPut
// Call = long 	enumGreek			PREMIUM        , DELTA        , etc
//	out_dVol			vol used to price.  Will depend on the
// volatility model
//

char *MidCurveCaller( // For fwds and vols
    long lToday, char *szYieldCurveName, char *szVolCurve, char *szRefRate,
    char *(*getCashVol)(char *szVolCurve, double start_date, double end_date,
                        double cash_strike, int zero, char *ref_rate_name,
                        double *vol, double *power),
    char *(*getDF)(char *szYieldCurve, double dStart, double dEnd, double *dDF),
    char *(*getSpread)(long start_date, long end_date, const char *szRefRate,
                       double *dSpread),
    long lOptionExpiry, double dStrike, long lFutureExpiry, long lFutureStart,
    double dFuturePrice, int iConvexityModel, int iVolatilityModel,
    SrtCallPutType enumCallPut, // Call = long
    SrtGreekType enumGreek, double dLongCashFwd, double dShortCashFwd,
    double dLongCashVol, double dShortCashVol, double *out_dPrice,
    double *out_dVol) {
  // Variable declaration
  double dLongSpread, dCashStrike;
  long lFutureEnd;
  SrtCurvePtr srtCurvePtr;
  SrtCcyParam *pCcyParam;
  SrtBasisCode basis;
  SrtCompounding compounding;
  /*	double dLongFwd;
          double dLongStrike;
          double dLongNormVol;
          double dLongExpiry;
          double dLongVar;
          double dUndStrike;
  */
  char *err;

  // Initialize return values
  *out_dPrice = 0.0;
  *out_dVol = 0.0;

  // Check that the volatility model is allowed
  if (iVolatilityModel != 0 && iVolatilityModel != 1)
    return "Unknown volatility model should be 0:norm or 1:lognoormal";

  // Get the yield curve pointer
  if (!(srtCurvePtr = swp_f_lookup_curve(szYieldCurveName)))
    return "Could not find yield curve";

  // Get the market conventions
  if (err = swp_f_get_CcyParam_from_CcyStr(srtCurvePtr->curve_ccy, &pCcyParam))
    return err;

  // Get the compounding from the refrate
  if (err = swp_f_get_ref_rate_details(szRefRate, &basis, &compounding))
    return err;

  // Calculate the end date of the fwd
  lFutureEnd = add_unit(lFutureStart, 12 / compounding, SRT_MONTH,
                        pCcyParam->swap_bus_day_conv);

  // Calculate the spread
  if (err = (*getSpread)(lFutureStart, lFutureEnd, szRefRate, &dLongSpread))
    return err;

  // Get the long forward if it is not passed in
  if (dLongCashFwd < 0.0) {
    double dRatioDF, dCvg;
    if (err = (*getDF)(szYieldCurveName, (double)lFutureStart,
                       (double)lFutureEnd, &dRatioDF))
      return err;
    dCvg = coverage(lFutureStart, lFutureEnd, basis);
    dLongCashFwd = (1.0 / dRatioDF - 1.0) / dCvg;
  }

  // Calculate the modified strike depending on the convexity model
  if (iConvexityModel == 0)
    dCashStrike = dStrike - dLongSpread;
  else if (iConvexityModel == 1)
    dCashStrike = dLongCashFwd + (dFuturePrice - dStrike) / 100.0;
  else
    return "Unknown Convexity Model";

  // Get the long vol if it is not passed in
  if (dLongCashVol < 0.0) {
    if (err =
            getCashVolwithType(szVolCurve, lToday, lFutureExpiry, lFutureStart,
                               lFutureEnd, dLongCashFwd, dCashStrike, szRefRate,
                               iVolatilityModel, getCashVol, &dLongCashVol))
      return err;
  }

  // Get the short vol if it is not passed in
  if (dShortCashVol < 0.0) {
    double dATMLongVol, dATMShortVol, dStdDev, dLongExpiry, dShortCashStrike,
        dShortExpiry;
    long lShortExpiry, lShortStart, lShortEnd;
    lShortExpiry = lFutureExpiry - lOptionExpiry + lToday;
    lShortStart = add_unit(lShortExpiry, pCcyParam->spot_lag, SRT_BDAY,
                           pCcyParam->swap_bus_day_conv);
    lShortEnd = add_unit(lShortStart, 12 / compounding, SRT_MONTH,
                         pCcyParam->swap_bus_day_conv);

    // Get the short forward if it is not passed in
    if (dShortCashFwd < 0.0) {
      double dRatioDF, dCvg;
      if (err = (*getDF)(szYieldCurveName, (double)lShortStart,
                         (double)lShortEnd, &dRatioDF))
        return err;
      dCvg = coverage(lShortStart, lShortEnd, basis);
      dShortCashFwd = (1.0 / dRatioDF - 1.0) / dCvg;
    }
    // Get the ATM vols
    if (err = getCashVolwithType(szVolCurve, lToday, lFutureExpiry,
                                 lFutureStart, lFutureEnd, dLongCashFwd,
                                 dLongCashFwd, szRefRate, iVolatilityModel,
                                 getCashVol, &dATMLongVol))
      return err;

    if (err = getCashVolwithType(szVolCurve, lToday, lShortExpiry, lShortStart,
                                 lShortEnd, dShortCashFwd, dShortCashFwd,
                                 szRefRate, iVolatilityModel, getCashVol,
                                 &dATMShortVol))
      return err;

    // Calculate the std deviations that the option strike is from the long cash
    // fwd        , and hence the strike
    dLongExpiry = (lFutureExpiry - lToday) / 365.0;
    dShortExpiry = (lShortExpiry - lToday) / 365.0;
    if (iVolatilityModel == 0) {
      dStdDev = (dCashStrike - dLongCashFwd) / dATMLongVol / sqrt(dLongExpiry);
      dShortCashStrike =
          dShortCashFwd + dStdDev * dATMShortVol * sqrt(dShortExpiry);
    } else {
      dStdDev = (log(dCashStrike / dLongCashFwd) +
                 0.5 * dATMLongVol * dATMLongVol * dLongExpiry) /
                dATMLongVol / sqrt(dLongExpiry);
      dShortCashStrike = dShortCashFwd *
                         exp(-0.5 * dATMShortVol * dATMShortVol * dShortExpiry +
                             dStdDev * dATMShortVol * sqrt(dShortExpiry));
    }
    // Calculate the vols
    if (err = getCashVolwithType(szVolCurve, lToday, lShortExpiry, lShortStart,
                                 lShortEnd, dShortCashFwd, dShortCashStrike,
                                 szRefRate, iVolatilityModel, getCashVol,
                                 &dShortCashVol))
      return err;
  }

  // Calculate the price
  return MidCurveOption(lToday, lOptionExpiry, dCashStrike, lFutureExpiry,
                        dFuturePrice, 0, iVolatilityModel, enumCallPut,
                        enumGreek, dLongCashFwd, dLongSpread, dLongCashVol,
                        dShortCashVol, out_dPrice, out_dVol);
}
// A MidCurve option wrapper
// ----------------------------------------------------------------------------------
// //

// ----------------------------------------------------------------------------------
// // A MidCurve pricer
//
// return value:  MidCurve Option Price
//
//	mkt:				for getting fwds and vols
//	lOptionStart		delivery date of the option.  exercise date
// calculated using lag 	dStrike				strike of the
// option 	lFutureStart		start date for the Libor underlying the midcurve
// option 	dFuturePrice		market price of the future
// iFutureTenorInMonths length of underlying Libor in months (generally 1 or 3)
// iConvexityModel 0:  future = 100 * (1 - forward)        , strike is yield
// strike 						1:	midcurve strike
// adjusted by forward-future difference
//
//	iVolatilityModel	0:  stationary normal volatility
//	enumCallPut			Call/Put on price ( i.e.        , opposite
//on
// yield ) 	enumGreek			PREMIUM        , DELTA        , etc
// out_dVol vol used to price.  Will depend on the volatility model
//
char *MidCurveOption(long lToday, long lOptionExpiry, double dStrike,
                     long lFutureExpiry, double dFuturePrice,
                     int iConvexityModel, int iVolatilityModel,
                     SrtCallPutType enumCallPut, // Call = long
                     SrtGreekType enumGreek, double dLongCashFwd,
                     double dLongSpread, double dLongCashVol,
                     double dShortCashVol, double *out_dPrice,
                     double *out_dVol) {
  // Variable Declarations
  double dUndStrike, dOptionExpiry;
  //	double dIntrinsic;
  double dLongVar, dLongExpiry;
  double dShortExpiry;
  double dShortVar;

  // Initialize the return quantities
  *out_dVol = -1.0;
  *out_dPrice = -1.0;

  // If the option date is negative        , return 0 value
  dOptionExpiry = (lOptionExpiry - lToday) / 365.0;
  if (dOptionExpiry < 0.0) {
    *out_dPrice = 0.0;
    *out_dVol = 0.0;
    return 0;
  }

  // Calculate the modified strike depending on the convexity model
  if (iConvexityModel == 0) {
    dUndStrike = dStrike;
    //		dIntrinsic =  dStrike - dLongCashFwd - dLongSpread;
  } else if (iConvexityModel == 1) {
    dUndStrike = dLongCashFwd + dLongSpread + (dFuturePrice - dStrike) / 100.0;
    //		dIntrinsic = ( dFuturePrice - dStrike ) / 100.0;
  }

  // Reverse the call put so that it is on the yield and calculate the intrinsic
  // value
  if (enumCallPut == SRT_CALL) {
    enumCallPut = SRT_PUT;
  } else {
    enumCallPut = SRT_CALL;
    //		dIntrinsic = -dIntrinsic;
  }
  //	dIntrinsic = dIntrinsic < 0.0 ? 0.0 : dIntrinsic;

  // If the option expiry is today        , return the intrinsic value
  /*	if ( dOptionExpiry == 0.0 )
          {
                  *out_dPrice = dIntrinsic;
                  *out_dVol = 0.0;
                  return 0;
          }
  */

  // Get the normal variance for the long option
  dLongExpiry = (lFutureExpiry - lToday) / 365.0;
  dLongVar = dLongCashVol * dLongCashVol * dLongExpiry;

  // Calculate the normal variance for the short option
  dShortExpiry = (lFutureExpiry - lOptionExpiry) / 365.0;
  dShortVar = dShortCashVol * dShortCashVol * dShortExpiry;

  // Calculate the normal volatility
  if (dLongVar < dShortVar)
    return "Variance decreasing over time:  no midcurve price for this model";
  *out_dVol = sqrt((dLongVar - dShortVar) / dOptionExpiry);

  // Adjust the strike by the spread
  dUndStrike -= dLongSpread;

  // Calculate the price
  if (iVolatilityModel == 0)
    *out_dPrice = srt_f_optblknrm(dLongCashFwd, dUndStrike, *out_dVol,
                                  dOptionExpiry, 1.0, enumCallPut, enumGreek);
  else if (iVolatilityModel == 1)
    *out_dPrice = srt_f_optblknrm(dLongCashFwd, dUndStrike, *out_dVol,
                                  dOptionExpiry, 1.0, enumCallPut, enumGreek);
  else
    return "Unknown volatility model";

  // Check that the the price is greater than the intrinsic
  switch (enumGreek) {
    //		case SRT_PREMIUM:
    //			*out_dPrice = *out_dPrice < dIntrinsic ? dIntrinsic :
    //*out_dPrice; 			break;
  case DELTA:
    *out_dPrice = -(*out_dPrice);
    break;
  default:
    break;
  }

  // Return the price
  return 0;
}
// A MidCurve option pricer
// ----------------------------------------------------------------------------------
// //

// ----------------------------------------------------------------------------------
// // Auxiliary functions for the correlation bootstrap calculator
static char *StationaryAutocorrel(int iLow, int iHigh, int iTime, int nDates,
                                  double *dvVol, double *dvDeltaT,
                                  double *out_dAutoCorrel) {
  /* Variable declaration */
  int i;
  double dSum;

  /* Check that the indices make sense */
  if (iLow > iHigh || iTime > iLow || iHigh > nDates || iLow < 0 || iTime < 0)
    return "Incorrect indices in StationaryAutocorrel";

  /* Loop over the relevant dates f*/
  dSum = 0.0;
  for (i = iLow - iTime; i < iLow; i++)
    dSum += dvVol[i] * dvVol[i + iHigh - iLow] * dvDeltaT[i];

  /* store the value and return */
  *out_dAutoCorrel = dSum;
  return 0;
}

static char *StationaryMidCurve(int iCap, int iTime, int nDates, double *dvVol,
                                double *dvDeltaT, double *out_dMidCurve) {
  /* Variable declaration */
  int i;
  double dSum;

  /* Check that the indices make sense */
  if (iTime > iCap || iCap > nDates || iCap < 0 || iTime < 0)
    return "Incorrect indices in StationaryMidCurve";

  /* Loop over the relevant dates f*/
  dSum = 0.0;
  for (i = iCap - iTime; i < iCap; i++)
    dSum += dvVol[i] * dvVol[i] * dvDeltaT[i];

  /* store the value and return */
  *out_dMidCurve = dSum;
  return 0;
}

char *sumInner(int iLow, int iHigh, int iTime, double **dmLvl, int nDates,
               double *dvVol, double *dvDeltaT, double **dmCorrel,
               double *out_dSum) {
  /* Variable declaration */
  int i;
  double dTemp = 0.0, dAC;
  char *err = 0;

  /* Check that the indices make sense */
  if (iLow > iHigh || iTime > iLow || iHigh > nDates || iLow < 0 || iTime < 0)
    return "Incorrect indices in StationaryAutocorrel";

  /* Loop over the relevant indices */
  for (i = iLow + 1; i <= iHigh; i++) {
    if (err =
            StationaryAutocorrel(iLow, i, iTime, nDates, dvVol, dvDeltaT, &dAC))
      return err;
    dTemp += dmLvl[i][i] * dmCorrel[iLow][i] * dAC;
  }

  /* store the value and return */
  *out_dSum = 2.0 * dTemp;
  return 0;
}

char *sumOuter(int iLow, int iHigh, int iTime, double **dmLvl, int nDates,
               double *dvVol, double *dvDeltaT, double **dmCorrel,
               double *out_dSum) {
  /* Variable declaration */
  int i;
  double dTemp = 0.0, dMC, dSumInner;
  char *err = 0;

  /* Check that the indices make sense */
  if (iLow > iHigh || iTime > iLow || iHigh > nDates || iLow < 0 || iTime < 0)
    return "Incorrect indices in StationaryAutocorrel";

  /* Loop over the relevant indices */
  for (i = iLow + 1; i <= iHigh; i++) {
    if (err = StationaryMidCurve(i, iTime, nDates, dvVol, dvDeltaT, &dMC))
      return err;

    if (err = sumInner(i, iHigh, iTime, dmLvl, nDates, dvVol, dvDeltaT,
                       dmCorrel, &dSumInner))
      return err;

    dTemp += dmLvl[i][i] * (dmLvl[i][i] * dMC + dSumInner);
  }

  /* store the value and return */
  *out_dSum = dTemp;
  return 0;
}

// ----------------------------------------------------------------------------------
// // A correlation bootstrap calculator
char *CorrelFromStationaryCaps(
    long lToday,
    /*								char* szVolCurve
             , char* getCashVol_NYK( char *szVolCurve        , double dStartDate
       , double dEndDate        , double dCashStrike        , int iZero        ,
       char *szRefRate        , double *dVol        , double *dPower )        ,
    */
    int nDates, long *lStartDates, long *lEndDates, long *lvExpiry,
    double **dmVol, double **dmLvl,
    /*							    long lStart        ,
                                                                    long lEnd ,
                                                                    char*
       szRefRate        ,
    */
    double ***out_dmCorrel) {
  /* Variable declaration */
  int i, iDiag, iCol, iRow;
  char *err;
  double dPrevVar, dPrevTime, dVar, dSum, dAC;
  double *dvTime, *dvDeltaT, *dvSigma;

  /* Allocate the memory */
  dvTime = malloc(nDates * sizeof(double));
  dvDeltaT = malloc(nDates * sizeof(double));
  dvSigma = malloc(nDates * sizeof(double));

  /* Generate the schedule */

  /* Generate the Vol        , Fwd and Level matrices matrix */
  /*	for ( i=0; i<nDates; i++ )
                  for ( j=i; j<nDates; j++ )
                  {
                          if( err = getCashVol_NYK( szVolCurve        , (double)
     lStartDates[i]        , (double) lEndDates[i]        , dvFwds[i]        , 0
     , szRefRate        , &dVol        , &dPower )  ) return err; dmVar[i][j] =
     dVol * dVol * dvTime[i];
                  }
  */

  /* Bootstrap the stationary vol function */
  dPrevVar = 0.0;
  dPrevTime = 0.0;
  for (i = 0; i < nDates; i++) {
    dvTime[i] = (lvExpiry[i] - lToday) / 365.0;
    dvDeltaT[i] = dvTime[i] - dPrevTime;
    /* Check that the variance is not decreasing */
    dVar = dmVol[i][i] * dmVol[i][i] * dvTime[i];
    if (dVar < dPrevVar)
      goto FREE_RETURN;
    /* calculate the piecewise constant vol */
    dvSigma[i] = sqrt((dVar - dPrevVar) / dvDeltaT[i]);
    /* set the values for next loop */
    dPrevVar = dVar;
    dPrevTime = dvTime[i];
  }

  /* Bootstrap the correlation matrix */
  for (iDiag = 0; iDiag < nDates; iDiag++) {
    for (iCol = iDiag; iCol < nDates - iDiag; iCol++) {
      iRow = iDiag + iCol;
      if (err = sumOuter(iRow, iCol, iRow, dmLvl, nDates, dvSigma, dvDeltaT,
                         *out_dmCorrel, &dSum))
        return err;

      if (err = StationaryAutocorrel(iRow + 1, iCol, iRow, nDates, dvSigma,
                                     dvDeltaT, &dAC))
        return err;

      (*out_dmCorrel)[iRow + 1][iCol] = 0.5 *
                                        (dmVol[iRow][iCol] * dmVol[iRow][iCol] *
                                             dvTime[iRow] * dmLvl[iRow][iCol] -
                                         dSum) /
                                        dmLvl[iRow + 1][iRow + 1] *
                                        dmLvl[iCol][iCol] / dAC;
    }
  }

FREE_RETURN:
  free(dvTime);
  free(dvDeltaT);
  free(dvSigma);

  return 0;
}
