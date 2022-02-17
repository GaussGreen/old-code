/*****************************************************************************************/
/*                                                                                       */
/* FUNCTION     	: */
/*														                                 */
/*                                                                                       */
/* PURPOSE      	:
 */
/*                                                                                       */
/* DESCRIPTION  	:
 */
/*																						 */
/*															                             */
/*                                                                                       */
/*		                                                                                 */
/* PARAMETERS */
/*	INPUT	    :
 */
/*              :
 */
/*              :
 */
/*              :
 */
/*				:
 */
/*												                                         */
/*												                                         */
/* RETURNS      : */
/*                                                                                       */
/*****************************************************************************************/

/* -------------------------------------------------------------------------- */
#include "srt_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_tme.h"
#include <opfnctns.h"
#include <swp_h_all.h"
#include <swp_h_cmsopt.h"
/* -------------------------------------------------------------------------- */

Err swp_f_tme(long TMEDate, double dTECSpread, char *szSwapYieldCurveName,
              char *szBondYieldCurveName, char *szBondVolCurveName,
              char *cSwapRefRateCode, char *cBondRefRateCode,
              int iNumStrikesInVol, double *dStrikes, int bAdjForSpread,
              Err (*GetVol)(double dStart, double dEnd, double dStrike,
                            SRT_Boolean bAdjForSpread, double dForward,
                            double dSpread, double *pdBsVol),
              double *dAnswer) {
  Err err = NULL;

  long lToday;
  long lFraSpot;
  long lEndDate;
  double dMaturity;
  double dNumperiods;
  long lSwapTheoricalEndDate;
  long lBond10YTheoricalEndDate;
  long lBond12YTheoricalEndDate;
  int iSpotLag;

  SrtBasisCode Basis;
  SrtCompounding Compounding;

  SwapDP sdpSwap;
  SwapDP sdpBond10Y;
  SwapDP sdpBond12Y;

  double dFwdSwap = 0.;
  double dFwdTEC;
  double dTECCms;
  double dBond10Y = 0.;
  double dBond12Y = 0.;
  double dBondSpread;
  double dSpread;

  SrtCrvPtr SwapYldCrv;
  SrtCrvPtr BondYldCrv;
  SrtDiffusionType CMSlognorm = SRT_LOGNORMAL;

  SwapYldCrv = lookup_curve(szSwapYieldCurveName);
  BondYldCrv = lookup_curve(szBondYieldCurveName);
  lToday = get_today_from_curve(SwapYldCrv);
  iSpotLag = get_spotlag_from_curve(BondYldCrv);

  /* Initialise the inputs of the swapDP */
  lFraSpot = add_unit(TMEDate, iSpotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
  lSwapTheoricalEndDate = add_unit(TMEDate, 10, SRT_YEAR, MODIFIED_SUCCEEDING);
  lBond10YTheoricalEndDate = lSwapTheoricalEndDate;
  lBond12YTheoricalEndDate =
      add_unit(TMEDate, 12, SRT_YEAR, MODIFIED_SUCCEEDING);
  lEndDate = add_unit(lFraSpot, 10, SRT_YEAR, MODIFIED_SUCCEEDING);
  Compounding = SRT_ANNUAL;
  Basis = BASIS_ACT_360;

  /* Initialise the swap and the two bonds*/
  err = swp_f_setSwapDP(lFraSpot, lSwapTheoricalEndDate, Compounding, Basis,
                        &sdpSwap);
  if (err)
    return err;
  err = swp_f_setSwapDP(lFraSpot, lBond10YTheoricalEndDate, Compounding, Basis,
                        &sdpBond10Y);
  if (err)
    return err;
  err = swp_f_setSwapDP(lFraSpot, lBond12YTheoricalEndDate, Compounding, Basis,
                        &sdpBond12Y);
  if (err)
    return err;

  /* Input the spot lag  */
  sdpSwap.spot_lag = iSpotLag;
  sdpBond10Y.spot_lag = iSpotLag;
  sdpBond12Y.spot_lag = iSpotLag;

  /* computation of the Forwards */
  err = swp_f_ForwardRate_SwapDP(&sdpSwap, szSwapYieldCurveName,
                                 cSwapRefRateCode, &dFwdSwap);
  dFwdTEC = dFwdSwap + dTECSpread;

  err = swp_f_ForwardRate_SwapDP(&sdpBond10Y, szBondYieldCurveName,
                                 cBondRefRateCode, &dBond10Y);

  err = swp_f_ForwardRate_SwapDP(&sdpBond12Y, szBondYieldCurveName,
                                 cBondRefRateCode, &dBond12Y);

  dBondSpread = dBond12Y - dBond10Y;

  /* maturity of the TME */
  dMaturity = coverage(lToday, TMEDate, BASIS_ACT_365);
  dNumperiods = (lBond10YTheoricalEndDate - lFraSpot) * Compounding / 365.0;
  dSpread = swp_f_spread(lFraSpot, lEndDate, cBondRefRateCode);

  if (err = swp_f_Cms_Rate(
          dFwdTEC, dMaturity, dNumperiods, Compounding, 0.0, 1.0, CMSlognorm,
          0.0, 1 /* Full Smile Approx */, lFraSpot, lBond10YTheoricalEndDate,
          (SRT_Boolean)bAdjForSpread, bAdjForSpread * dSpread,
          szBondVolCurveName, iNumStrikesInVol, dStrikes, &dTECCms))

  {
    smessage("Error in swp_f_cmsrate");
  }

  *dAnswer = dTECCms + dBondSpread;

  return err;
}