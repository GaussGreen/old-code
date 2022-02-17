/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CMT_CAPFLOOR                                    */
/*                                                                            */
/*      PURPOSE:  CAPS        , FLOORS        , IMPLIED_VOL CAP/FLOOR        ,
 * CAPLETS        , FLOORLETS     */
/*					on CMT MARKETS
 */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "num_h_allhdr.h"
#include "opfnctns.h"
#include "opsabrcalib.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_cmt_capfloor.h"
#include "swp_h_cmt_swap.h"
#include "swp_h_vol.h"

#define MAXITERATION 10
#define PREM_TOL_PCT2 0.01
#define PREM_TOL2 (double)1e-08

/* ------------------------------------------------------------------------ */
/* Function that fills in a GenFullSwap for a CMT cap/floor */

static Err set_CMT_capfloor(long start, long end_or_nfp, String cap_freq_str,
                            String cap_basis_str, double strike,
                            SrtCurvePtr crvptr, GenFullSwap *capfloor_swap) {
  Err err;
  GenFullSwap cmt_swap;
  SwapDP cmt_swapdp;
  String cmt_name, yc_name;
  SrtCrvPtr yc_crv;

  /* Gets yield curve from cmt_und */
  cmt_name = get_curve_name(crvptr);
  yc_name = get_ycname_from_cmtcrv(crvptr);
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find the %s yield curve", yc_name);

  /* Set up the underlying CMT swapdp */
  err = swp_f_initSwapDP(start, end_or_nfp, cap_freq_str, cap_basis_str,
                         &cmt_swapdp);
  if (err)
    return err;

  cmt_swapdp.spot_lag = get_spotlag_from_curve(yc_crv);

  /* Generates a CMT vs Fixed swap to generate dates        , cvg */
  err = swp_f_create_GenFullSwap(&cmt_swapdp, &cmt_swapdp, strike, 0.0, 0.0,
                                 CMT_MARGIN_FIXED_SWAP, cmt_name, &cmt_swap);

  /* Populates a CMT vs Fixed swap to generate strike        , vols... */
  err = swp_f_populate_GenFullSwap(0.0, crvptr, &cmt_swap);
  if (err)
    return err;

  /* Values a CMT vs Fixed swap to generate dfs and get intrinsic */
  err = swp_f_value_GenFullSwap(&cmt_swap, crvptr);
  if (err)
    return err;

  /* Sets the capfloor_swap to the newly created swap*/
  *capfloor_swap = cmt_swap;

  /* Returns a success message */
  return NULL;
}

/* ------------------------------------------------------------------------- */
/* Function that prices a full CMT capfloor from a GenFullSwap */

static Err price_CMT_capfloor(GenFullSwap *capfloor_swap,
                              SrtReceiverType rec_pay, SrtGreekType greek,
                              SrtCurvePtr crvptr, double *result) {
  Err err;
  int i;
  int list_length;
  long start_price_date;
  long fixing_date;
  DateList date_list;
  String yc_name, cmt_name;
  SrtCurvePtr yc_crv;
  SrtCallPutType call_put;
  CMT_Param *cmt_param;
  SRT_Boolean use_prop_vol;
  long spot, today;
  double caplet_grk;
  double cap_grk;
  double treas_rate, swap_rate;
  double fwd, strike, vol, mat, cvg, df;
  Date period_start, period_end;

  /* Checks the curve type */
  if (!ISCURVETYPE(crvptr, CMT_CURVE))
    return serror("Wrong underlying type: need a CMT_CURVE");

  /* Get cmt_crv and yc_crv names  */
  cmt_name = get_curve_name(crvptr);
  yc_name = get_ycname_from_cmtcrv(crvptr);

  /* Gets cmt parameters */
  cmt_param = get_cmtparam_from_cmtcrv(crvptr);
  use_prop_vol = cmt_param->use_prop_vol_flg;

  /* Get spot and clcn date*/
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find the %s yield curve", yc_name);
  today = get_today_from_curve(yc_crv);
  spot = get_spotdate_from_yldcrv(yc_crv);
  crvptr = lookup_curve(cmt_name);

  /* OVE: LET NOT USE START_PRICE_DATE */
  /* If cap has already started        , first fixing used in the price is after
   * spot
   */
  if (spot > capfloor_swap->leg[0].pay_date.date[0]) {
    start_price_date = spot;
  } else
  /* If cap starts today or after        , or is BROKEN        ,
          first period used is the first one of the swap STRICTLY after spot */
  {
    start_price_date = capfloor_swap->leg[0].pay_date.date[0];
  }

  /* Compute the value of each caplet independently (do not use last fixing)*/
  list_length = capfloor_swap->leg[0].leg_length;
  date_list = capfloor_swap->leg[0].pay_date;
  cap_grk = 0.0;
  for (i = 0; i < list_length - 1; i++) {
    /* We reject all the dates that are before OR EQUAL to spot */
    /* as if cap starts today (spot == date[0])         , the first fixing used
    /* is the next one in the swap*/
    if (date_list.date[i] > spot) {
      fwd = capfloor_swap->leg[0].cxxpn.d[i + 1]; /* convexity adjusted fwd */
      strike = capfloor_swap->leg[1].cxxpn.d[i];  /* strike (2nd leg) */

      if (cmt_param->flatvol != 0.) {
        vol = cmt_param->flatvol;
      } else {
        period_start = capfloor_swap->leg[0].start_date.date[i];
        period_end = add_unit((long)period_start, cmt_param->cmt_mat, SRT_YEAR,
                              SUCCEEDING);
        err = cmt_param->cmt_getvol(period_start, period_end, strike, &vol);
        if (err)
          return err;
      }

      treas_rate = capfloor_swap->leg[0].fwd.d[i];
      swap_rate = capfloor_swap->leg[0].temp_fwd.d[i];
      if (use_prop_vol == SRT_YES) {
        vol *= swap_rate / treas_rate;
      }

      /* OVE fix: a swap is described by theoretical payment dates: spot lag
      after fixing        , vol is effective only until fixing date */
      fixing_date = add_unit(capfloor_swap->leg[0].pay_date.date[i],
                             -cmt_param->spot_lag, SRT_BDAY, SUCCEEDING);
      mat = (double)(fixing_date - today) / 365.0;
      cvg = capfloor_swap->leg[0].cxxvg.d[i];
      /* Payment of the caplet is done on the next fixng date */
      df = capfloor_swap->leg[0].df.d[i + 1];
      if (rec_pay == SRT_RECEIVER)   /* Floor */
        call_put = SRT_PUT;          /* Put */
      else if (rec_pay == SRT_PAYER) /* Cap */
        call_put = SRT_CALL;         /* Call */

      caplet_grk = srt_f_optblksch(fwd, strike, vol, mat, df, call_put, greek);
      caplet_grk *= cvg;
      cap_grk += caplet_grk;

    } /* END if (date_list.date[i]>start_price_date) */
  }   /* END for i loop on fixing dates */

  /* Sets the result */
  *result = cap_grk;

  /* Returns a success message*/
  return NULL;
}

/* ------------------------------------------------------------------------ */
Err swp_f_CMT_capfloor(long start, long end_or_nfp, String cap_freq_str,
                       String cap_basis_str, double strike, double flatvol,
                       String vol_type_str, String rec_pay_str,
                       String greek_str, SrtCurvePtr crvptr, double *value) {
  Err err;
  CMT_Param *cmt_param;
  SrtDiffusionType vol_type = SRT_LOGNORMAL;
  SrtReceiverType rec_pay_type;
  SrtGreekType greek_type = PREMIUM;
  GenFullSwap capfloor_swap;
  double result;

  err = interp_rec_pay(rec_pay_str, &rec_pay_type);
  if (err)
    return err;

  if (greek_str) {
    err = interp_greeks(greek_str, &greek_type);
    if (err)
      return err;
  }

  /* Fill the flatvolatility        , we need the cmt parameters */
  cmt_param = get_cmtparam_from_cmtcrv(crvptr);
  if ((flatvol != 0.) && (vol_type_str)) {
    err = interp_diffusion_type(vol_type_str, &vol_type);
    if (err)
      return err;
  }

  cmt_param->voltype = vol_type;
  cmt_param->flatvol = flatvol;

  /* Creates the CMT_swap needed to generate cap/floor cashflows */
  err = set_CMT_capfloor(start, end_or_nfp, cap_freq_str, cap_basis_str, strike,
                         crvptr, &capfloor_swap);
  if (err)
    return err;

  /* Prices the CMT capfloor from the CMT swap */
  err = price_CMT_capfloor(&capfloor_swap, rec_pay_type, greek_type, crvptr,
                           &result);
  if (err)
    return err;

  /* Sets the result */
  *value = result;

  /* Free the content of the swap */
  swp_f_freein_GenFullSwap(&capfloor_swap);

  /* Returns a success message */
  return NULL;
}

/* ========================================================================= */
Err swp_f_CMT_capfloor_impvol(long start, long end_or_nfp, String cap_freq_str,
                              String cap_basis_str, double strike,
                              double premium, String rec_pay_str,
                              SrtCurvePtr crvptr, double *imp_vol) {
  Err err;
  SrtReceiverType rec_pay_type;
  CMT_Param *cmt_param;
  SrtDiffusionType vol_type = SRT_LOGNORMAL;
  GenFullSwap capfloor_swap;
  double result;
  double vol_guess;
  double shift;
  double a[3], b[3];
  int niter, count;
  double nstop, target;

  err = interp_rec_pay(rec_pay_str, &rec_pay_type);
  if (err)
    return err;

  /* Sets the initial vol to 0.0 : chack the intrinsic */
  cmt_param = get_cmtparam_from_cmtcrv(crvptr);
  cmt_param->voltype = vol_type;
  cmt_param->flatvol = NULL_VOL;

  /* Creates the CMT_swap needed to generate cap/floor cashflows */
  err = set_CMT_capfloor(start, end_or_nfp, cap_freq_str, cap_basis_str, strike,
                         crvptr, &capfloor_swap);
  if (err)
    return err;

  /* Prices the CMT capfloor from the CMT swap */
  err = price_CMT_capfloor(&capfloor_swap, rec_pay_type, PREMIUM, crvptr,
                           &result);
  if (err)
    return err;

  /* Free the content of the swap */
  swp_f_freein_GenFullSwap(&capfloor_swap);

  /* Checks premium >= price with vol 0.0 */
  if (premium < result)
    return serror("The premium input is smaller than intrinsic");

  /* If premium == price: return NULL vol */
  if (fabs(result / premium - 1) < DBL_EPSILON)
    return NULL;

  /* Sets Newton useful parameters */
  niter = 5;
  target = premium;
  count = 0;
  shift = 0.0001;
  nstop = 0.0;

  /* Sets Newton first values */
  vol_guess = 0.1;
  cmt_param->flatvol = vol_guess;

  err = set_CMT_capfloor(start, end_or_nfp, cap_freq_str, cap_basis_str, strike,
                         crvptr, &capfloor_swap);
  if (err)
    return err;

  err = price_CMT_capfloor(&capfloor_swap, rec_pay_type, PREMIUM, crvptr,
                           &result);
  if (err)
    return err;

  /* Free the content of the swap */
  swp_f_freein_GenFullSwap(&capfloor_swap);

  b[0] = result;

  vol_guess += shift;

  cmt_param->flatvol = vol_guess;

  /* Creates the CMT_swap needed to generate cap/floor cashflows */
  err = set_CMT_capfloor(start, end_or_nfp, cap_freq_str, cap_basis_str, strike,
                         crvptr, &capfloor_swap);
  if (err)
    return err;

  err = price_CMT_capfloor(&capfloor_swap, rec_pay_type, PREMIUM, crvptr,
                           &result);
  if (err)
    return err;

  /* Free the content of the swap */
  swp_f_freein_GenFullSwap(&capfloor_swap);

  b[1] = result;

  vol_guess += shift;

  a[2] = vol_guess;
  /* Starts Newton iterations */
  while ((count < MAXITERATION) && (nstop < 1.0)) {
    cmt_param->flatvol = vol_guess;
    /* Creates the CMT_swap needed to generate cap/floor cashflows */
    err = set_CMT_capfloor(start, end_or_nfp, cap_freq_str, cap_basis_str,
                           strike, crvptr, &capfloor_swap);
    if (err)
      return err;

    err = price_CMT_capfloor(&capfloor_swap, rec_pay_type, PREMIUM, crvptr,
                             &result);
    if (err)
      return err;
    /* Free the content of the swap */
    swp_f_freein_GenFullSwap(&capfloor_swap);

    b[2] = result;

    newton(target, niter, a, b, &nstop);
  }

  /* Sets the result */
  *imp_vol = a[2];

  /* Returns a success message */
  return NULL;
}

/* ========================================================================= */

Err swp_f_CMCapFloor(
    long TheoEnd, long numperiodcap, double Underlying, char *cCMTFrequency,
    char *cCMTBasis, char *cCMTRefRateCode, double Strike, char *cCapFrequency,
    char *cCapBasis, char *szCMTVolName, char *szSwaptionVolName,
    char *szYieldCurveName, char *cRefRateCode, char *cCapFloor,
    char *method, /*BS: BS on CMS Full Smile + Fwd Spread - FS: Full Smile on
               Fwd Swap + Fwd Spread - FVBS: BS on CMS Flat Vol + Fwd
               Spread - FV: Flat Vol on Fwd Swap + Fwd Spread*/
    SRT_Boolean AdjForSpread, int NumStrikesInVol, double *Strikes,
    double **result) {
  Err err = NULL;
  int i;
  double FwdSwapRate;
  double *CMS;
  double *CMT;
  double Cap = 0;
  double *Caplet;
  long *FixingDates;
  long *PayDates;
  long *StartDates;
  long numperiods;
  long TheoEndCMT;
  double dspread;
  double dFlatVol;
  double *CMTvol;
  double lognorm;
  long today, spotdate, spotlag;
  SrtCurvePtr Crv;
  SwapDP sdpSwap;
  SrtCompounding CapFreq, CMTFreq;
  SrtBasisCode CMTBasis;
  SrtBasisCode CapBasis;
  SrtCallPutType CapFloor;

  err = interp_call_put(cCapFloor, &CapFloor);
  err = interp_compounding(cCapFrequency, &CapFreq);
  err = interp_compounding(cCMTFrequency, &CMTFreq);
  err = interp_basis(cCMTBasis, &CMTBasis);
  err = interp_basis(cCapBasis, &CapBasis);

  Crv = lookup_curve(szYieldCurveName);
  today = get_clcndate_from_yldcrv(Crv);
  spotdate = get_spotdate_from_yldcrv(Crv);

  spotlag = get_spotlag_from_curve(Crv);

  numperiods = (long)Underlying * CMTFreq;

  /*Memory Allocation*/
  CMS = dvector(0, numperiodcap - 1);
  CMT = dvector(0, numperiodcap - 1);
  CMTvol = dvector(0, numperiodcap - 1);
  FixingDates = lngvector(0, numperiodcap - 1);
  StartDates = lngvector(0, numperiodcap - 1);
  PayDates = lngvector(0, numperiodcap - 1);
  Caplet = dvector(0, numperiodcap - 1);

  /*create Cap/Floor schedule*/
  for (i = 0; i < numperiodcap; i++) {
    PayDates[i] =
        add_unit(add_unit(TheoEnd, -(numperiodcap - i - 1) * 12 / CapFreq,
                          SRT_MONTH, NO_BUSDAY_CONVENTION),
                 0, SRT_BDAY, MODIFIED_SUCCEEDING);
    StartDates[i] =
        add_unit(add_unit(TheoEnd, -(numperiodcap - i) * 12 / CapFreq,
                          SRT_MONTH, NO_BUSDAY_CONVENTION),
                 0, SRT_BDAY, MODIFIED_SUCCEEDING);
    FixingDates[i] =
        add_unit(StartDates[i], -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);
  }

  for (i = 0; i < numperiodcap; i++) {

    /* first initialise the sdpFra */
    TheoEndCMT = add_unit(StartDates[i], numperiods * 12 / CMTFreq, SRT_MONTH,
                          NO_BUSDAY_CONVENTION);
    err =
        swp_f_setSwapDP(StartDates[i], TheoEndCMT, CMTFreq, CMTBasis, &sdpSwap);
    if (err)
      return err;

    /* input the spot lag  */
    sdpSwap.spot_lag = spotlag;

    /* computation of the Fra */
    err = swp_f_ForwardRate_SwapDP(&sdpSwap, szYieldCurveName, cRefRateCode,
                                   &FwdSwapRate);

    dspread = swp_f_spread(FixingDates[i], StartDates[i], cCMTRefRateCode);

    if (strcmp(method, "BS") == 0) {
      if (err = swp_f_Cms_Rate(
              FwdSwapRate, (FixingDates[i] - today) / 365.0, numperiods,
              (double)CMTFreq, (PayDates[i] - StartDates[i]) / 365.0, 1.0,
              SRT_LOGNORMAL, 0.0, 2 /* Full Smile Approx */, StartDates[i],
              TheoEndCMT, 0, 0.0, szSwaptionVolName, NumStrikesInVol, Strikes,
              &CMS[i]))

      {
        return err;
      }

      CMT[i] = CMS[i] + dspread;
      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, Strike, &(CMTvol[i]),
                &lognorm);
      Caplet[i] = srt_f_optblksch(CMT[i], Strike, CMTvol[i],
                                  (FixingDates[i] - today) / 365.25, 1,
                                  CapFloor, PREMIUM);

    }

    else if (strcmp(method, "FVBS") == 0) {
      swp_f_vol(szSwaptionVolName, StartDates[i], TheoEndCMT, FwdSwapRate,
                &dFlatVol, &lognorm);
      if (err = swp_f_Cms_Rate(
              FwdSwapRate, (FixingDates[i] - today) / 365.0, numperiods,
              (double)CMTFreq, (PayDates[i] - StartDates[i]) / 365.0, 1.0,
              SRT_LOGNORMAL, dFlatVol, 0 /* Flat Vol */, StartDates[i],
              TheoEndCMT, 0, 0.0, szSwaptionVolName, NumStrikesInVol, Strikes,
              &CMS[i]))

      {
        return err;
      }

      CMT[i] = CMS[i] + dspread;
      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, Strike, &(CMTvol[i]),
                &lognorm);
      Caplet[i] = srt_f_optblksch(CMT[i], Strike, CMTvol[i],
                                  (FixingDates[i] - today) / 365.25, 1,
                                  CapFloor, PREMIUM);

    }

    else if (strcmp(method, "FV") == 0) {

      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, FwdSwapRate, &dFlatVol,
                &lognorm);
      if (err = swp_f_Cms_Rate(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.0,
              numperiods, (double)CMTFreq,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL,
              dFlatVol, 0 /* Flat Vol */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &CMT[i]))

      {
        return err;
      }

      if (err = swp_f_Cms_Option(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.25,
              numperiods, Strike, (double)CMTFreq, CapFloor,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL,
              dFlatVol, 0 /* Flat vol */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &Caplet[i]))

      {
        return err;
      }

    }

    else {

      if (err = swp_f_Cms_Rate(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.0,
              numperiods, (double)CMTFreq,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL, 0.0,
              2 /* Full Smile Approx */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &CMT[i]))

      {
        return err;
      }

      if (err = swp_f_Cms_Option(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.25,
              numperiods, Strike, (double)CMTFreq, CapFloor,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL, 0.0,
              2 /* Full Smile Approx */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &Caplet[i]))

      {
        return err;
      }
    }

    Cap += Caplet[i] * coverage(StartDates[i], PayDates[i], CapBasis) *
           swp_f_df(today, PayDates[i], szYieldCurveName);
  }

  result[0][0] = Cap;

  for (i = 0; i < numperiodcap; i++) {
    result[i + 1][0] = FixingDates[i];
    result[i + 1][1] = StartDates[i];
    result[i + 1][2] = PayDates[i];
    result[i + 1][3] = CMT[i];
    result[i + 1][4] = Caplet[i] *
                       coverage(StartDates[i], PayDates[i], CapBasis) *
                       swp_f_df(today, PayDates[i], szYieldCurveName);
  }

  /*Free Memory*/
  free_dvector(CMS, 0, numperiodcap - 1);
  free_dvector(CMT, 0, numperiodcap - 1);
  free_dvector(CMTvol, 0, numperiodcap - 1);
  free_lngvector(FixingDates, 0, numperiodcap - 1);
  free_lngvector(StartDates, 0, numperiodcap - 1);
  free_lngvector(PayDates, 0, numperiodcap - 1);
  free_dvector(Caplet, 0, numperiodcap - 1);
  return err;
}

Err swp_f_CMCapFloorNY(
    long TheoEnd, long numperiodcap, double Underlying, char *cCMTFrequency,
    char *cCMTBasis, char *cCMTRefRateCode, double Strike, char *cCapFrequency,
    char *cCapBasis, char *szCMTVolName, char *szSwaptionVolName,
    char *szYieldCurveName, char *cRefRateCode, char *cCapFloor,
    char *method, /*BS: BS on CMS Full Smile + Fwd Spread - FS: Full Smile on
               Fwd Swap + Fwd Spread - FVBS: BS on CMS Flat Vol + Fwd
               Spread - FV: Flat Vol on Fwd Swap + Fwd Spread*/
    SRT_Boolean AdjForSpread, int NumStrikesInVol, double *Strikes,
    double **result) {
  Err err = NULL;
  int i;
  double FwdSwapRate;
  double *CMS;
  double *CMT;
  double Cap = 0;
  double *Caplet;
  long *FixingDates;
  long *PayDates;
  long *StartDates;
  long numperiods;
  long TheoEndCMT;
  double dspread;
  double dFlatVol;
  double *CMTvol;
  double lognorm;
  long today, spotdate, spotlag;
  SrtCurvePtr Crv;
  SwapDP sdpSwap;
  SrtCompounding CapFreq, CMTFreq;
  SrtBasisCode CMTBasis;
  SrtBasisCode CapBasis;
  SrtCallPutType CapFloor;

  err = interp_call_put(cCapFloor, &CapFloor);
  err = interp_compounding(cCapFrequency, &CapFreq);
  err = interp_compounding(cCMTFrequency, &CMTFreq);
  err = interp_basis(cCMTBasis, &CMTBasis);
  err = interp_basis(cCapBasis, &CapBasis);

  Crv = lookup_curve(szYieldCurveName);
  today = get_clcndate_from_yldcrv(Crv);
  spotdate = get_spotdate_from_yldcrv(Crv);

  spotlag = get_spotlag_from_curve(Crv);

  numperiods = (long)Underlying * CMTFreq;

  /*Memory Allocation*/
  CMS = dvector(0, numperiodcap - 1);
  CMT = dvector(0, numperiodcap - 1);
  CMTvol = dvector(0, numperiodcap - 1);
  FixingDates = lngvector(0, numperiodcap - 1);
  StartDates = lngvector(0, numperiodcap - 1);
  PayDates = lngvector(0, numperiodcap - 1);
  Caplet = dvector(0, numperiodcap - 1);

  /*create Cap/Floor schedule*/
  for (i = 0; i < numperiodcap; i++) {
    PayDates[i] =
        add_unit(add_unit(TheoEnd, -(numperiodcap - i - 1) * 12 / CapFreq,
                          SRT_MONTH, NO_BUSDAY_CONVENTION),
                 0, SRT_BDAY, MODIFIED_SUCCEEDING);
    StartDates[i] =
        add_unit(add_unit(TheoEnd, -(numperiodcap - i) * 12 / CapFreq,
                          SRT_MONTH, NO_BUSDAY_CONVENTION),
                 0, SRT_BDAY, MODIFIED_SUCCEEDING);
    FixingDates[i] =
        add_unit(StartDates[i], -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);
  }

  for (i = 0; i < numperiodcap; i++) {

    /* first initialise the sdpFra */
    TheoEndCMT = add_unit(StartDates[i], numperiods * 12 / CMTFreq, SRT_MONTH,
                          NO_BUSDAY_CONVENTION);
    err =
        swp_f_setSwapDP(StartDates[i], TheoEndCMT, CMTFreq, CMTBasis, &sdpSwap);
    if (err)
      return err;

    /* input the spot lag  */
    sdpSwap.spot_lag = spotlag;

    /* computation of the Fra */
    err = swp_f_ForwardRate_SwapDP(&sdpSwap, szYieldCurveName, cRefRateCode,
                                   &FwdSwapRate);

    dspread = swp_f_spread(FixingDates[i], StartDates[i], cCMTRefRateCode);

    if (strcmp(method, "BS") == 0) {
      if (err = swp_f_Cms_Rate(
              FwdSwapRate, (FixingDates[i] - today) / 365.0, numperiods,
              (double)CMTFreq, (PayDates[i] - StartDates[i]) / 365.0, 1.0,
              SRT_LOGNORMAL, 0.0, 2 /* Full Smile Approx */, StartDates[i],
              TheoEndCMT, 0, 0.0, szSwaptionVolName, NumStrikesInVol, Strikes,
              &CMS[i]))

      {
        return err;
      }

      CMT[i] = CMS[i] + dspread;
      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, Strike, &(CMTvol[i]),
                &lognorm);
      Caplet[i] = srt_f_optblksch(CMT[i], Strike, CMTvol[i],
                                  (FixingDates[i] - today) / 365.25, 1,
                                  CapFloor, PREMIUM);

    }

    else if (strcmp(method, "FVBS") == 0) {
      swp_f_vol(szSwaptionVolName, StartDates[i], TheoEndCMT, FwdSwapRate,
                &dFlatVol, &lognorm);
      if (err = swp_f_Cms_Rate(
              FwdSwapRate, (FixingDates[i] - today) / 365.0, numperiods,
              (double)CMTFreq, (PayDates[i] - StartDates[i]) / 365.0, 1.0,
              SRT_LOGNORMAL, dFlatVol, 0 /* Flat Vol */, StartDates[i],
              TheoEndCMT, 0, 0.0, szSwaptionVolName, NumStrikesInVol, Strikes,
              &CMS[i]))

      {
        return err;
      }

      CMT[i] = CMS[i] + dspread;
      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, Strike, &(CMTvol[i]),
                &lognorm);
      Caplet[i] = srt_f_optblksch(CMT[i], Strike, CMTvol[i],
                                  (FixingDates[i] - today) / 365.25, 1,
                                  CapFloor, PREMIUM);

    }

    else if (strcmp(method, "FV") == 0) {

      swp_f_vol(szCMTVolName, StartDates[i], TheoEndCMT, FwdSwapRate, &dFlatVol,
                &lognorm);
      if (err = swp_f_Cms_Rate(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.0,
              numperiods, (double)CMTFreq,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL,
              dFlatVol, 0 /* Flat Vol */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &CMT[i]))

      {
        return err;
      }

      if (err = swp_f_Cms_Option(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.25,
              numperiods, Strike, (double)CMTFreq, CapFloor,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL,
              dFlatVol, 0 /* Flat vol */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &Caplet[i]))

      {
        return err;
      }

    }

    else {

      if (err = swp_f_Cms_Rate(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.0,
              numperiods, (double)CMTFreq,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL, 0.0,
              2 /* Full Smile Approx */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &CMT[i]))

      {
        return err;
      }

      if (err = swp_f_Cms_OptionNY(
              FwdSwapRate + dspread, (FixingDates[i] - today) / 365.25,
              numperiods, Strike, (double)CMTFreq, CapFloor,
              (PayDates[i] - StartDates[i]) / 365.0, 1.0, SRT_LOGNORMAL, 0.0,
              2 /* Full Smile Approx */, StartDates[i], TheoEndCMT, 0, 0.0,
              szCMTVolName, NumStrikesInVol, Strikes, &(Caplet[i]),
              szYieldCurveName, cRefRateCode))

      {
        return err;
      }
    }

    Cap += Caplet[i] * coverage(StartDates[i], PayDates[i], CapBasis) *
           swp_f_df(today, PayDates[i], szYieldCurveName);
  }

  result[0][0] = Cap;

  for (i = 0; i < numperiodcap; i++) {
    result[i + 1][0] = FixingDates[i];
    result[i + 1][1] = StartDates[i];
    result[i + 1][2] = PayDates[i];
    result[i + 1][3] = CMT[i];
    result[i + 1][4] = Caplet[i] *
                       coverage(StartDates[i], PayDates[i], CapBasis) *
                       swp_f_df(today, PayDates[i], szYieldCurveName);
  }

  /*Free Memory*/
  free_dvector(CMS, 0, numperiodcap - 1);
  free_dvector(CMT, 0, numperiodcap - 1);
  free_dvector(CMTvol, 0, numperiodcap - 1);
  free_lngvector(FixingDates, 0, numperiodcap - 1);
  free_lngvector(StartDates, 0, numperiodcap - 1);
  free_lngvector(PayDates, 0, numperiodcap - 1);
  free_dvector(Caplet, 0, numperiodcap - 1);
  return err;
}

Err swp_f_getcmsvoladjusted(double FwdSwapRate, long StartSwap,
                            long TheoEndSwap, double Mat, double CMS,
                            double FwdSpread, double Strike,
                            char *szSwaptionVolName, double alpha, double beta,
                            double rho, double *CMSVol) {
  Err err = NULL;
  double ATMVolSwap, ATMVolCMS;
  double lognorm;

  swp_f_vol(szSwaptionVolName, StartSwap, TheoEndSwap, FwdSwapRate, &ATMVolSwap,
            &lognorm);
  ATMVolCMS = ATMVolSwap * FwdSwapRate / CMS;

  err = srt_f_optsarbvol(CMS, Strike - FwdSpread, Mat, ATMVolCMS, alpha, beta,
                         rho, SRT_LOGNORMAL, SRT_LOGNORMAL, CMSVol);

  return err;
}

Err swp_f_getcmtvol(double CMS, double FwdSpread, double CMSVol, double Strike,
                    double dLambda, double dSpreadLocVol, double dCorrel,
                    double Mat, double *CMTVol) {
  Err err = NULL;
  double premium;
  double CMSNormVol, CMTNormVol;
  double FwdSpreadNormVol;

  premium = srt_f_optblksch(CMS, Strike - FwdSpread, CMSVol, Mat, 1.0, SRT_CALL,
                            PREMIUM);
  err = srt_f_optimpvol(premium, CMS, Strike - FwdSpread, Mat, 1.0, SRT_CALL,
                        SRT_NORMAL, &CMSNormVol);

  FwdSpreadNormVol =
      dSpreadLocVol * sqrt((1 - exp(-2 * dLambda * Mat)) / (2 * dLambda * Mat));

  CMTNormVol = sqrt(CMSNormVol * CMSNormVol -
                    2 * dCorrel * CMSNormVol * FwdSpreadNormVol +
                    FwdSpreadNormVol * FwdSpreadNormVol);

  premium = srt_f_optblknrm(CMS + FwdSpread, Strike, CMTNormVol, Mat, 1.0,
                            SRT_CALL, PREMIUM);
  err = srt_f_optimpvol(premium, CMS + FwdSpread, Strike, Mat, 1.0, SRT_CALL,
                        SRT_LOGNORMAL, CMTVol);
  return err;
}

Err swp_f_getcmtvolMAD(double StartSwap, double TheoEndSwap,
                       char *szSwaptionVolName, double FwdSwapRate, double CMS,
                       double FwdSpread, int NumStrikes, double *Strikes,
                       double Mat, double alpha, double beta, double rho,
                       double *CMTVol) {
  Err err = NULL;
  int i;
  double SwapStrikeVol;
  double lognorm;

  for (i = 0; i < NumStrikes; i++) {
    swp_f_vol(szSwaptionVolName, StartSwap, TheoEndSwap, Strikes[i] - FwdSpread,
              &SwapStrikeVol, &lognorm);
    CMTVol[i] = SwapStrikeVol * FwdSwapRate / (FwdSwapRate + FwdSpread);
    /*		err = srt_f_optsarbvol(	CMS+FwdSpread        ,Strikes[i]        ,Mat
     * ,ATMVolCMT ,alpha        ,beta        ,rho        ,SRT_LOGNORMAL
     * ,SRT_LOGNORMAL        ,&(CMTVol[i]));*/
  }
  return err;
}

Err swp_f_getcmtvolCHI(long StartSwap, long TheoEndSwap,
                       char *szSwaptionVolName, double FwdSwapRate, double CMS,
                       double FwdSpread, double CMTFreq, double numperiods,
                       int NumStrikes, double *Strikes, double Mat,
                       double alpha, double beta, double rho, double CHI,
                       double *CMTVol, double *SabrComp) {
  Err err = NULL;
  int i;
  double ATMVolCMT, ATMVolCMS, ATMNormalVolCMT;
  double error;
  double *CMSOption = NULL;
  double *implied_vol = NULL;
  double CMSATMOption;
  double ATMNormalimpvol, ATMPrice;
  double alphacmt, betacmt, rhocmt;
  double *CMTStrikes = NULL;
  double *CMTVols = NULL;
  /*memory allocation*/
  CMSOption = dvector(0, NumStrikes - 1);
  implied_vol = dvector(0, NumStrikes - 1);
  CMTStrikes = dvector(0, NumStrikes);
  CMTVols = dvector(0, NumStrikes);

  alphacmt = alpha;
  betacmt = beta;
  rhocmt = rho;

  for (i = 0; i < NumStrikes; i++) {
    err = swp_f_Cms_Option(FwdSwapRate, Mat, numperiods, Strikes[i] - FwdSpread,
                           (double)CMTFreq, SRT_PAYER, 0, 1.0, SRT_LOGNORMAL,
                           0.0, 2 /* Full Smile Approx */, StartSwap,
                           TheoEndSwap, 0, 0.0, szSwaptionVolName, NumStrikes,
                           Strikes, &(CMSOption[i]));

    err = srt_f_optimpvol(CMSOption[i], CMS + FwdSpread, Strikes[i], Mat, 1.0,
                          SRT_PAYER, SRT_LOGNORMAL, &(implied_vol[i]));
  }

  err = swp_f_Cms_Option(
      FwdSwapRate, Mat, numperiods, CMS, (double)CMTFreq, SRT_PAYER, 0.0, 1.0,
      SRT_LOGNORMAL, 0.0, 2 /* Full Smile Approx */, StartSwap, TheoEndSwap, 0,
      0.0, szSwaptionVolName, NumStrikes, Strikes, &CMSATMOption);

  err = srt_f_optimpvol(CMSATMOption, CMS + FwdSpread, CMS + FwdSpread, Mat,
                        1.0, SRT_PAYER, SRT_NORMAL, &ATMNormalimpvol);

  err = srt_f_optimpvol(CMSATMOption, CMS + FwdSpread, CMS + FwdSpread, Mat,
                        1.0, SRT_PAYER, SRT_LOGNORMAL, &ATMVolCMS);

  i = 0;

  while (Strikes[i] < CMS + FwdSpread) {
    CMTStrikes[i] = Strikes[i];
    CMTVols[i] = implied_vol[i];
    i++;
  }

  CMTStrikes[i] = CMS + FwdSpread;
  CMTVols[i] = ATMVolCMS;
  i++;

  while (i <= NumStrikes) {
    CMTStrikes[i] = Strikes[i - 1];
    CMTVols[i] = implied_vol[i - 1];
    i++;
  }

  err = opsabrcalib(CMS + FwdSpread, Mat, NumStrikes, CMTStrikes, CMTVols,
                    &ATMVolCMS, &alphacmt, 1, &betacmt, 0, &rhocmt, 1, &error);

  ATMNormalVolCMT = CHI * ATMNormalimpvol;

  ATMPrice = srt_f_optblknrm(CMS + FwdSpread, CMS + FwdSpread, ATMNormalVolCMT,
                             Mat, 1.0, SRT_CALL, PREMIUM);

  err = srt_f_optimpvol(ATMPrice, CMS + FwdSpread, CMS + FwdSpread, Mat, 1.0,
                        SRT_PAYER, SRT_LOGNORMAL, &ATMVolCMT);

  SabrComp[0] = CMS + FwdSpread;
  SabrComp[1] = ATMVolCMT;
  SabrComp[2] = alphacmt;
  SabrComp[3] = betacmt;
  SabrComp[4] = rhocmt;

  for (i = 0; i < NumStrikes; i++) {

    err = srt_f_optsarbvol(CMS + FwdSpread, Strikes[i], Mat, ATMVolCMT,
                           alphacmt, betacmt, rhocmt, SRT_LOGNORMAL,
                           SRT_LOGNORMAL, &(CMTVol[i]));
  }

  /*free memory*/
  free_dvector(CMSOption, 0, NumStrikes - 1);
  free_dvector(implied_vol, 0, NumStrikes - 1);
  free_dvector(CMTStrikes, 0, NumStrikes);
  free_dvector(CMTVols, 0, NumStrikes);

  return err;
}

Err swp_f_GetCMVols(char *szYieldCurveName, char *szSwaptionVolName,
                    char *cRefRateCode, double Underlying, char *cCMTFrequency,
                    char *cCMTBasis, char *cCMTRefRateCode, int NumStrikes,
                    double *Strikes, int NumMats, double *Mats, double *alpha,
                    double *beta, double *rho,
                    char *Method, /*MAD: CMT ATM Normal Vol = CMS ATM Normal Vol
                               - CHI: CMT ATM Normal Vol = CHI * CMS ATM Normal
                               Vol (In fact dLambda is taken to be CHI - ELSE:
                               Use OU on the Spread with 3 parameters*/
                    double *chi, double dLambda, double dSpreadLocVol,
                    double dCorrel, double **CMTVols, double **SabrComponents) {
  Err err = NULL;
  int i, j;
  long today, spotdate, spotlag;
  SrtCurvePtr Crv;
  SwapDP sdpSwap;
  long StartSwap, TheoEndSwap;
  double FwdSwapRate, FwdSpread;
  double CMS;
  double CMSVol;
  SrtCompounding CMTFreq;
  SrtBasisCode CMTBasis;

  err = interp_compounding(cCMTFrequency, &CMTFreq);
  err = interp_basis(cCMTBasis, &CMTBasis);
  Crv = lookup_curve(szYieldCurveName);
  today = get_clcndate_from_yldcrv(Crv);
  spotdate = get_spotdate_from_yldcrv(Crv);
  spotlag = get_spotlag_from_curve(Crv);

  /*Loop on all the different Maturities*/
  for (i = 0; i < NumMats; i++) {
    /*First compute the CMS rate*/

    /* initialise the sdpFra */
    TheoEndSwap = add_unit(spotdate, (long)(Mats[i] + Underlying) * 12,
                           SRT_MONTH, NO_BUSDAY_CONVENTION);
    StartSwap = add_unit(add_unit(TheoEndSwap, (long)-Underlying * 12,
                                  SRT_MONTH, NO_BUSDAY_CONVENTION),
                         0, SRT_BDAY, MODIFIED_SUCCEEDING);

    err = swp_f_setSwapDP(StartSwap, TheoEndSwap, CMTFreq, CMTBasis, &sdpSwap);
    if (err)
      return err;

    /* input the spot lag  */
    sdpSwap.spot_lag = spotlag;

    /* computation of the Fra */
    err = swp_f_ForwardRate_SwapDP(&sdpSwap, szYieldCurveName, cRefRateCode,
                                   &FwdSwapRate);

    if (err = swp_f_Cms_Rate(
            FwdSwapRate, Mats[i], Underlying * CMTFreq, (double)CMTFreq, 0, 1.0,
            SRT_LOGNORMAL, 0.0, 2 /* Full Smile Approx */, StartSwap,
            TheoEndSwap, 0, 0.0, szSwaptionVolName, NumStrikes, Strikes, &CMS))

    {
      return err;
    }

    FwdSpread = swp_f_spread(
        add_unit(StartSwap, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING), StartSwap,
        cCMTRefRateCode);

    if (strcmp(Method, "MAD") == 0) {
      err = swp_f_getcmtvolMAD(StartSwap, TheoEndSwap, szSwaptionVolName,
                               FwdSwapRate, CMS, FwdSpread, NumStrikes, Strikes,
                               Mats[i], alpha[i], beta[i], rho[i], CMTVols[i]);

    }

    else if (strcmp(Method, "CHI") == 0) {
      err = swp_f_getcmtvolCHI(StartSwap, TheoEndSwap, szSwaptionVolName,
                               FwdSwapRate, CMS, FwdSpread, (double)CMTFreq,
                               Underlying * CMTFreq, NumStrikes, Strikes,
                               Mats[i], alpha[i], beta[i], rho[i], chi[i],
                               CMTVols[i], SabrComponents[i]);
    }

    else {
      for (j = 0; j < NumStrikes; j++) {
        err = swp_f_getcmsvoladjusted(
            FwdSwapRate, StartSwap, TheoEndSwap, Mats[i], CMS, FwdSpread,
            Strikes[j], szSwaptionVolName, alpha[i], beta[i], rho[i], &CMSVol);
        err = swp_f_getcmtvol(CMS, FwdSpread, CMSVol, Strikes[j], dLambda,
                              dSpreadLocVol, dCorrel, Mats[i], &CMTVols[i][j]);
      }
    }
  }
  return err;
}
