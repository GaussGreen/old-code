#include "OPFNCTNS.H>
#include "SRT_H_ALL.H>
#include "math.h"
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

/* Helper function to get the cash fra */
LGMErr CashRate(Date tNow, Date Start, Date End, String ycName,
                double *outFwd) {
  double dfSt, dfEnd;

  dfSt = swp_f_df(tNow, Start, ycName);
  dfEnd = swp_f_df(tNow, End, ycName);
  if (dfSt == SRT_DF_ERROR || dfEnd == SRT_DF_ERROR)
    return ("no discount factor");
  *outFwd = (dfSt - dfEnd) / coverage(Start, End, BASIS_ACT_360) / dfEnd;
  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Routine called by CompVal_ForwardsTimeSwaps in order to calculate the PV of
 * one coupon of the forward time swaps */
LGMErr CompValCouponTimeSwap(
    Date tNow,               /* Today */
    String ycName,           /* Yield Curve name */
    SrtCallTimeSwapPtr deal, /* Pointer to the deal */
    int eod_pay_flag,        /* if we can pay today */
    int noDaysAccrued,       /* number of days accrued */
    long index_coupon, /* the index of the coupon being calculated (actually it
                          is index -1) */
    double floatingFixing,         /* past fixing if the coupon is floating */
    double *tfras_subperiod,       /* output:  the fwd Libors */
    double *tcms_subperiod,        /* output:  the CMS corrected fwd Libors */
    double *tmaturities_subperiod, /* output:  the time until the Libor fixes */
    double **tvols_subperiod,      /* output:  the log vols for each barrier */
    double *
        *tvols_discr, /* output:  market log vols at discretization points */
    long *index_discret_obsfreq, /* output:  indices between vol discreization
                                    dates and Libor observation dates */
    double *forward_pvFixed,     /* output:  value of the accrual leg */
    double
        *forward_pvPeriod) /* output:  value of the accrual and funding legs */
{
  SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
  LGMMarkConv conventions;   /* standard market conventions */
  LGMErr error = NULL;
  Date tEx;
  long j, l, lag;
  double pvCpn;
  double dfSt, dfPay;
  SrtBasisCode basis;
  double fra_accrued, fra_subperiod;
  double coverage_Libor, level;
  double strikelow, strikelowshift, strikeup, strikeupshift;
  double vollow, vollowshift, volup, volupshift;
  double maturity;
  double proba_above, proba_below, proba_between; /* Probs for barriers */
  double proba_above_adj, proba_below_adj,
      proba_between_adj; /* Adjusted probs for floating component */
  double expected_accrued_ratio, expected_accrued_ratio_adj;
  double callup, callupshift, callup_adj, callupshift_adj;
  double putlow, putlowshift, putlow_adj, putlowshift_adj;
  double volLiborup, volLiborlow;
  double correl, adjustmentup, adjustmentlow;
  long Start, End;
  long obsfreq;
  long num_subdiscret_obsfreq;
  long index, m;
  double xx;
  double fra_for_adj_cms, vol_for_adj_cms;
  long index_for_adj_cms;
  long start_accrued, end_accrued;
  double df_start_accrued, df_end_accrued;

  /* Check for coupons */
  if (deal->nCpn < 1)
    return ("no coupons");

  /* Initializations */
  basis = BASIS_ACT_360; /* basis money market for LIBOR */
  obsfreq = deal->observation_freq;
  num_subdiscret_obsfreq = deal->num_subdiscret_obsfreq;

  /* Computes the Libor obsevation index corresponding to the vol discretisation
   */
  for (index = 0; index < num_subdiscret_obsfreq; index++)
    index_discret_obsfreq[index] =
        (int)(index * (obsfreq - 1.0) / (num_subdiscret_obsfreq - 1.0));

  /* Get market information */
  yldcrv = lookup_curve(ycName);

  /* Get discount factors to last pay and start dates */
  dfPay = swp_f_df(tNow, deal->tCpnPay[index_coupon], ycName);
  dfSt = swp_f_df(tNow, deal->tCpnStart[index_coupon], ycName);
  if (dfSt == SRT_DF_ERROR || dfPay == SRT_DF_ERROR)
    return ("no discount factor");

  /* computes the four strikes */
  strikelow = deal->barriers[index_coupon][0];
  strikelowshift = strikelow + deal->call_spread_low;
  strikeup = deal->barriers[index_coupon][1];
  strikeupshift = strikeup + deal->call_spread_up;

  /* if we have a floating digital component */
  if (deal->gear[index_coupon]) {
    start_accrued = deal->tCpnStart[index_coupon];
    end_accrued = add_unit(start_accrued, deal->undAccrued, SRT_MONTH,
                           MODIFIED_SUCCEEDING);
    df_start_accrued = swp_f_df(tNow, start_accrued, ycName);
    df_end_accrued = swp_f_df(tNow, end_accrued, ycName);
    if (df_start_accrued == SRT_DF_ERROR || df_end_accrued == SRT_DF_ERROR)
      return ("no discount factor");
    coverage_Libor = coverage(start_accrued, end_accrued, basis);
    level = coverage_Libor * df_end_accrued;
    fra_accrued = (df_start_accrued - df_end_accrued) / level;
    error = GetVol_FudgeTimeSwap(
        tNow, tNow, start_accrued, end_accrued,
        /*	fra_accrued  ,		 Changes the adjustment		*/
        strikeup, fra_accrued, deal, &volLiborup);
    error = GetVol_FudgeTimeSwap(tNow, tNow, start_accrued, end_accrued,
                                 /*			fra_accrued  ,	 Changes
                                    the adjustment			*/
                                 strikelow, fra_accrued, deal, &volLiborlow);
  }

  /* Initializations */
  expected_accrued_ratio = 0;
  expected_accrued_ratio_adj = 0;
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  lag = conventions.lag; /* get lag */

  /* computes the fras at the observation days */
  for (j = 1; j <= obsfreq; j++) {
    Start = deal->observationdays[index_coupon][j - 1]
                                 [0]; /* Start date of the fra */
    End =
        deal->observationdays[index_coupon][j - 1][1]; /* End date of the fra */

    /* Computes the fixing time for the fra */
    tEx = add_unit(Start, -lag, SRT_BDAY, SUCCEEDING);
    tmaturities_subperiod[j - 1] = (double)(tEx - tNow + 0.0) / 365.0;

    /* Computes the cash fra */
    if (CashRate(tNow, Start, End, ycName, &tfras_subperiod[j - 1]))
      return ("no discount factor");
  }

  /* computes the cms at the observation days */
  /* Get the ATM log volatility for the fra at the final observation date */
  index_for_adj_cms = obsfreq - 1;
  fra_for_adj_cms = tfras_subperiod[index_for_adj_cms];
  error = GetVol_FudgeTimeSwap(
      tNow, tNow, deal->observationdays[index_coupon][index_for_adj_cms][0],
      deal->observationdays[index_coupon][index_for_adj_cms][1],
      fra_for_adj_cms, fra_for_adj_cms, deal, &vol_for_adj_cms);
  /* For each FRA apply the correction CMS_j = FRA_j * exp[ sigma^2 FRA * tEx_j
   * * (tEnd_j - tPay)  ] */
  for (j = 0; j < obsfreq; j++)
    tcms_subperiod[j] =
        tfras_subperiod[j] * exp(vol_for_adj_cms * vol_for_adj_cms *
                                 fra_for_adj_cms * tmaturities_subperiod[j] *
                                 (deal->observationdays[index_coupon][j][1] -
                                  deal->tCpnPay[index_coupon] + 0.0) /
                                 365.0);

  /* computes the vols at the discretisation points */
  for (index = 0; index < num_subdiscret_obsfreq; index++) {
    l = index_discret_obsfreq[index];
    /* Computes Start and End for the observation days */
    Start = deal->observationdays[index_coupon][l][0];
    End = deal->observationdays[index_coupon][l][1];
    /* gets volup and volupshift */
    error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikeup,
                                 tfras_subperiod[l], deal,
                                 &(tvols_discr[index][1]));
    error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikeupshift,
                                 tfras_subperiod[l], deal,
                                 &(tvols_discr[index][3]));
    if (0 < strikelow) /* gets vollow and vollowshift */
    {
      error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikelow,
                                   tfras_subperiod[l], deal,
                                   &(tvols_discr[index][0]));
      error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikelowshift,
                                   tfras_subperiod[l], deal,
                                   &(tvols_discr[index][2]));
    }
  }

  /* Computes all vols for the period by interpolation of vols at discretisation
   * points */
  index = 0;
  for (l = 0; l < obsfreq; l++) /* Loop on the observation days */
  {
    if (l > index_discret_obsfreq[index +
                                  1]) /* if we have passed the vol observation
                                         point  , move to the next one */
      index += 1;

    if (tmaturities_subperiod[index_discret_obsfreq[index + 1]] ==
        tmaturities_subperiod[index_discret_obsfreq[index]])
      xx = 0.5; /* if the subperiods dates are the same.  This shouldn't really
                   happen. */
    else
      xx = (tmaturities_subperiod[l] -
            tmaturities_subperiod[index_discret_obsfreq[index]]) /
           (tmaturities_subperiod[index_discret_obsfreq[index + 1]] -
            tmaturities_subperiod[index_discret_obsfreq[index]]);
    /* ratio of time between discretization and previous observation point and
     * time between surrounding observation points */
    for (m = 0; m <= 3; m++) /* linear interpolation on vols */
      tvols_subperiod[l][m] =
          (1 - xx) * tvols_discr[index][m] + xx * tvols_discr[index + 1][m];
  }

  /* computes probabilities to be between barriers  , and corresponding ratios
   */
  for (j = 1; j <= obsfreq; j++) {
    maturity = tmaturities_subperiod[j - 1];
    fra_subperiod = tcms_subperiod[j - 1];
    volup = tvols_subperiod[j - 1][1];
    volupshift = tvols_subperiod[j - 1][3];

    /* call spread */
    callup = srt_f_optblksch(fra_subperiod, strikeup, volup, maturity, 1,
                             SRT_CALL, PREMIUM);
    callupshift = srt_f_optblksch(fra_subperiod, strikeupshift, volupshift,
                                  maturity, 1, SRT_CALL, PREMIUM);

    /* probability above upper barrier */
    proba_above = (callup - callupshift) / deal->call_spread_up;

    /* LOWER BARRIER */
    /* probability below lower barrier */
    if (deal->barriers[index_coupon][0] == 0)
      proba_below = 0;
    else {
      vollow = tvols_subperiod[j - 1][0];
      vollowshift = tvols_subperiod[j - 1][2];

      /* put spread */
      putlow = srt_f_optblksch(fra_subperiod, strikelow, vollow, maturity, 1,
                               SRT_PUT, PREMIUM);
      putlowshift = srt_f_optblksch(fra_subperiod, strikelowshift, vollowshift,
                                    maturity, 1, SRT_PUT, PREMIUM);
      proba_below = (putlowshift - putlow) / deal->call_spread_low;
    }

    if (deal->gear[index_coupon] == 0) {
      proba_below_adj = 0;
      proba_above_adj = 0;
    } else /* Floating digital */
    {
      /* computes the correl (Libor  ,index) */
      correl =
          (maturity - tmaturities_subperiod[0]) /
              (tmaturities_subperiod[obsfreq - 1] - tmaturities_subperiod[0]) *
              deal->CorrelEnd +
          (tmaturities_subperiod[obsfreq - 1] - maturity) /
              (tmaturities_subperiod[obsfreq - 1] - tmaturities_subperiod[0]) *
              deal->CorrelStart;
      adjustmentup = exp(correl * volLiborup * volup *
                         tmaturities_subperiod[0]); /* Changes the adjustment */
      callup_adj = srt_f_optblksch(fra_subperiod * adjustmentup, strikeup,
                                   volup, maturity, 1, SRT_CALL, PREMIUM);
      callupshift_adj =
          srt_f_optblksch(fra_subperiod * adjustmentup, strikeupshift,
                          volupshift, maturity, 1, SRT_CALL, PREMIUM);
      proba_above_adj = (callup_adj - callupshift_adj) / deal->call_spread_up;
      if (deal->barriers[index_coupon][0] == 0)
        proba_below_adj = 0;
      else {
        adjustmentlow =
            exp(correl * volLiborlow * vollow *
                tmaturities_subperiod[0]); /* Changes the adjustment */
        putlow_adj = srt_f_optblksch(fra_subperiod * adjustmentlow, strikelow,
                                     vollow, maturity, 1, SRT_PUT, PREMIUM);
        putlowshift_adj =
            srt_f_optblksch(fra_subperiod * adjustmentlow, strikelowshift,
                            vollowshift, maturity, 1, SRT_PUT, PREMIUM);
        proba_below_adj =
            (putlowshift_adj - putlow_adj) / deal->call_spread_low;
      }
    }
    proba_between = 1 - proba_below - proba_above;
    proba_between_adj = 1 - proba_below_adj - proba_above_adj;

    expected_accrued_ratio +=
        proba_between * (deal->ratiodays_for_subperiod)[index_coupon][j - 1];
    expected_accrued_ratio_adj +=
        proba_between_adj *
        (deal->ratiodays_for_subperiod)[index_coupon][j - 1];
  }

  /* Fix Coupon */
  pvCpn = dfPay * deal->tCvgCpn[index_coupon] * deal->tCpn[index_coupon] *
          expected_accrued_ratio;

  /* Floating Digital */
  pvCpn += dfPay * deal->tCvgCpn[index_coupon] * deal->gear[index_coupon] *
           fra_accrued * expected_accrued_ratio_adj;

  /* Past Fixings */
  if (deal->tCpnPay[index_coupon] < tNow + eod_pay_flag) {
    pvCpn +=
        dfPay * deal->tCvgCpn[index_coupon] * deal->tCpn[index_coupon] *
        (noDaysAccrued + 0.0) /
        (deal->tCpnPay[index_coupon] - deal->tCpnStart[index_coupon] + 0.0);
    pvCpn +=
        dfPay * deal->tCvgCpn[index_coupon] * deal->gear[index_coupon] *
        floatingFixing * (noDaysAccrued + 0.0) /
        (deal->tCpnPay[index_coupon] - deal->tCpnStart[index_coupon] + 0.0);
  }

  *forward_pvFixed = pvCpn;

  /* Funding Margin (only for option valuation) includes margin and Libor leg */
  pvCpn += dfPay * deal->tFundingPayment[index_coupon];
  dfSt = swp_f_df(tNow, deal->tCpnStart[index_coupon], ycName);
  if (deal->nofunding == 0)
    pvCpn += dfPay - dfSt;

  *forward_pvPeriod = pvCpn;
  return NULL;
}

/* -----------------------------------------------------------------------------------------------------------------------------------
 */
/* Routine called by LGMExtractInfoFromDeals to get out forward value fo the
 * time swaps for possible use in determining the equivalent */
/* strikes */
LGMErr CompVal_ForwardsTimeSwaps(
    SrtCallTimeSwapPtr deal, /* in:  Pointer to the CTS */
    Date tNow,               /* in:  calculation date */
    String ycName,           /* in:  yield curve name */
    double *ratios,          /* out: first index is 1 !!! */
    double *intVal, /* out: intrinsic value (of coupon  , or remaining swap?) */
    double *forwardsTimeSwaps, /* out: vector of forward time swap value PVs */
    double *pvfixedleg,        /* out: the PV of the fixed leg (accrual leg) */
    int nEx,                   /* number of exercise dates */
    int eod_pay_flag,          /* end of day pay flag (0) */
    int noDaysAccrued,         /* number of days accrued (0) */
    double floatingfixing)     /* (0) */
{
  /* Variable declarations */
  LGMErr err = NULL;
  long obsfreq, num_subdiscret_obsfreq;
  double *tfras_subperiod = NULL;
  double *tcms_subperiod = NULL;
  double *tmaturities_subperiod = NULL;
  double **tvols_subperiod = NULL;
  double **tvols_discr = NULL;
  long *index_discret_obsfreq = NULL;
  long index_exo;
  long index_coupon;
  double forward_pvFixed, forward_pvPeriod;
  double dfSt;

  /* Initializations */
  obsfreq = deal->observation_freq;
  num_subdiscret_obsfreq = deal->num_subdiscret_obsfreq;
  index_coupon = deal->nCpn;
  (*pvfixedleg) = 0;
  (*intVal) = 0;

  /* Memory allocation */
  tvols_discr = dmatrix(0, num_subdiscret_obsfreq - 1, 0, 3);
  tfras_subperiod = dvector(0, obsfreq - 1);
  tcms_subperiod = dvector(0, obsfreq - 1);
  tmaturities_subperiod = dvector(0, obsfreq - 1);
  tvols_subperiod = dmatrix(0, obsfreq - 1, 0, 3);
  index_discret_obsfreq = lngvector(0, num_subdiscret_obsfreq - 1);

  for (index_coupon = deal->nCpn - 1; index_coupon >= deal->iSet[0];
       index_coupon--) {
    /* compute a number of quantities for each coupon */
    err = CompValCouponTimeSwap(
        tNow,          /* in:  today */
        ycName,        /* in:  yield curve name */
        deal,          /* in:  pointer to the CTS */
        eod_pay_flag,  /* in:  can we pay today */
        noDaysAccrued, /* in:  number of days accrued so far (0) */
        index_coupon, /* in:  the index of the coupon that we are calculating */
        floatingfixing, /* in:  the fixing for the Libor component of the coupon
                         */
        tfras_subperiod, /* out: the fwd rates at the Libor observation dates*/
        tcms_subperiod,  /* out: the adjusted fwd rates at the Libor observation
                            dates */
        tmaturities_subperiod, /* out: the time (from today) to the fixing for
                                  each Libor observation dates */
        tvols_subperiod, /* out: the log volatilities (interpolated) for each
                            Libor observation date */
        tvols_discr,     /* out: the log volatilities at the volatility
                            discretization dates */
        index_discret_obsfreq, /* out: the map from the discretization dates to
                                  the observation dates */
        &forward_pvFixed,      /* out: the PV of the accrual leg */
        &forward_pvPeriod); /* out: the PV of the period (accrual - funding) */

    /* keep track of the accrual leg PV */
    (*pvfixedleg) += forward_pvFixed;

    /* loop over the exercise dates before the current coupon adding the pv to
     * the (earlier) forward time swap */
    for (index_exo = 0;
         (deal->iSet[index_exo] < index_coupon) && (index_exo < nEx);
         index_exo++)
      forwardsTimeSwaps[index_exo] += forward_pvPeriod;

    /* if the exercise index is equal to the index coupon  , then we have
     * completed the forward time swap */
    if (deal->iSet[index_exo] == index_coupon) {
      /* add in the pv.  The forwardTS is the sum of the accrual coupons  , the
       * float margins  , the final DF minus the start DF */
      forwardsTimeSwaps[index_exo] += forward_pvPeriod;

      /* calculate the discount factor at the start of the forward time swap */
      dfSt = swp_f_df(tNow, deal->tSet[index_exo], ycName);

      /* calculate the ratio (notice index is 1 higher) of the fixed PV to the
         floating PV (including the final but not the initial notional
         repayment). */
      ratios[index_exo + 1] = (forwardsTimeSwaps[index_exo] + dfSt) / dfSt;

      /* If the ratio is too high or too low  , set it to 0 (not sure what the
       * impact of this is) */
      if (ratios[index_exo + 1] > 5.0 || ratios[index_exo + 1] < 0.2)
        ratios[index_exo + 1] = 0.;

      /* Calculate what the intrinsic value is (0 if OTM  , otherwise the fwdTS
       * value).  */
      if (deal->PayRec == SRT_RECEIVER) {
        if ((*intVal) < forwardsTimeSwaps[index_exo])
          (*intVal) = forwardsTimeSwaps[index_exo];
      } else {
        if ((*intVal) < (-forwardsTimeSwaps[index_exo]))
          (*intVal) = (-forwardsTimeSwaps[index_exo]);
      }
    } /* end if statement */

  } /* end for loop on coupon index */

  /* Memory freeage */
  free_dmatrix(tvols_discr, 0, num_subdiscret_obsfreq - 1, 0, 3);
  free_dvector(tfras_subperiod, 0, obsfreq - 1);
  free_dvector(tcms_subperiod, 0, obsfreq - 1);
  free_dvector(tmaturities_subperiod, 0, obsfreq - 1);
  free_dmatrix(tvols_subperiod, 0, obsfreq - 1, 0, 3);
  free_lngvector(index_discret_obsfreq, 0, num_subdiscret_obsfreq - 1);

  return NULL;
}

LGMErr CompVal_TimeSwap(SrtCallTimeSwapPtr deal, Date tNow, String ycName,
                        double *pvfixedleg, int eod_pay_flag, int noDaysAccrued,
                        double floatingfixing) {
  LGMErr err = NULL;
  long obsfreq, num_subdiscret_obsfreq;
  double *tfras_subperiod = NULL;
  double *tcms_subperiod = NULL;
  double *tmaturities_subperiod = NULL;
  double **tvols_subperiod = NULL;
  double **tvols_discr = NULL;
  long *index_discret_obsfreq = NULL;
  long index_coupon;
  double forward_pvFixed, forward_pvPeriod;

  obsfreq = deal->observation_freq;
  num_subdiscret_obsfreq = deal->num_subdiscret_obsfreq;

  /* Memory allocation */
  tvols_discr = dmatrix(0, num_subdiscret_obsfreq - 1, 0, 3);
  tfras_subperiod = dvector(0, obsfreq - 1);
  tcms_subperiod = dvector(0, obsfreq - 1);
  tmaturities_subperiod = dvector(0, obsfreq - 1);
  tvols_subperiod = dmatrix(0, obsfreq - 1, 0, 3);
  index_discret_obsfreq = lngvector(0, num_subdiscret_obsfreq - 1);

  (*pvfixedleg) = 0;
  //	*pvSwap = 0.0;

  for (index_coupon = 0; index_coupon < deal->nCpn; index_coupon++) {
    // Calculate the PV of a single coupon
    if (err = CompValCouponTimeSwap(
            tNow, ycName, deal, eod_pay_flag, noDaysAccrued, index_coupon,
            floatingfixing, tfras_subperiod, tcms_subperiod,
            tmaturities_subperiod, tvols_subperiod, tvols_discr,
            index_discret_obsfreq, &forward_pvFixed, &forward_pvPeriod))
      return err;
    // Add it on to the fixed leg PV.
    (*pvfixedleg) += forward_pvFixed;
    //		*pvSwap += forward_pvPeriod;
  }

  /* Memory freeage */
  free_dmatrix(tvols_discr, 0, num_subdiscret_obsfreq - 1, 0, 3);
  free_dvector(tfras_subperiod, 0, obsfreq - 1);
  free_dvector(tcms_subperiod, 0, obsfreq - 1);
  free_dvector(tmaturities_subperiod, 0, obsfreq - 1);
  free_dmatrix(tvols_subperiod, 0, obsfreq - 1, 0, 3);
  free_lngvector(index_discret_obsfreq, 0, num_subdiscret_obsfreq - 1);

  return NULL;
}

/********************************************************************************/
/*	Returns intrinsic value and the ratio PV of coupon leg to PV of the
strike for a timeswap. If deal is strictly positive or strictly negative  ,
ratio is set to 0. */
LGMErr CompValTimeSwap(
    Date tNow,                        /* Value date of the deal */
    String ycName,                    /* yield curve */
    double *ratioPtr,                 /* ratio */
    long indexexo, double *intValPtr, /* intrinsic value */
    double Strike, /* fee to pay by the folder of the option */
    SrtCallTimeSwapPtr deal, int bullet_or_option, /* 0=bullet  ,1=option */
    int eod_pay_flag,
    int noDaysAccrued, /* for period which has already started */
    double
        floatingFixing /* fixing at the beginning of the period */) /* basis for
                                                                       rates */
{
  SrtCurvePtr yldcrv = NULL; /* yield curve pointer */
  LGMMarkConv conventions;   /* standard market conventions */
  LGMErr error = NULL;
  Date tEx;
  long i, j, l, lag;
  double pvCpn;
  double dfSt, dfEnd, dfPay;
  double sum, ratio, intVal;
  SrtBasisCode basis;
  double fra_accrued, fra_subperiod;
  double coverage_Libor, level;
  double strikelow, strikelowshift, strikeup, strikeupshift;
  double vollow, vollowshift, volup, volupshift;
  double maturity;
  double proba_above, proba_above_adj, proba_below, proba_below_adj,
      proba_between, proba_between_adj;
  double expected_accrued_ratio, expected_accrued_ratio_adj;
  double callup, callupshift, callup_adj, callupshift_adj;
  double putlow, putlowshift, putlow_adj, putlowshift_adj;
  double volLiborup, volLiborlow;
  double correl, adjustmentup, adjustmentlow;
  long Start, End;
  long obsfreq;
  double **tvols_discr = NULL; /* low  ,up  ,lowshift  ,upshift */
  double *tfras_subperiod = NULL;
  double *tcms_subperiod = NULL;
  double *tmaturities_subperiod = NULL;
  double **tvols_subperiod = NULL;
  long *index_discret_obsfreq = NULL;
  long num_subdiscret_obsfreq;
  long index, m;
  double xx;
  double fra_for_adj_cms, vol_for_adj_cms;
  long index_for_adj_cms;
  long start_accrued, end_accrued;
  double df_start_accrued, df_end_accrued;

  yldcrv = lookup_curve(ycName);
  error = LGMCcyDefaults(yldcrv, &conventions); /* get currency conventions */
  if (error != NULL)
    return (error);
  lag = conventions.lag; /* get lag */
  if (deal->nCpn < 1)
    return ("no coupons");

  basis = BASIS_ACT_360; /* basis money market for LIBOR */

  sum = 0.;
  obsfreq = deal->observation_freq;
  num_subdiscret_obsfreq = deal->num_subdiscret_obsfreq;

  /* Memory allocation */
  tvols_discr = dmatrix(0, num_subdiscret_obsfreq - 1, 0, 3);
  tfras_subperiod = dvector(0, obsfreq - 1);
  tcms_subperiod = dvector(0, obsfreq - 1);
  tmaturities_subperiod = dvector(0, obsfreq - 1);
  tvols_subperiod = dmatrix(0, obsfreq - 1, 0, 3);
  index_discret_obsfreq = lngvector(0, num_subdiscret_obsfreq - 1);

  /* Computes indexes of discretisation for the observation days */
  for (index = 0; index < num_subdiscret_obsfreq; index++) {
    index_discret_obsfreq[index] =
        (int)(index * (obsfreq - 1.0) / (num_subdiscret_obsfreq - 1.0));
  }

  for (i = deal->iSet[indexexo] + 1; i <= deal->nCpn;
       i++) /* for each period */ /* find PV of coupon leg */
  {

    dfPay = swp_f_df(tNow, deal->tCpnPay[i - 1],
                     ycName); /* calculate discount factor to pay date */
    /*	if (bullet_or_option==0 && deal->tCpnPay[i-1] < (today + eod_pay_flag))
            {
            */

    dfSt = swp_f_df(tNow, deal->tCpnStart[i - 1], ycName);
    if (dfSt == SRT_DF_ERROR || dfPay == SRT_DF_ERROR)
      return ("no discount factor");
    /*	coverage_Libor = coverage(deal->tCpnStart[i-1]  ,deal->tCpnPay[i-1]
       ,basis); level = coverage_Libor * dfPay; fra_period = (dfSt -
       dfPay)/level;   */

    /* computes the four strikes */
    strikelow = deal->barriers[i - 1][0];
    strikelowshift = strikelow + deal->call_spread_low;
    strikeup = deal->barriers[i - 1][1];
    strikeupshift = strikeup + deal->call_spread_up;

    if (!(0 == deal->gear[i - 1])) /* floating digital computations */
    {
      start_accrued = deal->tCpnStart[i - 1];
      end_accrued = add_unit(start_accrued, deal->undAccrued, SRT_MONTH,
                             MODIFIED_SUCCEEDING);
      df_start_accrued = swp_f_df(tNow, start_accrued, ycName);
      df_end_accrued = swp_f_df(tNow, end_accrued, ycName);
      if (df_start_accrued == SRT_DF_ERROR || df_end_accrued == SRT_DF_ERROR)
        return ("no discount factor");
      coverage_Libor = coverage(start_accrued, end_accrued, basis);
      level = coverage_Libor * df_end_accrued;
      fra_accrued = (df_start_accrued - df_end_accrued) / level;
      error = GetVol_FudgeTimeSwap(
          tNow, tNow, start_accrued, end_accrued,
          /*	fra_accrued  ,		 Changes the adjustment		*/
          strikeup, fra_accrued, deal, &volLiborup);
      error = GetVol_FudgeTimeSwap(tNow, tNow, start_accrued, end_accrued,
                                   /*			fra_accrued  ,	 Changes the
                                      adjustment			*/
                                   strikelow, fra_accrued, deal, &volLiborlow);
    }

    expected_accrued_ratio = 0;
    expected_accrued_ratio_adj = 0;

    /* computes the fras at the observation days */
    for (j = 1; j <= obsfreq; j++) {
      Start = deal->observationdays[i - 1][j - 1][0];
      End = deal->observationdays[i - 1][j - 1][1];

      /* Computes the maturity  */
      tEx = add_unit(Start, -lag, SRT_BDAY, SUCCEEDING);
      tmaturities_subperiod[j - 1] = (double)(tEx - tNow + 0.0) / 365.0;

      /* Computes the fra */
      dfSt = swp_f_df(tNow, Start, ycName);
      dfEnd = swp_f_df(tNow, End, ycName);
      if (dfSt == SRT_DF_ERROR || dfEnd == SRT_DF_ERROR)
        return ("no discount factor");
      coverage_Libor = coverage(Start, End, basis);
      level = coverage_Libor * dfEnd;
      tfras_subperiod[j - 1] = (dfSt - dfEnd) / level;
    }

    /* computes the cms at the observation days */
    index_for_adj_cms = obsfreq - 1;
    fra_for_adj_cms = tfras_subperiod[index_for_adj_cms];
    error = GetVol_FudgeTimeSwap(
        tNow, tNow, deal->observationdays[i - 1][index_for_adj_cms][0],
        deal->observationdays[i - 1][index_for_adj_cms][1], fra_for_adj_cms,
        fra_for_adj_cms, deal, &vol_for_adj_cms);
    for (j = 1; j <= obsfreq; j++) {
      tcms_subperiod[j - 1] =
          tfras_subperiod[j - 1] *
          exp(vol_for_adj_cms * vol_for_adj_cms * fra_for_adj_cms *
              tmaturities_subperiod[j - 1] *
              (deal->observationdays[i - 1][j - 1][1] - deal->tCpnPay[i - 1] +
               0.0) /
              365.0);
    }

    /* computes the vols at the discretisation points */
    for (index = 0; index < num_subdiscret_obsfreq; index++) {
      l = index_discret_obsfreq[index];
      /* Computes Start and End for the observation days */
      Start = deal->observationdays[i - 1][l][0];
      End = deal->observationdays[i - 1][l][1];
      /* gets volup and volupshift */
      error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikeup,
                                   tfras_subperiod[l], deal,
                                   &(tvols_discr[index][1]));
      error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikeupshift,
                                   tfras_subperiod[l], deal,
                                   &(tvols_discr[index][3]));
      if (0 < strikelow) { /* gets vollow and vollowshift */
        error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikelow,
                                     tfras_subperiod[l], deal,
                                     &(tvols_discr[index][0]));
        error = GetVol_FudgeTimeSwap(tNow, tNow, Start, End, strikelowshift,
                                     tfras_subperiod[l], deal,
                                     &(tvols_discr[index][2]));
      }
    }

    /* Computes all vols for the period by interpolation of vols at
     * discretisation points */
    index = 0;
    for (l = 0; l < obsfreq; l++) {
      if (l > index_discret_obsfreq[index + 1]) {
        index += 1;
      }
      if (tmaturities_subperiod[index_discret_obsfreq[index + 1]] ==
          tmaturities_subperiod[index_discret_obsfreq[index]]) {
        xx = 0.5;
      } else {
        xx = (tmaturities_subperiod[l] -
              tmaturities_subperiod[index_discret_obsfreq[index]]) /
             (tmaturities_subperiod[index_discret_obsfreq[index + 1]] -
              tmaturities_subperiod[index_discret_obsfreq[index]]);
      }
      for (m = 0; m <= 3; m++) {
        /* TEMPORARY */
        tvols_subperiod[l][m] =
            (1 - xx) * tvols_discr[index][m] + xx * tvols_discr[index + 1][m];
      }
    } /* Loop on the observation days */

    /* computes probabilities to be between barriers  , and corresponding ratios
     */
    for (j = 1; j <= obsfreq; j++) {
      maturity = tmaturities_subperiod[j - 1];
      fra_subperiod = tcms_subperiod[j - 1];
      volup = tvols_subperiod[j - 1][1];
      volupshift = tvols_subperiod[j - 1][3];

      /* call spread */
      callup = srt_f_optblksch(fra_subperiod, strikeup, volup, maturity, 1,
                               SRT_CALL, PREMIUM);
      callupshift = srt_f_optblksch(fra_subperiod, strikeupshift, volupshift,
                                    maturity, 1, SRT_CALL, PREMIUM);

      /* probability above upper barrier */
      proba_above = (callup - callupshift) / deal->call_spread_up;
      /*				*/

      /* LOWER BARRIER */
      /* probability below lower barrier */
      if (deal->barriers[i - 1][0] == 0) {
        proba_below = 0;
      } else {
        /*		error = GetVol_FudgeTimeSwap(
                                                tNow  ,
                                                tNow  ,
                                                Start  ,
                                                End  ,
                                                strikelow  ,
                                                fra_subperiod  ,
                                                deal  ,
                                                &vollow);
                        error = GetVol_FudgeTimeSwap(
                                                tNow  ,
                                                tNow  ,
                                                Start  ,
                                                End  ,
                                                strikelowshift  ,
                                                fra_subperiod  ,
                                                deal  ,
                                                &vollowshift);
                                                */

        vollow = tvols_subperiod[j - 1][0];
        vollowshift = tvols_subperiod[j - 1][2];

        /* put spread */
        putlow = srt_f_optblksch(fra_subperiod, strikelow, vollow, maturity, 1,
                                 SRT_PUT, PREMIUM);
        putlowshift =
            srt_f_optblksch(fra_subperiod, strikelowshift, vollowshift,
                            maturity, 1, SRT_PUT, PREMIUM);

        proba_below = (putlowshift - putlow) / deal->call_spread_low;
      }

      if (deal->gear[i - 1] == 0) {
        proba_below_adj = 0;
        proba_above_adj = 0;
      } else /* Floating digital */
      {
        /* computes the correl (Libor  ,index) */
        correl = (maturity - tmaturities_subperiod[0]) /
                     (tmaturities_subperiod[obsfreq - 1] -
                      tmaturities_subperiod[0]) *
                     deal->CorrelEnd +
                 (tmaturities_subperiod[obsfreq - 1] - maturity) /
                     (tmaturities_subperiod[obsfreq - 1] -
                      tmaturities_subperiod[0]) *
                     deal->CorrelStart;
        adjustmentup =
            exp(correl * volLiborup * volup *
                tmaturities_subperiod[0]); /* Changes the adjustment */
        callup_adj = srt_f_optblksch(fra_subperiod * adjustmentup, strikeup,
                                     volup, maturity, 1, SRT_CALL, PREMIUM);
        callupshift_adj =
            srt_f_optblksch(fra_subperiod * adjustmentup, strikeupshift,
                            volupshift, maturity, 1, SRT_CALL, PREMIUM);
        proba_above_adj = (callup_adj - callupshift_adj) / deal->call_spread_up;
        if (deal->barriers[i - 1][0] == 0) {
          proba_below_adj = 0;
        } else {
          adjustmentlow =
              exp(correl * volLiborlow * vollow *
                  tmaturities_subperiod[0]); /* Changes the adjustment */
          putlow_adj = srt_f_optblksch(fra_subperiod * adjustmentlow, strikelow,
                                       vollow, maturity, 1, SRT_PUT, PREMIUM);
          putlowshift_adj =
              srt_f_optblksch(fra_subperiod * adjustmentlow, strikelowshift,
                              vollowshift, maturity, 1, SRT_PUT, PREMIUM);
          proba_below_adj =
              (putlowshift_adj - putlow_adj) / deal->call_spread_low;
        }
      }
      proba_between = 1 - proba_below - proba_above;
      proba_between_adj = 1 - proba_below_adj - proba_above_adj;

      expected_accrued_ratio +=
          proba_between * (deal->ratiodays_for_subperiod)[i - 1][j - 1];
      expected_accrued_ratio_adj +=
          proba_between_adj * (deal->ratiodays_for_subperiod)[i - 1][j - 1];
    }

    /* TO BE ERASED...JUST FOR TESTING
expected_accrued_ratio = 1;
expected_accrued_ratio_adj = 1;   */

    /* Fix Coupon */
    pvCpn = dfPay * deal->tCvgCpn[i - 1] * deal->tCpn[i - 1] *
            expected_accrued_ratio;
    /* Floating Digital */
    pvCpn += dfPay * deal->tCvgCpn[i - 1] * deal->gear[i - 1] * fra_accrued *
             expected_accrued_ratio_adj;
    /* Past Fixings */
    if (bullet_or_option == 0 && deal->tCpnPay[i - 1] < tNow + eod_pay_flag) {
      pvCpn += dfPay * deal->tCvgCpn[i - 1] * deal->tCpn[i - 1] *
               (noDaysAccrued + 0.0) /
               (deal->tCpnPay[i - 1] - deal->tCpnStart[i - 1] + 0.0);
      pvCpn += dfPay * deal->tCvgCpn[i - 1] * deal->gear[i - 1] *
               floatingFixing * (noDaysAccrued + 0.0) /
               (deal->tCpnPay[i - 1] - deal->tCpnStart[i - 1] + 0.0);
    }

    /* Funding Margin (only for opion valuation)*/
    if (bullet_or_option == 1) {
      pvCpn += dfPay * deal->tFundingPayment[i - 1];
    }
    sum = sum + pvCpn;
    if (deal->nCpn == i &&
        bullet_or_option == 1) /* final notional (only for option valuation)  */
    {
      sum += dfPay;
    }
  }

  /*
          dfSt = swp_f_df(tNow  ,deal->tCpnStart[indexexo]  ,ycName);
  */
  dfSt = swp_f_df(tNow, deal->tSet[indexexo], ycName);

  intVal = sum;
  if (bullet_or_option == 1) /* for bullet  , no floating leg */
    intVal -= dfSt * Strike; /* intrinsic value */

  /* David's trade*/
  if (deal->nofunding == 1) {
    intVal += dfSt * Strike - dfPay;
  }

  if (deal->PayRec == SRT_PAYER)
    intVal = -intVal;

  if (Strike == 0.)
    ratio = 0.0;
  else
    ratio = sum / (dfSt * Strike); /* ratio of PV's */

  if (ratio > 5.0 || ratio < 0.2) /* clearly  , one sign deal */
    ratio = 0.;

  /* Memory freeage */
  free_dmatrix(tvols_discr, 0, num_subdiscret_obsfreq - 1, 0, 3);
  free_dvector(tfras_subperiod, 0, obsfreq - 1);
  free_dvector(tcms_subperiod, 0, obsfreq - 1);
  free_dvector(tmaturities_subperiod, 0, obsfreq - 1);
  free_dmatrix(tvols_subperiod, 0, obsfreq - 1, 0, 3);
  free_lngvector(index_discret_obsfreq, 0, num_subdiscret_obsfreq - 1);

  *ratioPtr = ratio;
  *intValPtr = intVal;
  return (NULL);
}
