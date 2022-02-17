
#include "LGM2FMC.h"
#include "LGMSVUtil.h"
#include "RainbowOpt.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define BIG_BIG_NUMBER 1E40

static int nCheckCall = 1;

/*	Functions for the funding leg */

Err ccf_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */
    double fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    CF_FUND_LEG fund_leg) {
  int i, j;
  SrtBasisCode bas;
  CF_FUND_CPN cpn;
  Err err = NULL;

  /*	Notional */
  fund_leg->notional = fund_not;

  /*	Initialise pointers to NULL */
  fund_leg->cpn = NULL;

  /*	Skip coupons fixed before today */
  i = 0;
  while (i < fund_ncpn && fund_fix[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one coupon is left */
  if (i == fund_ncpn) {
    /*	err = "All funding coupons start before today in ccf_fill_fund_leg"; */
    fund_leg->num_cpn = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  fund_leg->num_cpn = fund_ncpn - i;
  fund_leg->cpn = (cf_fund_cpn *)calloc(fund_leg->num_cpn, sizeof(cf_fund_cpn));
  if (!fund_leg->cpn) {
    err = "Allocation error in ccf_fill_fund_leg";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;
  while (i < fund_ncpn) {
    cpn = fund_leg->cpn + j;

    /*	Dates */
    cpn->start_date = fund_start[i];
    cpn->pay_date = fund_pay[i];

    /*	Times */
    cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
    cpn->pay_time = (cpn->pay_date - today) * YEARS_IN_DAY;

    /*	Coupon */
    err = interp_basis(fund_basis[i], &bas);
    if (err) {
      goto FREE_RETURN;
    }
    cpn->cvg = coverage(fund_start[i], fund_pay[i], bas);
    cpn->cpn = fund_not * cpn->cvg * (fund_spr[i] + fund_mrg[i]);

    i++;
    j++;
  }

  err = ccf_check_fund_leg(fund_leg);

FREE_RETURN:

  if (err) {
    ccf_free_fund_leg(fund_leg);
  }

  return err;
}

/*	Check dates consistency */
Err ccf_check_fund_leg(CF_FUND_LEG fund_leg) {
  int i;

  /*	Check that start and pay dates are increasing */
  for (i = 1; i < fund_leg->num_cpn; i++) {
    if (fund_leg->cpn[i].start_date < fund_leg->cpn[i - 1].start_date) {
      return "Start dates should be increasing in funding leg";
    }

    if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i - 1].pay_date) {
      return "Pay dates should be increasing in funding leg";
    }
  }

  /*	Check that pay dates are after start dates */
  for (i = 0; i < fund_leg->num_cpn; i++) {
    if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i].start_date) {
      return "Pay dates should be after start dates in funding leg";
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err ccf_free_fund_leg(CF_FUND_LEG fund_leg) {
  if (fund_leg->cpn) {
    free(fund_leg->cpn);
    fund_leg->cpn = NULL;
  }

  return NULL;
}

#define ONE_MONTH 0.083333333
/*	Functions for the exotic leg */
/*	ATTENTION: cms lambdas are not filled and must be handled separetely */
Err ccf_fill_exo_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */
    double cf_not, int cf_ncpn, long *cf_fix, long *cf_start, long *cf_pay,
    char **cf_basis, char **cf_cms_tenor1, char *cf_cms_freq1,
    char *cf_cms_basis1, double *cf_cms_spread1, char **cf_cms_tenor2,
    char *cf_cms_freq2, char *cf_cms_basis2, double *cf_cms_spread2,
    double *cf_alpha, double *cf_beta, double *cf_gamma, int *cf_floored,
    double *cf_floor, int *cf_capped, double *cf_cap,

    long spread_vol_n, double *spread_vol_time,
    double *spread_slbeta1, /* Shifted log beta on the CMS1 */
    double *spread_slbeta2, /* Shifted log beta on the CMS2 */

    int cf_nopt,           /* Number of spread options */
    double **cf_notopt,    /* Notional of the spread options */
    double **cf_strikeopt, /* spread option strikes */
    int **cf_typeopt,      /* spread option type 0 call 1 put */

    int use_SL,   /*  1: use Shifted Log
                      0: don't use				*/
    int calib_SL, /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                                  0: use the Shifted Log beta given */
    double NbStdcalib_SL, /*  Nb of std for the calibration of to the skew */
    int calib_correl_SL,  /*	1: Calibrate the correlation between the two SL
                             to get the same ATM  normal spread vol  0: use the
                             normal spread correl for the sl correl */
    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                               the floor and cap in the cf
                          coupons
                                       0: take into account the floor and cap in
                          the cf coupons */

    int cms_adj, int cms_smile, int cms_vol_adj, /*	1: adjust for CMS vol
                                                    effect 0: don't */
    double cms_beta1, double cms_beta2,
    int num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double *strikes_in_vol,
    SrtDiffusionType
        vol_type,     /*	Type of vol in matrix  , SRT_NORMAL or SRT_LOGNORMAL */
    int cash_vol,     /*	1: matrix is a cash vol
/*	Needed to calculate ATM (converging) normal volatility associated to
    CMS coupons */
    char *yc,         /*	yc */
    char *vc,         /*	vc */
    char *ref,        /*	ref rate */
    char *swap_freq,  /*	swap freq */
    char *swap_basis, /*	swap basis */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    CF_EXO_LEG exo_leg) {
  double temp, strike, libor;
  int i, j, k, l;
  SrtBasisCode bas;
  SwapDP Swap;
  SrtDateList SwapDates;
  CF_EXO_CPN cpn;
  CF_CMS_DESC cms;
  int use_cms2;
  char *cmstenor, *cmsfreq, *cmsbasis;
  double cmsspr, cmsbeta, rconv, forward;
  int cmscpnd;
  SrtBasisCode cmsbas;
  int free_dates = 0;
  Err err = NULL;

  double *spread_slbeta = NULL;
  int iNumOpt;
  double dStrike1, dStrike2, dLogVol1, dLogVol2;

  /*	Notional */
  exo_leg->notional = cf_not;

  /*	Initialise pointers to NULL */
  exo_leg->cpn = NULL;

  /*	Skip coupons fixed before today */
  i = 0;
  while (i < cf_ncpn && cf_fix[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one coupon is left */
  if (i == cf_ncpn) {
    /* err = "All funding coupons start before today in ccf_fill_exo_leg"; */
    exo_leg->num_cpn = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  exo_leg->num_cpn = cf_ncpn - i;
  exo_leg->cpn = (cf_exo_cpn *)calloc(exo_leg->num_cpn, sizeof(cf_exo_cpn));
  if (!exo_leg->cpn) {
    err = "Allocation error in ccf_fill_exo_leg";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;
  while (i < cf_ncpn) {
    cpn = exo_leg->cpn + j;

    /*	Dates */
    cpn->start_date = cf_start[i];
    cpn->pay_date = cf_pay[i];
    cpn->cms_fix_date = cf_fix[i];

    /*	Times */
    cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
    cpn->pay_time = (cpn->pay_date - today) * YEARS_IN_DAY;
    cpn->cms_fix_time = (cpn->cms_fix_date - today) * YEARS_IN_DAY;

    /*	Forward Libor and coverage */
    err = interp_basis(cf_basis[i], &bas);
    if (err) {
      goto FREE_RETURN;
    }
    cpn->cvg = coverage(cpn->start_date, cpn->pay_date, bas);

    /* Equivalent libor calculation */
    libor =
        (1.0 / swp_f_df(max(cpn->start_date, today), cpn->pay_date, yc) - 1.0) /
        cpn->cvg;

    cpn->cvg *= exo_leg->notional;

    /*	Definition */
    cpn->alphabeta[0] = cf_alpha[i];
    cpn->alphabeta[1] = cf_beta[i];
    cpn->gamma = cf_gamma[i];
    cpn->floored = cf_floored[i];
    cpn->capped = cf_capped[i];
    cpn->floor = cf_floor[i];
    cpn->cap = cf_cap[i];

    cpn->use_SL = use_SL;
    cpn->use_cfoptions = use_cfoptions;
    if (use_cfoptions) {
      cpn->nopt = cf_nopt;
      for (iNumOpt = 0; iNumOpt < cf_nopt; iNumOpt++) {
        cpn->notopt[iNumOpt] = cf_notopt[i][iNumOpt];
        cpn->strikeopt[iNumOpt] = cf_strikeopt[i][iNumOpt];
        cpn->typeopt[iNumOpt] = cf_typeopt[i][iNumOpt];
      }
    } else {
      cpn->nopt = 0;
    }

    /*	Type */
    use_cms2 = 0;
    if (use_cfoptions) {
      cpn->floored = cpn->capped = 0;
    }

    if (cpn->floored && cpn->capped && fabs(cpn->cap - cpn->floor) < 1.0e-08) {
      cpn->type = 0;
      cpn->ncms = 0;
      cpn->alphabeta[0] = cpn->alphabeta[1] = 0.0;
      cpn->gamma = cpn->floor;
      cpn->floored = cpn->capped = 0;
    } else if (fabs(cpn->alphabeta[0]) < 1.0e-08 &&
               fabs(cpn->alphabeta[1]) < 1.0e-08) {
      cpn->type = 0;
      cpn->ncms = 0;
      cpn->alphabeta[0] = cpn->alphabeta[1] = 0.0;
      cpn->floored = cpn->capped = 0;
    } else if (fabs(cpn->alphabeta[0]) < 1.0e-08 ||
               fabs(cpn->alphabeta[1]) < 1.0e-08) {
      cpn->type = 2;
      cpn->ncms = 1;
      if (fabs(cpn->alphabeta[0]) < 1.0e-08) {
        cpn->alphabeta[0] = cpn->alphabeta[1];
        cpn->alphabeta[1] = 0.0;
        use_cms2 = 1;
      } else {
        cpn->alphabeta[1] = 0.0;
      }
    } else {
      cpn->type = 3;
      cpn->ncms = 2;
    }

    /*	CMS */
    for (l = 0; l < cpn->ncms; l++) {
      if (l == 0 && use_cms2 == 0) {
        cmstenor = cf_cms_tenor1[i];
        cmsfreq = cf_cms_freq1;
        cmsbasis = cf_cms_basis1;
        cmsspr = cf_cms_spread1[i];
        cmsbeta = cms_beta1;
        spread_slbeta = spread_slbeta1;
      } else {
        cmstenor = cf_cms_tenor2[i];
        cmsfreq = cf_cms_freq2;
        cmsbasis = cf_cms_basis2;
        cmsspr = cf_cms_spread2[i];
        cmsbeta = cms_beta2;
        spread_slbeta = spread_slbeta2;
      }

      /*	Interpret frequencies and basis */
      err = interp_basis(cmsbasis, &cmsbas);

      if (cmsbas == BASIS_ACT_360) {
        rconv = 365.0 / 360.0;
      } else {
        rconv = 1.0;
      }

      err = interp_compounding(cmsfreq, &cmscpnd);

      cms = &(cpn->cms[l]);

      cms->lambda = 0.0;
      cms->beta = cmsbeta;
      cms->start_date =
          add_unit(cpn->cms_fix_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
      cms->start_time = (cms->start_date - today) * YEARS_IN_DAY;

      /* Default value for shifted log model */
      cms->slbeta = interp(spread_vol_time, spread_slbeta, spread_vol_n,
                           cms->start_time, 0, &temp);

      cms->slshift = 0;
      cms->slvol = 0;

      err = add_tenor(cms->start_date, cmstenor, NO_BUSDAY_CONVENTION,
                      &(cms->end_date));
      if (err) {
        goto FREE_RETURN;
      }

      err = swp_f_initSwapDP(cms->start_date, cms->end_date, cmsfreq, cmsbasis,
                             &Swap);
      if (err) {
        goto FREE_RETURN;
      }

      SwapDates = SwapDP_to_DateList(&Swap, MODIFIED_SUCCEEDING);
      free_dates = 1;

      cms->num_cpn = SwapDates.len - 1;
      for (k = 0; k < cms->num_cpn; k++) {
        cms->cpn_pay_date[k] = SwapDates.date[k + 1];
        cms->cpn_pay_time[k] = (SwapDates.date[k + 1] - today) * YEARS_IN_DAY;
        cms->cvg[k] =
            coverage(SwapDates.date[k], SwapDates.date[k + 1], cmsbas);
      }

      cms->end_date = SwapDates.date[SwapDates.len - 1];
      cms->end_time = (cms->end_date - today) * YEARS_IN_DAY;
      cms->swap_spread = cmsspr;
      cms->cpnd = cmscpnd;
      cms->nfp =
          ((int)((cms->end_date - cms->start_date) * YEARS_IN_DAY * cms->cpnd +
                 0.1)) *
          1.0;
      cms->delay = cpn->pay_time - cms->start_time;
      cms->fix_time = cpn->cms_fix_time;

      swp_f_free_in_DateList(SwapDates);
      free_dates = 0;

      /*	Calculate forward and volatility */

      err = swp_f_ForwardRate(cms->start_date, cms->end_date, cmsfreq, cmsbasis,
                              yc, "CASH", &(cms->fwd));

      if (err) {
        goto FREE_RETURN;
      }

      cms->lib_spread = cms->fwd - libor;

      if (cpn->cms_fix_date > today) {
        if (cms_smile && cms_adj) {
          /* calculate the Fwd CMS */

          forward = cms->fwd + cms->swap_spread;

          err = swp_f_Cms_Rate(forward, cpn->cms_fix_time, cms->nfp, cms->cpnd,
                               cms->delay, rconv, vol_type, 0.0, 1,
                               cms->start_date, cms->end_date, cash_vol,
                               cms->swap_spread, vc, num_strikes_in_vol,
                               strikes_in_vol, &(cms->fwd_cms));

          if (err) {
            goto FREE_RETURN;
          }

          /* now price the ATM option */
          if (cms_vol_adj) {
            /* use the smile */
            err = swp_f_Cms_Option(
                forward, cpn->cms_fix_time, cms->nfp, cms->fwd_cms, cms->cpnd,
                SRT_RECEIVER, cms->delay, rconv, vol_type, 0.0, 1,
                cms->start_date, cms->end_date, cash_vol, cms->swap_spread, vc,
                num_strikes_in_vol, strikes_in_vol, &temp);

            err =
                srt_f_optimpvol(temp, cms->fwd_cms, cms->fwd_cms, cms->fix_time,
                                1.0, SRT_PUT, SRT_NORMAL, &(cms->atmvar));
          } else {
            /* disregard the smile */
            err = get_cash_vol(vc, cms->start_date, cms->end_date, cms->fwd, 0,
                               ref, &(cms->atmvar), &temp);

            if (err) {
              goto FREE_RETURN;
            }

            if (temp) {
              if (cms->fwd < 1e-10)
                cms->atmvar = 1e-6;
              else {
                temp = srt_f_optblksch(cms->fwd, cms->fwd, cms->atmvar,
                                       cms->fix_time, 1.0, SRT_CALL, PREMIUM);

                err =
                    srt_f_optimpvol(temp, cms->fwd, cms->fwd, cms->fix_time,
                                    1.0, SRT_CALL, SRT_NORMAL, &(cms->atmvar));
              }
            }
          }

          /* Case of Shifted Log model on CMS */
          if (use_SL) {
            if (!calib_SL) {
              /* Calibration of the shifted log parameter from the slbeta and
               * the atm norm vol */
              if (cms_vol_adj) {
                /* Calibration to the ATM cmsoption */
                cms->slshift = cms->fwd_cms * (1.0 / cms->slbeta - 1.0);
                cms->slvol = cms->atmvar * cms->slbeta / cms->fwd_cms;
              } else {
                /* Calibration to the ATM swaption */
                cms->slshift = forward * (1.0 / cms->slbeta - 1.0);
                cms->slvol = cms->atmvar * cms->slbeta / forward;
              }

            } else {
              /* Calibration of the shifted log model  */

              if (cms_vol_adj) {
                /* Calibration to cmsoptions */

                /* calibration of the slope */
                dStrike1 = cms->fwd_cms -
                           NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);
                dStrike2 = cms->fwd_cms +
                           NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);

                err = swp_f_Cms_Option(
                    forward, cpn->cms_fix_time, cms->nfp, dStrike1, cms->cpnd,
                    SRT_RECEIVER, cms->delay, rconv, vol_type, 0.0, 1,
                    cms->start_date, cms->end_date, cash_vol, cms->swap_spread,
                    vc, num_strikes_in_vol, strikes_in_vol, &temp);

                if (err)
                  goto FREE_RETURN;

                err =
                    srt_f_optimpvol(temp, cms->fwd_cms, dStrike1, cms->fix_time,
                                    1.0, SRT_PUT, SRT_LOGNORMAL, &dLogVol1);

                if (err)
                  goto FREE_RETURN;

                err = swp_f_Cms_Option(
                    forward, cpn->cms_fix_time, cms->nfp, dStrike2, cms->cpnd,
                    SRT_RECEIVER, cms->delay, rconv, vol_type, 0.0, 1,
                    cms->start_date, cms->end_date, cash_vol, cms->swap_spread,
                    vc, num_strikes_in_vol, strikes_in_vol, &temp);

                if (err)
                  goto FREE_RETURN;

                err =
                    srt_f_optimpvol(temp, cms->fwd_cms, dStrike2, cms->fix_time,
                                    1.0, SRT_PUT, SRT_LOGNORMAL, &(dLogVol2));

                if (err)
                  goto FREE_RETURN;

                /* Find the shifted log params */
                err = find_shifted_log_params(
                    cms->fix_time, cms->fwd_cms, dStrike1, dLogVol1, dStrike2,
                    dLogVol2, &cms->slshift, &cms->slvol, &temp);

                if (err)
                  goto FREE_RETURN;

                /* Recalibrate the ATM vol */
                temp = srt_f_optblknrm(cms->fwd_cms, cms->fwd_cms, cms->atmvar,
                                       cms->fix_time, 1.0, SRT_CALL, PREMIUM);

                if (cms->fwd_cms + cms->slshift > 0) {
                  err = srt_f_optimpvol(temp, cms->fwd_cms + cms->slshift,
                                        cms->fwd_cms + cms->slshift,
                                        cms->fix_time, 1.0, SRT_CALL,
                                        SRT_LOGNORMAL, &cms->slvol);

                  if (err)
                    goto FREE_RETURN;
                } else {
                  err = srt_f_optimpvol(temp, -(cms->fwd_cms + cms->slshift),
                                        -(cms->fwd_cms + cms->slshift),
                                        cms->fix_time, 1.0, SRT_CALL,
                                        SRT_LOGNORMAL, &cms->slvol);

                  if (err)
                    goto FREE_RETURN;
                }
              } else {
                /* Calibration to swaptions */

                /* calibration of the slope */
                dStrike1 =
                    forward - NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);
                dStrike2 =
                    forward + NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);

                /* Option 1 */
                err = get_cash_vol(vc, cms->start_date, cms->end_date,
                                   dStrike1 - cms->swap_spread, 0, ref,
                                   &(dLogVol1), &temp);

                if (err) {
                  goto FREE_RETURN;
                }

                /* Transform into LOG VOL */
                if (temp < 0.5) {
                  temp = srt_f_optblknrm(cms->fwd, dStrike1 - cms->swap_spread,
                                         dLogVol1, cms->fix_time, 1.0, SRT_CALL,
                                         PREMIUM);

                  err = srt_f_optimpvol(
                      temp, cms->fwd, dStrike1 - cms->swap_spread,
                      cms->fix_time, 1.0, SRT_CALL, SRT_LOGNORMAL, &(dLogVol1));
                }

                /* Option 2 */
                err = get_cash_vol(vc, cms->start_date, cms->end_date,
                                   dStrike2 - cms->swap_spread, 0, ref,
                                   &(dLogVol2), &temp);

                if (err) {
                  goto FREE_RETURN;
                }

                /* Transform into LOG VOL */
                if (temp < 0.5) {
                  temp = srt_f_optblknrm(cms->fwd, dStrike2 - cms->swap_spread,
                                         dLogVol2, cms->fix_time, 1.0, SRT_CALL,
                                         PREMIUM);

                  err = srt_f_optimpvol(
                      temp, cms->fwd, dStrike2 - cms->swap_spread,
                      cms->fix_time, 1.0, SRT_CALL, SRT_LOGNORMAL, &(dLogVol2));
                }

                /* Find the shifted log params */
                err = find_shifted_log_params(
                    cms->fix_time, forward, dStrike1, dLogVol1, dStrike2,
                    dLogVol2, &cms->slshift, &cms->slvol, &temp);

                if (err)
                  goto FREE_RETURN;

                /* Recalibrate the ATM vol */
                temp = srt_f_optblknrm(forward, forward, cms->atmvar,
                                       cms->fix_time, 1.0, SRT_CALL, PREMIUM);

                if (forward + cms->slshift > 0) {
                  err = srt_f_optimpvol(
                      temp, forward + cms->slshift, forward + cms->slshift,
                      cms->fix_time, 1.0, SRT_CALL, SRT_LOGNORMAL, &cms->slvol);

                  if (err)
                    goto FREE_RETURN;
                } else {
                  err = srt_f_optimpvol(temp, -(forward + cms->slshift),
                                        -(forward + cms->slshift),
                                        cms->fix_time, 1.0, SRT_CALL,
                                        SRT_LOGNORMAL, &cms->slvol);

                  if (err)
                    goto FREE_RETURN;
                }
              }
            }

            /* Calculation of the lambda for adjusting the fwd cms */
            cms->lambda = (cms->fwd_cms - forward) / (forward + cms->slshift) /
                          (forward + cms->slshift) / cpn->cms_fix_time;

          } else {
            /* Case no SL selected */

            /* Calculation of the lambda for adjusting the fwd cms */
            if (cms->beta > 1.0e-08) {
              cms->lambda = (cms->fwd_cms - forward) /
                            pow(fabs(forward), 2.0 * cms->beta) /
                            cpn->cms_fix_time;
            } else {
              cms->lambda = (cms->fwd_cms - forward) / cpn->cms_fix_time;
            }
          }

        } else {
          /* No smile adjustment */
          err = get_cash_vol(vc, cms->start_date, cms->end_date, cms->fwd, 0,
                             ref, &(cms->atmvar), &temp);

          if (err) {
            goto FREE_RETURN;
          }
          if (temp) {
            if (cms->fwd < 1e-10)
              cms->atmvar = 1e-6;
            else {
              temp = srt_f_optblksch(cms->fwd, cms->fwd, cms->atmvar,
                                     cms->fix_time, 1.0, SRT_CALL, PREMIUM);

              err = srt_f_optimpvol(temp, cms->fwd, cms->fwd, cms->fix_time,
                                    1.0, SRT_CALL, SRT_NORMAL, &(cms->atmvar));
            }
          }

          /* Case of Shifted Log model on CMS */
          if (use_SL) {
            /* calculate the Fwd CMS */
            forward = cms->fwd + cms->swap_spread;

            if (!calib_SL) {
              /* Calibration of the shifted log parameter from the slbeta and
               * the atm norm vol */
              /* Calibration to the ATM swaption */
              cms->slshift = forward * (1.0 / cms->slbeta - 1.0);
              cms->slvol = cms->atmvar * cms->slbeta / forward;
            } else {
              /* Calibration to swaptions */

              /* calibration of the slope */
              dStrike1 =
                  forward - NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);
              dStrike2 =
                  forward + NbStdcalib_SL * cms->atmvar * sqrt(cms->fix_time);

              /* Option 1 */
              err = get_cash_vol(vc, cms->start_date, cms->end_date,
                                 dStrike1 - cms->swap_spread, 0, ref,
                                 &(dLogVol1), &temp);

              if (err) {
                goto FREE_RETURN;
              }

              /* Transform into LOG VOL */
              if (temp < 0.5) {
                temp = srt_f_optblknrm(cms->fwd, dStrike1 - cms->swap_spread,
                                       dLogVol1, cms->fix_time, 1.0, SRT_CALL,
                                       PREMIUM);

                err = srt_f_optimpvol(
                    temp, cms->fwd, dStrike1 - cms->swap_spread, cms->fix_time,
                    1.0, SRT_CALL, SRT_LOGNORMAL, &(dLogVol1));
              }

              /* Option 2 */
              err = get_cash_vol(vc, cms->start_date, cms->end_date,
                                 dStrike2 - cms->swap_spread, 0, ref,
                                 &(dLogVol2), &temp);

              if (err) {
                goto FREE_RETURN;
              }

              /* Transform into LOG VOL */
              if (temp < 0.5) {
                temp = srt_f_optblknrm(cms->fwd, dStrike2 - cms->swap_spread,
                                       dLogVol2, cms->fix_time, 1.0, SRT_CALL,
                                       PREMIUM);

                err = srt_f_optimpvol(
                    temp, cms->fwd, dStrike2 - cms->swap_spread, cms->fix_time,
                    1.0, SRT_CALL, SRT_LOGNORMAL, &(dLogVol2));
              }

              /* Find the shifted log params */
              err = find_shifted_log_params(cms->fix_time, forward, dStrike1,
                                            dLogVol1, dStrike2, dLogVol2,
                                            &cms->slshift, &cms->slvol, &temp);

              if (err)
                goto FREE_RETURN;

              /* Recalibrate the ATM vol */
              temp = srt_f_optblknrm(forward, forward, cms->atmvar,
                                     cms->fix_time, 1.0, SRT_CALL, PREMIUM);

              if (forward + cms->slshift > 0) {
                err = srt_f_optimpvol(
                    temp, forward + cms->slshift, forward + cms->slshift,
                    cms->fix_time, 1.0, SRT_CALL, SRT_LOGNORMAL, &cms->slvol);

                if (err)
                  goto FREE_RETURN;
              } else {
                err = srt_f_optimpvol(
                    temp, -(forward + cms->slshift), -(forward + cms->slshift),
                    cms->fix_time, 1.0, SRT_CALL, SRT_LOGNORMAL, &cms->slvol);

                if (err)
                  goto FREE_RETURN;
              }
            }
          }
        }

        cms->atmvar *= cms->atmvar;
      } else {
        cms->atmvar = 0.0;
      }

      cms->floorvar = cms->capvar = cms->atmvar;

      /* Fill the option normal variance with the ATM normal vol by default for
       * string of options */
      if (use_cfoptions) {
        for (iNumOpt = 0; iNumOpt < cpn->nopt; iNumOpt++) {
          cms->optvar[iNumOpt] = cms->atmvar;
        }
      }

      if (cms_adj) {
        if (cms_smile) {
          cms->adj_type = 2;
        } else {
          cms->adj_type = 1;
        }
      } else {
        cms->adj_type = 0;
      }
    }

    /* end for on l */

    if (cpn->type == 2) {
      cms = &(cpn->cms[0]);
      if (fabs(cms->nfp - 1.0) < 1.0e-08 &&
          cms->end_time < cpn->pay_time + ONE_MONTH) {
        cpn->type = 1;
        cms->adj_type = 0;
      }

      if (cpn->cms_fix_date > today) {
        if (cpn->floored && !use_cfoptions) {
          strike =
              (cpn->floor - cpn->gamma) / cpn->alphabeta[0] - cms->swap_spread;

          if (cms_smile && cms_adj && cms_vol_adj) {
            if (strike + cms->swap_spread < 1.0e-08) {
              err = "Error in CIF/CMS definition: negative strike";
              goto FREE_RETURN;
            }

            /* now price the option */
            err = swp_f_Cms_Option(
                cms->fwd + cms->swap_spread, cpn->cms_fix_time, cms->nfp,
                strike + cms->swap_spread, cms->cpnd, SRT_RECEIVER, cms->delay,
                rconv, vol_type, 0.0, 1, cms->start_date, cms->end_date,
                cash_vol, cms->swap_spread, vc, num_strikes_in_vol,
                strikes_in_vol, &temp);

            err = srt_f_optimpvol(temp, cms->fwd_cms, strike + cms->swap_spread,
                                  cms->fix_time, 1.0, SRT_PUT, SRT_NORMAL,
                                  &(cms->floorvar));

          } else {
            if (strike < 1.0e-08) {
              err = "Error in CIF/CMS definition: negative embdedded strike";
              goto FREE_RETURN;
            }

            err = get_cash_vol(vc, cms->start_date, cms->end_date, strike, 0,
                               ref, &(cms->floorvar), &temp);
            if (err) {
              goto FREE_RETURN;
            }
            if (temp) {
              if (cms->fwd < 1e-10)
                cms->floorvar = 1e-6;
              else {
                temp = srt_f_optblksch(cms->fwd, strike, cms->floorvar,
                                       cms->fix_time, 1.0, SRT_CALL, PREMIUM);
                err =
                    srt_f_optimpvol(temp, cms->fwd, strike, cms->fix_time, 1.0,
                                    SRT_CALL, SRT_NORMAL, &(cms->floorvar));
              }
            }
          }

          cms->floorvar *= cms->floorvar;
        }

        if (cpn->capped && !use_cfoptions) {
          strike =
              (cpn->cap - cpn->gamma) / cpn->alphabeta[0] - cms->swap_spread;

          if (cms_smile && cms_adj && cms_vol_adj) {
            if (strike + cms->swap_spread < 1.0e-08) {
              err = "Error in CIF/CMS definition: negative strike";
              goto FREE_RETURN;
            }

            /* now price the option */
            err = swp_f_Cms_Option(
                cms->fwd + cms->swap_spread, cpn->cms_fix_time, cms->nfp,
                strike + cms->swap_spread, cms->cpnd, SRT_RECEIVER, cms->delay,
                rconv, vol_type, 0.0, 1, cms->start_date, cms->end_date,
                cash_vol, cms->swap_spread, vc, num_strikes_in_vol,
                strikes_in_vol, &temp);

            err = srt_f_optimpvol(temp, cms->fwd_cms, strike + cms->swap_spread,
                                  cms->fix_time, 1.0, SRT_PUT, SRT_NORMAL,
                                  &(cms->capvar));

          } else {
            if (strike < 1.0e-08) {
              err = "Error in CIF/CMS definition: negative embdedded strike";
              goto FREE_RETURN;
            }

            err = get_cash_vol(vc, cms->start_date, cms->end_date, strike, 0,
                               ref, &(cms->capvar), &temp);
            if (err) {
              goto FREE_RETURN;
            }
            if (temp) {
              if (cms->fwd < 1e-10)
                cms->capvar = 1e-6;
              else {
                temp = srt_f_optblksch(cms->fwd, strike, cms->capvar,
                                       cms->fix_time, 1.0, SRT_CALL, PREMIUM);
                err =
                    srt_f_optimpvol(temp, cms->fwd, strike, cms->fix_time, 1.0,
                                    SRT_CALL, SRT_NORMAL, &(cms->capvar));
              }
            }
          }

          cms->capvar *= cms->capvar;
        }

        /* string of options case */
        if (use_cfoptions) {
          for (iNumOpt = 0; iNumOpt < cpn->nopt; iNumOpt++) {
            /* Cash strike */
            strike =
                cpn->strikeopt[iNumOpt] / cpn->alphabeta[0] - cms->swap_spread;

            if (cms_smile && cms_adj && cms_vol_adj) {
              if (strike + cms->swap_spread < 1.0e-08) {
                err = "Error in CIF/CMS definition: negative strike";
                goto FREE_RETURN;
              }

              /* now price the option */
              err = swp_f_Cms_Option(
                  cms->fwd + cms->swap_spread, cpn->cms_fix_time, cms->nfp,
                  strike + cms->swap_spread, cms->cpnd, cpn->typeopt[iNumOpt],
                  cms->delay, rconv, vol_type, 0.0, 1, cms->start_date,
                  cms->end_date, cash_vol, cms->swap_spread, vc,
                  num_strikes_in_vol, strikes_in_vol, &temp);

              err =
                  srt_f_optimpvol(temp, cms->fwd_cms, strike + cms->swap_spread,
                                  cms->fix_time, 1.0, cpn->typeopt[iNumOpt],
                                  SRT_NORMAL, &(cms->optvar[iNumOpt]));

            } else {
              if (strike < 1.0e-08) {
                err = "Error in CIF/CMS definition: negative embdedded strike";
                goto FREE_RETURN;
              }

              err = get_cash_vol(vc, cms->start_date, cms->end_date, strike, 0,
                                 ref, &(cms->optvar[iNumOpt]), &temp);
              if (err) {
                goto FREE_RETURN;
              }

              /* Convert to normal vol if getcash vol return logvol */
              if (temp) {
                if (cms->fwd < 1e-10)
                  cms->optvar[iNumOpt] = 1e-6;
                else {
                  temp = srt_f_optblksch(cms->fwd, strike, cms->optvar[iNumOpt],
                                         cms->fix_time, 1.0, SRT_CALL, PREMIUM);
                  err = srt_f_optimpvol(temp, cms->fwd, strike, cms->fix_time,
                                        1.0, SRT_CALL, SRT_NORMAL,
                                        &(cms->optvar[iNumOpt]));
                }
              }
            }

            cms->optvar[iNumOpt] *= cms->optvar[iNumOpt];
          }
        }
      }
    }

    i++;
    j++;
  }

  exo_leg->type = 0;
  for (j = 0; j < exo_leg->num_cpn; j++) {
    cpn = exo_leg->cpn + j;
    if (exo_leg->type < cpn->type) {
      exo_leg->type = cpn->type;
    }
  }

  err = ccf_check_exo_leg(exo_leg);

FREE_RETURN:

  if (free_dates) {
    swp_f_free_in_DateList(SwapDates);
  }

  if (err) {
    ccf_free_exo_leg(exo_leg);
  }

  return err;
}

/*	Check dates consistency */
Err ccf_check_exo_leg(CF_EXO_LEG exo_leg) {
  int i, j, k;
  CF_EXO_CPN cpn;
  CF_CMS_DESC cms;

  /*	Notional has to be different from 0 */
  if (exo_leg->notional <= 0.0) {
    return "Domestic coupon notional has to greater than 0";
  }

  /*	Check that start  , pay  , fix and val dates are increasing */
  for (i = 1; i < exo_leg->num_cpn; i++) {
    if (exo_leg->cpn[i].start_date < exo_leg->cpn[i - 1].start_date) {
      return "Start dates should be increasing in exotic leg";
    }

    if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i - 1].pay_date) {
      return "Pay dates should be increasing in exotic leg";
    }

    if (exo_leg->cpn[i].cms_fix_date < exo_leg->cpn[i - 1].cms_fix_date) {
      return "Fixing dates should be increasing in exotic leg";
    }
  }

  /*	Check that pay dates are after start dates and fix dates */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i].start_date) {
      return "Pay dates should be after start dates in exotic leg";
    }

    if (exo_leg->cpn[i].pay_date < exo_leg->cpn[i].cms_fix_date) {
      return "Pay dates should be after fixing dates in exotic leg";
    }
  }
  /*	Check that floor is less than cap */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    if (exo_leg->cpn[i].floored && exo_leg->cpn[i].capped &&
        exo_leg->cpn[i].floor > exo_leg->cpn[i].cap) {
      return "Floor should be lower than cap";
    }
  }

  /*	CMS reference swaps */

  /*	Check that start dates are before end dates */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    cpn = exo_leg->cpn + i;
    for (k = 0; k < cpn->ncms; k++) {
      cms = &(cpn->cms[k]);
      if (cms->end_date < cms->start_date) {
        return "End dates should be after start dates in CMS reference swaps";
      }
    }
  }

  /*	Check that coupon pay dates are increasing */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    cpn = exo_leg->cpn + i;
    for (k = 0; k < cpn->ncms; k++) {
      cms = &(cpn->cms[k]);

      for (j = 1; j < cms->num_cpn; j++) {
        if (cms->cpn_pay_date[j] < cms->cpn_pay_date[j - 1]) {
          return "Pay dates should be increasing in CMS reference swaps";
        }
      }
    }
  }

  /*	Check that pay dates are between start and end dates */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    cpn = exo_leg->cpn + i;
    for (k = 0; k < cpn->ncms; k++) {
      cms = &(cpn->cms[k]);

      if (cms->end_date < cms->cpn_pay_date[cms->num_cpn - 1] ||
          cms->start_date > cms->cpn_pay_date[0]) {
        return "Pay dates should be between start and end dates in CMS "
               "reference swaps";
      }
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err ccf_free_exo_leg(CF_EXO_LEG exo_leg) {
  if (exo_leg->cpn) {
    free(exo_leg->cpn);
    exo_leg->cpn = NULL;
  }

  return NULL;
}

/*	Functions for the calls */

Err ccf_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag,           /*	0: I  , 1: E */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee, CCF_STR ccf) {
  int i, j, k;

  CF_FUND_LEG fund_leg;
  CF_EXO_LEG exo_leg;
  CF_CALL call;
  Err err = NULL;

  /*	Initialise pointers */
  ccf->call = NULL;
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  /*	Skip calls to be exercised before today */
  i = 0;
  while (i < ncall && ex_date[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one call is left */
  if (i == ncall) {
    /* err = "All calls are to be exercised before today in ccf_fill_calls"; */
    ccf->num_calls = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  ccf->num_calls = ncall - i;
  ccf->call = (cf_call *)calloc(ccf->num_calls, sizeof(cf_call));
  if (!ccf->num_calls) {
    err = "Allocation error in ccf_fill_calls";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;

  while (i < ncall) {
    call = ccf->call + j;

    /*	Dates */
    call->ex_date = ex_date[i];
    call->set_date = set_date[i];

    /*	Times */
    call->ex_time = (ccf->call[j].ex_date - today) * YEARS_IN_DAY;
    call->set_time = (ccf->call[j].set_date - today) * YEARS_IN_DAY;

    /*	Call on funding leg */
    /*	k = index of the first coupon to be called on funding leg  ,
                    i.e. first coupon with a start date >= ex date */

    k = 0;
    while (k < fund_leg->num_cpn &&
           fund_leg->cpn[k].start_date < call->ex_date) {
      k++;
    }
    if (k == fund_leg->num_cpn) {
      err = serror("Call number %d does not control any coupon in funding leg",
                   i);
      goto FREE_RETURN;
    }
    call->fund_idx = k;
    call->num_fund_cpn = fund_leg->num_cpn - k;

    /*	Call on exotic leg */
    /*	k = index of the first coupon to be called on exotic leg  ,
                    i.e. first coupon with a start date >= ex date */
    k = 0;
    while (k < exo_leg->num_cpn && exo_leg->cpn[k].start_date < call->ex_date) {
      k++;
    }
    if (k == exo_leg->num_cpn) {
      err =
          serror("Call number %d does not control any coupon in exotic leg", i);
      goto FREE_RETURN;
    }
    call->cf_idx = k;
    call->num_cf_cpn = exo_leg->num_cpn - k;

    /*	Payer or receiver */
    call->pay_rec = pay_rec;

    /*	Fee */
    call->fee = fee[i];

    i++;
    j++;
  }

  err = ccf_check_calls(ccf);

FREE_RETURN:

  if (err) {
    ccf_free_calls(ccf);
  }

  return err;
}

void SetCheckCall(int nCall) { nCheckCall = nCall; }

/*	Check dates consistency */
Err ccf_check_calls(CCF_STR ccf) {

  if (nCheckCall) {
    int i;

    /*	Check that ex and set dates are strictly increasing
                    Also check that funding and pd indices are strictly
       increasing  , i.e. there is no redundant calls */
    for (i = 1; i < ccf->num_calls; i++) {
      if (ccf->call[i].ex_date <= ccf->call[i - 1].ex_date) {
        return "Exercise dates should be increasing";
      }

      if (ccf->call[i].set_date <= ccf->call[i - 1].set_date) {
        return "Settlement dates should be increasing";
      }

      if (ccf->call[i].fund_idx < ccf->call[i - 1].fund_idx) {
        return "Number of funding coupons controlled by calls should be "
               "decreasing";
      }

      if (ccf->call[i].cf_idx < ccf->call[i - 1].cf_idx) {
        return "Number of exotic coupons controlled by calls should be "
               "decreasing";
      }

      if (ccf->call[i].fund_idx <= ccf->call[i - 1].fund_idx &&
          ccf->call[i].cf_idx <= ccf->call[i - 1].cf_idx) {
        return serror("Calls %d and %d -indexed after today- are redundant",
                      i - 1, i);
      }
    }

    /*	Check that set dates are after ex dates
                    Also check that the call date is before the start  , end and
       fixing dates of the coupons it controls */
    for (i = 0; i < ccf->num_calls; i++) {
      if (ccf->call[i].set_date < ccf->call[i].ex_date) {
        return "Settlement dates should be after exercise dates";
      }

      if (ccf->fund_leg->cpn[ccf->call[i].fund_idx].start_date <
              ccf->call[i].ex_date ||
          ccf->fund_leg->cpn[ccf->call[i].fund_idx].pay_date <
              ccf->call[i].ex_date) {
        return "A funding coupon starts before its exercise date";
      }

      if (ccf->cf_leg->cpn[ccf->call[i].cf_idx].start_date <
              ccf->call[i].ex_date ||
          ccf->cf_leg->cpn[ccf->call[i].cf_idx].pay_date <
              ccf->call[i].ex_date ||
          ccf->cf_leg->cpn[ccf->call[i].cf_idx].cms_fix_date <
              ccf->call[i].ex_date) {
        return "An exotic coupon starts or fixes before its exercise date";
      }
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err ccf_free_calls(CCF_STR ccf) {
  if (ccf->call) {
    free(ccf->call);
    ccf->call = NULL;
  }

  return NULL;
}

/*	Functions for the underlying */

/*	Fill underlying structure from a predefined underlying */
Err ccf_fill_und(char *lgm2dund, char *vc, char *ref, char *swap_freq,
                 char *swap_basis, long spread_vol_n, double *spread_vol_time,
                 double *spread_vol_floor, double *spread_vol_cap, int cvg_sv,
                 int is_corr, CCF_UND und) {
  double *sig_time = NULL, *sig = NULL;

  SrtUndPtr srtund;

  int i;
  Err err = NULL;

  /*	Initialise */
  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;
  und->lambda_date = NULL;
  und->lambda_time = NULL;
  und->lambda = NULL;
  und->spread_vol_time = NULL;
  und->spread_vol_floor = NULL;
  und->spread_vol_cap = NULL;

  und->cvg_sv = cvg_sv;
  und->is_corr = is_corr;

  und->spread_vol_n = spread_vol_n;
  if (spread_vol_n > 0) {
    und->spread_vol_time = (double *)calloc(spread_vol_n, sizeof(double));
    und->spread_vol_floor = (double *)calloc(spread_vol_n, sizeof(double));
    und->spread_vol_cap = (double *)calloc(spread_vol_n, sizeof(double));

    if (!und->spread_vol_time || !und->spread_vol_floor ||
        !und->spread_vol_cap) {
      err = "Allocation error (1) in ccf_fill_und";
      goto FREE_RETURN;
    }

    memcpy(und->spread_vol_time, spread_vol_time,
           spread_vol_n * sizeof(double));
    memcpy(und->spread_vol_floor, spread_vol_floor,
           spread_vol_n * sizeof(double));
    memcpy(und->spread_vol_cap, spread_vol_cap, spread_vol_n * sizeof(double));
  }

  strcpy(und->name, lgm2dund);

  /*	Get term structures */
  err = Get_LGM2F_TermStructure2(lgm2dund, &(und->sigma), &(und->sigma_time),
                                 &(und->sigma_n), &(und->lambda),
                                 &(und->lambda_time), &(und->lambda_n),
                                 &(und->alpha), &(und->gamma), &(und->rho));
  if (err) {
    goto FREE_RETURN;
  }

  /*	Fill dates */
  srtund = lookup_und(lgm2dund);
  und->today = get_today_from_underlying(srtund);
  strcpy(und->yc, get_ycname_from_irund(srtund));
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));
  und->lambda_date = (double *)calloc(und->lambda_n, sizeof(double));

  if (!und->sigma_date || !und->lambda_date) {
    err = "Allocation error in ccf_fill_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] =
        und->today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

  for (i = 0; i < und->lambda_n; i++) {
    und->lambda_date[i] =
        und->today + und->lambda_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));

FREE_RETURN:

  if (err) {
    ccf_free_und(und);
  }

  return err;
}

/*	Find calibration strikes */
static Err find_strikes(
    /*	Market parameters */
    long today, long end_date, char *yc, /*	yc */
    char *vc,                            /*	vc (only if calib) */
    char *ref,                           /*	ref rate (only if calib) */
    char *swap_freq,                     /*	swap freq (only if calib) */
    char *swap_basis,                    /*	swap basis (only if calib) */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	The structure */
    CCF_STR ccf, /*	structure */
    /*	The results */
    double *short_strikes, double *long_strikes, int *strike_type) {
  Err err = NULL;
  int i, j, k;
  int cpnn, fundn;
  double *cpnki = NULL, *bi = NULL, *fundki = NULL;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  CF_CALL call;
  double pvk1, pvk2, pvlib, pvfee;
  double dfra, dswap;
  double cvg;
  double fl_lvl, fra_lvl;
  double fwd, vol, lvl, std, temp;

  fund_leg = ccf->fund_leg;
  fundn = fund_leg->num_cpn;
  fundki = (double *)calloc(fundn, sizeof(double));

  exo_leg = ccf->cf_leg;
  cpnn = exo_leg->num_cpn;
  cpnki = (double *)calloc(cpnn, sizeof(double));
  bi = (double *)calloc(cpnn, sizeof(double));

  if (!fundki || !bi || !cpnki) {
    err = "Allocation error in find_strikes";
    goto FREE_RETURN;
  }

  /*	1.)	Convert both legs into coupons of the form ki + bi * Libor */

  /*	11.) Fund leg */
  /*	-- includes notional but excludes coverage -- */
  for (i = 0; i < fund_leg->num_cpn; i++) {
    fund_cpn = fund_leg->cpn + i;
    fundki[i] = -fund_cpn->cpn / fund_cpn->cvg;
  }

  /*	11.) Exotic leg */
  /*	-- excludes notional and coverage -- */
  for (i = 0; i < exo_leg->num_cpn; i++) {
    exo_cpn = exo_leg->cpn + i;
    cpnki[i] = exo_cpn->gamma;
    bi[i] = fund_leg->notional;
    for (j = 0; j < exo_cpn->ncms; j++) {
      cms = &(exo_cpn->cms[j]);
      cpnki[i] += exo_cpn->alphabeta[j] * (cms->lib_spread + cms->swap_spread);
      bi[i] -= exo_cpn->alphabeta[j] * exo_leg->notional;
    }
  }

  /*	2.) Loop on exercise dates */
  for (i = 0; i < ccf->num_calls; i++) {
    call = ccf->call + i;
    /*	21.) Compute the FRA move necessary to put the structure at par */
    pvk1 = 0.0;
    for (j = 0; j < call->num_cf_cpn; j++) {
      k = call->cf_idx + j;
      exo_cpn = exo_leg->cpn + k;
      pvk1 += cpnki[k] * exo_cpn->cvg * swp_f_df(today, exo_cpn->pay_date, yc);
      /*	cvg includes notional */
    }
    pvk2 = 0.0;
    for (j = 0; j < call->num_fund_cpn; j++) {
      k = call->fund_idx + j;
      fund_cpn = fund_leg->cpn + k;
      pvk2 +=
          fundki[k] * fund_cpn->cvg * swp_f_df(today, fund_cpn->pay_date, yc);
      /*	Coupons are negative and include notional */
    }
    pvlib = 0.0;
    fl_lvl = 0.0;
    fra_lvl = 0.0;
    for (j = 0; j < call->num_cf_cpn; j++) {
      k = call->cf_idx + j;
      exo_cpn = exo_leg->cpn + k;
      cvg = exo_cpn->cvg / exo_leg->notional;
      pvlib += bi[k] * (swp_f_df(today, exo_cpn->start_date, yc) -
                        swp_f_df(today, exo_cpn->pay_date, yc));
      fl_lvl += bi[k] * cvg * swp_f_df(today, exo_cpn->pay_date, yc);
      fra_lvl += cvg * swp_f_df(today, exo_cpn->pay_date, yc);
      /*	cvg includes notional */
    }
    if (call->pay_rec == 0) {
      pvfee = -call->fee * swp_f_df(today, call->set_date, yc);
    } else {
      pvfee = call->fee * swp_f_df(today, call->set_date, yc);
    }
    if (fabs(fl_lvl) > 1.0e-08) {
      dfra = (pvk1 + pvk2 + pvfee - pvlib) / fl_lvl;
      /*	Special case  , no sensitivity to level of rates  , do atm */
    } else {
      dfra = 0.0;
    }

    /* 22.) Calculate forward swap  , level and normal vol */
    err = swp_f_ForwardRate(
        add_unit(call->ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING), end_date,
        swap_freq, swap_basis, yc, "CASH", &fwd);
    if (err) {
      goto FREE_RETURN;
    }
    err = swp_f_LevelPayment(
        add_unit(call->ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING), end_date,
        swap_freq, swap_basis, yc, "CASH", &lvl);
    if (err) {
      goto FREE_RETURN;
    }
    err = get_cash_vol(
        vc, add_unit(call->ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING), end_date,
        fwd, 0, ref, &vol, &temp);
    if (err) {
      goto FREE_RETURN;
    }
    if (temp) {
      if (call->ex_time > 1.0E-10) {
        if (fwd < 1e-10)
          vol = 1e-6;
        else {
          temp = srt_f_optblksch(fwd, fwd, vol, call->ex_time, 1.0, SRT_CALL,
                                 PREMIUM);
          err = srt_f_optimpvol(temp, fwd, fwd, call->ex_time, 1.0, SRT_CALL,
                                SRT_NORMAL, &vol);
        }
      } else {
        vol *= fwd;
      }

      if (err) {
        goto FREE_RETURN;
      }
    }
    std = vol * sqrt(call->ex_time);

    /* 23.) Compute the equivalent swap move */
    dswap = dfra * fra_lvl / lvl;

    /* 24.) Compute the equivalent number of standard deviations */
    if (std > 1.0e-04) {
      short_strikes[i] = long_strikes[i] = dswap / std;
    } else {
      short_strikes[i] = long_strikes[i] = 0.0;
    }
  }

  *strike_type = 3;

FREE_RETURN:

  if (fundki)
    free(fundki);
  if (bi)
    free(bi);
  if (cpnki)
    free(cpnki);

  return err;
}

/*	Fill underlying structure from calibration instruments */
Err ccf_calib_und(
    long today,
    /*	EOD Flag */
    int eod_flag,      /*	0: I  , 1: E */
    char *yc,          /*	yc */
    char *vc,          /*	vc (only if calib) */
    char *ref,         /*	ref rate (only if calib) */
    char *swap_freq,   /*	swap freq (only if calib) */
    char *swap_basis,  /*	swap basis (only if calib) */
    int lam_ts,        /*	0: use unique lambda  , 1: use ts */
    double lambda,     /*	lambda if unique */
    int tsnlam,        /*	number of lambdas if ts */
    double *tslamtime, /*	lambda times i.e. (date - today) / 365 if ts */
    double *tslam,     /*	corresponding lambdas if ts */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	Calib params */
    int force_atm, /*	force atm calib */
    double max_std_long, double max_std_short,
    int fix_lambda,          /*	0: calib lambda to cap  , 1: fix lambda calib
                                                             to diagonal */
    int cal_vol_shift_type,  /*	vol shift type for volatility */
    double cal_vol_shift,    /*	vol shift */
    double cal_lambda_shift, /*	shift on lambda after calibration */
    int one_f_equi,          /*	1F equivalent flag:
                                                             if set to 1  , then 2F
                                lambda will calibrate          to the cap priced within calibrated
                                1F          with the given lambda */
    int skip_last,     /*	If 1  , the last option is disregarded and the forward
                          volatility is flat from option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int proba_weight, int use_exe_bound, double *proba, double *exe_bound,
    /*	End of calib params */
    CCF_STR ccf,        /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    CCF_UND und, long spread_vol_n, double *spread_vol_time,
    double *spread_vol_floor, double *spread_vol_cap, int cvg_sv, int is_corr) {
  int i, nex;
  long *lex = NULL;

  long end_struct, end_swap;
  double *long_strikes = NULL, *short_strikes = NULL;

  double *copy_tslam = NULL;

  int strike_type;

  CF_FUND_LEG fund_leg;
  CF_EXO_LEG exo_leg;

  Err err = NULL;

  /*	Check ts of lambda if needed */
  if (lam_ts) {
    for (i = 0; i < tsnlam; i++) {
      if (tslamtime[i] <= 0.0) {
        err = "Past lambda times are forbidden - please check your tau term "
              "struct";
        goto FREE_RETURN;
      }
    }

    copy_tslam = calloc(tsnlam, sizeof(double));

    if (!copy_tslam) {
      err = "Memory allocation faillure in ccf_calib_und";
      goto FREE_RETURN;
    }

    memcpy(copy_tslam, tslam, tsnlam * sizeof(double));
  }

  /*	Eliminate zero lambdas */
  if (lam_ts) {
    for (i = 0; i < tsnlam; i++) {
      if (fabs(copy_tslam[i]) < 1.0e-08) {
        copy_tslam[i] = 1.0e-08;
      }
    }
  } else {
    if (fabs(lambda) < 1.0e-08) {
      lambda = 1.0e-08;
    }
  }

  /*	Initialise */
  exo_leg = ccf->cf_leg;
  fund_leg = ccf->fund_leg;

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;
  und->lambda_date = NULL;
  und->lambda_time = NULL;
  und->lambda = NULL;
  und->spread_vol_time = NULL;
  und->spread_vol_floor = NULL;
  und->spread_vol_cap = NULL;
  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));
  und->alpha = alpha;
  und->gamma = gamma;
  und->rho = rho;
  und->today = today;
  und->cvg_sv = cvg_sv;
  und->is_corr = is_corr;
  strcpy(und->yc, yc);
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);
  strcpy(und->name, "CALIB");

  /*	Spread vol */
  und->spread_vol_n = spread_vol_n;
  if (spread_vol_n > 0) {
    und->spread_vol_time = (double *)calloc(spread_vol_n, sizeof(double));
    und->spread_vol_floor = (double *)calloc(spread_vol_n, sizeof(double));
    und->spread_vol_cap = (double *)calloc(spread_vol_n, sizeof(double));

    if (!und->spread_vol_time || !und->spread_vol_floor ||
        !und->spread_vol_cap) {
      err = "Allocation error (1) in ccf_calib_und";
      goto FREE_RETURN;
    }

    memcpy(und->spread_vol_time, spread_vol_time,
           spread_vol_n * sizeof(double));
    memcpy(und->spread_vol_floor, spread_vol_floor,
           spread_vol_n * sizeof(double));
    memcpy(und->spread_vol_cap, spread_vol_cap, spread_vol_n * sizeof(double));
  }

  /*	Exercise dates for calibration */

  end_struct = fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
  end_swap = end_struct;
  if (exo_leg->cpn[exo_leg->num_cpn - 1].pay_date > end_swap) {
    end_swap = exo_leg->cpn[exo_leg->num_cpn - 1].pay_date;
  }

  if (ccf->num_calls > 0 &&
      !(ccf->num_calls == 1 && ccf->call[0].ex_date <= today + eod_flag)) {
    /*	If call dates  , choose call dates as option expiries for calibration */
    nex = ccf->num_calls;
    lex = (long *)calloc(nex, sizeof(long));
    if (!lex) {
      err = "Allocation error (2) in ccf_calib_und";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      lex[i] = ccf->call[i].ex_date;
    }
  } else {
    /*	If no call dates  , exit */
    ccf->num_calls = 0;
    goto FREE_RETURN;
  }

  /*	Implement force atm */
  if (force_atm) {
    strike_type = 0;
    short_strikes = NULL;
    long_strikes = NULL;
    smessage("FORCE ATM flag detected - calibrating ATM");
  } else {
    /*	Find strikes */
    short_strikes = (double *)calloc(nex, sizeof(double));
    long_strikes = (double *)calloc(nex, sizeof(double));
    if (!short_strikes || !long_strikes) {
      err = "Allocation error (3) in ccf_calib_und";
      goto FREE_RETURN;
    }

    err = find_strikes(today, end_swap, yc, vc, ref, swap_freq, swap_basis,
                       get_cash_vol, ccf, short_strikes, long_strikes,
                       &strike_type);
    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Apply lambda shift if no calibration required */
  if (fabs(cal_lambda_shift) > 1.0E-08 && fix_lambda) {
    if (!lam_ts) {
      /* shift the lambda */
      lambda += cal_lambda_shift;
    } else {
      for (i = 0; i < tsnlam; i++) {
        copy_tslam[i] += cal_lambda_shift;
      }
    }
  }

  /*	LGM 2D underlying */
  if (lam_ts) {
    err = cpd_calib_diagonal_tauts(
        yc, vc, ref, get_cash_vol, 0.0, 0, nex, lex, end_swap, long_strikes,
        strike_type, max_std_long, swap_freq, swap_basis, skip_last, min_fact,
        max_fact, use_jumps, tsnlam, tslamtime, copy_tslam, 2, alpha, gamma,
        rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma),
        &(und->inst_data));
  } else {
    err = cpd_calib_diagonal(
        yc, vc, ref, get_cash_vol, 0.0, 0, nex, lex, end_swap, long_strikes,
        short_strikes, strike_type, max_std_long, max_std_short, swap_freq,
        swap_basis, fix_lambda, one_f_equi, skip_last, long_prec, short_prec,
        min_fact, max_fact, use_jumps, proba_weight, proba, &lambda, 2, alpha,
        gamma, rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma),
        &(und->inst_data));
  }
  if (err) {
    goto FREE_RETURN;
  }
  und->has_inst_data = 1;

  if (fabs(cal_lambda_shift) > 1.0E-08 && !fix_lambda) {
    und->sigma_n = 0;
    if (und->sigma_time)
      free(und->sigma_time);
    und->sigma_time = NULL;
    if (und->sigma)
      free(und->sigma);
    und->sigma = NULL;

    if (!lam_ts) {
      /* shift the lambda */
      lambda += cal_lambda_shift;

      /* recalibrate to swaptions */
      err = cpd_calib_diagonal(
          yc, vc, ref, get_cash_vol, 0.0, 0, nex, lex, end_swap, long_strikes,
          short_strikes, strike_type, max_std_long, max_std_short, swap_freq,
          swap_basis, 1, one_f_equi, skip_last, long_prec, short_prec, min_fact,
          max_fact, use_jumps, proba_weight, proba, &lambda, 2, alpha, gamma,
          rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma), NULL);

      if (err) {
        goto FREE_RETURN;
      }
    } else {
      for (i = 0; i < tsnlam; i++) {
        copy_tslam[i] += cal_lambda_shift;
      }

      err = cpd_calib_diagonal_tauts(
          yc, vc, ref, get_cash_vol, 0.0, 0, nex, lex, end_swap, long_strikes,
          strike_type, max_std_long, swap_freq, swap_basis, skip_last, min_fact,
          max_fact, use_jumps, tsnlam, tslamtime, copy_tslam, 2, alpha, gamma,
          rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma),
          &(und->inst_data));
    }
  }

  if (lam_ts) {
    und->lambda_n = tsnlam;
  } else {
    und->lambda_n = 1;
  }

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));
  und->lambda_time = (double *)calloc(und->lambda_n, sizeof(double));
  und->lambda_date = (double *)calloc(und->lambda_n, sizeof(double));
  und->lambda = (double *)calloc(und->lambda_n, sizeof(double));

  if (!und->sigma_date || !und->lambda_time || !und->lambda_date) {
    err = "Allocation error (6) in ccf_calib_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] = today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

  if (lam_ts) {
    for (i = 0; i < und->lambda_n; i++) {
      und->lambda_time[i] = tslamtime[i];
      und->lambda_date[i] =
          today + und->lambda_time[i] * DAYS_IN_YEAR + 1.0e-08;
      und->lambda[i] = copy_tslam[i];
    }
  } else {
    und->lambda_time[0] = und->sigma_time[0];
    und->lambda_date[0] = und->sigma_date[0];
    und->lambda[0] = lambda;
  }

  if (fabs(cal_vol_shift) > 1.0E-08) {
    /* Shift LGM vols */
    for (i = 0; i < und->sigma_n; i++) {
      if (cal_vol_shift_type == 0) {
        und->sigma[i] += cal_vol_shift;
      } else {
        und->sigma[i] += und->sigma[i] * cal_vol_shift;
      }
    }
  }

FREE_RETURN:

  if (err)
    ccf_free_und(und);

  if (lex)
    free(lex);
  if (long_strikes)
    free(long_strikes);
  if (short_strikes)
    free(short_strikes);

  if (copy_tslam)
    free(copy_tslam);

  return err;
}

void ccf_copy_und(CCF_UND src, CCF_UND dest) {
  strcpy(dest->name, src->name);
  dest->today = src->today;
  strcpy(dest->yc, src->yc);
  strcpy(dest->vc, src->vc);
  strcpy(dest->ref, src->ref);
  strcpy(dest->swap_freq, src->swap_freq);
  strcpy(dest->swap_basis, src->swap_basis);

  dest->sigma_n = src->sigma_n;
  dest->lambda_n = src->lambda_n;
  dest->spread_vol_n = src->spread_vol_n;

  if (src->sigma_date) {
    dest->sigma_date = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_date, src->sigma_date, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_date = NULL;
  }

  if (src->sigma_time) {
    dest->sigma_time = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_time, src->sigma_time, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_time = NULL;
  }

  if (src->sigma) {
    dest->sigma = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma, src->sigma, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma = NULL;
  }

  if (src->lambda_date) {
    dest->lambda_date = (double *)calloc(dest->lambda_n, sizeof(double));
    memcpy(dest->lambda_date, src->lambda_date,
           dest->lambda_n * sizeof(double));
  } else {
    dest->lambda_date = NULL;
  }

  if (src->lambda_time) {
    dest->lambda_time = (double *)calloc(dest->lambda_n, sizeof(double));
    memcpy(dest->lambda_time, src->lambda_time,
           dest->lambda_n * sizeof(double));
  } else {
    dest->lambda_time = NULL;
  }

  if (src->lambda) {
    dest->lambda = (double *)calloc(dest->lambda_n, sizeof(double));
    memcpy(dest->lambda, src->lambda, dest->lambda_n * sizeof(double));
  } else {
    dest->lambda = NULL;
  }

  if (src->spread_vol_time) {
    dest->spread_vol_time =
        (double *)calloc(dest->spread_vol_n, sizeof(double));
    memcpy(dest->spread_vol_time, src->spread_vol_time,
           dest->spread_vol_n * sizeof(double));
  } else {
    dest->spread_vol_time = NULL;
  }

  if (src->spread_vol_floor) {
    dest->spread_vol_floor =
        (double *)calloc(dest->spread_vol_n, sizeof(double));
    memcpy(dest->spread_vol_floor, src->spread_vol_floor,
           dest->spread_vol_n * sizeof(double));
  } else {
    dest->spread_vol_floor = NULL;
  }

  if (src->spread_vol_cap) {
    dest->spread_vol_cap = (double *)calloc(dest->spread_vol_n, sizeof(double));
    memcpy(dest->spread_vol_cap, src->spread_vol_cap,
           dest->spread_vol_n * sizeof(double));
  } else {
    dest->spread_vol_cap = NULL;
  }

  if (src->has_inst_data) {
    cpd_copy_calib_inst_data(&(dest->inst_data), &(src->inst_data));
    dest->has_inst_data = 1;
  } else {
    dest->has_inst_data = 0;
  }

  if (src->has_fwd_iv) {
    dest->has_fwd_iv = 1;
    dest->nb_fwdiv = src->nb_fwdiv;

    dest->exercise_date = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->exercise_date, src->exercise_date,
           dest->nb_fwdiv * sizeof(double));

    dest->market_fwdiv = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->market_fwdiv, src->market_fwdiv,
           dest->nb_fwdiv * sizeof(double));

    dest->model_fwdiv = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->model_fwdiv, src->model_fwdiv,
           dest->nb_fwdiv * sizeof(double));

    dest->extra_fees = (double *)calloc(dest->nb_fwdiv, sizeof(double));
    memcpy(dest->extra_fees, src->extra_fees, dest->nb_fwdiv * sizeof(double));
  } else {
    dest->has_fwd_iv = 0;
    dest->nb_fwdiv = 0;
    dest->exercise_date = NULL;
    dest->market_fwdiv = NULL;
    dest->model_fwdiv = NULL;
    dest->extra_fees = NULL;
  }

  dest->alpha = src->alpha;
  dest->gamma = src->gamma;
  dest->rho = src->rho;

  dest->cvg_sv = src->cvg_sv;
  dest->is_corr = src->is_corr;
}

Err ccf_free_und(CCF_UND und) {
  if (und->sigma_date)
    free(und->sigma_date);
  if (und->sigma_time)
    free(und->sigma_time);
  if (und->sigma)
    free(und->sigma);
  if (und->lambda_date)
    free(und->lambda_date);
  if (und->lambda_time)
    free(und->lambda_time);
  if (und->lambda)
    free(und->lambda);
  if (und->spread_vol_time)
    free(und->spread_vol_time);
  if (und->spread_vol_floor)
    free(und->spread_vol_floor);
  if (und->spread_vol_cap)
    free(und->spread_vol_cap);

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;
  und->lambda_date = NULL;
  und->lambda_time = NULL;
  und->lambda = NULL;
  und->spread_vol_time = NULL;
  und->spread_vol_floor = NULL;
  und->spread_vol_cap = NULL;

  if (und->has_inst_data) {
    cpd_free_calib_inst_data(&und->inst_data);
    und->has_inst_data = 0;
  }

  if (und->has_fwd_iv) {
    if (und->market_fwdiv)
      free(und->market_fwdiv);
    if (und->model_fwdiv)
      free(und->model_fwdiv);
    if (und->extra_fees)
      free(und->extra_fees);
    if (und->exercise_date)
      free(und->exercise_date);

    und->has_fwd_iv = 0;
    und->nb_fwdiv = 0;
  }

  return NULL;
}

/*	Function for constants used for reconstruction and evaluation  */

/*	Just round a double to the closest long */
static long strnd(double x) {
  static long y;
  y = (long)x;
  if (fabs(y + 1 - x) < fabs(y - x) && fabs(y + 1 - x) < fabs(y - 1 - x)) {
    return y + 1;
  } else if (fabs(y - x) < fabs(y - 1 - x)) {
    return y;
  } else {
    return y - 1;
  }
}

/*	THE function */
Err ccf_fill_eval_const(CCF_UND und, CCF_STR ccf,
                        /*	Index of the current call */
                        int call_idx, CCF_EVAL_CONST eval_const) {
  CF_CALL call;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  long evt_date;
  double evt_time;
  double lam1, lam2, phi1, phi2, phi12;
  long today;
  long temp_date;
  double temp_time;
  int tmpidx, tmpidx2;
  char *yc;
  int i, j, k, l;
  int iNumOpt;
  double temp, spread_vol;
  double sv_time;
  double *ts_time;
  double ta, tb, t1, t2, dt;
  int nb_ts, nb_ts_minus1;
  double gam1, gam2, gam12, fact1, fact1temp, fact2, fact2temp;
  Err err = NULL;

  /*	Initialise
          ----------	*/

  /*	Extract info from call */
  call = ccf->call + call_idx;
  evt_date = call->ex_date;
  evt_time = call->ex_time;

  /*	Extract info for structure */
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  /*	Extract info from und */
  today = und->today;
  yc = (char *)(und->yc);

  err = LGM2FDtsPhi(evt_time, und->sigma_time, und->sigma_n, und->sigma,
                    und->lambda_time, und->lambda_n, und->lambda, und->alpha,
                    und->gamma, und->rho, &phi1, &phi2, &phi12);

  /*	Get the maturities of the required discount factors
          --------------------------------------------------- */
  eval_const->num_df = 0;

  /*	Funding leg */
  if (call->num_fund_cpn > 0) {
    eval_const->do_fund = 1;
    eval_const->num_fund_cpn = call->num_fund_cpn;
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = fund_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }
  } else {
    eval_const->do_fund = 0;
    eval_const->num_fund_cpn = 0;
    err = "No funding coupons left at call date";
    goto FREE_RETURN;
  }
  eval_const->df_mat[eval_const->num_df] =
      fund_leg->cpn[call->fund_idx].start_time - evt_time;
  (eval_const->num_df)++;

  /*	CMS Leg */
  if (call->num_cf_cpn > 0) {
    eval_const->do_cf_disc = 1;
    eval_const->do_cf_fwd = 1;
    eval_const->num_cf_cpn = call->num_cf_cpn;

    /*	Discounting */
    for (i = call->cf_idx; i < exo_leg->num_cpn; i++) {
      exo_cpn = exo_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = exo_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }

    /*	CMS Coupons	*/
    for (i = 0; i < call->num_cf_cpn; i++) {
      exo_cpn = exo_leg->cpn + call->cf_idx + i;

      for (l = 0; l < exo_cpn->ncms; l++) {
        cms = &(exo_cpn->cms[l]);

        eval_const->df_mat[eval_const->num_df] = cms->start_time - evt_time;
        (eval_const->num_df)++;
        for (k = 0; k < cms->num_cpn; k++) {
          eval_const->df_mat[eval_const->num_df] =
              cms->cpn_pay_time[k] - evt_time;
          (eval_const->num_df)++;
        }
      }
    }
  } else {
    eval_const->do_cf_disc = 0;
    eval_const->do_cf_fwd = 0;
    eval_const->num_cf_cpn = 0;
    err = "No exotic coupons left at call date";
    goto FREE_RETURN;
  }

  /*	Fee	*/
  eval_const->df_mat[eval_const->num_df] = call->set_time - evt_time;
  (eval_const->num_df)++;

  /*	Sort and unique */
  num_f_sort_vector(eval_const->num_df, eval_const->df_mat);
  num_f_unique_vector(&(eval_const->num_df), eval_const->df_mat);

  /*	Precalculate DFF  , gamma  , 0.5 * gamma * gamma
          --------------------------------------------	*/

  /* Find first used lambda */
  j = 0;
  nb_ts = und->lambda_n;
  nb_ts_minus1 = nb_ts - 1;
  ts_time = und->lambda_time;

  while (j < nb_ts && ts_time[j] < evt_time) {
    j++;
  }
  j = min(j, nb_ts - 1);

  gam1 = gam2 = gam12 = 0.0;
  lam1 = und->lambda[j];
  lam2 = lam1 + und->gamma;

  fact1 = 0.0;
  fact2 = 0.0;

  t1 = evt_time;

  for (i = 0; i < eval_const->num_df; i++) {
    temp_time = eval_const->df_mat[i];
    temp_date = evt_date + strnd(temp_time * 365.0);
    t2 = evt_time + temp_time;

    ta = t1;
    tb = ts_time[j];

    gam1 = 0.0;
    gam2 = 0.0;
    fact1temp = 0.0;
    fact2temp = 0.0;

    while (tb < t2 && j < nb_ts_minus1) {
      dt = (tb - ta);

      gam1 += exp(-fact1temp) * (1.0 - exp(-lam1 * dt)) / lam1;
      gam2 += exp(-fact2temp) * (1.0 - exp(-lam2 * dt)) / lam2;

      fact1temp += lam1 * dt;
      fact2temp += lam2 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = und->lambda[j];
      lam2 = lam1 + und->gamma;
    }

    dt = (t2 - ta);

    gam1 += exp(-fact1temp) * (1.0 - exp(-lam1 * dt)) / lam1;
    gam2 += exp(-fact2temp) * (1.0 - exp(-lam2 * dt)) / lam2;

    fact1temp += lam1 * dt;
    fact2temp += lam2 * dt;

    if (i > 0) {
      eval_const->df_beta[i] = eval_const->df_beta[i - 1] + exp(-fact1) * gam1;
      eval_const->df_gamma[i] =
          eval_const->df_gamma[i - 1] + exp(-fact2) * gam2;
    } else {
      eval_const->df_beta[i] = gam1;
      eval_const->df_gamma[i] = gam2;
    }

    eval_const->df_alpha[i] =
        swp_f_zr(evt_date, temp_date, yc) * temp_time +
        0.5 * (eval_const->df_beta[i] * eval_const->df_beta[i] * phi1 +
               eval_const->df_gamma[i] * eval_const->df_gamma[i] * phi2) +
        und->rho * eval_const->df_beta[i] * eval_const->df_gamma[i] * phi12;

    fact1 += fact1temp;
    fact2 += fact2temp;

    t1 = t2;
  }

  /*	Find all dfs indices
          --------------------	*/

  /*	Funding Leg */
  tmpidx = 0;
  temp_time = fund_leg->cpn[call->fund_idx].start_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->start_idx = tmpidx;
  if (eval_const->do_fund) {
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      temp_time = fund_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->fund_idx[i - call->fund_idx] = tmpidx;
    }
  }

  /*	CMS Leg */
  if (eval_const->do_cf_disc && eval_const->do_cf_fwd) {
    /*	Discounting */
    tmpidx = 0;
    for (i = call->cf_idx; i < exo_leg->num_cpn; i++) {
      exo_cpn = exo_leg->cpn + i;
      temp_time = exo_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->cf_disc_idx[i - call->cf_idx] = tmpidx;
    }

    /*	CMS Coupons	*/
    tmpidx = 0;
    for (i = 0; i < call->num_cf_cpn; i++) {
      exo_cpn = exo_leg->cpn + call->cf_idx + i;

      for (l = 0; l < exo_cpn->ncms; l++) {
        cms = &(exo_cpn->cms[l]);

        temp_time = cms->start_time - evt_time;
        while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
          tmpidx++;
        eval_const->cf_fwd_start_idx[l][i] = tmpidx;
        tmpidx2 = tmpidx;
        for (k = 0; k < cms->num_cpn; k++) {
          temp_time = cms->cpn_pay_time[k] - evt_time;
          while (eval_const->df_mat[tmpidx2] < temp_time - 1.0e-08)
            tmpidx2++;
          eval_const->cf_fwd_idx[l][i][k] = tmpidx2;
        }
      }
    }
  }

  /*	Fee	*/
  tmpidx = 0;
  temp_time = call->set_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->fee_idx = tmpidx;

  /*	Precalculate vols
          -----------------	*/
  if (eval_const->do_cf_disc && eval_const->do_cf_fwd) {
    for (i = 0; i < call->num_cf_cpn; i++) {
      exo_cpn = ccf->cf_leg->cpn + call->cf_idx + i;
      temp_time = (exo_cpn->cms_fix_time - evt_time);

      /* Switch if use string of options or shifted log model */
      if (!exo_cpn->use_cfoptions && !exo_cpn->use_SL) {
        /* Classical Case */

        switch (exo_cpn->type) {
        case 0:
        default:
          /*	Midat */
          break;

        case 1:
        case 2:
          /*	CIF or CMS */
          if (temp_time > 1.0E-08) {
            /* CMS Vol */
            cms = &(exo_cpn->cms[0]);
            eval_const->cf_cms_vol[0][i] = cms->atmvar * temp_time;
            eval_const->cf_spr_floor_std[i] =
                fabs(exo_cpn->alphabeta[0]) * sqrt(cms->floorvar * temp_time);
            eval_const->cf_spr_cap_std[i] =
                fabs(exo_cpn->alphabeta[0]) * sqrt(cms->capvar * temp_time);
          } else {
            eval_const->cf_cms_vol[0][i] = eval_const->cf_spr_floor_std[i] =
                eval_const->cf_spr_cap_std[i] = 0.0;
          }

          if (exo_cpn->floored) {
            if (eval_const->cf_spr_floor_std[i] > 1.0E-08) {
              eval_const->cf_spr_floor_type[i] = 4; /*	Put */
            } else {
              eval_const->cf_spr_floor_type[i] = 2; /*	Put IV */
            }

            eval_const->cf_spr_floor_str[i] = exo_cpn->floor - exo_cpn->gamma;
          } else {
            eval_const->cf_spr_floor_type[i] = 0; /*	None */
            eval_const->cf_spr_floor_str[i] = 0.0;
          }

          if (exo_cpn->capped) {
            if (eval_const->cf_spr_cap_std[i] > 1.0E-08) {
              eval_const->cf_spr_cap_type[i] = 3; /*	Call */
            } else {
              eval_const->cf_spr_cap_type[i] = 1; /*	Call IV */
            }

            eval_const->cf_spr_cap_str[i] = exo_cpn->cap - exo_cpn->gamma;
          } else {
            eval_const->cf_spr_cap_type[i] = 0; /*	None */
            eval_const->cf_spr_cap_str[i] = 0.0;
          }
          break;

        case 3:
          /* CMS Spread */

          if (temp_time > 1.0E-08) {
            /* CMS Vol */
            for (l = 0; l < 2; l++) {
              cms = &(exo_cpn->cms[l]);
              eval_const->cf_cms_vol[l][i] = cms->atmvar * temp_time;
            }

            /*	Converging/Sliding model on spread vol/correl */
            if (und->cvg_sv) {
              sv_time = exo_cpn->cms_fix_time;
            } else {
              sv_time = temp_time;
            }

            /*	Spread Vol in Floor */
            spread_vol = interp(und->spread_vol_time, und->spread_vol_floor,
                                und->spread_vol_n, sv_time, 0, &temp);

            if (und->is_corr) {
              eval_const->cf_spr_floor_std[i] =
                  sqrt(exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] *
                           exo_cpn->cms[0].floorvar +
                       exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] *
                           exo_cpn->cms[1].floorvar +
                       2.0 * spread_vol * exo_cpn->alphabeta[0] *
                           exo_cpn->alphabeta[1] *
                           sqrt(exo_cpn->cms[0].floorvar *
                                exo_cpn->cms[1].floorvar));
            } else {
              eval_const->cf_spr_floor_std[i] = spread_vol;
            }
            eval_const->cf_spr_floor_std[i] *= sqrt(temp_time);

            /*	Spread Vol in Cap */
            spread_vol = interp(und->spread_vol_time, und->spread_vol_cap,
                                und->spread_vol_n, sv_time, 0, &temp);

            if (und->is_corr) {
              eval_const->cf_spr_cap_std[i] = sqrt(
                  exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] *
                      exo_cpn->cms[0].capvar +
                  exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] *
                      exo_cpn->cms[1].capvar +
                  2.0 * spread_vol * exo_cpn->alphabeta[0] *
                      exo_cpn->alphabeta[1] *
                      sqrt(exo_cpn->cms[0].capvar * exo_cpn->cms[1].capvar));
            } else {
              eval_const->cf_spr_cap_std[i] = spread_vol;
            }
            eval_const->cf_spr_cap_std[i] *= sqrt(temp_time);
          } else {
            eval_const->cf_cms_vol[0][i] = eval_const->cf_cms_vol[1][i] =
                eval_const->cf_spr_floor_std[i] =
                    eval_const->cf_spr_cap_std[i] = 0.0;
          }

          if (exo_cpn->floored) {
            if (eval_const->cf_spr_floor_std[i] > 1.0E-08) {
              eval_const->cf_spr_floor_type[i] = 4; /*	Put */
            } else {
              eval_const->cf_spr_floor_type[i] = 2; /*	Put IV */
            }

            eval_const->cf_spr_floor_str[i] = exo_cpn->floor - exo_cpn->gamma;
          } else {
            eval_const->cf_spr_floor_type[i] = 0; /*	None */
            eval_const->cf_spr_floor_str[i] = 0.0;
          }

          if (exo_cpn->capped) {
            if (eval_const->cf_spr_cap_std[i] > 1.0E-08) {
              eval_const->cf_spr_cap_type[i] = 3; /*	Call */
            } else {
              eval_const->cf_spr_cap_type[i] = 1; /*	Call IV */
            }

            eval_const->cf_spr_cap_str[i] = exo_cpn->cap - exo_cpn->gamma;
          } else {
            eval_const->cf_spr_cap_type[i] = 0; /*	None */
            eval_const->cf_spr_cap_str[i] = 0.0;
          }
          break;
        }
      } else {
        if (exo_cpn->use_SL) {
          if (!exo_cpn->use_cfoptions) {
            /* Shifted log case wo string of options */

            if (temp_time > 1.0E-08) {
              /* CMS Vol */
              for (l = 0; l < 2; l++) {
                cms = &(exo_cpn->cms[l]);
                eval_const->cf_cms_vol[l][i] = cms->atmvar * temp_time;
              }

              /*	Converging/Sliding model on spread vol/correl */
              if (und->cvg_sv) {
                sv_time = exo_cpn->cms_fix_time;
              } else {
                sv_time = temp_time;
              }

              /*	Use the Spread information in Floor as the correlation
               */
              eval_const->cf_spr_slcorr[i] =
                  interp(und->spread_vol_time, und->spread_vol_floor,
                         und->spread_vol_n, sv_time, 0, &temp);
            } else {
              for (l = 0; l < 2; l++) {
                eval_const->cf_cms_vol[l][i] = 0.0;
              }

              eval_const->cf_spr_slcorr[i] = 0.0;
            }

            if (exo_cpn->floored) {
              eval_const->cf_spr_floor_type[i] = 4; /*	Put */
              eval_const->cf_spr_floor_str[i] = exo_cpn->floor - exo_cpn->gamma;
            } else {
              eval_const->cf_spr_floor_type[i] = 0; /*	None */
              eval_const->cf_spr_floor_str[i] = 0.0;
            }

            if (exo_cpn->capped) {
              eval_const->cf_spr_cap_type[i] = 3; /*	Call */
              eval_const->cf_spr_cap_str[i] = exo_cpn->cap - exo_cpn->gamma;
            } else {
              eval_const->cf_spr_cap_type[i] = 0; /*	None */
              eval_const->cf_spr_cap_str[i] = 0.0;
            }
          } else {
            /* Shifted log case with string of options */
            if (temp_time > 1.0E-08) {
              /* CMS Vol */
              for (l = 0; l < 2; l++) {
                cms = &(exo_cpn->cms[l]);
                eval_const->cf_cms_vol[l][i] = cms->atmvar * temp_time;
              }

              /*	Converging/Sliding model on spread vol/correl */
              if (und->cvg_sv) {
                sv_time = exo_cpn->cms_fix_time;
              } else {
                sv_time = temp_time;
              }

              /*	Use the Spread information in Floor as the correlation
               */
              eval_const->cf_spr_slcorr[i] =
                  interp(und->spread_vol_time, und->spread_vol_floor,
                         und->spread_vol_n, sv_time, 0, &temp);
            } else {
              for (l = 0; l < 2; l++) {
                eval_const->cf_cms_vol[l][i] = 0.0;
              }

              eval_const->cf_spr_slcorr[i] = 0.0;
            }
          }

        } else {
          /* string of options wo sl */
          if (temp_time > 1.0E-08) {
            /* CMS Vol */
            for (l = 0; l < 2; l++) {
              cms = &(exo_cpn->cms[l]);
              eval_const->cf_cms_vol[l][i] = cms->atmvar * temp_time;
            }

            /*	Converging/Sliding model on spread vol/correl */
            if (und->cvg_sv) {
              sv_time = exo_cpn->cms_fix_time;
            } else {
              sv_time = temp_time;
            }

            /*	Use the Spread information in Floor as the correlation */
            spread_vol = interp(und->spread_vol_time, und->spread_vol_floor,
                                und->spread_vol_n, sv_time, 0, &temp);

            /* Calculation of the normal std */
            for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
              eval_const->cf_spr_opt_std[i][iNumOpt] = sqrt(
                  temp_time * (exo_cpn->alphabeta[0] * exo_cpn->alphabeta[0] *
                                   exo_cpn->cms[0].optvar[iNumOpt] +
                               exo_cpn->alphabeta[1] * exo_cpn->alphabeta[1] *
                                   exo_cpn->cms[1].optvar[iNumOpt] +
                               2.0 * spread_vol * exo_cpn->alphabeta[0] *
                                   exo_cpn->alphabeta[1] *
                                   sqrt(exo_cpn->cms[0].optvar[iNumOpt] *
                                        exo_cpn->cms[1].optvar[iNumOpt])));
            }
          } else {
            for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
              eval_const->cf_spr_opt_std[i][iNumOpt] = 0.0;
            }
          }
        }
      }
    }
  }

  /*	End of precalculations
          ---------------------- */

FREE_RETURN:

  return err;
}

Err ccf_fill_adi_arg(
    CCF_UND und, CCF_STR ccf,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, CCF_ADI_ARG adi_arg) {
  Err err = NULL;
  CCF_PAY_ARG ccf_prm;
  int i, j;

  /*	Initialise */
  adi_arg->time = NULL;
  adi_arg->date = NULL;

  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Compute time steps */

  /*	Copy event dates */

  adi_arg->nstp = ccf->num_calls;
  adi_arg->nstpx = req_stpx;

  adi_arg->time = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->time) {
    err = "Memory allocation error (1) in ccf_fill_adi_arg";
    goto FREE_RETURN;
  }
  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->time[i] = ccf->call[i].ex_time;
  }

  /*	Fill time vector */

  /*	Add today if required */
  if (adi_arg->time[0] < -EPS) {
    err = "Past event date in ccf_fill_adi_arg";
    goto FREE_RETURN;
  }
  if (adi_arg->time[0] > EPS) {
    num_f_add_number(&(adi_arg->nstp), &(adi_arg->time), 0.0);
    num_f_sort_vector(adi_arg->nstp, adi_arg->time);
    num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);
  }

  /*	If only one event today  , add empty event */
  if (adi_arg->nstp == 1) {
    num_f_add_number(&(adi_arg->nstp), (&adi_arg->time), 1.0);
  }

  /*	Fill the vector */
  num_f_fill_vector_newalgo(&(adi_arg->nstp), &(adi_arg->time), req_stp);

  /*	Make dates */
  adi_arg->date = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->date) {
    err = "Memory allocation error (2) in ccf_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->date[i] = und->today + DAYS_IN_YEAR * adi_arg->time[i];

    if (i > 0 && adi_arg->date[i] - adi_arg->date[i - 1] >= 1) {
      adi_arg->date[i] = (long)(adi_arg->date[i] + 1.0e-08);
      adi_arg->time[i] = YEARS_IN_DAY * (adi_arg->date[i] - und->today);
    }
  }

  /*	Sigma */

  adi_arg->sig_time = und->sigma_time;
  adi_arg->sig1 = und->sigma;
  adi_arg->nb_sig = und->sigma_n;

  /*	Lambdas and ratio */

  adi_arg->lam_time = und->lambda_time;
  adi_arg->lam = und->lambda;
  adi_arg->nb_lam = und->lambda_n;

  adi_arg->alpha = und->alpha;
  adi_arg->gamma = und->gamma;
  adi_arg->rho = und->rho;

  /*	Spot fx and yield curves */

  strcpy(adi_arg->yc, und->yc);

  /*	Fill distributions */

  adi_arg->ifr = (double *)calloc(adi_arg->nstp, sizeof(double));

  if (!adi_arg->ifr) {
    err = "Memory allocation error (3) in ccf_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp - 1; i++) {
    adi_arg->ifr[i] =
        swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], adi_arg->yc);
  }

  /*	Fill limit conditions (product) */

  adi_arg->is_event = (int *)calloc(adi_arg->nstp, sizeof(int));
  adi_arg->void_prm = (void **)calloc(adi_arg->nstp, sizeof(void *));

  if (!adi_arg->is_event || !adi_arg->void_prm) {
    err = "Memory allocation error (4) in ccf_fill_adi_arg";
    goto FREE_RETURN;
  }

  j = ccf->num_calls - 1;

  for (i = adi_arg->nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(adi_arg->date[i] - ccf->call[j].ex_date) < 1.0e-08) {
      ccf_prm = malloc(sizeof(ccf_pay_arg));
      ccf_prm->und = und;
      ccf_prm->ccf = ccf;

      ccf_prm->call_idx = j;

      err = ccf_fill_eval_const(und, ccf, j, &(ccf_prm->eval_const));

      adi_arg->is_event[i] = 1;
      adi_arg->void_prm[i] = (void *)ccf_prm;

      j--;
    } else {
      adi_arg->is_event[i] = 0;
      adi_arg->void_prm[i] = NULL;
    }
  }

FREE_RETURN:

  if (err) {
    ccf_free_adi_arg(adi_arg);
  }

  return err;
}

Err ccf_free_adi_arg(CCF_ADI_ARG adi_arg) {
  int i;
  CCF_PAY_ARG ccf_prm;

  if (adi_arg->time)
    free(adi_arg->time);
  if (adi_arg->date)
    free(adi_arg->date);
  if (adi_arg->ifr)
    free(adi_arg->ifr);
  if (adi_arg->is_event)
    free(adi_arg->is_event);

  if (adi_arg->void_prm) {
    for (i = 0; i < adi_arg->nstp; i++) {
      ccf_prm = (CCF_PAY_ARG)(adi_arg->void_prm[i]);
      free(ccf_prm);
    }

    free(adi_arg->void_prm);
  }

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->sig1 = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;

  return NULL;
}

Err ccf_fill_mc_arg(
    CCF_UND und, CCF_STR ccf,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    long req_paths, int jumping_num, CCF_MC_ARG mc_arg) {
  Err err = NULL;
  CCF_PAY_ARG ccf_prm;
  int i;
  int flag;

  // New variables for lambda TS
  double *new_time = NULL, *new_lambda = NULL, *new_sigma = NULL,
         *new_lambda2 = NULL;

  int n_new_time;

  /*	Initialise */
  mc_arg->time = NULL;
  mc_arg->date = NULL;

  mc_arg->dom_fwd1 = NULL;
  mc_arg->dom_fwd2 = NULL;
  mc_arg->dom_exp1 = NULL;
  mc_arg->dom_exp2 = NULL;
  mc_arg->dom_phi1 = NULL;
  mc_arg->dom_phi2 = NULL;
  mc_arg->dom_phi12 = NULL;
  mc_arg->dom_gam1_fwd = NULL;
  mc_arg->dom_gam2_fwd = NULL;
  mc_arg->dom_bond_pay = NULL;
  mc_arg->dom_gam1_pay = NULL;
  mc_arg->dom_gam2_pay = NULL;
  mc_arg->covar = NULL;

  mc_arg->void_prm = NULL;

  /*	Copy event dates */

  mc_arg->nb_dates = ccf->num_calls;
  mc_arg->npaths = req_paths;
  mc_arg->jumping_num = jumping_num;

  mc_arg->time = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  if (!mc_arg->time) {
    err = "Memory allocation error (1) in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }
  for (i = 0; i < mc_arg->nb_dates; i++) {
    mc_arg->time[i] = ccf->call[i].ex_time;
  }

  /*	Fill time vector */

  /*	Add today if required */
  if (mc_arg->time[0] < -EPS) {
    err = "Past event date in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }

  flag = 0;

  if (mc_arg->time[0] > EPS) {
    num_f_add_number(&(mc_arg->nb_dates), &(mc_arg->time), 0.0);
    num_f_sort_vector(mc_arg->nb_dates, mc_arg->time);
    num_f_unique_vector(&(mc_arg->nb_dates), mc_arg->time);
    flag = 1;
  }

  /*	Make dates */
  mc_arg->date = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  if (!mc_arg->date) {
    err = "Memory allocation error (2) in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < mc_arg->nb_dates; i++) {
    mc_arg->date[i] = und->today + DAYS_IN_YEAR * mc_arg->time[i];

    if (i > 0 && mc_arg->date[i] - mc_arg->date[i - 1] >= 1) {
      mc_arg->date[i] = (long)(mc_arg->date[i] + 1.0e-08);
      mc_arg->time[i] = YEARS_IN_DAY * (mc_arg->date[i] - und->today);
    }
  }

  // Merge the term structures of lambda and sigmas
  err = merge_lambda_sigma_ts(und->sigma, und->sigma_time, und->sigma_n,
                              und->lambda, und->lambda_time, und->lambda_n,
                              &new_time, &new_lambda, &new_sigma, &n_new_time);
  if (err) {
    goto FREE_RETURN;
  }

  // Create the array of new_lambda2[i] = new_lambda[i] + gamma
  new_lambda2 = (double *)calloc(n_new_time, sizeof(double));

  if (!new_lambda2) {
    err = "Memory allocation failure in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < n_new_time; i++) {
    new_lambda2[i] = new_lambda[i] + und->gamma;
  }

  /*	Distributions */

  mc_arg->dom_fwd1 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_fwd2 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_exp1 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_exp2 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_phi1 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_phi2 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_phi12 = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_gam1_fwd = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_gam2_fwd = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_bond_pay = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_gam1_pay = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->dom_gam2_pay = (double *)calloc(mc_arg->nb_dates, sizeof(double));
  mc_arg->covar = f3tensor(0, mc_arg->nb_dates - 1, 0, 1, 0, 1);

  if (!mc_arg->dom_fwd1 || !mc_arg->dom_fwd2 || !mc_arg->dom_phi1 ||
      !mc_arg->dom_phi2 || !mc_arg->dom_phi12 || !mc_arg->dom_exp1 ||
      !mc_arg->dom_exp2 || !mc_arg->dom_gam1_fwd || !mc_arg->dom_gam2_fwd ||
      !mc_arg->dom_bond_pay || !mc_arg->dom_gam1_pay || !mc_arg->dom_gam2_pay ||
      !mc_arg->covar) {
    err = "Memory allocation error (3) in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }

  if (mc_arg->jumping_num) {
    mc_arg->pay_date = ccf->call[0].ex_date;
    mc_arg->pay_time = ccf->call[0].ex_time;
  } else {
    mc_arg->pay_date = ccf->call[ccf->num_calls - 1].ex_date;
    mc_arg->pay_time = ccf->call[ccf->num_calls - 1].ex_time;
  }

  fill_mc_init_lgm2f_lambda(
      mc_arg->jumping_num, mc_arg->pay_date, mc_arg->pay_time, mc_arg->date,
      mc_arg->time, mc_arg->nb_dates, new_time, n_new_time, new_sigma,
      new_lambda, new_lambda2, und->alpha, und->gamma, und->rho, und->yc,
      mc_arg->dom_fwd1, mc_arg->dom_fwd2, mc_arg->dom_exp1, mc_arg->dom_exp2,
      mc_arg->dom_phi1, mc_arg->dom_phi2, mc_arg->dom_phi12,
      mc_arg->dom_gam1_fwd, mc_arg->dom_gam2_fwd, mc_arg->dom_bond_pay,
      mc_arg->dom_gam1_pay, mc_arg->dom_gam2_pay, mc_arg->covar);

  /*	Fill limit conditions (product) */
  mc_arg->void_prm = (void **)calloc(mc_arg->nb_dates, sizeof(void *));

  if (!mc_arg->void_prm) {
    err = "Memory allocation error (4) in ccf_fill_mc_arg";
    goto FREE_RETURN;
  }

  for (i = mc_arg->nb_dates - 1; i >= flag; i--) {
    ccf_prm = malloc(sizeof(ccf_pay_arg));
    ccf_prm->und = und;
    ccf_prm->ccf = ccf;

    ccf_prm->call_idx = i - flag;

    err = ccf_fill_eval_const(und, ccf, ccf_prm->call_idx,
                              &(ccf_prm->eval_const));

    mc_arg->void_prm[i] = (void *)ccf_prm;
  }

  if (flag) {
    mc_arg->void_prm[0] = NULL;
  }

FREE_RETURN:

  if (err) {
    ccf_free_mc_arg(mc_arg);
  }

  // Free new structures
  if (new_time)
    free(new_time);
  if (new_lambda)
    free(new_lambda);
  if (new_lambda2)
    free(new_lambda2);
  if (new_sigma)
    free(new_sigma);

  return err;
}

Err ccf_free_mc_arg(CCF_MC_ARG mc_arg) {
  int i;
  CCF_PAY_ARG ccf_prm;

  if (mc_arg) {
    if (mc_arg->time)
      free(mc_arg->time);
    if (mc_arg->date)
      free(mc_arg->date);

    if (mc_arg->void_prm) {
      for (i = 0; i < mc_arg->nb_dates; i++) {
        ccf_prm = (CCF_PAY_ARG)(mc_arg->void_prm[i]);
        free(ccf_prm);
      }

      free(mc_arg->void_prm);
    }

    if (mc_arg->dom_fwd1)
      free(mc_arg->dom_fwd1);
    if (mc_arg->dom_fwd2)
      free(mc_arg->dom_fwd2);
    if (mc_arg->dom_exp1)
      free(mc_arg->dom_exp1);
    if (mc_arg->dom_exp2)
      free(mc_arg->dom_exp2);
    if (mc_arg->dom_phi1)
      free(mc_arg->dom_phi1);
    if (mc_arg->dom_phi2)
      free(mc_arg->dom_phi2);
    if (mc_arg->dom_phi12)
      free(mc_arg->dom_phi12);
    if (mc_arg->dom_gam1_fwd)
      free(mc_arg->dom_gam1_fwd);
    if (mc_arg->dom_gam2_fwd)
      free(mc_arg->dom_gam2_fwd);
    if (mc_arg->dom_bond_pay)
      free(mc_arg->dom_bond_pay);
    if (mc_arg->dom_gam1_pay)
      free(mc_arg->dom_gam1_pay);
    if (mc_arg->dom_gam2_pay)
      free(mc_arg->dom_gam2_pay);

    if (mc_arg->covar)
      free_f3tensor(mc_arg->covar, 0, mc_arg->nb_dates - 1, 0, 1, 0, 1);

    mc_arg->time = NULL;
    mc_arg->date = NULL;
    mc_arg->void_prm = NULL;

    mc_arg->dom_fwd1 = NULL;
    mc_arg->dom_fwd2 = NULL;
    mc_arg->dom_exp1 = NULL;
    mc_arg->dom_exp2 = NULL;
    mc_arg->dom_phi1 = NULL;
    mc_arg->dom_phi2 = NULL;
    mc_arg->dom_phi12 = NULL;
    mc_arg->dom_gam1_fwd = NULL;
    mc_arg->dom_gam2_fwd = NULL;
    mc_arg->dom_bond_pay = NULL;
    mc_arg->dom_gam1_pay = NULL;
    mc_arg->dom_gam2_pay = NULL;
    mc_arg->covar = NULL;

    mc_arg->void_prm = NULL;
  }

  return NULL;
}

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

Err ccf_fill_check_all_struct(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund  , 1: calibrate */
    /*		if calib */
    char *yc,          /*	yc */
    char *vc,          /*	vc (only if calib) */
    char *ref,         /*	ref rate (only if calib) */
    char *swap_freq,   /*	swap freq (only if calib) */
    char *swap_basis,  /*	swap basis (only if calib) */
    int lam_ts,        /*	0: use unique lambda  , 1: use ts */
    double lambda,     /*	lambda if unique */
    int tsnlam,        /*	number of lambdas if ts */
    double *tslamtime, /*	lambda times i.e. (date - today) / 365 if ts */
    double *tslam,     /*	corresponding lambdas if ts */
    double alpha,      /*	alpha */
    double gamma,      /*	gamma */
    double rho,        /*	rho */
    /*	End of calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*		if no calilb */
    char *lgm2dund,
    /*	The structure */
    /*		funding */
    double fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_pay, char **fund_basis, double *fund_spr, double *fund_mrg,
    /*		cf */
    double cf_not, int cf_ncpn, long *cf_fix, long *cf_start, long *cf_pay,
    char **cf_basis, char **cf_cms_tenor1, char *cf_cms_freq1,
    char *cf_cms_basis1, double *cf_cms_spread1, char **cf_cms_tenor2,
    char *cf_cms_freq2, char *cf_cms_basis2, double *cf_cms_spread2,
    long spread_vol_n, double *spread_vol_time, double *spread_vol_floor,
    double *spread_vol_cap,

    double *spread_slbeta1, /* Shifted log beta on the CMS1 */
    double *spread_slbeta2, /* Shifted log beta on the CMS2 */

    int cvg_sv, int is_corr, double *cf_alpha, double *cf_beta,
    double *cf_gamma, int *cf_floored, double *cf_floor, int *cf_capped,
    double *cf_cap,

    int cf_nopt,           /* Number of spread options */
    double **cf_notopt,    /* Notional of the spread options */
    double **cf_strikeopt, /* spread option strikes */
    int **cf_typeopt,      /* spread option type 0 call 1 put */

    int use_SL,   /*  1: use Shifted Log
                      0: don't use				*/
    int calib_SL, /*  1: Calibrate the Shifted Log models on CMS1 and on CMS2 *
                                  0: use the Shifted Log beta given */
    double NbStdcalib_SL, /*  Nb of std for the calibration of to the skew */
    int calib_correl_SL,  /*	1: Calibrate the correlation between the two SL
                             to get the same ATM  normal spread vol  0: use the
                             normal spread correl for the sl correl */
    int use_cfoptions, /*  1: Use the spread options and don't take into account
                                               the floor and cap in the cf
                          coupons
                                       0: take into account the floor and cap in
                          the cf coupons */

    int cms_adj, int cms_smile, int cms_vol_adj, /*	1: adjust for CMS vol
                                                    effect 0: don't */
    double cms_beta1, double cms_beta2,
    int num_strikes_in_vol, /*	Array of strikes in vol matrix */
    double *strikes_in_vol,
    SrtDiffusionType
        vol_type, /*	Type of vol in matrix  , SRT_NORMAL or SRT_LOGNORMAL */
    int cash_vol, /*	1: matrix is a cash vol
/*		calls */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee,
    /*	Numerical params */
    int req_stp, int req_stpx, long req_paths,
    /*	Calib params */
    int force_atm, double max_std_long, double max_std_short,
    int fix_lambda,          /*	0: calib lambda to cap  , 1: fix lambda calib
                                                                     to diagonal */
    int cal_vol_shift_type,  /*	vol shift type for volatility */
    double cal_vol_shift,    /*	vol shift */
    double cal_lambda_shift, /*	shift on the calibrated lambda */
    int one_f_equi,          /*	1F equivalent flag:
                                                             if set to 1  , then 2F
                                lambda will calibrate          to the cap priced within calibrated
                                1F          with the given lambda */
    int skip_last,     /*	If 1  , the last option is disregarded and the forward
                          volatility is flat from option n-1 */
    double long_prec,  /*	Precision on primary instruments */
    double short_prec, /*	Precision on secondary instruments */
    double min_fact,   /*	Maximum down jump on variance */
    double max_fact,   /*	Maximum up jump on variance */
    int use_jumps,     /*	Allow vol term structure to jump */
    int proba_weight,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */
    /*	Results */
    CCF_STR ccf, CCF_UND und,
    int *call_feat, /*	0: No callable feature to be valued
                            1: Callable feature to be valued through adi */
    CCF_ADI_ARG adi_arg) {
  Err err = NULL;
  CCF_MC_ARG mc_arg = NULL;
  double *probas = NULL, *exe_bound = NULL;
  double call;

  /*	Initialisation */
  ccf->fund_leg = NULL;
  ccf->cf_leg = NULL;
  ccf->call = NULL;

  und->sigma_n = 0;
  und->spread_vol_n = 0;
  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;
  und->lambda_n = 0;
  und->lambda_date = NULL;
  und->lambda_time = NULL;
  und->lambda = NULL;
  und->spread_vol_time = NULL;
  und->spread_vol_floor = NULL;
  und->spread_vol_cap = NULL;
  cpd_init_calib_inst_data(&(und->inst_data));
  und->has_inst_data = 0;
  und->today = today;

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->sig1 = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Funding leg */
  ccf->fund_leg = (CF_FUND_LEG)malloc(sizeof(cf_fund_leg));
  if (!ccf->fund_leg) {
    err = "Memory allocation error (1) in ccf_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = ccf_fill_fund_leg(und->today, eod_fix_flag, fund_not, fund_ncpn,
                          fund_fix, fund_start, fund_pay, fund_basis, fund_spr,
                          fund_mrg, ccf->fund_leg);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Exotic leg */

  ccf->cf_leg = (CF_EXO_LEG)malloc(sizeof(cf_exo_leg));
  if (!ccf->cf_leg) {
    err = "Memory allocation error (2) in ccf_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = ccf_fill_exo_leg(
      und->today, eod_fix_flag, cf_not, cf_ncpn, cf_fix, cf_start, cf_pay,
      cf_basis, cf_cms_tenor1, cf_cms_freq1, cf_cms_basis1, cf_cms_spread1,
      cf_cms_tenor2, cf_cms_freq2, cf_cms_basis2, cf_cms_spread2, cf_alpha,
      cf_beta, cf_gamma, cf_floored, cf_floor, cf_capped, cf_cap,

      spread_vol_n, spread_vol_time,
      spread_slbeta1, /* Shifted log beta on the CMS1 */
      spread_slbeta2, /* Shifted log beta on the CMS2 */

      cf_nopt, cf_notopt, cf_strikeopt, cf_typeopt, use_SL, calib_SL,
      NbStdcalib_SL, calib_correl_SL, use_cfoptions,

      cms_adj, cms_smile, cms_vol_adj, cms_beta1, cms_beta2, num_strikes_in_vol,
      strikes_in_vol, vol_type, cash_vol, yc, vc, ref, swap_freq, swap_basis,
      get_cash_vol, ccf->cf_leg);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calls */

  if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag) {
    err = ccf_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date,
                         set_date, fee, ccf);
    if (err) {
      goto FREE_RETURN;
    }
  } else {
    ccf->num_calls = 0;
    ccf->call = NULL;
  }

  /*	Underlying */
  if (ccf->num_calls > 0 && ccf->call) {
    if (use_calib) {
      err = ccf_calib_und(
          today, eod_ex_flag, yc, vc, ref, swap_freq, swap_basis, lam_ts,
          lambda, tsnlam, tslamtime, tslam, alpha, gamma, rho, force_atm,
          max_std_long, max_std_short, fix_lambda, cal_vol_shift_type,
          cal_vol_shift, cal_lambda_shift, one_f_equi, skip_last, long_prec,
          short_prec, min_fact, max_fact, use_jumps, 0, 0, probas, exe_bound,
          ccf, get_cash_vol, und, spread_vol_n, spread_vol_time,
          spread_vol_floor, spread_vol_cap, cvg_sv, is_corr);
    } else {
      err = ccf_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, spread_vol_n,
                         spread_vol_time, spread_vol_floor, spread_vol_cap,
                         cvg_sv, is_corr, und);
    }
  }

  if (err) {
    goto FREE_RETURN;
  }

  if (!fix_lambda && proba_weight && ccf->num_calls > 0 && ccf->call &&
      use_calib) {
    mc_arg = calloc(1, sizeof(ccf_mc_arg));

    if (!mc_arg) {
      err = "Memory allocation faillure in ccf_fill_check_all_struct";
      goto FREE_RETURN;
    }

    /* launch MC to get the probas */
    err = ccf_fill_mc_arg(und, ccf, get_cash_vol, req_paths, 1, mc_arg);
    if (err) {
      goto FREE_RETURN;
    }

    probas = calloc(ccf->num_calls, sizeof(double));
    exe_bound = calloc(ccf->num_calls, sizeof(double));

    if (!probas || !exe_bound) {
      err = "Memory allocation faillure in fill_adi_arg";
      goto FREE_RETURN;
    }

    err = ccf_launch_mc(ccf, und, mc_arg, &call, probas, exe_bound);

    if (err) {
      goto FREE_RETURN;
    }

    /* update the first guess */
    lambda = und->lambda[0];

    ccf_free_und(und);

    err = ccf_calib_und(today, eod_ex_flag, yc, vc, ref, swap_freq, swap_basis,
                        lam_ts, lambda, tsnlam, tslamtime, tslam, alpha, gamma,
                        rho, force_atm, max_std_long, max_std_short, fix_lambda,
                        cal_vol_shift_type, cal_vol_shift, cal_lambda_shift,
                        one_f_equi, skip_last, long_prec, short_prec, min_fact,
                        max_fact, use_jumps, proba_weight, 0, probas, exe_bound,
                        ccf, get_cash_vol, und, spread_vol_n, spread_vol_time,
                        spread_vol_floor, spread_vol_cap, cvg_sv, is_corr);
  }

  /*	Tree */

  if (ccf->num_calls > 0 && ccf->call) {
    err = ccf_fill_adi_arg(und, ccf, get_cash_vol, req_stp, req_stpx, adi_arg);

    if (err) {
      goto FREE_RETURN;
    }

    *call_feat = 1;
  } else {
    *call_feat = 0;
  }

FREE_RETURN:

  if (err) {
    ccf_free_all_struct(ccf, und, *call_feat, adi_arg);
  }

  if (probas)
    free(probas);
  if (exe_bound)
    free(exe_bound);

  if (mc_arg) {
    ccf_free_mc_arg(mc_arg);
    free(mc_arg);
  }

  return err;
}

/*	Free all structures */
Err ccf_free_all_struct(CCF_STR ccf, CCF_UND und, int call_feat,
                        CCF_ADI_ARG adi_arg) {
  ccf_free_und(und);
  ccf_free_calls(ccf);

  if (ccf->fund_leg) {
    ccf_free_fund_leg(ccf->fund_leg);
    free(ccf->fund_leg);
  }

  if (ccf->cf_leg) {
    ccf_free_exo_leg(ccf->cf_leg);
    free(ccf->cf_leg);
  }

  ccf_free_adi_arg(adi_arg);

  return NULL;
}

static Err LGM2FCalcPhi(int nstept, double *time, double *sigma,
                        double *sigma_time, int nb_sigma, double *lambda,
                        double *lambda_time, int nb_lambda, double alpha,
                        double gamma, double rho, double *phi1, double *phi2,
                        double *phi12) {
  double t1, t2, ta, tb, dt;
  double alpha2, gamma2;
  double lamtemp, sig1, lam1;

  LGMSVSolFunc Phi1_, Phi2_, Phi12_;

  LGMSVSolFunc *Phi1 = &Phi1_, *Phi2 = &Phi2_, *Phi12 = &Phi12_;

  double vector[20];
  int i, j, nb_sig_minus1;

  double *ts_time = NULL, *sig = NULL, *lam = NULL;

  int nb_ts;

  Err err = NULL;

  memset(vector, 0, 20);

  /* First merge the two term structure */

  ts_time = calloc(nb_sigma, sizeof(double));

  if (!ts_time) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  memcpy(ts_time, sigma_time, nb_sigma * sizeof(double));
  nb_ts = nb_sigma;

  for (i = 0; i < nb_lambda; i++) {
    num_f_add_number(&nb_ts, &ts_time, lambda_time[i]);
  }

  num_f_sort_vector(nb_ts, ts_time);
  num_f_unique_vector(&nb_ts, ts_time);

  sig = calloc(nb_ts, sizeof(double));
  lam = calloc(nb_ts, sizeof(double));

  if (!sig || !lam) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_ts; i++) {
    sig[i] = sigma[Get_Index(ts_time[i], sigma_time, nb_sigma)];
    lam[i] = lambda[Get_Index(ts_time[i], lambda_time, nb_lambda)];
  }

  nb_sig_minus1 = nb_ts - 1;
  alpha2 = alpha * alpha;
  gamma2 = 2.0 * gamma;

  /* Initialisation */
  phi1[0] = 0.0;
  phi2[0] = 0.0;
  phi12[0] = 0.0;

  Phi1->bIsft1 = 1;
  Phi1->b = 0.0;
  Phi1->bIsgt1 = 1;
  Phi1->c = 0.0;
  Phi1->bIsht1 = 1;

  Phi2->bIsft1 = 1;
  Phi2->b = 0.0;
  Phi2->bIsgt1 = 1;
  Phi2->c = 0.0;
  Phi2->bIsht1 = 1;

  Phi12->bIsft1 = 1;
  Phi12->b = 0.0;
  Phi12->bIsgt1 = 1;
  Phi12->c = 0.0;
  Phi12->bIsht1 = 1;

  lam1 = lam[0];
  sig1 = sig[0];

  Phi1->a = sig1 * sig1;
  Phi1->dLambda = 2.0 * lam1;
  Phi2->a = alpha2 * Phi1->a;
  Phi2->dLambda = Phi1->dLambda + gamma2;
  Phi12->a = alpha * Phi1->a;
  Phi12->dLambda = Phi1->dLambda + gamma;

  t1 = 0.0;
  j = 0;

  for (i = 1; i < nstept; i++) {
    t2 = time[i];
    dt = t2 - t1;

    ta = t1;
    tb = ts_time[j];

    lamtemp = 0.0;

    Phi1->dXt1 = phi1[i - 1];
    Phi2->dXt1 = phi2[i - 1];
    Phi12->dXt1 = phi12[i - 1];

    while (tb < t2 && j < nb_sig_minus1) {
      dt = (tb - ta);

      /* Calculation at intermedary time tb */
      LGMSVFuncValue2(Phi1, dt, vector, 0, &Phi1->dXt1);
      LGMSVFuncValue2(Phi2, dt, vector, 0, &Phi2->dXt1);
      LGMSVFuncValue2(Phi12, dt, vector, 0, &Phi12->dXt1);

      lamtemp += lam1 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = lam[j];
      sig1 = sig[j];

      /* Update the model parameters */
      Phi1->a = sig1 * sig1;
      Phi1->dLambda = 2.0 * lam1;
      Phi2->a = alpha2 * Phi1->a;
      Phi2->dLambda = Phi1->dLambda + gamma2;
      Phi12->a = alpha * Phi1->a;
      Phi12->dLambda = Phi1->dLambda + gamma;
    }

    dt = (t2 - ta);

    LGMSVFuncValue2(Phi1, dt, vector, 0, &phi1[i]);
    LGMSVFuncValue2(Phi2, dt, vector, 0, &phi2[i]);
    LGMSVFuncValue2(Phi12, dt, vector, 0, &phi12[i]);

    lamtemp += lam1 * dt;

    t1 = t2;
  }

FREE_RETURN:

  if (sig)
    free(sig);
  if (lam)
    free(lam);
  if (ts_time)
    free(ts_time);

  return err;
}

Err ccf_calc_mdl_iv_fwd(CCF_STR ccf, CCF_UND und, CCF_ADI_ARG adi_arg,

                        int num_hermite,
                        /*	Result */
                        double *premium) {

  int i, j, k, index;
  double *time = NULL, *x = NULL, *w = NULL, *phi1 = NULL, *phi2 = NULL,
         *phi12 = NULL, *r1 = NULL, **r2 = NULL, *r3 = NULL, ***pv = NULL;

  double std1, std2, std3;
  double coef, sum, sum_part, df, fee;

  CF_CALL call;

  Err err = NULL;

  /* memory allocation */

  time = (double *)calloc(ccf->num_calls + 1, sizeof(double));

  x = dvector(1, num_hermite);
  w = dvector(1, num_hermite);

  phi1 = (double *)calloc(ccf->num_calls + 1, sizeof(double));
  phi2 = (double *)calloc(ccf->num_calls + 1, sizeof(double));
  phi12 = (double *)calloc(ccf->num_calls + 1, sizeof(double));

  r1 = dvector(0, num_hermite - 1);
  r2 = dmatrix(0, num_hermite - 1, 0, num_hermite - 1);
  r3 = dvector(0, num_hermite - 1);

  pv = f3tensor(0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

  if (!time || !x || !w || !phi1 || !phi2 || !phi12 || !r1 || !r2 || !r3 ||
      !pv) {
    err = "Memory allocation failure in ccf_calc_mdl_iv_fwd";
    goto FREE_RETURN;
  }

  /* Hermite calculation */
  err = HermiteStandard(x, w, num_hermite);

  if (err) {
    goto FREE_RETURN;
  }

  time[0] = 0.0;
  for (i = 0; i < ccf->num_calls; i++) {
    call = ccf->call + i;
    time[i + 1] = call->ex_time;
  }

  /* Phi Calculation */
  err = LGM2FCalcPhi(ccf->num_calls + 1, time, und->sigma, und->sigma_time,
                     und->sigma_n, und->lambda, und->lambda_time, und->lambda_n,
                     und->alpha, und->gamma, und->rho, phi1, phi2, phi12);

  if (err) {
    goto FREE_RETURN;
  }

  index = 0;

  for (i = 0; i < ccf->num_calls; i++) {
    /* first remove the fee */
    call = ccf->call + i;
    fee = call->fee;
    call->fee = 0.0;

    std1 = sqrt(phi1[i + 1]);
    std2 = sqrt(phi2[i + 1]);

    if (fabs(time[i + 1]) < 1.0E-08) {
      coef = 0.0;
    } else {
      coef = -und->rho * phi12[i + 1] / phi1[i + 1];
    }

    /* r3 = r2 + coef * r1 */
    std3 = sqrt(phi2[i + 1] + coef * coef * phi1[i + 1] +
                2.0 * und->rho * coef * std1 * std2);

    /* construct the grid */
    for (j = 0; j < num_hermite; j++) {
      r1[j] = std1 * x[j + 1];
      r3[j] = std3 * x[j + 1];
    }

    /* reconstruction */
    for (j = 0; j < num_hermite; j++) {
      for (k = 0; k < num_hermite; k++) {
        r2[j][k] = (r3[k] - coef * r1[j]);
        pv[j][k][0] = -BIG_BIG_NUMBER;
      }
    }

    /* launch evaluation */

    while (!((adi_arg->void_prm)[index]) && index < 10000) {
      index++;
    }

    if (ccf->cf_leg->cpn->use_SL || ccf->cf_leg->cpn->use_cfoptions) {
      err = ccf_payoff_4_3dfx_adi_sl(
          call->ex_date, call->ex_time, (adi_arg->void_prm)[index], und->yc,
          und->lambda, und->lambda_time, und->lambda_n, und->gamma, und->rho,
          phi1[i + 1], phi2[i + 1], phi12[i + 1], 0, num_hermite - 1, 0,
          num_hermite - 1, r1, r2, 1, pv);
    } else {
      err = ccf_payoff_4_3dfx_adi(
          call->ex_date, call->ex_time, (adi_arg->void_prm)[index], und->yc,
          und->lambda, und->lambda_time, und->lambda_n, und->gamma, und->rho,
          phi1[i + 1], phi2[i + 1], phi12[i + 1], 0, num_hermite - 1, 0,
          num_hermite - 1, r1, r2, 1, pv);
    }

    index++;

    /* put the fee back */
    call->fee = fee;

    /* convolution */

    sum = 0.0;
    sum_part = 0.0;

    for (j = 0; j < num_hermite; j++) {
      sum_part = 0.0;
      for (k = 0; k < num_hermite; k++) {
        sum_part += w[k + 1] * pv[j][k][0];
      }

      sum = sum + w[j + 1] * sum_part;
    }

    df = swp_f_df(und->today, call->ex_date, und->yc);

    sum *= df;

    premium[i] = sum;
  }

FREE_RETURN:

  if (time)
    free(time);

  if (x)
    free_dvector(x, 1, num_hermite);
  if (w)
    free_dvector(w, 1, num_hermite);

  if (phi1)
    free(phi1);
  if (phi2)
    free(phi2);
  if (phi12)
    free(phi12);

  if (r1)
    free_dvector(r1, 0, num_hermite - 1);
  if (r2)
    free_dmatrix(r2, 0, num_hermite - 1, 0, num_hermite - 1);
  if (r3)
    free_dvector(r3, 0, num_hermite - 1);

  if (pv)
    free_f3tensor(pv, 0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

  return err;
}

/*	Payoff function */
/*	---------------	*/

static double cms_rate(double delay, double nfp, double compd, double forward,
                       double nvar) {
  static double lvl, dlvl, cms_adj, del_adj, tot_adj;

  if (fabs(forward) > 1.0E-08) {
    lvl = (1.0 - pow(1.0 + forward / compd, -nfp)) / (forward / compd);
    dlvl =
        (nfp / compd) * pow(1.0 + forward / compd, -nfp - 1) * forward / compd -
        (1.0 - pow(1.0 + forward / compd, -nfp)) / compd;
    dlvl /= (forward / compd) * (forward / compd);
  } else {
    lvl = nfp;
    dlvl = -0.5 * nfp * (nfp + 1.0) / compd;
  }

  /*	Adjustments	*/
  cms_adj = -dlvl / lvl;
  del_adj = -delay / (1.0 + forward * delay);
  tot_adj = cms_adj + del_adj;

  return (forward + tot_adj * nvar);
}

static double dens(double x) { return INV_SQRT_TWO_PI * exp(-x * x / 2.0); }

#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL_N(FWD, STRIKE, STD, D)                                        \
  (((FWD) - (STRIKE)) * norm((D)) + (STD)*dens((D)))

#define PUT_VAL_N(FWD, STRIKE, STD, D)                                         \
  (((STRIKE) - (FWD)) * norm(-(D)) + (STD)*dens((D)))

#define OPT_VAL_MACRO_N(TYPE, FWD, STRIKE, STD)                                \
  ((TYPE) == 0                                                                 \
       ? 0.0                                                                   \
       : ((TYPE) == 1                                                          \
              ? POS_VAL((FWD) - (STRIKE))                                      \
              : ((TYPE) == 2                                                   \
                     ? POS_VAL((STRIKE) - (FWD))                               \
                     : ((TYPE) == 3                                            \
                            ? CALL_VAL_N((FWD), (STRIKE), (STD),               \
                                         ((FWD) - (STRIKE)) / (STD))           \
                            : PUT_VAL_N((FWD), (STRIKE), (STD),                \
                                        ((FWD) - (STRIKE)) / (STD))))))

Err ccf_payoff_4_3dfx_adi(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double *lam, double *ts_time, int nb_ts, double gamma, double rho,
    double phi1, double phi2, double phi12,
    /* Nodes data */
    int l1, int u1, int l2, int u2, double *r1, double **r2, int nprod,
    /* Vector of results to be updated */
    double ***prod_val) {
  CCF_PAY_ARG ccf_arg;
  CCF_STR ccf;
  CF_CALL call;
  CCF_UND und;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  CCF_EVAL_CONST eval_const;
  int call_idx;

  int i, j, k, l, m;

  double R1, R2, *r2i;

  double fund_leg_pv, df, df_start, df_end, level, forward, cmsrate[2], spread,
      coupon, floor, cap, exo_leg_pv, iv;

  int num_fund_cpn, num_cf_cpn;

  double mat;

  double fee;

  Err err = NULL;

  /*	Get the event */

  ccf_arg = (CCF_PAY_ARG)func_parm;
  eval_const = (CCF_EVAL_CONST)(&(ccf_arg->eval_const));
  ccf = ccf_arg->ccf;
  call_idx = ccf_arg->call_idx;
  call = ccf->call + call_idx;
  und = (CCF_UND)(ccf_arg->und);
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  num_fund_cpn = call->num_fund_cpn;
  num_cf_cpn = call->num_cf_cpn;

  /*	Eval payoff
          ----------- */

  for (i = l1; i <= u1; i++) {
    R1 = r1[i];
    r2i = r2[i];
    for (j = l2; j <= u2; j++) {
      R2 = r2i[j];

      /*	Calculate discount factors	*/
      for (l = 0; l < eval_const->num_df; l++) {
        eval_const->df[l] =
            exp(-eval_const->df_alpha[l] - eval_const->df_beta[l] * R1 -
                eval_const->df_gamma[l] * R2);
      }

      /*	PV of funding leg */
      fund_leg_pv = 0.0;

      /*	DF initial */

      fund_leg_pv += eval_const->df[eval_const->start_idx] * fund_leg->notional;

      /*	Fund Spread Coupons */

      for (l = 0; l < num_fund_cpn; l++) {
        fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
        fund_leg_pv += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cpn;
      }

      /*	Notional */
      fund_leg_pv -= eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]] *
                     fund_leg->notional;

      /*	PV of CMS leg */
      exo_leg_pv = 0.0;

      /*	Coupons */
      for (l = 0; l < num_cf_cpn; l++) {
        /*	Coupon access */
        exo_cpn = exo_leg->cpn + call->cf_idx + l;

        /*	Discount */
        df = eval_const->df[eval_const->cf_disc_idx[l]];

        /*	CMS evaluation */
        for (m = 0; m < exo_cpn->ncms; m++) {
          cms = &(exo_cpn->cms[m]);

          /*	First we evaluate the Forward: */
          /*	DF initial */
          df_start = eval_const->df[eval_const->cf_fwd_start_idx[m][l]];

          /*	Level minus the last one */
          level = 0.0;
          for (k = 0; k < cms->num_cpn - 1; k++) {
            level +=
                eval_const->df[eval_const->cf_fwd_idx[m][l][k]] * cms->cvg[k];
          }

          /*	DF final */
          df_end =
              eval_const->df[eval_const->cf_fwd_idx[m][l][cms->num_cpn - 1]];
          level += df_end * cms->cvg[cms->num_cpn - 1];
          forward = (df_start - df_end) / level + cms->swap_spread;

          /*	CMS adjustment */
          if (eval_const->cf_cms_vol[m][l] > 1.0E-08) {
            switch (cms->adj_type) {
            case 0:
            default:
              cmsrate[m] = forward;
              break;

            case 1:
              cmsrate[m] = cms_rate(cms->delay, cms->nfp, cms->cpnd, forward,
                                    eval_const->cf_cms_vol[m][l]);
              break;

            case 2:
              mat = exo_cpn->cms_fix_time - call->ex_time;
              if (cms->beta > 1.0e-08) {
                cmsrate[m] = forward + cms->lambda *
                                           pow(fabs(forward), 2.0 * cms->beta) *
                                           mat;
              } else {
                cmsrate[m] = forward + cms->lambda * mat;
              }
              break;
            }
          } else {
            cmsrate[m] = forward;
          }
        }

        switch (exo_cpn->type) {
        case 0:
        default:
          /*	Midat */
          spread = 0.0;
          coupon = exo_cpn->gamma;
          floor = cap = 0.0;
          break;

        case 1:
        case 2:
          /*	CIF / CMS */
          spread = exo_cpn->alphabeta[0] * cmsrate[0];
          coupon = spread + exo_cpn->gamma;
          floor = OPT_VAL_MACRO_N(eval_const->cf_spr_floor_type[l], spread,
                                  eval_const->cf_spr_floor_str[l],
                                  eval_const->cf_spr_floor_std[l]);
          cap = OPT_VAL_MACRO_N(eval_const->cf_spr_cap_type[l], spread,
                                eval_const->cf_spr_cap_str[l],
                                eval_const->cf_spr_cap_std[l]);
          break;

        case 3:
          /*	CMS Spread */
          spread = exo_cpn->alphabeta[0] * cmsrate[0] +
                   exo_cpn->alphabeta[1] * cmsrate[1];
          coupon = spread + exo_cpn->gamma;
          floor = OPT_VAL_MACRO_N(eval_const->cf_spr_floor_type[l], spread,
                                  eval_const->cf_spr_floor_str[l],
                                  eval_const->cf_spr_floor_std[l]);
          cap = OPT_VAL_MACRO_N(eval_const->cf_spr_cap_type[l], spread,
                                eval_const->cf_spr_cap_str[l],
                                eval_const->cf_spr_cap_std[l]);
          break;
        }

        /*	Coupon pv */
        exo_leg_pv += df * (coupon + floor - cap) * exo_cpn->cvg;
      }

      /*	Intrinsic value */
      if (call->pay_rec == 0) {
        iv = exo_leg_pv - fund_leg_pv;
      } else {
        iv = fund_leg_pv - exo_leg_pv;
      }

      /*	Fee */
      if (fabs(call->fee) > 1.0e-08) {
        fee = eval_const->df[eval_const->fee_idx] * call->fee;
      } else {
        fee = 0.0;
      }

      /*	Process max */

      if (prod_val[i][j][0] < iv - fee) {
        prod_val[i][j][0] = iv - fee;
      }
    }
  }

  /*	End of payoff valuation
          ----------------------- */

  return err;
}

Err ccf_payoff_4_3dfx_adi_sl(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double *lam, double *ts_time, int nb_ts, double gamma, double rho,
    double phi1, double phi2, double phi12,
    /* Nodes data */
    int l1, int u1, int l2, int u2, double *r1, double **r2, int nprod,
    /* Vector of results to be updated */
    double ***prod_val) {
  CCF_PAY_ARG ccf_arg;
  CCF_STR ccf;
  CF_CALL call;
  CCF_UND und;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  CCF_EVAL_CONST eval_const;
  int call_idx;

  int i, j, k, l, m;

  double R1, R2, *r2i;

  double fund_leg_pv, df, df_start, df_end, level, forward, cmsrate[2], spread,
      coupon, floor, cap, exo_leg_pv, iv;

  int num_fund_cpn, num_cf_cpn, iNumOpt;

  double optionvalue;

  double mat;

  double fee;

  Err err = NULL;

  /*	Get the event */

  ccf_arg = (CCF_PAY_ARG)func_parm;
  eval_const = (CCF_EVAL_CONST)(&(ccf_arg->eval_const));
  ccf = ccf_arg->ccf;
  call_idx = ccf_arg->call_idx;
  call = ccf->call + call_idx;
  und = (CCF_UND)(ccf_arg->und);
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  num_fund_cpn = call->num_fund_cpn;
  num_cf_cpn = call->num_cf_cpn;

  /*	Eval payoff
          ----------- */

  for (i = l1; i <= u1; i++) {
    R1 = r1[i];
    r2i = r2[i];
    for (j = l2; j <= u2; j++) {
      R2 = r2i[j];

      /*	Calculate discount factors	*/
      for (l = 0; l < eval_const->num_df; l++) {
        eval_const->df[l] =
            exp(-eval_const->df_alpha[l] - eval_const->df_beta[l] * R1 -
                eval_const->df_gamma[l] * R2);
      }

      /*	PV of funding leg */
      fund_leg_pv = 0.0;

      /*	DF initial */

      fund_leg_pv += eval_const->df[eval_const->start_idx] * fund_leg->notional;

      /*	Fund Spread Coupons */

      for (l = 0; l < num_fund_cpn; l++) {
        fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
        fund_leg_pv += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cpn;
      }

      /*	Notional */
      fund_leg_pv -= eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]] *
                     fund_leg->notional;

      /*	PV of CMS leg */
      exo_leg_pv = 0.0;

      /*	Coupons */
      for (l = 0; l < num_cf_cpn; l++) {
        /*	Coupon access */
        exo_cpn = exo_leg->cpn + call->cf_idx + l;

        /*	Discount */
        df = eval_const->df[eval_const->cf_disc_idx[l]];

        /*	CMS evaluation */
        for (m = 0; m < exo_cpn->ncms; m++) {
          cms = &(exo_cpn->cms[m]);

          /*	First we evaluate the Forward: */
          /*	DF initial */
          df_start = eval_const->df[eval_const->cf_fwd_start_idx[m][l]];

          /*	Level minus the last one */
          level = 0.0;
          for (k = 0; k < cms->num_cpn - 1; k++) {
            level +=
                eval_const->df[eval_const->cf_fwd_idx[m][l][k]] * cms->cvg[k];
          }

          /*	DF final */
          df_end =
              eval_const->df[eval_const->cf_fwd_idx[m][l][cms->num_cpn - 1]];
          level += df_end * cms->cvg[cms->num_cpn - 1];
          forward = (df_start - df_end) / level + cms->swap_spread;

          /*	CMS adjustment */
          if (eval_const->cf_cms_vol[m][l] > 1.0E-08) {
            switch (cms->adj_type) {
            case 0:
            default:
              cmsrate[m] = forward;
              break;

            case 1:
              cmsrate[m] = cms_rate(cms->delay, cms->nfp, cms->cpnd, forward,
                                    eval_const->cf_cms_vol[m][l]);
              break;

            case 2:
              mat = exo_cpn->cms_fix_time - call->ex_time;

              /* SL adjustement */
              if (exo_cpn->use_SL) {
                cmsrate[m] = forward + cms->lambda * (forward + cms->slshift) *
                                           (forward + cms->slshift) * mat;
              } else {
                /* Beta adjustement */
                if (cms->beta > 1.0e-08) {
                  cmsrate[m] =
                      forward +
                      cms->lambda * pow(fabs(forward), 2.0 * cms->beta) * mat;
                } else {
                  cmsrate[m] = forward + cms->lambda * mat;
                }
              }
              break;
            }
          } else {
            cmsrate[m] = forward;
          }
        }

        /*	CMS Spread */
        if (!exo_cpn->use_cfoptions) {
          /* Use the Shifted log model */
          spread = exo_cpn->alphabeta[0] * cmsrate[0] +
                   exo_cpn->alphabeta[1] * cmsrate[1];
          coupon = spread + exo_cpn->gamma;
          mat = exo_cpn->cms_fix_time - call->ex_time;

          /* Evaluation of the floor */
          if (eval_const->cf_spr_floor_type[l] == 0) {
            floor = 0.0;
          } else {
            err = OptSpreadNew(
                cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
                exo_cpn->alphabeta[0],                /*  nx  ,  */
                (eval_const->cf_spr_floor_type[l] > 2 ? exo_cpn->cms[0].slvol
                                                      : 0.0), /* sigx  , */
                cmsrate[1] + exo_cpn->cms[1].slshift,         /* fwdy  , */
                exo_cpn->alphabeta[1],                        /* ny  ,   */
                (eval_const->cf_spr_floor_type[l] > 2 ? exo_cpn->cms[1].slvol
                                                      : 0.0), /* sigy  , */
                eval_const->cf_spr_floor_str[l] +
                    exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                    exo_cpn->cms[1].slshift *
                        exo_cpn->alphabeta[1], /* K  ,	 */
                mat,                           /* mat  ,  */
                eval_const->cf_spr_slcorr[l],  /* rho  ,  */
                SRT_PUT, &floor);
          }

          /* Evaluation of the cap */
          if (eval_const->cf_spr_cap_type[l] == 0) {
            cap = 0.0;
          } else {
            err = OptSpreadNew(
                cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
                exo_cpn->alphabeta[0],                /*  nx  ,  */
                (eval_const->cf_spr_cap_type[l] > 2 ? exo_cpn->cms[0].slvol
                                                    : 0.0), /* sigx  , */
                cmsrate[1] + exo_cpn->cms[1].slshift,       /* fwdy  , */
                exo_cpn->alphabeta[1],                      /* ny  ,   */
                (eval_const->cf_spr_cap_type[l] > 2 ? exo_cpn->cms[1].slvol
                                                    : 0.0), /* sigy  , */
                eval_const->cf_spr_cap_str[l] +
                    exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                    exo_cpn->cms[1].slshift *
                        exo_cpn->alphabeta[1], /* K  ,	 */

                mat,                          /* mat  ,  */
                eval_const->cf_spr_slcorr[l], /* rho  ,  */
                SRT_CALL, &cap);
          }

        } else {
          /* use of a string of options */
          spread = exo_cpn->alphabeta[0] * cmsrate[0] +
                   exo_cpn->alphabeta[1] * cmsrate[1];
          mat = exo_cpn->cms_fix_time - call->ex_time;

          /* Put all the value of the options in the floor variable */
          coupon = exo_cpn->gamma;
          floor = 0;
          cap = 0;

          if (!exo_cpn->use_SL) {
            for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
              if (exo_cpn->notopt[iNumOpt] != 0) {
                floor +=
                    exo_cpn->notopt[iNumOpt] *
                    OPT_VAL_MACRO_N(
                        (exo_cpn->typeopt[iNumOpt] == SRT_CALL
                             ? (eval_const->cf_spr_opt_std[l][iNumOpt] == 0 ? 1
                                                                            : 3)
                             : (eval_const->cf_spr_opt_std[l][iNumOpt] == 0
                                    ? 2
                                    : 4)),
                        spread, exo_cpn->strikeopt[iNumOpt],
                        eval_const->cf_spr_opt_std[l][iNumOpt]);
              }
            }
          } else {

            /* Use the Shifted log model */
            for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
              if (exo_cpn->notopt[iNumOpt] != 0) {
                err = OptSpreadNew(
                    cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
                    exo_cpn->alphabeta[0],                /*  nx  ,  */
                    exo_cpn->cms[0].slvol,                /* sigx  , */
                    cmsrate[1] + exo_cpn->cms[1].slshift, /* fwdy  , */
                    exo_cpn->alphabeta[1],                /* ny  ,   */
                    exo_cpn->cms[1].slvol,                /* sigy  , */
                    exo_cpn->strikeopt[iNumOpt] +
                        exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                        exo_cpn->cms[1].slshift *
                            exo_cpn->alphabeta[1], /* K  ,	 */
                    mat,                           /* mat  ,  */
                    eval_const->cf_spr_slcorr[l],  /* rho  ,  */
                    exo_cpn->typeopt[iNumOpt], &optionvalue);

                floor += exo_cpn->notopt[iNumOpt] * optionvalue;
              }
            }
          }
        }

        /*	Coupon pv */
        exo_leg_pv += df * (coupon + floor - cap) * exo_cpn->cvg;
      }

      /*	Intrinsic value */
      if (call->pay_rec == 0) {
        iv = exo_leg_pv - fund_leg_pv;
      } else {
        iv = fund_leg_pv - exo_leg_pv;
      }

      /*	Fee */
      if (fabs(call->fee) > 1.0e-08) {
        fee = eval_const->df[eval_const->fee_idx] * call->fee;
      } else {
        fee = 0.0;
      }

      /*	Process max */

      if (prod_val[i][j][0] < iv - fee) {
        prod_val[i][j][0] = iv - fee;
      }
    }
  }

  /*	End of payoff valuation
          ----------------------- */

  return err;
}

Err ccf_payoff_4_3dfx_mc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double R1, double R2, int nprod,
    /* Vector of results to be updated */
    double *res, int *stop_path) {
  CCF_PAY_ARG ccf_arg;
  CCF_STR ccf;
  CF_CALL call;
  CCF_UND und;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  CCF_EVAL_CONST eval_const;
  int call_idx;

  int k, l, m;

  double fund_leg_pv, df, df_start, df_end, level, forward, cmsrate[2], spread,
      coupon, floor, cap, exo_leg_pv, iv;

  int num_fund_cpn, num_cf_cpn;

  double mat;

  double fee;

  Err err = NULL;

  /*	Get the event */

  ccf_arg = (CCF_PAY_ARG)func_parm;
  eval_const = (CCF_EVAL_CONST)(&(ccf_arg->eval_const));
  ccf = ccf_arg->ccf;
  call_idx = ccf_arg->call_idx;
  call = ccf->call + call_idx;
  und = (CCF_UND)(ccf_arg->und);
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  num_fund_cpn = call->num_fund_cpn;
  num_cf_cpn = call->num_cf_cpn;

  /*	Eval payoff
          ----------- */

  /*	Calculate discount factors	*/
  for (l = 0; l < eval_const->num_df; l++) {
    eval_const->df[l] =
        exp(-eval_const->df_alpha[l] - eval_const->df_beta[l] * R1 -
            eval_const->df_gamma[l] * R2);
  }

  /*	PV of funding leg */
  fund_leg_pv = 0.0;

  /*	DF initial */

  fund_leg_pv += eval_const->df[eval_const->start_idx] * fund_leg->notional;

  /*	Fund Spread Coupons */

  for (l = 0; l < num_fund_cpn; l++) {
    fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
    fund_leg_pv += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cpn;
  }

  /*	Notional */
  fund_leg_pv -= eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]] *
                 fund_leg->notional;

  /*	PV of CMS leg */
  exo_leg_pv = 0.0;

  /*	Coupons */
  for (l = 0; l < num_cf_cpn; l++) {
    /*	Coupon access */
    exo_cpn = exo_leg->cpn + call->cf_idx + l;

    /*	Discount */
    df = eval_const->df[eval_const->cf_disc_idx[l]];

    /*	CMS evaluation */
    for (m = 0; m < exo_cpn->ncms; m++) {
      cms = &(exo_cpn->cms[m]);

      /*	First we evaluate the Forward: */
      /*	DF initial */
      df_start = eval_const->df[eval_const->cf_fwd_start_idx[m][l]];

      /*	Level minus the last one */
      level = 0.0;
      for (k = 0; k < cms->num_cpn - 1; k++) {
        level += eval_const->df[eval_const->cf_fwd_idx[m][l][k]] * cms->cvg[k];
      }

      /*	DF final */
      df_end = eval_const->df[eval_const->cf_fwd_idx[m][l][cms->num_cpn - 1]];
      level += df_end * cms->cvg[cms->num_cpn - 1];
      forward = (df_start - df_end) / level + cms->swap_spread;

      /*	CMS adjustment */
      if (eval_const->cf_cms_vol[m][l] > 1.0E-08) {
        switch (cms->adj_type) {
        case 0:
        default:
          cmsrate[m] = forward;
          break;

        case 1:
          cmsrate[m] = cms_rate(cms->delay, cms->nfp, cms->cpnd, forward,
                                eval_const->cf_cms_vol[m][l]);
          break;

        case 2:
          mat = exo_cpn->cms_fix_time - call->ex_time;
          if (cms->beta > 1.0e-08) {
            cmsrate[m] = forward + cms->lambda *
                                       pow(fabs(forward), 2.0 * cms->beta) *
                                       mat;
          } else {
            cmsrate[m] = forward + cms->lambda * mat;
          }
          break;
        }
      } else {
        cmsrate[m] = forward;
      }
    }

    switch (exo_cpn->type) {
    case 0:
    default:
      /*	Midat */
      spread = 0.0;
      coupon = exo_cpn->gamma;
      floor = cap = 0.0;
      break;

    case 1:
    case 2:
      /*	CIF / CMS */
      spread = exo_cpn->alphabeta[0] * cmsrate[0];
      coupon = spread + exo_cpn->gamma;
      floor = OPT_VAL_MACRO_N(eval_const->cf_spr_floor_type[l], spread,
                              eval_const->cf_spr_floor_str[l],
                              eval_const->cf_spr_floor_std[l]);
      cap = OPT_VAL_MACRO_N(eval_const->cf_spr_cap_type[l], spread,
                            eval_const->cf_spr_cap_str[l],
                            eval_const->cf_spr_cap_std[l]);
      break;

    case 3:
      /*	CMS Spread */
      spread = exo_cpn->alphabeta[0] * cmsrate[0] +
               exo_cpn->alphabeta[1] * cmsrate[1];
      coupon = spread + exo_cpn->gamma;
      floor = OPT_VAL_MACRO_N(eval_const->cf_spr_floor_type[l], spread,
                              eval_const->cf_spr_floor_str[l],
                              eval_const->cf_spr_floor_std[l]);
      cap = OPT_VAL_MACRO_N(eval_const->cf_spr_cap_type[l], spread,
                            eval_const->cf_spr_cap_str[l],
                            eval_const->cf_spr_cap_std[l]);
      break;
    }

    /*	Coupon pv */
    exo_leg_pv += df * (coupon + floor - cap) * exo_cpn->cvg;
  }

  /*	Intrinsic value */
  if (call->pay_rec == 0) {
    iv = exo_leg_pv - fund_leg_pv;
  } else {
    iv = fund_leg_pv - exo_leg_pv;
  }

  /*	Fee */
  if (fabs(call->fee) > 1.0e-08) {
    fee = eval_const->df[eval_const->fee_idx] * call->fee;
  } else {
    fee = 0.0;
  }

  /*	Fill Result */
  res[0] = iv - fee;

  /* Calculate the Co Terminal Swap */
  /*	DF initial */
  df_start = eval_const->df[eval_const->start_idx];

  /*	Level minus the last one */
  level = 0.0;
  for (l = 0; l < num_fund_cpn; l++) {
    fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
    level += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cvg;
  }

  /*	DF final */
  df_end = eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]];
  forward = (df_start - df_end) / level;

  if (call->pay_rec == 0) {
    res[1] = -forward;
  } else {
    res[1] = forward;
  }

  /*	End of payoff valuation
          ----------------------- */

  return err;
}

Err ccf_payoff_4_3dfx_mc_sl(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double R1, double R2, int nprod,
    /* Vector of results to be updated */
    double *res, int *stop_path) {
  CCF_PAY_ARG ccf_arg;
  CCF_STR ccf;
  CF_CALL call;
  CCF_UND und;
  CF_FUND_LEG fund_leg;
  CF_FUND_CPN fund_cpn;
  CF_EXO_LEG exo_leg;
  CF_EXO_CPN exo_cpn;
  CF_CMS_DESC cms;
  CCF_EVAL_CONST eval_const;
  int call_idx;

  int k, l, m;

  double fund_leg_pv, df, df_start, df_end, level, forward, cmsrate[2], spread,
      coupon, floor, cap, exo_leg_pv, iv;

  int num_fund_cpn, num_cf_cpn, iNumOpt;

  double optionvalue;

  double mat;

  double fee;

  Err err = NULL;

  /*	Get the event */

  ccf_arg = (CCF_PAY_ARG)func_parm;
  eval_const = (CCF_EVAL_CONST)(&(ccf_arg->eval_const));
  ccf = ccf_arg->ccf;
  call_idx = ccf_arg->call_idx;
  call = ccf->call + call_idx;
  und = (CCF_UND)(ccf_arg->und);
  fund_leg = ccf->fund_leg;
  exo_leg = ccf->cf_leg;

  num_fund_cpn = call->num_fund_cpn;
  num_cf_cpn = call->num_cf_cpn;

  /*	Eval payoff
          ----------- */

  /*	Calculate discount factors	*/
  for (l = 0; l < eval_const->num_df; l++) {
    eval_const->df[l] =
        exp(-eval_const->df_alpha[l] - eval_const->df_beta[l] * R1 -
            eval_const->df_gamma[l] * R2);
  }

  /*	PV of funding leg */
  fund_leg_pv = 0.0;

  /*	DF initial */

  fund_leg_pv += eval_const->df[eval_const->start_idx] * fund_leg->notional;

  /*	Fund Spread Coupons */

  for (l = 0; l < num_fund_cpn; l++) {
    fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
    fund_leg_pv += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cpn;
  }

  /*	Notional */
  fund_leg_pv -= eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]] *
                 fund_leg->notional;

  /*	PV of CMS leg */
  exo_leg_pv = 0.0;

  /*	Coupons */
  for (l = 0; l < num_cf_cpn; l++) {
    /*	Coupon access */
    exo_cpn = exo_leg->cpn + call->cf_idx + l;

    /*	Discount */
    df = eval_const->df[eval_const->cf_disc_idx[l]];

    /*	CMS evaluation */
    for (m = 0; m < exo_cpn->ncms; m++) {
      cms = &(exo_cpn->cms[m]);

      /*	First we evaluate the Forward: */
      /*	DF initial */
      df_start = eval_const->df[eval_const->cf_fwd_start_idx[m][l]];

      /*	Level minus the last one */
      level = 0.0;
      for (k = 0; k < cms->num_cpn - 1; k++) {
        level += eval_const->df[eval_const->cf_fwd_idx[m][l][k]] * cms->cvg[k];
      }

      /*	DF final */
      df_end = eval_const->df[eval_const->cf_fwd_idx[m][l][cms->num_cpn - 1]];
      level += df_end * cms->cvg[cms->num_cpn - 1];
      forward = (df_start - df_end) / level + cms->swap_spread;

      /*	CMS adjustment */
      if (eval_const->cf_cms_vol[m][l] > 1.0E-08) {
        switch (cms->adj_type) {
        case 0:
        default:
          cmsrate[m] = forward;
          break;

        case 1:
          cmsrate[m] = cms_rate(cms->delay, cms->nfp, cms->cpnd, forward,
                                eval_const->cf_cms_vol[m][l]);
          break;

        case 2:
          mat = exo_cpn->cms_fix_time - call->ex_time;

          /* SL adjustement */
          if (exo_cpn->use_SL) {
            cmsrate[m] = forward + cms->lambda * (forward + cms->slshift) *
                                       (forward + cms->slshift) * mat;
          } else {
            /* Beta adjustement */
            if (cms->beta > 1.0e-08) {
              cmsrate[m] = forward + cms->lambda *
                                         pow(fabs(forward), 2.0 * cms->beta) *
                                         mat;
            } else {
              cmsrate[m] = forward + cms->lambda * mat;
            }
          }
          break;
        }
      } else {
        cmsrate[m] = forward;
      }
    }

    /*	CMS Spread */
    if (!exo_cpn->use_cfoptions) {
      /* Use the Shifted log model */
      spread = exo_cpn->alphabeta[0] * cmsrate[0] +
               exo_cpn->alphabeta[1] * cmsrate[1];
      coupon = spread + exo_cpn->gamma;
      mat = exo_cpn->cms_fix_time - call->ex_time;

      /* Evaluation of the floor */
      if (eval_const->cf_spr_floor_type[l] == 0) {
        floor = 0.0;
      } else {
        err = OptSpreadNew(
            cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
            exo_cpn->alphabeta[0],                /*  nx  ,  */
            (eval_const->cf_spr_floor_type[l] > 2 ? exo_cpn->cms[0].slvol
                                                  : 0.0), /* sigx  , */
            cmsrate[1] + exo_cpn->cms[1].slshift,         /* fwdy  , */
            exo_cpn->alphabeta[1],                        /* ny  ,   */
            (eval_const->cf_spr_floor_type[l] > 2 ? exo_cpn->cms[1].slvol
                                                  : 0.0), /* sigy  , */
            eval_const->cf_spr_floor_str[l] +
                exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                exo_cpn->cms[1].slshift *
                    exo_cpn->alphabeta[1], /* K  ,	 */
            mat,                           /* mat  ,  */
            eval_const->cf_spr_slcorr[l],  /* rho  ,  */
            SRT_PUT, &floor);
      }

      /* Evaluation of the cap */
      if (eval_const->cf_spr_cap_type[l] == 0) {
        cap = 0.0;
      } else {
        err = OptSpreadNew(
            cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
            exo_cpn->alphabeta[0],                /*  nx  ,  */
            (eval_const->cf_spr_cap_type[l] > 2 ? exo_cpn->cms[0].slvol
                                                : 0.0), /* sigx  , */
            cmsrate[1] + exo_cpn->cms[1].slshift,       /* fwdy  , */
            exo_cpn->alphabeta[1],                      /* ny  ,   */
            (eval_const->cf_spr_cap_type[l] > 2 ? exo_cpn->cms[1].slvol
                                                : 0.0), /* sigy  , */
            eval_const->cf_spr_cap_str[l] +
                exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                exo_cpn->cms[1].slshift *
                    exo_cpn->alphabeta[1], /* K  ,	 */

            mat,                          /* mat  ,  */
            eval_const->cf_spr_slcorr[l], /* rho  ,  */
            SRT_CALL, &cap);
      }

    } else {
      /* use of a string of options */
      spread = exo_cpn->alphabeta[0] * cmsrate[0] +
               exo_cpn->alphabeta[1] * cmsrate[1];
      mat = exo_cpn->cms_fix_time - call->ex_time;

      /* Put all the value of the options in the floor variable */
      coupon = exo_cpn->gamma;
      floor = 0;
      cap = 0;

      if (!exo_cpn->use_SL) {
        for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
          if (exo_cpn->notopt[iNumOpt] != 0) {
            floor +=
                exo_cpn->notopt[iNumOpt] *
                OPT_VAL_MACRO_N(
                    (exo_cpn->typeopt[iNumOpt] == SRT_CALL
                         ? (eval_const->cf_spr_opt_std[l][iNumOpt] == 0 ? 1 : 3)
                         : (eval_const->cf_spr_opt_std[l][iNumOpt] == 0 ? 2
                                                                        : 4)),
                    spread, exo_cpn->strikeopt[iNumOpt],
                    eval_const->cf_spr_opt_std[l][iNumOpt]);
          }
        }
      } else {

        /* Use the Shifted log model */
        for (iNumOpt = 0; iNumOpt < exo_cpn->nopt; iNumOpt++) {
          if (exo_cpn->notopt[iNumOpt] != 0) {
            err = OptSpreadNew(
                cmsrate[0] + exo_cpn->cms[0].slshift, /* fwd x */
                exo_cpn->alphabeta[0],                /*  nx  ,  */
                exo_cpn->cms[0].slvol,                /* sigx  , */
                cmsrate[1] + exo_cpn->cms[1].slshift, /* fwdy  , */
                exo_cpn->alphabeta[1],                /* ny  ,   */
                exo_cpn->cms[1].slvol,                /* sigy  , */
                exo_cpn->strikeopt[iNumOpt] +
                    exo_cpn->cms[0].slshift * exo_cpn->alphabeta[0] +
                    exo_cpn->cms[1].slshift *
                        exo_cpn->alphabeta[1], /* K  ,	 */
                mat,                           /* mat  ,  */
                eval_const->cf_spr_slcorr[l],  /* rho  ,  */
                exo_cpn->typeopt[iNumOpt], &optionvalue);

            floor += exo_cpn->notopt[iNumOpt] * optionvalue;
          }
        }
      }
    }

    /*	Coupon pv */
    exo_leg_pv += df * (coupon + floor - cap) * exo_cpn->cvg;
  }

  /*	Intrinsic value */
  if (call->pay_rec == 0) {
    iv = exo_leg_pv - fund_leg_pv;
  } else {
    iv = fund_leg_pv - exo_leg_pv;
  }

  /*	Fee */
  if (fabs(call->fee) > 1.0e-08) {
    fee = eval_const->df[eval_const->fee_idx] * call->fee;
  } else {
    fee = 0.0;
  }

  /*	Fill Result */
  res[0] = iv - fee;

  /* Calculate the Co Terminal Swap */
  /*	DF initial */
  df_start = eval_const->df[eval_const->start_idx];

  /*	Level minus the last one */
  level = 0.0;
  for (l = 0; l < num_fund_cpn; l++) {
    fund_cpn = ccf->fund_leg->cpn + call->fund_idx + l;
    level += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cvg;
  }

  /*	DF final */
  df_end = eval_const->df[eval_const->fund_idx[num_fund_cpn - 1]];
  forward = (df_start - df_end) / level;

  if (call->pay_rec == 0) {
    res[1] = -forward;
  } else {
    res[1] = forward;
  }

  /*	End of payoff valuation
          ----------------------- */

  return err;
}

/*	Main pricing function */

/*	Launch the tree */
Err ccf_launch_adi(CCF_STR ccf, CCF_UND und, CCF_ADI_ARG adi_arg,
                   /*	Result */
                   double *prem) {
  double temp_val[2];
  Err err = NULL;
  int one_lam = 1;
  double fixed_lam;
  int i;

  /* launch the corresponding ADI */

  fixed_lam = (adi_arg->lam)[0];

  for (i = 1; i < adi_arg->nb_lam; i++) {
    if (fabs((adi_arg->lam)[i] - fixed_lam) > 1.0E-10) {
      one_lam = 0;
      i = 1000000;
    }
  }

  if (one_lam && fixed_lam > 0.0) {
    if (ccf->cf_leg->cpn->use_SL || ccf->cf_leg->cpn->use_cfoptions) {
      err = lgm2f_adi(adi_arg->nstp, adi_arg->time, adi_arg->date,
                      adi_arg->nstpx, fixed_lam, adi_arg->sig_time,
                      adi_arg->sig1, adi_arg->nb_sig, adi_arg->alpha,
                      adi_arg->gamma, adi_arg->rho, adi_arg->void_prm,
                      adi_arg->is_event, adi_arg->ifr, adi_arg->yc,
                      ccf_payoff_4_3dfx_adi_sl, 1, (double *)temp_val);
    } else {
      err = lgm2f_adi(adi_arg->nstp, adi_arg->time, adi_arg->date,
                      adi_arg->nstpx, fixed_lam, adi_arg->sig_time,
                      adi_arg->sig1, adi_arg->nb_sig, adi_arg->alpha,
                      adi_arg->gamma, adi_arg->rho, adi_arg->void_prm,
                      adi_arg->is_event, adi_arg->ifr, adi_arg->yc,
                      ccf_payoff_4_3dfx_adi, 1, (double *)temp_val);
    }
  } else {
    if (ccf->cf_leg->cpn->use_SL || ccf->cf_leg->cpn->use_cfoptions) {
      err = lgm2fTau2_adi(
          adi_arg->nstp, adi_arg->time, adi_arg->date, adi_arg->nstpx, 0,
          adi_arg->sig1, adi_arg->sig_time, adi_arg->nb_sig, adi_arg->lam,
          adi_arg->lam_time, adi_arg->nb_lam, adi_arg->alpha, adi_arg->gamma,
          adi_arg->rho, adi_arg->void_prm, adi_arg->is_event, adi_arg->ifr,
          adi_arg->yc, ccf_payoff_4_3dfx_adi_sl, 1, (double *)temp_val);
    } else {
      err = lgm2fTau2_adi(
          adi_arg->nstp, adi_arg->time, adi_arg->date, adi_arg->nstpx, 0,
          adi_arg->sig1, adi_arg->sig_time, adi_arg->nb_sig, adi_arg->lam,
          adi_arg->lam_time, adi_arg->nb_lam, adi_arg->alpha, adi_arg->gamma,
          adi_arg->rho, adi_arg->void_prm, adi_arg->is_event, adi_arg->ifr,
          adi_arg->yc, ccf_payoff_4_3dfx_adi, 1, (double *)temp_val);
    }
  }

  *prem = temp_val[0];

  return err;
}

/*	Launch the MC */
Err ccf_launch_mc(CCF_STR ccf, CCF_UND und, CCF_MC_ARG mc_arg,
                  /*	Result */
                  double *prem, double *probas, double *exe_bound) {
  Err err = NULL;
  int i, k;
  MCEBParams params;
  int *optimise = NULL;
  double **prod_val = NULL;
  int flag;
  int nb_row;

  nb_row = max(mc_arg->nb_dates, 3);

  optimise = calloc(mc_arg->nb_dates, sizeof(int));
  prod_val = dmatrix(0, nb_row - 1, 0, 3);

  if (!optimise || !prod_val) {
    err = "Memory allocation faillure in ccf_launch_mc";
    goto FREE_RETURN;
  }

  for (i = 0; i < mc_arg->nb_dates; i++) {
    if (mc_arg->void_prm) {
      optimise[i] = 1;
    } else {
      optimise[i] = 0;
    }
  }

  mceb_set_default_params(&params);

  params.iCallCurrent = 1;
  params.iColBound = 1;
  params.iColPay = 0;
  params.iDoInfos = 1;
  params.iIsKO = 0;
  params.iKnockInCol = 0;
  params.iMultiIndex = 0;
  params.iNbIndex = 1;

  err = mceb_allocate_params(&params, mc_arg->nb_dates);
  if (err)
    goto FREE_RETURN;

  /* launch the corresponding MC */
  if (ccf->cf_leg->cpn->use_SL || ccf->cf_leg->cpn->use_cfoptions) {
    err = mc_main_lgm2f(
        /*	Time data */
        mc_arg->npaths, 2, mc_arg->time, mc_arg->date, mc_arg->nb_dates,
        mc_arg->jumping_num, mc_arg->dom_fwd1, mc_arg->dom_fwd2,
        mc_arg->dom_exp1, mc_arg->dom_exp2, mc_arg->dom_phi1, mc_arg->dom_phi2,
        mc_arg->dom_phi12, mc_arg->dom_gam1_fwd, mc_arg->dom_gam2_fwd,
        mc_arg->dom_bond_pay, mc_arg->dom_gam1_pay, mc_arg->dom_gam2_pay,
        mc_arg->covar, mc_arg->void_prm, 0, 1, optimise, &params, NULL,
        /*	Payoff function */
        ccf_payoff_4_3dfx_mc_sl, /*	Result */
        prod_val);

    if (err)
      goto FREE_RETURN;
  } else {
    err = mc_main_lgm2f(
        /*	Time data */
        mc_arg->npaths, 2, mc_arg->time, mc_arg->date, mc_arg->nb_dates,
        mc_arg->jumping_num, mc_arg->dom_fwd1, mc_arg->dom_fwd2,
        mc_arg->dom_exp1, mc_arg->dom_exp2, mc_arg->dom_phi1, mc_arg->dom_phi2,
        mc_arg->dom_phi12, mc_arg->dom_gam1_fwd, mc_arg->dom_gam2_fwd,
        mc_arg->dom_bond_pay, mc_arg->dom_gam1_pay, mc_arg->dom_gam2_pay,
        mc_arg->covar, mc_arg->void_prm, 0, 1, optimise, &params, NULL,
        /*	Payoff function */
        ccf_payoff_4_3dfx_mc, /*	Result */
        prod_val);

    if (err)
      goto FREE_RETURN;
  }

  /* Recopy Barrier / CoefLin for the moment */
  for (i = 0; i < mc_arg->nb_dates; i++) {
    prod_val[i][2] = params.dBarrier[i];

    for (k = 0; k < params.iNbIndex; k++) {
      prod_val[i][3 + k] = params.dCoefLin[i][k + 1];
    }
  }

  *prem = prod_val[2][0] * swp_f_df(und->today, mc_arg->pay_date, und->yc);

  flag = mc_arg->nb_dates - ccf->num_calls;

  for (i = 0; i < ccf->num_calls; i++) {
    probas[i] = prod_val[i + flag][3];
    exe_bound[i] = prod_val[i + flag][2];
  }

FREE_RETURN:

  if (optimise)
    free(optimise);
  if (prod_val)
    free_dmatrix(prod_val, 0, nb_row - 1, 0, 3);
  mceb_free_params(&params);

  return err;
}
