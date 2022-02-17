//--------------------------------------------------------------------
//----------------From D.M.'s srt_f_rangeaccrual.cxx--------------------
//--------------------------------------------------------------------
#include "RangeAccrualProdStruct.h"
#include "OPFNCTNS.H"
#include "math.h"
#include "srt_h_all.h"
#include "swp_h_vol.h"

#define K_TINY 1.0e-06
#define MAXCPN 1000

Err ra_free(RangeAccrualStruct *RA) {
  RangeAccrualObservStep *stp, *tmp;
  //	RAVolBumpData *bump_data        , *bump_tmp;

  if (RA) {
    if (RA->dates)
      free(RA->dates);
    if (RA->float_cpn_dates)
      free(RA->float_cpn_dates);
    if (RA->dates_f)
      free(RA->dates_f);
    if (RA->cvg_f)
      free(RA->cvg_f);
    if (RA->cvg_float_cpn)
      free(RA->cvg_float_cpn);
    if (RA->idx)
      free(RA->idx);
    if (RA->obs_steps) {
      stp = tmp = RA->obs_steps[0];
      while (stp) {
        stp = stp->next;
        free(tmp);
        tmp = stp;
      }
      free(RA->obs_steps);
    }
    /*		bump_data = bump_tmp = RA->bumps;
                    while (bump_data)
                    {
                            bump_data = bump_data->next;
                            free(bump_tmp);
                            bump_tmp = bump_data;
                    }*/
    //		free(RA);
  }

  return NULL;
}

Err ra_init_struct(
    long today, char *yc_d, char *yc_f, char *volc_d, char *volc_f,
    char *refrate_d, char *refrate_f, double notional, int ra_cpn_type,
    int n_fxvol_dates, long *fxvol_dates, double *fxvol, double *qtocorrel,
    int typeVol, int n_periods, long *dates, /* n_periods + 1 */
    double *cpns, char *basis, int *ra_nfixings, long **ra_fixingdates,
    double **ra_fixings, /*	Past coupon fixing if relevant */

    // RA floating coupon
    int ra_float_refrate_is_dom_for, char *ra_float_refrate,
    long *ra_float_startl, double *ra_float_past_fixings,
    double *ra_float_gearings,

    double *upper_barr, double *lower_barr, char *recpay, char *ra_refrate,
    double c_spread, int obs_freq, double rho_df,

    // Params for floating coupon
    double correl_start, double correl_end, int float_adj_strike,

    int eod_flag, /*	0: I        , 1: E */
    RangeAccrualStruct *RA) {
  Err err = NULL;
  RangeAccrualObservStep *stp;
  int i, j, index;
  long d, d1, d2, d_cpn, d1_cpn, d2_cpn;
  SrtBasisCode basis_d, basis_f, ibasis, basis_float_cpn;
  SrtCompounding freq_d, freq_f, freq_float_cpn;
  double vol_fx, rho_ffx, coef, fra_f, vol_f, fra_d, cvg_d, vol_d, power,
      vol_fra_cpn, vol_fx_cpn, rho_ffx_cpn, coef_cpn, fra_d_cpn, cvg_d_cpn,
      vol_d_cpn, cpn_float_spread, vol_f_lower_barr, vol_f_upper_barr;
  double quanto_corr, DRS_corr, quanto_corr_cpn, DRS_corr_cpn;
  double fra_float_cpn, cpn_float_upper_barrier_vol,
      cpn_float_lower_barrier_vol, n_std_dev, std_strike, correl;
  SrtReceiverType irecpay;
  SrtCurvePtr ycd_ptr;
  long ltoday;

  char *yc_cpn_float, *volc_cpn_float;

  RA->typeVol = typeVol;

  RA->dates_f = NULL;
  RA->float_cpn_dates = NULL;
  RA->cvg_f = NULL;
  RA->cvg_float_cpn = NULL;
  RA->idx = NULL;
  RA->dates = NULL;
  RA->obs_steps = NULL;

  RA->float_cpn_dates = NULL;

  ycd_ptr = lookup_curve(yc_d);
  if (!ycd_ptr) {
    err = "cannot find domestic yield curve";
    goto RETURN_ERR;
  }

  ltoday = get_today_from_curve(ycd_ptr);

  err = srt_f_get_spot_lag_from_refrate(refrate_f, &RA->spotlag_f);
  if (err)
    goto RETURN_ERR;

  // for ( index=0; index < n_periods &&
  //		   eod_flag + ltoday > add_unit( dates[index]        ,
  //-RA->spotlag_f        , SRT_BDAY        , MODIFIED_SUCCEEDING );
  // index++ );
  // Changed on 12 Jan 2004 following Laurent's request
  for (index = 0;
       index < n_periods && eod_flag + ltoday > ra_fixingdates[index][0];
       index++)
    ;

  /* Allocate memory and copy data to RA */
  RA->notional = notional;

  RA->n_periods = n_periods - index;
  RA->dates = (long *)calloc(n_periods - index + 1, sizeof(long));
  if (!RA->dates) {
    err = "memory allocation error in ra_init";
    goto RETURN_ERR;
  }

  RA->obs_steps = (RangeAccrualObservStep **)calloc(
      n_periods - index + 1, sizeof(RangeAccrualObservStep *));

  if (!RA->obs_steps) {
    err = "memory allocation error in ra_init";
    goto RETURN_ERR;
  }

  RA->obs_steps[0] = NULL;

  memcpy(RA->dates, dates + index, (n_periods - index + 1) * sizeof(long));
  RA->quanto = (fxvol != NULL);

  err = swp_f_get_ref_rate_details(ra_refrate, &basis_f, &freq_f);
  if (err)
    goto RETURN_ERR;
  RA->tenor_f = 12 / freq_f;
  RA->basis_f = basis_f;

  sprintf(RA->tenor_f_char, "%dM", RA->tenor_f);

  err = swp_f_get_ref_rate_details(refrate_d, &basis_d, &freq_d);
  if (err)
    goto RETURN_ERR;

  err = interp_basis(basis, &ibasis);
  if (err)
    goto RETURN_ERR;

  err = interp_rec_pay(recpay, &irecpay);
  if (err)
    goto RETURN_ERR;

  /* For RA floating coupons */
  RA->cpn_type = ra_cpn_type;

  if (RA->cpn_type != 0) {
    err = swp_f_get_ref_rate_details(ra_float_refrate, &basis_float_cpn,
                                     &freq_float_cpn);
    if (err)
      goto RETURN_ERR;
    RA->float_cpn_is_dom_for = ra_float_refrate_is_dom_for;
    RA->tenor_float_cpn = 12 / freq_float_cpn;
    RA->basis_float_cpn = basis_float_cpn;
    sprintf(RA->tenor_float_cpn_char, "%dM", RA->tenor_float_cpn);
    RA->float_cpn_dates = (long *)calloc(2 * (n_periods - index), sizeof(long));
    if (!RA->float_cpn_dates) {
      err = "memory allocation error in ra_init";
      goto RETURN_ERR;
    }
    RA->nb_float_cpn_dates = n_periods - index;

    // Stores the start and end dates of the floating coupons
    // j= index
    j = 0;
    // for (i = index; i < n_periods - index; i++)
    for (i = index; i < n_periods; i++) {
      RA->float_cpn_dates[j] = ra_float_startl[i];
      RA->float_cpn_dates[j + 1] =
          add_unit(RA->float_cpn_dates[j], RA->tenor_float_cpn, SRT_MONTH,
                   MODIFIED_SUCCEEDING);
      j += 2;
    }
  }

  /* Initialize observation schedule */

  stp = RA->obs_steps[0] =
      (RangeAccrualObservStep *)malloc(sizeof(RangeAccrualObservStep));
  if (!stp) {
    err = "memory allocation error in ra_init";
    goto RETURN_ERR;
  }

  stp->next = NULL;
  i = 0;

  /*	stp->date = d = RA->dates[0];
          while (i < n_periods - index)
          {
                  stp->next = (RangeAccrualObservStep *)
     malloc(sizeof(RangeAccrualObservStep)); stp = stp->next; if (!stp) goto
     RETURN_ERR; stp->next = NULL;

                  d = add_unit(d        , obs_freq        , SRT_DAY        ,
     SUCCEEDING); if (d >= RA->dates[i+1])
                  {
                          i++;
                          RA->obs_steps[i] = stp;
                          d = RA->dates[i];
                  }
                  stp->date = d;
          }
  */

  stp->date = d = add_unit(ra_fixingdates[index][0], RA->spotlag_f, SRT_BDAY,
                           MODIFIED_SUCCEEDING);
  j = obs_freq;
  while (i < n_periods - index) {
    stp->next =
        (RangeAccrualObservStep *)malloc(sizeof(RangeAccrualObservStep));
    stp = stp->next;
    if (!stp) {
      err = "memory allocation error in ra_init";
      goto RETURN_ERR;
    }

    stp->next = NULL;

    if (j < ra_nfixings[i + index]) {
      d = add_unit(ra_fixingdates[i + index][j], RA->spotlag_f, SRT_BDAY,
                   MODIFIED_SUCCEEDING);
    } else {
      i++;
      j = 0;
      RA->obs_steps[i] = stp;

      if (i < n_periods - index) {
        d = add_unit(ra_fixingdates[i + index][j], RA->spotlag_f, SRT_BDAY,
                     MODIFIED_SUCCEEDING);
      } else {
        d = RA->dates[i];
      }

      j = 0;
    }

    j = j + obs_freq;
    stp->date = d;
  }

  /* Add today with zero bump to the vol bump list */
  d = get_today_from_curve(ycd_ptr);
  /*	RA->bumps = (RAVolBumpData *) malloc(sizeof(RAVolBumpData));
          if (!RA->bumps) goto RETURN_ERR;
          RA->bumps->date = d;
          RA->bumps->bump = 0.0;
          RA->bumps->next = NULL;
          RA->bump_data = RA->bumps;
  */
  RA->last_date = -1;

  /* Fill in observation steps */

  j = 0;
  for (i = index, stp = RA->obs_steps[0]; stp->next; stp = stp->next) {

    // stp->date = start date of the index and finds the corresponding i
    if (stp == RA->obs_steps[i + 1 - index] || stp == RA->obs_steps[0]) {
      // Increment the coupon index i if we are not at the beginnning
      if (stp != RA->obs_steps[0])
        i++;

      // Calculates the RA floating coupon stuff
      if (RA->cpn_type != 0) {

        d_cpn = add_unit(ra_float_startl[i], RA->tenor_float_cpn, SRT_MONTH,
                         MODIFIED_SUCCEEDING);

        // Get the FRA of the floating coupon
        if (RA->float_cpn_is_dom_for == 0) {
          yc_cpn_float = yc_d;
          volc_cpn_float = volc_d;
        } else {
          yc_cpn_float = yc_f;
          volc_cpn_float = volc_f;
        }

        fra_float_cpn =
            swp_f_fra(ra_float_startl[i], d_cpn, RA->basis_float_cpn,
                      yc_cpn_float, ra_float_refrate);
        cpn_float_spread =
            swp_f_spread(ra_float_startl[i], d_cpn, ra_float_refrate);
        err = swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn, 0.05,
                            &vol_fra_cpn, &power, SABR_ATMLOG);
        if (err)
          goto RETURN_ERR;

        // Quanto adjustment of the floating coupon
        if (RA->float_cpn_is_dom_for == 1) {
          // Get the FX ATM vol & FX correlations (using the same ones as for
          // the quanto index
          for (;
               j < n_fxvol_dates - 1 && fxvol_dates[j + 1] < ra_float_startl[i];
               j++)
            ;
          if (fxvol_dates[j] > d_cpn || j == n_fxvol_dates - 1) {
            vol_fx_cpn = fxvol[j];
            rho_ffx_cpn = qtocorrel[j];
          } else {
            coef_cpn = ((double)(d_cpn - fxvol_dates[j])) /
                       (fxvol_dates[j + 1] - fxvol_dates[j]);
            vol_fx_cpn = coef_cpn * fxvol[j + 1] + (1.0 - coef_cpn) * fxvol[j];
            rho_ffx_cpn =
                coef_cpn * qtocorrel[j + 1] + (1.0 - coef_cpn) * qtocorrel[j];
          }

          quanto_corr_cpn = -rho_ffx_cpn * vol_fx_cpn * vol_fra_cpn;
        } else {
          quanto_corr_cpn = 0.0;
        }

        // DRS adjustment of the floating coupon
        if (ra_float_startl[i] == RA->dates[i - index] &&
            d_cpn == RA->dates[i + 1 - index])
          DRS_corr_cpn = 0;
        else {
          d1_cpn = (RA->dates[i + 1 - index] < d_cpn ? RA->dates[i + 1 - index]
                                                     : d_cpn);
          d2_cpn = (RA->dates[i + 1 - index] > d_cpn ? RA->dates[i + 1 - index]
                                                     : d_cpn);

          if (fabs(d2_cpn - d1_cpn) < 10) {
            d2_cpn = d1_cpn + 10;
          }
          fra_d_cpn = swp_f_fra(d1_cpn, d2_cpn, basis_d, yc_d, refrate_d);
          cvg_d_cpn = coverage(d1_cpn, d2_cpn, basis_d);
          err = swp_f_SABRvol(volc_d, d1_cpn, d2_cpn, fra_d_cpn, &vol_d_cpn,
                              &power, SABR_LOGVOL);
          if (err)
            return err;

          DRS_corr_cpn = rho_df * vol_d_cpn * vol_fra_cpn * cvg_d_cpn *
                         fra_d_cpn / (1.0 + cvg_d_cpn * fra_d_cpn);
          if (RA->dates[i + 1 - index] > d_cpn)
            DRS_corr_cpn = -DRS_corr_cpn;
        }

        // Initially was getting the smile of the floating coupon
        // in order to get the correct LOGVOL when the LGM diffusion is done
        // and calculate the correct numeraire change. However        , since
        // Bertrand was already freezing the logvol today for the index
        // adjustment        , we will do the same.

        // We now get the logvol of the floating coupon for the numeraire
        // adjustment

        // Floating coupon volatility for the numeraire change if ATM or ATS
        if (float_adj_strike == 0) {
          // ATM
          err =
              swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn, 0.05,
                            &cpn_float_lower_barrier_vol, &power, SABR_ATMLOG);
          if (err)
            goto RETURN_ERR;
          cpn_float_upper_barrier_vol = cpn_float_lower_barrier_vol;

        } else if (float_adj_strike == 1) {
          // ATS
          err = swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn,
                              lower_barr[i], &cpn_float_lower_barrier_vol,
                              &power, SABR_LOGVOL);
          if (err)
            goto RETURN_ERR;
          err = swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn,
                              upper_barr[i], &cpn_float_upper_barrier_vol,
                              &power, SABR_LOGVOL);
          if (err)
            goto RETURN_ERR;
        }
      }
    }

    // Writes all params for floating coupon
    if (RA->cpn_type != 0) {
      // Params for floating coupon
      stp->cpn_float_spread = cpn_float_spread;
      stp->cpn_float_correction = DRS_corr_cpn + quanto_corr_cpn;
      stp->cpn_float_gearing = ra_float_gearings[i];
      stp->cvg = ((double)(stp->next->date - stp->date)) /
                 ((double)(RA->obs_steps[i + 1 - index]->date -
                           RA->obs_steps[i - index]->date)) *
                 coverage(RA->dates[i], RA->dates[i + 1], ibasis);
      stp->cpn_float_past_fixing = ra_float_past_fixings[i];

      // Params for numeraire adjustment (needs to add at the barrier strikes)
      correl = (correl_end - correl_start) /
                   (RA->dates[i + 1 - index] - RA->dates[i - index]) *
                   (stp->date - RA->dates[i - index]) +
               correl_start;
      stp->cpn_float_lower_barrier_adj = correl * cpn_float_lower_barrier_vol;
      stp->cpn_float_upper_barrier_adj = correl * cpn_float_upper_barrier_vol;
    }

    // d = end date of the index
    // d = add_unit(stp->date        , RA->tenor_f        , SRT_MONTH        ,
    // NO_BUSDAY_CONVENTION);
    d = add_unit(stp->date, RA->tenor_f, SRT_MONTH, MODIFIED_SUCCEEDING);
    stp->spread = swp_f_spread(stp->date, d, ra_refrate);

    // RA coupon (or margin)
    stp->cpn =
        cpns[i] * (stp->next->date - stp->date) /
        (RA->obs_steps[i + 1 - index]->date - RA->obs_steps[i - index]->date) *
        coverage(RA->dates[i - index], RA->dates[i + 1 - index], ibasis);

    // Get the ATM vol of the index
    fra_f = swp_f_fra(stp->date, d, basis_f, yc_f, ra_refrate);
    err =
        swp_f_SABRvol(volc_f, stp->date, d, fra_f, &vol_f, &power, SABR_LOGVOL);
    if (err)
      goto RETURN_ERR;

    // If floating coupon AND numeraire adjustment method is std dev logvol of
    // the index stdev depends on the index so has to be done for every index
    if (RA->cpn_type != 0 && float_adj_strike == 2) {
      if (vol_f > 0.0) {
        // Lower barrier stdev logvol
        n_std_dev = log(lower_barr[i] / fra_f) /
                    (vol_f * sqrt((stp->date - today) / 365.0));
        std_strike =
            fra_float_cpn * exp(correl * vol_fra_cpn * n_std_dev *
                                sqrt((ra_float_startl[i] - today) / 365.0));
        err =
            swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn, std_strike,
                          &cpn_float_lower_barrier_vol, &power, SABR_LOGVOL);
        if (err)
          goto RETURN_ERR;

        // Upper barrier stdev logvol
        n_std_dev = log(upper_barr[i] / fra_f) /
                    (vol_f * sqrt((stp->date - today) / 365.0));
        std_strike =
            fra_float_cpn * exp(correl * vol_fra_cpn * n_std_dev *
                                sqrt((ra_float_startl[i] - today) / 365.0));
        err =
            swp_f_SABRvol(volc_cpn_float, ra_float_startl[i], d_cpn, std_strike,
                          &cpn_float_upper_barrier_vol, &power, SABR_LOGVOL);
        if (err)
          goto RETURN_ERR;

        // Storage but still needs to add at the barrier vols
        stp->cpn_float_lower_barrier_adj = correl * cpn_float_lower_barrier_vol;
        stp->cpn_float_upper_barrier_adj = correl * cpn_float_upper_barrier_vol;
      } else {
        stp->cpn_float_lower_barrier_adj = 0.0;
        stp->cpn_float_upper_barrier_adj = 0.0;
      }
    }

    // Finally multiply by the ATBarrier log vols of the index
    if (RA->cpn_type != 0) {
      err = swp_f_SABRvol(volc_f, stp->date, d, lower_barr[i],
                          &vol_f_lower_barr, &power, SABR_LOGVOL);
      if (err)
        goto RETURN_ERR;
      stp->cpn_float_lower_barrier_adj *= vol_f_lower_barr;

      err = swp_f_SABRvol(volc_f, stp->date, d, upper_barr[i],
                          &vol_f_upper_barr, &power, SABR_LOGVOL);
      if (err)
        goto RETURN_ERR;
      stp->cpn_float_upper_barrier_adj *= vol_f_upper_barr;
    }

    if (RA->quanto) /* Calculate quanto adjustment */
    {
      /* Interpolate FX vol and rho_ffx at d */
      for (; j < n_fxvol_dates - 1 && fxvol_dates[j + 1] < d; j++)
        ;
      if (fxvol_dates[j] > d || j == n_fxvol_dates - 1) {
        vol_fx = fxvol[j];
        rho_ffx = qtocorrel[j];
      } else {
        coef = ((double)(d - fxvol_dates[j])) /
               (fxvol_dates[j + 1] - fxvol_dates[j]);
        vol_fx = coef * fxvol[j + 1] + (1.0 - coef) * fxvol[j];
        rho_ffx = coef * qtocorrel[j + 1] + (1.0 - coef) * qtocorrel[j];
      }
      quanto_corr = -rho_ffx * vol_fx * vol_f;
    } else
      quanto_corr = 0.0;

    /* Calculate DRS adjustment */
    if (RA->dates[i + 1 - index] == d)
      DRS_corr = 0;
    else {
      d1 = (RA->dates[i + 1 - index] < d ? RA->dates[i + 1 - index] : d);
      d2 = (RA->dates[i + 1 - index] > d ? RA->dates[i + 1 - index] : d);

      if (fabs(d2 - d1) < 10) {
        d2 = d1 + 10;
      }
      fra_d = swp_f_fra(d1, d2, basis_d, yc_d, refrate_d);
      cvg_d = coverage(d1, d2, basis_d);
      err = swp_f_SABRvol(volc_d, d1, d2, fra_d, &vol_d, &power, SABR_LOGVOL);
      if (err) {
        return err;
      }

      DRS_corr = rho_df * vol_d * vol_f * cvg_d * fra_d / (1.0 + cvg_d * fra_d);
      if (RA->dates[i + 1 - index] > d)
        DRS_corr = -DRS_corr;
    }

    /* Correction for index = quanto correction and DRS correction */
    stp->correction = quanto_corr + DRS_corr;

    /* Calculate vols at barriers */
    if (irecpay == SRT_RECEIVER) {
      if (lower_barr[i] == 0) {
        lower_barr[i] = -999;
      }
      stp->K[0] = lower_barr[i];
      stp->K[1] = lower_barr[i] + c_spread;
      stp->K[2] = upper_barr[i] - c_spread;
      stp->K[3] = upper_barr[i];

    } else {
      stp->K[0] = lower_barr[i] - c_spread;
      stp->K[1] = lower_barr[i];
      stp->K[2] = upper_barr[i];
      stp->K[3] = upper_barr[i] + c_spread;
    }

    /* Get SABR parameters for the step */
    err = swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->betavol, &power,
                        SABR_BETAVOL);
    if (err) {
      goto RETURN_ERR;
    }
    err = swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->normvol, &power,
                        SABR_ATMNORM);
    if (err) {
      goto RETURN_ERR;
    }
    err = swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->lognvol, &power,
                        SABR_ATMLOG);
    if (err) {
      goto RETURN_ERR;
    }
    err = swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->alpha, &power,
                        SABR_ALPHA);
    if (err) {
      goto RETURN_ERR;
    }
    err = swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->beta, &power,
                        SABR_BETA);
    if (err) {
      goto RETURN_ERR;
    }
    err =
        swp_f_SABRvol(volc_f, stp->date, d, 0.05, &stp->rho, &power, SABR_RHO);
    if (err) {
      goto RETURN_ERR;
    }
  }

RETURN_ERR:
  if (err) {
    ra_free(RA);
  }
  return err;
}

Err RA_RequestDfDates(RangeAccrualStruct *RA, long date, int *pn_unds,
                      int **pn_dfs, long ***pdates) {
  int i, j, ndfs_d, ndfs_f, ndfs_float_cpn, cpn_col;
  long d, d1, dates_d[MAXCPN], *dates_f;
  long idx[MAXCPN];
  Err err = NULL;

  // Need to add 1 more column for floating coupon
  // *pn_unds = (RA->quanto ? 2 : 1);

  if (RA->cpn_type != 0)
    // For floating coupons
    *pn_unds = (RA->quanto ? 3 : 2);
  else
    // For fix coupons
    *pn_unds = (RA->quanto ? 2 : 1);

  *pn_dfs = (int *)calloc(*pn_unds, sizeof(int));
  if (!(*pn_dfs)) {
    err = "Memory allocation failed in RA_RequestDfDates";
    goto FREE_RETURN;
  }
  *pdates = (long **)calloc(*pn_unds, sizeof(long *));
  if (!(*pdates)) {
    err = "Memory allocation failed in RA_RequestDfDates";
    goto FREE_RETURN;
  }

  // i stores the index of the coupon corresponding to the current event date
  for (i = 0;
       i < RA->n_periods && date > add_unit(RA->dates[i], -RA->spotlag_f,
                                            SRT_BDAY, MODIFIED_SUCCEEDING);
       i++)
    ;

  if (i == RA->n_periods) {
    for (j = 0; j < *pn_unds; j++) {
      (*pn_dfs)[j] = 0;
      (*pdates)[j] = NULL;
    }
    return NULL;
  }

  // For the funding        , we only need RA->n_periods - i
  ndfs_d = RA->n_periods - i;

  // Number of discount factors for floating coupons
  ndfs_float_cpn = ndfs_d;

  dates_f = dates_d + ndfs_d;

  dates_f[j = 0] = d = RA->dates[i];
  d1 = d + 15 * RA->tenor_f;
  while (1) {
    dates_f[++j] = bus_date_method(d1, MODIFIED_SUCCEEDING);
    d = add_unit(d, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
    d1 = add_unit(d1, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
    dates_f[++j] = bus_date_method(d, MODIFIED_SUCCEEDING);
    if (dates_f[j - 2] >= RA->dates[RA->n_periods])
      break;
  }
  ndfs_f = j + 1;

  memcpy(dates_d, RA->dates + i + 1, ndfs_d * sizeof(long));

  if (RA->quanto) {
    (*pn_dfs)[0] = ndfs_d;
    (*pn_dfs)[1] = ndfs_f;
    (*pdates)[0] = (long *)calloc(ndfs_d, sizeof(long));
    if (!((*pdates)[0])) {
      err = "Memory allocation failed in RA_RequestDfDates";
      goto FREE_RETURN;
    }
    (*pdates)[1] = (long *)calloc(ndfs_f, sizeof(long));
    if (!((*pdates)[1])) {
      err = "Memory allocation failed in RA_RequestDfDates";
      goto FREE_RETURN;
    }
    memcpy((*pdates)[0], dates_d, ndfs_d * sizeof(long));
    memcpy((*pdates)[1], dates_f, ndfs_f * sizeof(long));
  } else {
    (*pn_dfs)[0] = ndfs_d + ndfs_f;
    (*pdates)[0] = (long *)calloc(ndfs_d + ndfs_f, sizeof(long));
    if (!((*pdates)[0])) {
      err = "Memory allocation failed in RA_RequestDfDates";
      goto FREE_RETURN;
    }
    indexx_ll(dates_d, idx, ndfs_d + ndfs_f);
    for (j = 0; j < ndfs_d + ndfs_f; j++)
      (*pdates)[0][j] = dates_d[idx[j]];
  }

  // If floating coupons
  if (RA->cpn_type != 0) {
    // Calculates correct index column (1 if non quanto and 2 if quanto)
    if (RA->quanto)
      cpn_col = 2; // 0 and 1 are taken
    else
      cpn_col = 1; // 0 is taken

    (*pn_dfs)[cpn_col] = ndfs_float_cpn;
    (*pdates)[cpn_col] = (long *)calloc(2 * ndfs_float_cpn, sizeof(long));
    if (!((*pdates)[cpn_col])) {
      err = "Memory allocation failed in RA_RequestDfDates";
      goto FREE_RETURN;
    }

    // A checker si on veut + 2*i ou + 2*i + 2 !!!
    memcpy((*pdates)[cpn_col], RA->float_cpn_dates + 2 * i,
           2 * ndfs_float_cpn * sizeof(long));
  }

FREE_RETURN:

  return err;
}

Err RA_Payoff(RangeAccrualStruct *RA, long today, long date, double **dfs,
              double *payoff) {
  Err err = NULL;
  int i, j, k, nfwds;
  long d, d1;
  double dfs_d[MAXCPN], dfs_cpn_float[MAXCPN], *dfs_f;
  double fwds[MAXCPN], fwds_cpn_float[MAXCPN];
  long dates_d[MAXCPN], dates_cpn_float[MAXCPN], *dates_f;
  double coef, fwd, fixing, mat, df, sum, vol, bs[4], bump_dir[4],
      bs_float_cpn[4];
  double fwd_cpn_float, mat_cpn_float, fwd_lower_barr, fwd_upper_barr,
      fixing_cpn_float;
  double normvol, lognvol, betavol;
  RangeAccrualObservStep *stp;
  //	RAVolBumpData *bump_tmp;

  *payoff = 0.0;

  /* If stored deterministic values are no more valid recalculate them (when
   * event date changes)*/
  if (date != RA->last_date) {
    if (RA->idx) {
      free(RA->idx);
      RA->idx = NULL;
    }

    if (RA->dates_f) {
      free(RA->dates_f);
      RA->dates_f = NULL;
    }

    if (RA->cvg_f) {
      free(RA->cvg_f);
      RA->cvg_f = NULL;
    }

    if (RA->cvg_float_cpn) {
      free(RA->cvg_float_cpn);
      RA->cvg_float_cpn = NULL;
    }

    for (i = 0;
         i < RA->n_periods && date > add_unit(RA->dates[i], -RA->spotlag_f,
                                              SRT_BDAY, MODIFIED_SUCCEEDING);
         i++)
      ;
    RA->last_date = date;
    RA->period_idx = i;
    if (i == RA->n_periods)
      return NULL;

    RA->ndfs_d = RA->n_periods - i;

    // Number of discount factors for floating coupons
    if (RA->cpn_type != 0) {

      RA->ndfs_cpn_float = RA->nb_float_cpn_dates - i;

      memcpy(dates_cpn_float, RA->float_cpn_dates + 2 * i,
             2 * RA->ndfs_cpn_float * sizeof(long));

      RA->cvg_float_cpn = (double *)calloc(RA->ndfs_cpn_float, sizeof(double));
      if (!RA->cvg_float_cpn) {
        smessage("Memory allocation failed in RA_Payoff");
        err = "Memory allocation failed in RA_Payoff";
        goto FREE_RETURN;
      }

      k = 0;
      for (j = 0; j < 2 * RA->ndfs_cpn_float; j += 2) {
        RA->cvg_float_cpn[k] = coverage(
            dates_cpn_float[j], dates_cpn_float[j + 1], RA->basis_float_cpn);
        k++;
      }
    }

    dates_f = dates_d + RA->ndfs_d;

    dates_f[j = 0] = d = RA->dates[i];
    d1 = d + 15 * RA->tenor_f;
    while (1) {
      dates_f[++j] = bus_date_method(d1, MODIFIED_SUCCEEDING);
      d = add_unit(d, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
      d1 = add_unit(d1, RA->tenor_f, SRT_MONTH, NO_BUSDAY_CONVENTION);
      dates_f[++j] = bus_date_method(d, MODIFIED_SUCCEEDING);
      if (dates_f[j - 2] >= RA->dates[RA->n_periods])
        break;
    }
    RA->ndfs_f = j + 1;

    if (!RA->quanto) {
      memcpy(dates_d, RA->dates + i + 1, RA->ndfs_d * sizeof(long));
      if (RA->idx) {
        free(RA->idx);
        RA->idx = NULL;
      }

      RA->idx = (long *)calloc((RA->ndfs_d + RA->ndfs_f), sizeof(long));
      if (!RA->idx) {
        smessage("Memory allocation failed in RA_Payoff");
        err = "Memory allocation failed in RA_Payoff";
        goto FREE_RETURN;
      }

      //			RA->idx = (long *) realloc( RA->idx        ,
      //				(RA->ndfs_d + RA->ndfs_f) * sizeof(long)
      //);

      indexx_ll(dates_d, RA->idx, RA->ndfs_d + RA->ndfs_f);
    }

    if (RA->dates_f) {
      free(RA->dates_f);
      RA->dates_f = NULL;
    }
    RA->dates_f = (long *)calloc((RA->ndfs_f - 2), sizeof(long));
    if (!RA->dates_f) {
      smessage("Memory allocation failed in RA_Payoff");
      err = "Memory allocation failed in RA_Payoff";
      goto FREE_RETURN;
    }

    if (RA->cvg_f) {
      free(RA->cvg_f);
      RA->cvg_f = NULL;
    }
    RA->cvg_f = (double *)calloc((RA->ndfs_f - 2), sizeof(double));
    if (!RA->cvg_f) {
      smessage("Memory allocation failed in RA_Payoff");
      err = "Memory allocation failed in RA_Payoff";
      goto FREE_RETURN;
    }

    memcpy(RA->dates_f, dates_f, (RA->ndfs_f - 2) * sizeof(long));
    for (j = 0; j < RA->ndfs_f - 2; j++)
      RA->cvg_f[j] = coverage(dates_f[j], dates_f[j + 2], RA->basis_f);
  }
  if (RA->period_idx == RA->n_periods)
    return NULL;

  dfs_f = dfs_d + RA->ndfs_d;

  if (RA->quanto) {
    memcpy(dfs_d, dfs[0], RA->ndfs_d * sizeof(double));
    memcpy(dfs_f, dfs[1], RA->ndfs_f * sizeof(double));
  } else
    for (j = 0; j < RA->ndfs_d + RA->ndfs_f; j++)
      dfs_d[RA->idx[j]] = dfs[0][j];

  /* Calculate forwards cash */
  nfwds = RA->ndfs_f - 2;
  for (j = 0; j < nfwds; j++)
    fwds[j] = (dfs_f[j] / dfs_f[j + 2] - 1.0) / RA->cvg_f[j];

  /* In case of floating coupons */
  if (RA->cpn_type != 0) {
    memcpy(dfs_cpn_float, dfs[2], 2 * RA->ndfs_cpn_float * sizeof(double));

    /* Calculate forwards cash for the floating coupon */
    k = 0;
    for (j = 0; j < 2 * RA->ndfs_cpn_float; j += 2) {
      fwds_cpn_float[k] = (dfs_cpn_float[j] / dfs_cpn_float[j + 1] - 1.0) /
                          RA->cvg_float_cpn[k];
      k++;
    }
  }

  /* Find vol bump corresponding to today */
  /*	if (RA->bump_data->date != today)
          {
                  RA->bump_data = RA->bumps;
                  while (RA->bump_data && RA->bump_data->date != today)
                  {
                          bump_tmp = RA->bump_data;
                          RA->bump_data = RA->bump_data->next;
                  }
                  if (!RA->bump_data)		// today not in the list
                  {
                          bump_tmp->next = (RAVolBumpData *)
     malloc(sizeof(RAVolBumpData)); RA->bump_data = bump_tmp->next;
                          RA->bump_data->date = today;
                          RA->bump_data->bump = 0.0;
                          RA->bump_data->next = NULL;
                  }
          }
  */

  // If floating coupon        , setup the floating coupon params for the first
  // period
  if (RA->cpn_type != 0) {
    i = RA->period_idx;
    stp = RA->obs_steps[i];
    fwd_cpn_float = fwds_cpn_float[i - RA->period_idx] + stp->cpn_float_spread;

    // Do not need to do RA->nb - i since RA->float_cpn_dates is the entire
    // array (not a partial array)
    fixing_cpn_float = add_unit(RA->float_cpn_dates[2 * i], -RA->spotlag_f,
                                SRT_BDAY, MODIFIED_SUCCEEDING);
    mat_cpn_float = (fixing_cpn_float - today) * YEARS_IN_DAY;

    // If floating coupon has already fixed...
    if (mat_cpn_float <= 0.0)
      fwd_cpn_float = stp->cpn_float_past_fixing;
    else
      fwd_cpn_float *= exp(stp->cpn_float_correction * mat_cpn_float);
  }

  /* Calculate the sum of call spreads interpolating forwards */
  for (k = 0, stp = RA->obs_steps[i = RA->period_idx]; stp->next;
       stp = stp->next) {
    if (stp == RA->obs_steps[i + 1]) {
      i++;

      if (RA->cpn_type != 0) {
        // Calculates the fwd corresponding to the floating coupon
        fwd_cpn_float =
            fwds_cpn_float[i - RA->period_idx] + stp->cpn_float_spread;

        // Do not need to do RA->nb - i since RA->float_cpn_dates is the entire
        // array (not a partial array)
        fixing_cpn_float = add_unit(RA->float_cpn_dates[2 * i], -RA->spotlag_f,
                                    SRT_BDAY, MODIFIED_SUCCEEDING);
        mat_cpn_float = (fixing_cpn_float - today) * YEARS_IN_DAY;

        // If floating coupon has already fixed...
        if (mat_cpn_float <= 0.0)
          fwd_cpn_float = stp->cpn_float_past_fixing;
        else
          fwd_cpn_float *= exp(stp->cpn_float_correction * mat_cpn_float);
      }
    }

    for (; k < nfwds - 1 && RA->dates_f[k + 1] < stp->date; k++)
      ;

    coef = ((double)(stp->date - RA->dates_f[k])) /
           (RA->dates_f[k + 1] - RA->dates_f[k]);
    fwd = coef * fwds[k + 1] + (1.0 - coef) * fwds[k];

    fixing = add_unit(stp->date, -RA->spotlag_f, SRT_BDAY, MODIFIED_SUCCEEDING);
    mat = (fixing - today) * YEARS_IN_DAY;
    df = dfs_d[i - RA->period_idx]; /* df(dates[i+1]) */

    /* Adjust forward */
    fwd += stp->spread;
    fwd *= exp(stp->correction * mat);

    /* Calculate bump_dir */
    if (fwd > stp->K[3])
      bump_dir[3] = bump_dir[2] = -1.0;
    else if (fwd < stp->K[2])
      bump_dir[3] = bump_dir[2] = 1.0;
    else
      bump_dir[3] = bump_dir[2] = 0.0;

    if (fwd < stp->K[0])
      bump_dir[0] = bump_dir[1] = -1.0;
    else if (fwd > stp->K[1])
      bump_dir[0] = bump_dir[1] = 1.0;
    else
      bump_dir[0] = bump_dir[1] = 0.0;

    /* Calculate puts */
    for (j = 0; j < 4; j++) {
      if (stp->K[j] < K_TINY) {
        bs[j] = 0.0;
        if (RA->cpn_type != 0)
          bs_float_cpn[j] = 0.0;
      } else if (fwd < K_TINY) {
        bs[j] = df * (stp->K[j] - fwd);
        if (RA->cpn_type != 0) {
          if (j == 0 || j == 1) {
            fwd_lower_barr = fwd * exp(stp->cpn_float_lower_barrier_adj *
                                       min(mat, mat_cpn_float));
            bs_float_cpn[j] = df * (stp->K[j] - fwd_lower_barr);
          } else {
            fwd_upper_barr = fwd * exp(stp->cpn_float_upper_barrier_adj *
                                       min(mat, mat_cpn_float));
            bs_float_cpn[j] = df * (stp->K[j] - fwd_upper_barr);
          }
        }
      } else {
        if (RA->typeVol == 1) // Normal
        {
          normvol = stp->normvol;
          if (normvol < K_TINY)
            normvol = K_TINY;
          err = srt_f_optsarbvol(fwd, stp->K[j], mat, normvol, stp->alpha,
                                 stp->beta, stp->rho, SRT_NORMAL, SRT_LOGNORMAL,
                                 &vol);
          if (err)
            return err;
        } else if (RA->typeVol == 2) // Lognormal
        {
          lognvol = stp->lognvol;
          if (lognvol < K_TINY)
            lognvol = K_TINY;
          err = srt_f_optsarbvol(fwd, stp->K[j], mat, lognvol, stp->alpha,
                                 stp->beta, stp->rho, SRT_LOGNORMAL,
                                 SRT_LOGNORMAL, &vol);
          if (err)
            return err;
        } else // Beta
        {
          betavol = stp->betavol;
          if (betavol < K_TINY)
            betavol = K_TINY;
          err = srt_f_optsarbvol(fwd, stp->K[j], mat, betavol, stp->alpha,
                                 stp->beta, stp->rho, SRT_BETAVOL,
                                 SRT_LOGNORMAL, &vol);
          if (err)
            return err;
        }

        // FIX coupon
        bs[j] = srt_f_optblksch(fwd, stp->K[j], vol, mat, df, SRT_PUT, PREMIUM);
        if (RA->cpn_type != 0) {
          // FLOAT coupon
          if (j == 0 || j == 1) {
            fwd_lower_barr = fwd * exp(stp->cpn_float_lower_barrier_adj *
                                       min(mat, mat_cpn_float));
            bs_float_cpn[j] = srt_f_optblksch(fwd_lower_barr, stp->K[j], vol,
                                              mat, df, SRT_PUT, PREMIUM);
          } else {
            fwd_upper_barr = fwd * exp(stp->cpn_float_upper_barrier_adj *
                                       min(mat, mat_cpn_float));
            bs_float_cpn[j] = srt_f_optblksch(fwd_upper_barr, stp->K[j], vol,
                                              mat, df, SRT_PUT, PREMIUM);
          }
        }
      }
    }

    sum = (bs[3] - bs[2]) / (stp->K[3] - stp->K[2]) -
          (bs[1] - bs[0]) / (stp->K[1] - stp->K[0]);

    *payoff += sum * stp->cpn;

    if (RA->cpn_type != 0) {
      sum = (bs_float_cpn[3] - bs_float_cpn[2]) / (stp->K[3] - stp->K[2]) -
            (bs_float_cpn[1] - bs_float_cpn[0]) / (stp->K[1] - stp->K[0]);

      *payoff += sum * stp->cpn_float_gearing * fwd_cpn_float * stp->cvg;
    }
  }

  *payoff = *payoff * RA->notional;

FREE_RETURN:

  if (err) {
    if (RA->idx) {
      free(RA->idx);
      RA->idx = NULL;
    }
    if (RA->dates_f) {
      free(RA->dates_f);
      RA->dates_f = NULL;
    }
    if (RA->cvg_f) {
      free(RA->cvg_f);
      RA->cvg_f = NULL;
    }
    if (RA->cvg_float_cpn) {
      free(RA->cvg_float_cpn);
      RA->cvg_float_cpn = NULL;
    }
  }

  return NULL;
}

Err RA_FwdPV(char *yc_d, char *yc_f, RangeAccrualStruct *RA, long today,
             int nbDates, long *Dates, double *fwdPV) {
  Err err = NULL;
  int i, j;
  int n_cpns;

  int pn_unds;
  int *pn_dfs = NULL;
  long **pdates = NULL;

  double **ra_dfs = NULL;

  long date;

  double payoff;

  n_cpns = RA->n_periods;

  for (i = 0; i < nbDates; ++i) {
    //		date = add_unit( RA->dates[i]        , -RA->spotlag_f - 1        ,
    //SRT_BDAY        , MODIFIED_SUCCEEDING );
    date = Dates[i];
    err = RA_RequestDfDates(RA, date, &pn_unds, &pn_dfs, &pdates);
    if (err)
      goto FREE_RETURN;

    ra_dfs = (double **)calloc(pn_unds, sizeof(double *));
    if (!ra_dfs) {
      err = "Memory allocation failed";
      goto FREE_RETURN;
    }
    ra_dfs[0] = (double *)calloc(pn_dfs[0], sizeof(double));
    if (!(ra_dfs[0])) {
      err = "Memory allocation failed";
      goto FREE_RETURN;
    }
    ra_dfs[1] = (double *)calloc(pn_dfs[1], sizeof(double));
    if (!(ra_dfs[1])) {
      err = "Memory allocation failed";
      goto FREE_RETURN;
    }

    // In case of floating coupons
    if (RA->cpn_type != 0) {
      ra_dfs[2] = (double *)calloc(2 * pn_dfs[2], sizeof(double));
      if (!ra_dfs[2]) {
        err = "Memory allocation error failed";
        goto FREE_RETURN;
      }
    }

    for (j = 0; j < pn_dfs[0]; ++j) {
      ra_dfs[0][j] = swp_f_df(today, pdates[0][j], yc_d);
    }
    for (j = 0; j < pn_dfs[1]; ++j) {
      ra_dfs[1][j] = swp_f_df(today, pdates[1][j], yc_f);
    }

    // DFs in order to calculate the floating coupons of the RA leg
    if (RA->cpn_type != 0) {
      if (RA->float_cpn_is_dom_for == 0) {
        // j++ ou ++j ????
        for (j = 0; j < 2 * pn_dfs[2]; ++j) {
          ra_dfs[2][j] = swp_f_df(today, pdates[2][j], yc_d);
        }
      } else {
        for (j = 0; j < 2 * pn_dfs[2]; ++j) {
          ra_dfs[2][j] = swp_f_df(today, pdates[2][j], yc_f);
        }
      }
    }

    err = RA_Payoff(RA, today, date, ra_dfs, &payoff);
    if (err)
      goto FREE_RETURN;

    fwdPV[i] = payoff;

    if (ra_dfs[0])
      free(ra_dfs[0]);
    ra_dfs[0] = NULL;
    if (ra_dfs[1])
      free(ra_dfs[1]);
    ra_dfs[1] = NULL;
    if (RA->cpn_type != 0) {
      if (ra_dfs[2])
        free(ra_dfs[2]);
      ra_dfs[2] = NULL;
    }
    if (ra_dfs)
      free(ra_dfs);
    ra_dfs = NULL;

    if (pn_dfs)
      free(pn_dfs);
    pn_dfs = NULL;
    if (pdates[0])
      free(pdates[0]);
    pdates[0] = NULL;
    if (pdates[1])
      free(pdates[1]);
    pdates[1] = NULL;
    if (RA->cpn_type != 0) {
      if (pdates[2])
        free(pdates[2]);
      pdates[2] = NULL;
    }
    if (pdates)
      free(pdates);
    pdates = NULL;
  }

FREE_RETURN:

  /*	if(ra_dfs[0]) free(ra_dfs[0]);
          ra_dfs[0] = NULL;
          if(ra_dfs[1]) free(ra_dfs[1]);
          ra_dfs[1] = NULL;
          if(ra_dfs) free(ra_dfs);
          ra_dfs = NULL;

          if(pn_dfs) free(pn_dfs);
          pn_dfs = NULL;
          if(pdates[0]) free(pdates[0]);
          pdates[0] = NULL;
          if(pdates[1]) free(pdates[1]);
          pdates[1] = NULL;
          if(pdates) free(pdates);
          pdates = NULL;
  */
  return err;
}

//-------------------------------------------------------------------------
//----------------------New Structures-------------------------------------
//-------------------------------------------------------------------------
/*
Err init_ra_leg(
                        char			*yc_d        ,
                        char			*yc_f        ,
                        char			*volc_d        ,
                        char			*volc_f        ,
                        char			*refrate_d        ,
                        char			*refrate_f        ,
                        double			notional        ,
                        int				n_fxvol_dates        ,
                        long			*fxvol_dates        ,
                        double			*fxvol        ,
                        double			*qtocorrel        ,
                        int				num_dates        ,
                        long			*start_dates        ,
                        long			*end_dates        ,
                        double			*cpns        ,
                        char			*basis        ,
                        int				*ra_nfixings        ,
                        long			**ra_fixingdates        ,
                        double			**ra_fixings        ,
//	Past coupon fixing if relevant double *upper_barr        , double
*lower_barr        , char			*recpay        , double
c_spread        , int				obs_freq        , double
rho_df        , int				eod_flag        ,	// 0: I
, 1: E ra_leg *RA)
{
        Err err = NULL;
        RA_OBS	*stp;
        RA_CPN	*cpn;
        int i        , j        , index;
        long d        , d1        , d2;
        SrtBasisCode basis_d        , basis_f        , ibasis;
        SrtCompounding freq_d        , freq_f;
        double vol_fx        , rho_ffx        , coef        , fra_f        ,
vol_f        , fra_d        , cvg_d        , vol_d        , power; double
quanto_corr        , DRS_corr; SrtReceiverType irecpay; SrtCurvePtr ycd_ptr;
        long ltoday;

        double weight;

        ycd_ptr = lookup_curve(yc_d);
        if (!ycd_ptr) return NULL;

        ltoday = get_today_from_curve(ycd_ptr);

        err = srt_f_get_spot_lag_from_refrate(refrate_f        ,
&RA->spotlag_f); if (err) goto RETURN_ERR;

        for ( index=0; index < num_dates &&
                           eod_flag + ltoday > add_unit( start_dates[index] ,
-RA->spotlag_f        , SRT_BDAY        , MODIFIED_SUCCEEDING ); index++ );

        if(index >= num_dates)
        {
                err = "All start dates are in the past in init_ra_leg";

                goto FREE_RETURN;
        }


        // Allocate memory and copy data to RA

        RA->notional = notional;

        RA->num_cpn = num_dates - index;

        err = swp_f_get_ref_rate_details(refrate_f        , &basis_f        ,
&freq_f); if (err) goto RETURN_ERR; RA->tenor_f = 12 / freq_f; RA->basis_f =
basis_f;

        sprintf(RA->tenor_f_char        , "%dM"        , RA->tenor_f);

        err = swp_f_get_ref_rate_details(refrate_d        , &basis_d        ,
&freq_d); if (err) goto RETURN_ERR;

        err = interp_basis(basis        , &ibasis);
        if (err) goto RETURN_ERR;

        err = interp_rec_pay(recpay        , &irecpay);
        if (err) goto RETURN_ERR;

        RA->cpn = (ra_cpn*) calloc (RA->num_cpn        , sizeof (ra_cpn) );

        for(i=0;i < num_dates - index;++i)
        {
                RA->cpn[i].start_date = start_dates[i];
                RA->cpn[i].start_time = (start_dates[i] - ltoday) *
YEARS_IN_DAY; RA->cpn[i].pay_date = end_dates[i];

                RA->cpn[i].num_obs = 1;
                RA->cpn[i].obs = (ra_obs*) calloc (1        , sizeof (ra_obs));
                RA->cpn[i].obs[0].fixing_date = ra_fixingdates[i+index][0];
                RA->cpn[i].obs[0].start_date = ra_fixingdates[i+index][0];
                RA->cpn[i].obs[0].end_date = add_unit(
RA->cpn[i].obs[0].start_date        , RA->tenor_f        , SRT_MONTH        ,
                                                                NO_BUSDAY_CONVENTION);
                RA->cpn[i].obs[0].spread =
swp_f_spread(RA->cpn[i].obs[0].start_date        , RA->cpn[i].obs[0].end_date ,
                                                                                                refrate_f);

                if(obs_freq < ra_nfixings[i+index])
                {
                        weight = ra_fixingdates[i+index][obs_freq] -
ra_fixingdates[i+index][0];
                }
                else
                {
                        weight = end_dates[i] - start_dates[i];
                }

                RA->cpn[i].obs[0].cxxpn = cpns[i]
                                * weight
                                / (end_dates[i] - start_dates[i])
                                * coverage(end_dates[i]        , start_dates[i]
, ibasis);

                for(j=obs_freq ; j < ra_nfixings[i+index]; j = j + obs_freq)
                {
                        RA->cpn[i].num_obs += 1;
                        RA->cpn[i].obs = (ra_obs*) realloc (RA->cpn[i].obs ,
RA->cpn[i].num_obs * sizeof (ra_obs));
                        RA->cpn[i].obs[(int)(j/obs_freq)].fixing_date =
ra_fixingdates[i+index][j]; RA->cpn[i].obs[(int)(j/obs_freq)].start_date =
add_unit( RA->cpn[i].obs[(int)(j/obs_freq)].fixing_date        , RA->spotlag_f ,
                                                                SRT_BDAY ,
                                                                MODIFIED_SUCCEEDING
); RA->cpn[i].obs[(int)(j/obs_freq)].end_date = add_unit(
RA->cpn[i].obs[(int)(j/obs_freq)].start_date        , RA->tenor_f        ,
SRT_MONTH        , NO_BUSDAY_CONVENTION);

                        if(j + obs_freq < ra_nfixings[i+index])
                        {
                                weight = ra_fixingdates[i+index][j + obs_freq ]
- ra_fixingdates[i+index][j];
                        }
                        else
                        {
                                weight = end_dates[i] -
ra_fixingdates[i+index][j];
                        }

                        RA->cpn[i].obs[0].cxxpn = cpns[i]
                                        * weight
                                        / (end_dates[i] - start_dates[i])
                                        * coverage(end_dates[i]        ,
start_dates[i]        , ibasis);
                }
        }

        // Fill in observation steps

        j = 0;
        for (i=0        , stp = RA->obs_steps[0]; stp->next; stp = stp->next)
        {
                if (add_unit( stp->date        , -RA->spotlag_f        ,
SRT_BDAY        , MODIFIED_SUCCEEDING ) >= RA->dates[i+1]) i++; d =
add_unit(stp->date        , RA->tenor_f        , SRT_MONTH        ,
NO_BUSDAY_CONVENTION); stp->spread = swp_f_spread(stp->date        , d        ,
refrate_f); stp->cpn = cpns[i] * (stp->next->date
- stp->date) / (RA->dates[i+1] - RA->dates[i])
                                                        * coverage(RA->dates[i]
      , RA->dates[i+1]        , ibasis);

                fra_f = swp_f_fra(stp->date        , d        , basis_f        ,
yc_f        , refrate_f); err = swp_f_SABRvol(volc_f        , stp->date        ,
d        , fra_f        , &vol_f        , &power        , SABR_LOGVOL); if (err)
                {
                        goto RETURN_ERR;
                }

                if (RA->quanto)		// Calculate quanto adjustment
                {
                        // Interpolate FX vol and rho_ffx at d
                        for (; j < n_fxvol_dates-1 && fxvol_dates[j+1] < d;
j++); if (fxvol_dates[j] > d || j == n_fxvol_dates-1)
                        {
                                vol_fx = fxvol[j];
                                rho_ffx = qtocorrel[j];
                        }
                        else
                        {
                                coef = ((double)(d - fxvol_dates[j])) /
                                        (fxvol_dates[j+1] - fxvol_dates[j]);
                                vol_fx = coef * fxvol[j+1] + (1.0-coef) *
fxvol[j]; rho_ffx = coef * qtocorrel[j+1] + (1.0-coef) * qtocorrel[j];
                        }
                        quanto_corr = - rho_ffx * vol_fx * vol_f;
                }
                else quanto_corr = 0.0;

                // Calculate DRS adjustment
                if (RA->dates[i+1] == d) DRS_corr = 0;
                else
                {
                        d1 = (RA->dates[i+1] < d ? RA->dates[i+1] : d);
                        d2 = (RA->dates[i+1] > d ? RA->dates[i+1] : d);

                        fra_d = swp_f_fra(d1        , d2        , basis_d , yc_d
, refrate_d); cvg_d = coverage(d1        , d2        , basis_d); err =
swp_f_SABRvol(volc_d        , d1        , d2        , fra_d        , &vol_d ,
&power        , SABR_LOGVOL); if (err)
                        {
                                goto RETURN_ERR;
                        }

                        DRS_corr = rho_df * vol_d * vol_f * cvg_d * fra_d / (1.0
+ cvg_d * fra_d); if (RA->dates[i+1] > d) DRS_corr = -DRS_corr;
                }
                stp->correction = quanto_corr + DRS_corr;

                // Calculate vols at barriers
                if (irecpay == SRT_RECEIVER)
                {
                        stp->K[0] = lower_barr[i];
                        stp->K[1] = lower_barr[i] + c_spread;
                        stp->K[2] = upper_barr[i] - c_spread;
                        stp->K[3] = upper_barr[i];
                }
                else
                {
                        stp->K[0] = lower_barr[i] - c_spread;
                        stp->K[1] = lower_barr[i];
                        stp->K[2] = upper_barr[i];
                        stp->K[3] = upper_barr[i] + c_spread;
                }

                // Get SABR parameters for the step
                err = swp_f_SABRvol( volc_f        , stp->date        , d , 0.0
, &stp->betavol        , &power        , SABR_BETAVOL ); if (err)
                {
                        goto RETURN_ERR;
                }
                err = swp_f_SABRvol( volc_f        , stp->date        , d , 0.0
, &stp->alpha        , &power        , SABR_ALPHA ); if (err)
                {
                        goto RETURN_ERR;
                }
                err = swp_f_SABRvol( volc_f        , stp->date        , d , 0.0
, &stp->beta        , &power        , SABR_BETA ); if (err)
                {
                        goto RETURN_ERR;
                }
                err = swp_f_SABRvol( volc_f        , stp->date        , d , 0.0
, &stp->rho        , &power        , SABR_RHO ); if (err)
                {
                        goto RETURN_ERR;
                }
        }


RETURN_ERR:
        if(err)
        {
                ra_free(RA);
        }
        return NULL;
}



*/