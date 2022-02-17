/* ------------------------------------------------------------------------
        FILE NAME				SRT_F_FWDVOLFNC.C

    PURPOSE					UYTILITY FUNCTIONS USED IN
   SRT_F_FWDVOL.C
   ------------------------------------------------------------------------ */

/*	SRT LIBRARY */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmresetcap.h"
/*	FOR SPREAD NUMERAIRE FUNCTION */
#include "opfnctns.h"
/*	FOR GRFN TABLEAU EVAL */
#include "grf_h_all.h"

#define SWAP(a, b)                                                             \
  {                                                                            \
    long tempr;                                                                \
    tempr = (a);                                                               \
    (a) = (b);                                                                 \
    (b) = tempr;                                                               \
  }

Err cmsfwd_vol_set_calibration_instruments(
    long num_period, int fix_tau_flag, long *vol_start, long *vol_end,

    char *und_tenor, SrtCompounding swap_compd, SrtBasisCode swap_basis,
    String ref_rate,

    String yc_id, String bs_vol_type,
    Err(*GetVol) /* GET VOL FUNCTION */
    (Ddate start, Ddate end, double strike, double dForward, double dSpread,
     double *bs_vol),
    long *cal_numinst, long **cal_start, long **cal_end, String **cal_freq,
    String **cal_basis, double **cal_str, double **cal_bndstr,
    String **cal_type, String **cal_recpay, String **cal_refrate,
    double **cal_price, double **cal_weights) {
  long i, j;
  double swap;
  Date spot_lag, spot_date;
  Date *short_swap_start_dates = NULL;
  Date *long_swap_start_dates = NULL;
  SrtCurvePtr yldcrv;
  Err err;
  SwapDP swapdp;
  SrtDiffusionType srt_vol_type;
  double ATMVol;
  int num_long_swaptions;

  yldcrv = lookup_curve(yc_id);
  spot_date = get_spotdate_from_yldcrv(yldcrv);
  spot_lag = get_spotlag_from_curve(yldcrv);

  long_swap_start_dates = lngvector(0, 2 * num_period - 1);

  if (fix_tau_flag == 0) /*---FREEZETAU = NO---*/
  {
    short_swap_start_dates = lngvector(0, num_period - 1);
  }

  if (fix_tau_flag == 0) /*---FREEZETAU = NO---*/
  {
    *cal_start = lngvector(0, 3 * num_period - 1);
  }

  else if (fix_tau_flag == 1) /*---FREEZETAU = YES---*/
  {
    *cal_start = lngvector(0, 2 * num_period - 1);
  }
  /* STORE ALL START DATES */

  for (i = 0; i < num_period; i++) {

    long_swap_start_dates[i] = DTOL(vol_start[i]);
    long_swap_start_dates[num_period + i] = DTOL(vol_end[i]);

    if (fix_tau_flag == 0) /*---FREEZETAU = NO---*/
    {
      short_swap_start_dates[i] = DTOL(vol_start[i]);
    }
  }

  /* LONG SWAPTION FIRSTt */

  /* SORT THE LONG SWAPTIONS START DATES */
  for (i = 0; i < 2 * num_period; i++) {
    for (j = i + 1; j < 2 * num_period; j++) {
      if (long_swap_start_dates[j] < long_swap_start_dates[i]) {
        SWAP(long_swap_start_dates[j], long_swap_start_dates[i]);
      }
    }
  }

  j = 0;
  while (long_swap_start_dates[j] <= spot_date) {
    j++;
  }

  *cal_start[0] = long_swap_start_dates[j];
  *cal_numinst = 1;

  for (i = j + 1; i < 2 * num_period; i++) {
    if (long_swap_start_dates[i] != (*cal_start)[*cal_numinst - 1]) {
      (*cal_numinst)++;
      (*cal_start)[*cal_numinst - 1] = long_swap_start_dates[i];
    }
  }

  num_long_swaptions = (*cal_numinst);
  free_lngvector(long_swap_start_dates, 0, 2 * num_period - 1);
  long_swap_start_dates = NULL;

  if (fix_tau_flag == 0) {
    /* SHORT SWAPTION CASE */

    /* SORT SHORT SWAP DATES  */
    for (i = 0; i < num_period; i++) {
      for (j = i + 1; j < num_period; j++) {
        if (short_swap_start_dates[j] < short_swap_start_dates[i]) {
          SWAP(short_swap_start_dates[j], short_swap_start_dates[i]);
        }
      }
    }

    j = 0;
    while (short_swap_start_dates[j] <= spot_date) {
      j++;
    }

    (*cal_start)[*cal_numinst] = short_swap_start_dates[j];
    (*cal_numinst) += 1;

    for (i = j + 1; i < num_period; i++) {
      if (short_swap_start_dates[i] != (*cal_start)[*cal_numinst - 1]) {
        (*cal_numinst)++;
        (*cal_start)[*cal_numinst - 1] = short_swap_start_dates[i];
      }
    }
  }

  *cal_end = lngvector(0, *cal_numinst - 1);
  *cal_freq = (String *)malloc((*cal_numinst) * sizeof(String));
  *cal_basis = (String *)malloc((*cal_numinst) * sizeof(String));
  *cal_type = svector_size(0, *cal_numinst - 1, 32);
  *cal_recpay = svector_size(0, *cal_numinst - 1, 32);
  *cal_refrate = svector_size(0, *cal_numinst - 1, 32);
  *cal_str = dvector(0, *cal_numinst - 1);
  *cal_bndstr = dvector(0, *cal_numinst - 1);
  *cal_price = dvector(0, *cal_numinst - 1);
  *cal_weights = dvector(0, *cal_numinst - 1);

  /* LONG SWAPTION CASE */

  for (i = 0; i < num_long_swaptions; i++) {
    (*cal_start)[i] =
        add_unit((*cal_start)[i], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

    add_tenor(DTOL((*cal_start)[i]), und_tenor, NO_BUSDAY_CONVENTION,
              &(*cal_end)[i]);

    err = translate_compounding(*cal_freq + i, swap_compd);
    if (err)
      return err;

    err = translate_basis(*cal_basis + i, swap_basis);
    if (err)
      return err;

    strcpy((*cal_type)[i], "swaption");
    strcpy((*cal_recpay)[i], "rec");
    strcpy((*cal_refrate)[i], ref_rate);

    err = swp_f_ForwardRate((*cal_start)[i], (*cal_end)[i], (*cal_freq)[i],
                            (*cal_basis)[i], yc_id, (*cal_refrate)[i], &swap);
    if (err)
      return err;

    (*cal_str)[i] = swap;
    (*cal_bndstr)[i] = 1.00;

    err = swp_f_initSwapDP((*cal_start)[i], (*cal_end)[i], (*cal_freq)[i],
                           (*cal_basis)[i], &swapdp);

    if (err)
      return err;

    err = GetVol((*cal_start)[i], (*cal_end)[i], swap, swap, 0, &ATMVol);
    if (err)
      return err;

    err = interp_diffusion_type(bs_vol_type, &srt_vol_type);
    if (err)
      return err;

    err = swp_f_Swaption_SwapDP(&swapdp, ATMVol, (*cal_str)[i], SRT_RECEIVER,
                                (*cal_refrate)[i], yc_id, PREMIUM, srt_vol_type,
                                &((*cal_price)[i]));
    if (err)
      return err;

    (*cal_weights)[i] = 1.0;
  }

  /* SHORT SWAPTION CASE - NBR OF SHORT SWAPTIONS = NBR OF PERIOD */

  for (i = num_long_swaptions; i < (*cal_numinst); i++) {
    (*cal_start)[i] = add_unit(short_swap_start_dates[i - num_long_swaptions],
                               spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);
    (*cal_end)[i] = DTOL(vol_end[i - num_long_swaptions]);

    err = translate_compounding(*cal_freq + i, swap_compd);
    if (err)
      return err;

    err = translate_basis(*cal_basis + i, swap_basis);
    if (err)
      return err;

    strcpy((*cal_type)[i], "swaption");
    strcpy((*cal_recpay)[i], "rec");
    strcpy((*cal_refrate)[i], ref_rate);

    err = swp_f_ForwardRate((*cal_start)[i], (*cal_end)[i], (*cal_freq)[i],
                            (*cal_basis)[i], yc_id, (*cal_refrate)[i], &swap);
    if (err) {
      return err;
    }
    (*cal_str)[i] = swap;
    (*cal_bndstr)[i] = 1.00;

    err = swp_f_initSwapDP((*cal_start)[i], (*cal_end)[i], (*cal_freq)[i],
                           (*cal_basis)[i], &swapdp);

    err = GetVol((*cal_start)[i], (*cal_end)[i], swap, swap, 0, &ATMVol);

    err = interp_diffusion_type(bs_vol_type, &srt_vol_type);

    err = swp_f_Swaption_SwapDP(&swapdp, ATMVol, (*cal_str)[i], SRT_RECEIVER,
                                (*cal_refrate)[i], yc_id, PREMIUM, srt_vol_type,
                                &((*cal_price)[i]));

    if (err)
      return err;

    (*cal_weights)[i] = 1.0;
  }

  if (short_swap_start_dates)
    free_lngvector(short_swap_start_dates, 0, num_period - 1);
  short_swap_start_dates = NULL;
  return NULL;
}

Err capfwd_vol_set_calibration_instruments(
    long num_period, long *vol_start, long *vol_end, String ref_rate,
    SrtCompounding compd, SrtBasisCode basis, String yc_id, String bs_vol_type,
    Err (*GetVol)(Ddate start, Ddate end, double strike, double dForward,
                  double dSpread, double *bs_vol),
    long *cal_numinst, Date **cal_start, Date **cal_end, String **cal_freq,
    String **cal_basis, double **cal_str, double **cal_bndstr,
    String **cal_type, String **cal_recpay, String **cal_refrate,
    double **cal_price, double **cal_weights)

{
  long i, j;
  double fra;
  Date spot_lag, spot_date, *caplet_start_dates;
  SrtCurvePtr yldcrv;
  Err err;

  /*	Determine start dates (No memory error assumed) */
  caplet_start_dates = lngvector(0, 2 * num_period - 1);
  *cal_start = lngvector(0, 2 * num_period - 1);

  /*	Fill large array of caplet start dates = vol start dates + vol end dates
   */
  for (i = 0; i < num_period; i++) {
    caplet_start_dates[i] = DTOL(vol_start[i]);
    caplet_start_dates[num_period + i] = DTOL(vol_end[i]);
  }

  /* Sorts dates: from smallest (i=0) to biggest (i=2*num_period-1) */
  for (i = 0; i < 2 * num_period; i++) {
    for (j = i + 1; j < 2 * num_period; j++) {
      if (caplet_start_dates[j] < caplet_start_dates[i]) {
        SWAP(caplet_start_dates[j], caplet_start_dates[i]);
      }
    }
  }

  /* Skip dates < spot and store number of skipped dates in j */
  yldcrv = lookup_curve(yc_id);
  spot_date = get_spotdate_from_yldcrv(yldcrv);
  spot_lag = get_spotlag_from_curve(yldcrv);

  j = 0;
  while (caplet_start_dates[j] <= spot_date) {
    j++;
  }
  (*cal_start)[0] = caplet_start_dates[j];
  *cal_numinst = 1;
  for (i = j + 1; i < 2 * num_period; i++) {
    if (caplet_start_dates[i] != (*cal_start)[*cal_numinst - 1]) {
      (*cal_numinst)++;
      (*cal_start)[*cal_numinst - 1] = caplet_start_dates[i];
    }
  }
  free_lngvector(caplet_start_dates, 0, 2 * num_period - 1);

  /*	No memory error assumed */
  *cal_end = lngvector(0, *cal_numinst - 1);
  *cal_freq = (String *)malloc((*cal_numinst) * sizeof(String));
  *cal_basis = (String *)malloc((*cal_numinst) * sizeof(String));
  *cal_type = svector_size(0, *cal_numinst - 1, 32);
  *cal_recpay = svector_size(0, *cal_numinst - 1, 32);
  *cal_refrate = svector_size(0, *cal_numinst - 1, 32);
  *cal_str = dvector(0, *cal_numinst - 1);
  *cal_bndstr = dvector(0, *cal_numinst - 1);
  *cal_price = dvector(0, *cal_numinst - 1);
  *cal_weights = dvector(0, *cal_numinst - 1);

  for (i = 0; i < *cal_numinst; i++) {
    (*cal_start)[i] =
        add_unit((*cal_start)[i], spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

    (*cal_end)[i] = add_unit(DTOL((*cal_start)[i]), 12 / (int)compd, SRT_MONTH,
                             NO_BUSDAY_CONVENTION);

    err = translate_compounding(*cal_freq + i, compd);
    if (err)
      return err;

    err = translate_basis(*cal_basis + i, basis);
    if (err)
      return err;

    strcpy((*cal_type)[i], "capfloor");
    strcpy((*cal_recpay)[i], "cap");
    strcpy((*cal_refrate)[i], ref_rate);

    err = swp_f_ForwardRate((*cal_start)[i], (*cal_end)[i], (*cal_freq)[i],
                            (*cal_basis)[i], yc_id, (*cal_refrate)[i], &fra);
    if (err) {
      return err;
    }
    (*cal_str)[i] = fra;
    (*cal_bndstr)[i] = 1.00;

    err = swp_f_CapFloor((*cal_start)[i], (*cal_end)[i], (*cal_str)[i], GetVol,
                         (*cal_recpay)[i], (*cal_refrate)[i], yc_id, "premium",
                         bs_vol_type, &((*cal_price)[i]));

    if (err)
      return err;

    if (err)
      return err;
    (*cal_weights)[i] = 1.0;
  }

  return NULL;
}

/* ----------------------------------------------------------------------------
 */

void fwd_vol_set_grfn_options(SrtMdlType mdl_type, String **grfn_param,
                              String **grfn_value, long *num_grfn_param) {
  *num_grfn_param = 7;

  /*	No memory error assumed */
  *grfn_param = svector_size(0, *num_grfn_param - 1, 32);
  *grfn_value = svector_size(0, *num_grfn_param - 1, 32);
  strcpy((*grfn_param)[0], "FORCEMC");
  strcpy((*grfn_value)[0], "NO");
  strcpy((*grfn_param)[1], "MAXTIME");
  strcpy((*grfn_value)[1], "0.25");
  strcpy((*grfn_param)[2], "MINNODE");
  strcpy((*grfn_value)[2], "25");
  strcpy((*grfn_param)[3], "SAMPLETYPE");
  strcpy((*grfn_value)[3], "BALANTISAM");
  strcpy((*grfn_param)[4], "NUMPATH");
  strcpy((*grfn_value)[4], "1000");
  strcpy((*grfn_param)[5], "RENORM");
  strcpy((*grfn_value)[5], "YES");
  strcpy((*grfn_param)[6], "RANDSEED");
  strcpy((*grfn_value)[6], "-123456789");
}

/* ----------------------------------------------------------------------------
 */

void fwd_vol_set_calib_options(SrtMdlType mdl_type, int option_type,
                               int fix_tau_flag, String **calib_param,
                               String **calib_value, long *num_calib_param) {
  *num_calib_param = 7;

  /*	No memory error assumed */
  *calib_param = svector_size(0, *num_calib_param - 1, 32);
  *calib_value = svector_size(0, *num_calib_param - 1, 32);
  strcpy((*calib_param)[0], "CALIBTYPE");
  strcpy((*calib_value)[0], "FIXED");
  strcpy((*calib_param)[1], "CALIBALGO");
  strcpy((*calib_value)[1], "LEVENBERG");
  strcpy((*calib_param)[2], "NITER");
  strcpy((*calib_value)[2], "3");
  strcpy((*calib_param)[3], "ONETAU");
  strcpy((*calib_value)[3], "YES");
  strcpy((*calib_param)[4], "FREEZETAU");

  if (option_type == RESETCAPFLOOR) {
    strcpy((*calib_value)[4], "NO");
  }
  if ((option_type == RESETCMSOPTION) & (fix_tau_flag == 0)) {
    strcpy((*calib_value)[4], "NO");
  }

  else if ((option_type == RESETCMSOPTION) & (fix_tau_flag == 1)) {
    strcpy((*calib_value)[4], "YES");
  }

  strcpy((*calib_param)[5], "LGMSTARTINGPOINT");
  strcpy((*calib_value)[5], "YES");
  strcpy((*calib_param)[6], "LGMITER");
  strcpy((*calib_value)[6], "3");
}

/* ----------------------------------------------------------------------------
 */

/*	Forward Reset Caplet is the following profile:
        At T4        , pay: MAX (F(T2        ,T3        ,T4) - F(T1        ,T3
   ,T4)        , 0) * cvg (T3        ,
   T4) Where	T1 = vol start date T2 = vol end date T3 = fra start date T4 =
   fra end date */
Err srt_f_fwd_resetcaplet(SrtUndPtr undptr, SrtGrfnParam *grfnparam,
                          String yc_id, Date today, String ref_rate,
                          SrtMdlType mdl_type, SrtReceiverType rec_pay,
                          Ddate str_fix, Ddate spot_fix, Ddate start_act,
                          Ddate end_th, String compdStr, String basisStr,
                          double *premium) {
  double start_end[2];
  double cvg, level;
  SrtBasisCode basis;
  Date *eventdates;
  GrfnCell **sprdsht;
  SrtIOStruct *iolist;
  char *und_name = get_underlying_name(undptr);
  Err err;
  long ncol;
  long nrow;

  start_end[0] = start_act;
  start_end[1] = bus_date_method(DTOL(end_th), MODIFIED_SUCCEEDING);
  err = interp_basis(basisStr, &basis);
  if (err) {
    return err;
  }
  cvg = coverage(DTOL(start_end[0]), DTOL(start_end[1]), basis);

  if (mdl_type == LGM)
  /* LGM case: call srt_f_lgm_resetcaplet */
  {
    err = swp_f_LevelPayment(DTOL(start_act), DTOL(end_th), compdStr, basisStr,
                             yc_id, ref_rate, &level);
    if (err) {
      return err;
    }
    undptr = lookup_und(und_name);
    err = srt_f_lgm_resetcaplet(undptr, str_fix, start_end, cvg, spot_fix,
                                start_end, cvg, start_end[1], level, rec_pay,
                                premium, 0.0);
    if (err) {
      return err;
    }
  } else
  /* CHEY case: make grfn tableau */
  {
    eventdates = (Date *)srt_calloc(3, sizeof(Date));
    eventdates[0] = DTOL(str_fix);
    eventdates[1] = DTOL(spot_fix);
    eventdates[2] = DTOL(start_end[1]);
    sprdsht = GrfnCellmatrix(3, 2, GRFN_DEF_ARGBUFSZ);

    /*---------------------------C[0]--------------------------------*/

    sprintf(sprdsht[0][0].sval,
            "fra(%d        ,%d        ,\"%s\"        ,\"%s\"        ,\"%s\")",
            (long)start_act, (long)end_th, basisStr, und_name, ref_rate);
    sprdsht[0][0].type = GRFNSCELL;

    sprintf(sprdsht[1][0].sval,
            "fra(%d        ,%d        ,\"%s\"        ,\"%s\"        ,\"%s\")",
            (long)start_act, (long)end_th, basisStr, und_name, ref_rate);
    sprdsht[1][0].type = GRFNSCELL;

    sprintf(sprdsht[2][0].sval, "0.0");
    sprdsht[2][0].type = GRFNSCELL;

    /*---------------------------C[1]--------------------------------*/

    sprintf(sprdsht[0][1].sval, "0.0");
    sprdsht[0][1].type = GRFNSCELL;

    sprintf(sprdsht[1][1].sval, "0.0");
    sprdsht[1][1].type = GRFNSCELL;

    sprintf(sprdsht[2][1].sval,
            "max(c[0        ,1]-c[0        ,0]        ,0)*%.12f", cvg);
    sprdsht[2][1].type = GRFNSCELL;

    err = srt_f_IOstructcreate(&iolist, "xxx");

    if (!err) {
      ncol = 2;
      nrow = 3;
      err = srt_f_grfn(undptr, grfnparam, 3, &eventdates, &nrow, &ncol,
                       &sprdsht, 0, 0, 0, 0, 0, iolist, 0, 0);
    }

    if (!err) {
      err = srt_f_IOstructgetpremiumval(*iolist, premium);
    }
    if (!err) {
      err = srt_f_IOstructfree(&iolist);
    }
    if (eventdates) {
      srt_free(eventdates);
    }
    if (sprdsht) {
      grfn_free_GrfnCellmatrix(sprdsht, 3, 2);
    }
  }

  return err;
}

/*	Forward Reset CMS Option is the following profile:
        At T3        , pay: MAX (CMS(T2        ,T3        ,T4) - CMS(T1 ,T3 ,T4)
   , 0) Where	T1 = vol start date = strike_fix T2 = vol end date = CMS_fix T3
   = CMS start date = start_act T4 = CMS end date = end_theo */

Err srt_f_fwd_resetcms(SrtUndPtr irundptr, SrtGrfnParam *grfnparam,
                       Date strike_fix, Date cms_fix, Date start_act,
                       Date end_theo, String basisStr, String compStr,
                       String ref_rate, String yc_id,
                       Err (*GetVol)(Ddate, Ddate, double, double, double,
                                     double *),
                       String bs_vol_type, double *premium)

{

  Date *eventdates;
  GrfnCell **sprdsht;
  String ir_und_name;
  double cmsvol;
  double swap;
  Err err;
  SrtIOStruct *iolist;
  long ncol;
  long nrow;
  double SwpSpread;

  err = swp_f_ForwardRate(start_act, end_theo, compStr, basisStr, yc_id,
                          ref_rate, &swap);

  SwpSpread = swp_f_spread(start_act, end_theo, ref_rate);
  err = GetVol(start_act, end_theo, swap, swap, SwpSpread, &cmsvol);

  ir_und_name = get_underlying_name(irundptr);

  eventdates = (Date *)srt_calloc(3, sizeof(Date));
  eventdates[0] = DTOL(strike_fix);
  eventdates[1] = DTOL(cms_fix);
  eventdates[2] = add_unit(DTOL(start_act), 0, SRT_DAY, MODIFIED_SUCCEEDING);

  sprdsht = GrfnCellmatrix(3, 2, GRFN_DEF_ARGBUFSZ);

  /*---------------------------C[0]--------------------------------*/

  sprintf(sprdsht[0][0].sval,
          "CMS(%d        ,%d        ,\"%s\"        ,\"%s\"        ,%lf        "
          ",%d        ,\"%s\"        ,\"%s\"        ,\"%s\")",
          start_act, end_theo, compStr, basisStr, cmsvol, start_act,
          bs_vol_type, ir_und_name, ref_rate);
  sprdsht[0][0].type = GRFNSCELL;

  sprintf(sprdsht[1][0].sval,
          "CMS(%d        ,%d        ,\"%s\"        ,\"%s\"        ,%lf        "
          ",%d        ,\"%s\"        ,\"%s\"        ,\"%s\")",
          start_act, end_theo, compStr, basisStr, cmsvol, start_act,
          bs_vol_type, ir_und_name, ref_rate);
  sprdsht[1][0].type = GRFNSCELL;

  sprintf(sprdsht[2][0].sval, "0.0");
  sprdsht[2][0].type = GRFNSCELL;

  /*---------------------------C[1]--------------------------------*/

  sprintf(sprdsht[0][1].sval, "0.0");
  sprdsht[0][1].type = GRFNSCELL;

  sprintf(sprdsht[1][1].sval, "0.0");
  sprdsht[1][1].type = GRFNSCELL;

  sprintf(sprdsht[2][1].sval, "max(c[0        ,1]-c[0        ,0]        ,0)");
  sprdsht[2][1].type = GRFNSCELL;

  err = srt_f_IOstructcreate(&iolist, "xxx");

  if (!err) {
    ncol = 2;
    nrow = 3;
    err = srt_f_grfn(irundptr, grfnparam, 3, &eventdates, &nrow, &ncol,
                     &sprdsht, 0, 0, 0, 0, 0, iolist, 0, 0);
  }

  if (!err) {
    err = srt_f_IOstructgetpremiumval(*iolist, premium);
  }
  if (!err) {
    err = srt_f_IOstructfree(&iolist);
  }
  if (eventdates) {
    srt_free(eventdates);
  }
  if (sprdsht) {
    grfn_free_GrfnCellmatrix(sprdsht, 3, 2);
  }

  return err;
}

/* ----------------------------------------------------------------------------
 */

/*	Prices half atm parabola with smile */
Err parabola(Date today, Date start, SrtCompounding compd, double fwd_cms,
             SrtCallPutType call_put,
             Err (*GetVol)(Ddate, Ddate, double, double, double, double *),
             char *vol_type, double *answer, long max_strike, long max_vol,
             double delta_strike, double tol, long nvol) {
  long MyWay = (call_put == SRT_CALL ? 1 : -1);
  long j, k, l;
  Date end;
  SrtDiffusionType srt_vol_type;
  double strike = fwd_cms;
  double maturity = (start - today) * YEARS_IN_DAY;
  double opt, vol;
  double amt = 2 * delta_strike;
  Err err;

  *answer = 0.0;

  end = add_unit(start, 12 / (int)compd, SRT_MONTH, MODIFIED_SUCCEEDING);

  if (err = GetVol(start, end, strike, 1.0, 0.0, &vol)) {
    return err;
  }

  if (err = interp_diffusion_type(vol_type, &srt_vol_type)) {
    return err;
  }

  if (srt_vol_type == SRT_NORMAL) {
    opt =
        srt_f_optblknrm(fwd_cms, strike, vol, maturity, 1.0, call_put, PREMIUM);
  } else {
    opt =
        srt_f_optblksch(fwd_cms, strike, vol, maturity, 1.0, call_put, PREMIUM);
  }

  opt *= delta_strike;

  *answer += opt;

  strike += MyWay * delta_strike;

  k = 0;
  l = 1;
  for (j = 2; j <= max_strike; j++) {

    k = (k + 1) % nvol;

    if ((k == 0) && (l <= max_vol)) {
      if (err = GetVol(start, end, strike, 1.0, 0.0, &vol)) {
        return err;
      }

      l++;
    }

    if (srt_vol_type == SRT_NORMAL) {
      opt = srt_f_optblknrm(fwd_cms, strike, vol, maturity, 1.0, call_put,
                            PREMIUM);
    } else {
      opt = srt_f_optblksch(fwd_cms, strike, vol, maturity, 1.0, call_put,
                            PREMIUM);
    }

    *answer += amt * opt;

    strike += MyWay * delta_strike;

    if (fabs(opt) < tol)
      break;
  }

  return NULL;
}

/* ----------------------------------------------------------------------------
 */
