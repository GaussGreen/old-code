/* ===========================================================================

         FILENAME:     swp_f_swap_utils.cxx

     PURPOSE:      A few useful functions when dealing with a swap structure

   ===========================================================================
 */
#include "math.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_swap_compute.h"
#include "swp_h_swap_utils.h"
#include "swp_h_vol_interpol.h"

/* ------------- Functions operating on a Dlist
 * --------------------------------- */

Err new_Dlist(int len, Dlist *x) {

  x->d = (double *)(srt_calloc(len, sizeof(double)));
  if (x->d == NULL)
    return ("Fatal: calloc failed in new_Dlist");
  x->len = len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */
int free_inDlist(Dlist *dl) {
  srt_free(dl->d);
  dl->len = 0;

  return (0);
}

/* ---------------------------------------------------------------------------
 */
void Dlist_oper(Dlist l, double x, char oper) {
  int i;
  switch (oper) {
  case '+':
    for (i = 0; i < l.len; i++)
      l.d[i] += x;
    break;
  case '-':
    for (i = 0; i < l.len; i++)
      l.d[i] -= x;
    break;
  case '*':
    for (i = 0; i < l.len; i++)
      l.d[i] *= x;
    break;
  case '/':
    for (i = 0; i < l.len; i++)
      l.d[i] /= x;
    break;
  }
}

/* ---------------------------------------------------------------------------
 */
Err Dlist_func(Dlist l1, Dlist l2, char oper, Dlist *l3) {
  int i, l = IMIN(l1.len, l2.len);
  Err err = NULL;

  err = new_Dlist(l, l3);
  if (err)
    return (err);

  switch (oper) {
  case '+':
    for (i = 0; i < l; i++)
      l3->d[i] = l1.d[i] + l2.d[i];
    break;
  case '-':
    for (i = 0; i < l; i++)
      l3->d[i] = l1.d[i] - l2.d[i];
    break;
  case '*':
    for (i = 0; i < l; i++)
      l3->d[i] = l1.d[i] * l2.d[i];
    break;
  case '/':
    for (i = 0; i < l; i++)
      l3->d[i] = l1.d[i] / l2.d[i];
    break;
  }

  l3->len = l1.len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */
Err Dlist_copy(Dlist l1, Dlist *l2) {
  Err err = NULL;

  err = new_Dlist(l1.len, l2);
  if (err)
    return (err);

  memcpy(l2->d, l1.d, l1.len * sizeof(double));

  l2->len = l1.len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The time  correspond to time from today to swap dates in the date list
   Allocation is done inside.
   ------------------------------------------------------------------------- */
Err time_list(DateList dl, Date today, Dlist *t) {
  int i;
  Err err = NULL;

  err = new_Dlist(dl.len, t);
  if (err)
    return (err);

  for (i = 0; i < dl.len; i++)
    t->d[i] = (dl.date[i] - today) * YEARS_IN_DAY;

  t->len = dl.len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */
/* -------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The maturity  correspond to time from today to the fixing date as reported
   in the list (spot_lag bd before start date)
   (no maturity for last date: it is not a fixing )
   Allocation is done inside.
   ------------------------------------------------------------------------- */
Err maturity_list(DateList fixing_date, Date today, Dlist *t) {
  int i;
  Err err = NULL;

  err = new_Dlist(fixing_date.len, t);
  if (err)
    return (err);

  for (i = 0; i < fixing_date.len; i++) {
    t->d[i] = (double)(fixing_date.date[i] - today) * YEARS_IN_DAY;
  }

  t->len = fixing_date.len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */
/* -------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The fixing date corresponds to SpotLag business days before the start date
   Allocation is done inside.
   ------------------------------------------------------------------------- */
Err fixing_list(DateList start_date, int spot_lag, DateList *fixing) {
  int i;
  Err err = NULL;

  *fixing = new_DateList(start_date.len);
  if (err)
    return (err);

  for (i = 0; i < start_date.len; i++)
    fixing->date[i] =
        add_unit(start_date.date[i], -spot_lag, SRT_BDAY, MODIFIED_SUCCEEDING);

  fixing->len = start_date.len;

  return NULL;
}

/* ---------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The df are stored on each date of the datelist
   Allocation is done inside
   ------------------------------------------------------------------------- */
Err df_list(DateList pay_date, String ycname, Dlist *x) {
  int i;
  Err err = NULL;

  err = new_Dlist(pay_date.len, x);
  if (err)
    return (err);

  for (i = 0; i < pay_date.len; i++) {
    x->d[i] = swp_f_df(0, pay_date.date[i], ycname);
  }

  x->len = pay_date.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The cvg correspond to the coverage from start date to end date of each period
   Allocation is done inside
   ---------------------------------------------------------------------------
 */
Err cvg_list(DateList start_date, DateList end_date, BasisCode bc, Dlist *x) {
  int i;
  Err err = NULL;

  if (start_date.len != end_date.len)
    return serror("Not the same number of start and end dates in cvg_list");

  err = new_Dlist(start_date.len, x);
  if (err)
    return (err);

  for (i = 0; i < start_date.len; i++)
    x->d[i] = coverage(start_date.date[i], end_date.date[i], bc);

  x->len = start_date.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */

Err const_list(DateList coupon_date, double constant, Dlist *l) {
  int i;
  Err err = NULL;

  /* Allocates space for the list */
  err = new_Dlist(coupon_date.len, l);
  if (err)
    return (err);

  /* Sets all the constant values in the list */
  for (i = 0; i < coupon_date.len; i++) {
    l->d[i] = constant;
  }

  /* Returns the list */
  l->len = coupon_date.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
   Each coupon corresponds to the strike * accrual period in the right basis
   cpn->d[0] corresponds to the initial exchange only
   cpn->d[i] corresponds to the strike * accrual from Ti-1 to Ti
   cpn->d[len-1] also takes the final exchange into account
   (each coupon is paid at the dates stored in pay_date)
   ---------------------------------------------------------------------------
 */

Err cpn_list(BasisCode bc, SrtCompounding compd, DateList start_date,
             DateList end_date, DateList pay_date, double strike, double ini,
             double fin, StructType type, Dlist *cpn) {
  int i;
  Dlist cvg;
  double k;
  Err err = NULL;

  /* Check on the number of dates */
  if (start_date.len != end_date.len)
    return serror("Start and End dates number do not match");

  if (fabs(ini) > EPS)
    k = ini;
  else
    k = 1.0;

  strike *= k;

  err = new_Dlist(pay_date.len, cpn);
  if (err)
    return (err);

  switch (type) {
  case BOND:
  case BOND_OPTION:

    cpn->d[0] = -ini;
    for (i = 1; i < pay_date.len; i++) {
      cpn->d[i] = strike / ((double)compd);
    }
    cpn->d[pay_date.len - 1] += fin;
    break;
  case RESETCAPFLOOR:
    err = cvg_list(start_date, end_date, bc, &cvg);
    for (i = 1; i < pay_date.len; i++) {
      cpn->d[i] = cvg.d[i];
    }
    free_inDlist(&cvg);
    break;

  case SWAP:
  case SWAPTION:
    cpn->d[pay_date.len - 1] = fin;
    cpn->d[0] = -ini;
  case CAPFLOOR:
  default:
    err = cvg_list(start_date, end_date, bc, &cvg);
    for (i = 1; i < pay_date.len; i++) {
      cpn->d[i] += cvg.d[i] * strike;
    }
    free_inDlist(&cvg);
    break;
  }

  cpn->len = pay_date.len;

  return NULL;
}
/* ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The FRAs (cash) are computed in between each start and end date
   ( spot lag after fixing date )
   the payment will be on the corresponding end date
   (assumed to be three months apart)
   Allocation is done inside
   ---------------------------------------------------------------------------
 */
Err fwd_cash_list(DateList start_date, DateList end_date,
                  SrtBasisCode basiscode, String ycname, Dlist *fwd_cash) {
  int i;
  Err err = NULL;

  /* Check on the number of dates */
  if (start_date.len != end_date.len)
    return serror("Start and End dates number do not match");

  /* Allocates space for the list */
  err = new_Dlist(start_date.len, fwd_cash);
  if (err)
    return (err);

  /* Sets all the fra values in the list (not the last one: no more points) */
  for (i = 0; i < start_date.len; i++) {
    fwd_cash->d[i] =
        swp_f_fwdcash(start_date.date[i], end_date.date[i], basiscode, ycname);
    if (fwd_cash->d[i] == SRT_DF_ERROR)
      return ("Error: fwd_cash_list found a SRT_DF_ERROR");
  }

  /* Returns the list */
  fwd_cash->len = start_date.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
   DateList is the list of modified accrual period start and end dates
   The spreads (RefRate vs Cash fra) are computed for a period running from
   start date to the corresponding end date
   Allocation is done inside
   ---------------------------------------------------------------------------
 */
Err spread_list(DateList start_date, DateList end_date, String ref_rate_code,
                Dlist *spread) {
  int i;
  Err err = NULL;

  /* Allocates space for the list */
  err = new_Dlist(start_date.len, spread);
  if (err)
    return (err);

  /* Sets all the fra values in the list (not the first one: no more points) */
  for (i = 0; i < start_date.len; i++) {
    spread->d[i] =
        swp_f_spread(start_date.date[i], end_date.date[i], ref_rate_code);
    if (spread->d[i] == SRT_SPREAD_ERROR) {
      free_inDlist(spread);
      return ("Error: swp_f_spread");
    }
  }

  /* Returns the list of spreads */
  spread->len = start_date.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */

/* ------------------ These functions are mainly used for CMT        , CMS
 * ---------
 */

/* ------------------------------------------------------------------------- */
/* The treasury forwards are stored on the date where they are fixed:
   the payment will be on the next date in the list for a Libor leg
   Allocation is done inside
*/

Err fwd_treas_list(DateList dl, SrtCurvePtr cmt_crv, Dlist *fwd_treas) {
  int i;
  Err err = NULL;

  /* Allocates space for the fwd_treas*/
  err = new_Dlist(dl.len, fwd_treas);
  if (err)
    return (err);
  /* Computes all the forward        , though the very last one is useless */
  for (i = 0; i < dl.len; i++) {
    err = fwd_treas_rate(dl.date[i], cmt_crv, &(fwd_treas->d[i]));
    if (err)
      return (err);
  }

  /* Returns the list */
  fwd_treas->len = dl.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* The swap forwards are stored on the date where they are fixed:
   Allocation is done inside
*/

Err fwd_swap_list(DateList dl, SrtCurvePtr cmt_crv, Dlist *fwd_swap) {
  int i;
  CMT_Param *cmt_param;
  String cmt_name;
  String yc_name;
  SrtCurvePtr yc_crv;
  SwapDP swapdp;
  Err err = NULL;

  /* Allocates space for the fwd_swap*/
  err = new_Dlist(dl.len, fwd_swap);
  if (err)
    return (err);

  /* Gets cmt parameters */
  cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);
  cmt_name = get_curve_name(cmt_crv);
  yc_name = get_ycname_from_cmtcrv(cmt_crv);
  yc_crv = lookup_curve(yc_name);

  /* Sets in the swapdp the elements that are not fixing dependent */
  swapdp.nfp = (int)cmt_param->cmt_mat * cmt_param->swap_compd;
  swapdp.end = swapdp.nfp;
  swapdp.direction = FWD;
  swapdp.cxxompd = cmt_param->swap_compd;
  swapdp.basis_code = cmt_param->swap_basis_code;
  swapdp.spot_lag = cmt_param->spot_lag;

  /* Computes all the forward swaps        , though the very last one is useless
   */
  for (i = 0; i < dl.len; i++) {
    swapdp.start = dl.date[i];
    swapdp.first_full_fixing = swapdp.start;

    swp_f_ForwardRate_SwapDP(&swapdp, yc_name, cmt_param->swap_ref_rate,
                             &(fwd_swap->d[i]));
  }

  /* Make sure the crv points where it should */
  cmt_crv = lookup_curve(cmt_name);

  /* Returns the list */
  fwd_swap->len = dl.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* The swap volatilities are stored on the date where the underlying
        swaption expires
   Allocation is done inside */

Err swap_vol_list(DateList dl, Dlist swap_rate, SrtCurvePtr cmt_crv,
                  Dlist *swap_vol) {
  int i;
  long period_start, period_end;
  double strike;
  /*	double tgt_exp_date;	*/
  CMT_Param *cmt_param;
  Err err = NULL;

  /* Allocates space for the swap_vol*/
  err = new_Dlist(dl.len, swap_vol);
  if (err)
    return (err);
  /* Gets cmt parameters */
  cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);

  for (i = 0; i < dl.len; i++) {
    period_start = dl.date[i];
    period_end =
        add_unit((long)period_start, cmt_param->cmt_mat, SRT_YEAR, SUCCEEDING);
    strike = swap_rate.d[i];
    /* OVE fix: a swap is described by theoretical payment dates: spot lag after
       fixing        , vol is effective only until fixing date tgt_exp_date =
       add_unit((long)period_start        , -cmt_param->spot_lag        ,
       SRT_BDAY        , SUCCEEDING);
    */
    /* use the getvol function */
    if (cmt_param->flatvol != 0.) {
      swap_vol->d[i] = cmt_param->flatvol;
    } else {
      err = cmt_param->cmt_getvol((Date)period_start, (Date)period_end, strike,
                                  &swap_vol->d[i]);
      if (err)
        return err;
    }
  }

  /* Returns the list */
  swap_vol->len = dl.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* The treasury volatilities are stored on the date where the underlying
   option expires (i.e. same date as the swap vols
   Allocation is done inside */

Err treas_vol_list(DateList dl, Dlist swap_rate, Dlist treas_rate,
                   SrtCurvePtr cmt_crv, Dlist *treas_vol) {
  int i;
  long period_start, period_end;
  double strike;
  /*	double tgt_exp_date;	*/
  CMT_Param *cmt_param;
  SRT_Boolean use_prop_vol;
  Err err = NULL;

  /* Allocates space for the treas_vol*/
  err = new_Dlist(dl.len, treas_vol);
  if (err)
    return (err);
  /* Gets cmt parameters */
  cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);
  use_prop_vol = cmt_param->use_prop_vol_flg;

  /* Interpolates all the vols        , though the very last one is useless */
  for (i = 0; i < treas_rate.len; i++) {
    period_start = dl.date[i];
    period_end =
        add_unit((long)period_start, cmt_param->cmt_mat, SRT_YEAR, SUCCEEDING);
    if (use_prop_vol == SRT_YES) {
      strike = swap_rate.d[i];
    } else {
      strike = treas_rate.d[i];
    }

    /* OVE fix: a swap is described by theoretical payment dates: spot lag after
fixing        , vol is effective only until fixing date tgt_exp_date =
add_unit((long)period_start        , -cmt_param->spot_lag        , SRT_BDAY ,
SUCCEEDING);
    */

    /* use the getvol function */
    if (cmt_param->flatvol != 0.) {
      treas_vol->d[i] = cmt_param->flatvol;
    } else {
      err = cmt_param->cmt_getvol((Date)period_start, (Date)period_end, strike,
                                  &treas_vol->d[i]);
      if (err)
        return err;
    }

    if (use_prop_vol == SRT_YES) {
      treas_vol->d[i] *= swap_rate.d[i] / treas_rate.d[i];
    }
  }

  /* Returns the list */
  treas_vol->len = dl.len;

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* The convexity forward are stored on the date where they are paid        ,
   (because of the delay between fixing and payment in the cms function):
   the payment i s assumed to be on the next date in the list after fixing
   Allocation is done inside
*/

Err conv_fwd_list(long today, DateList dl, Dlist fwd_rate, Dlist vol,
                  SrtCurvePtr cmt_crv, CMRateType rate_type,
                  Dlist *conv_fwd_treas) {
  int i;
  CMT_Param *cmt_param;
  SrtCompounding compd;
  BasisCode basis;
  Err err;
  int nfp;
  long fixing_date;
  double delay;
  double maturity;
  double cm_rate;

  /* Allocates space for the conv_fwd_treas*/
  err = new_Dlist(dl.len, conv_fwd_treas);
  if (err)
    return (err);

  conv_fwd_treas->d[0] = 0.0;

  /* Gets cmt parameters  */
  cmt_param = get_cmtparam_from_cmtcrv(cmt_crv);

  /* Sets the parameters for the underlying constant maturity rate */
  if (rate_type == SWAP_RATE) {
    basis = cmt_param->swap_basis_code;
    compd = cmt_param->swap_compd;
  } else if (rate_type == TREASURY_RATE) {
    basis = cmt_param->bond_basis_code;
    compd = cmt_param->bond_compd;
  }

  /* A few dates and coverage calculations (assume dl.date[0] is today) */
  nfp = (int)((cmt_param->cmt_mat) * compd);

  /* Computes all the convexity adjusted forwards        , storing themon the
   * next one
   */
  for (i = 0; i < dl.len - 1; i++) {
    /* OVE fix: a swap is described by theoretical payment dates: spot lag after
       fixing        , vol is effective only until fixing date */
    fixing_date =
        add_unit(dl.date[i], -cmt_param->spot_lag, SRT_BDAY, SUCCEEDING);
    delay = (double)coverage(dl.date[i], dl.date[i + 1], basis);
    maturity = (double)(fixing_date - today) * YEARS_IN_DAY;
    err = swp_f_cmsrate(fwd_rate.d[i], nfp, cmt_param->bond_compd, vol.d[i],
                        maturity, delay, DEFAULT_CMS_DELTA, MAX_CMS_SWAPS,
                        SRT_LOGNORMAL, &cm_rate);
    if (err)
      cm_rate = -1.0e30;

    if (basis == BASIS_ACT_360)
      cm_rate = cm_rate * 360.0 / 365.0;

    conv_fwd_treas->d[i + 1] = cm_rate;
  }

  /* Returns the list */
  conv_fwd_treas->len = dl.len;

  return NULL;
}
