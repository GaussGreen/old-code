
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"

/*	Functions used to deal with the Fx smile using a GRFN tableau generation
 */

/*	A reduced version of cpd_fill_check_all_struct that does only the cpd
 * (product) structure */
Err cpd_fill_check_prod_struct(
    /*	Today's date */
    long today,
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I        , 1: E */
    int eod_ex_flag,  /*	0: I        , 1: E */
    /*	The structure */
    /*		funding */
    double fund_not, int fund_ccy, /*	0: domestic        , 1: foreign */
    int fund_ncpn, long *fund_fix, long *fund_start, long *fund_pay,
    char **fund_basis, double *fund_spr, double *fund_mrg,
    /*		pd */
    double pd_not, int pd_ncpn, long *pd_fix, long *pd_start, long *pd_pay,
    char **pd_basis, double *pd_alpha, double *pd_beta, int *pd_floored,
    double *pd_floor, int *pd_capped, double *pd_cap,
    /*		pd not refund */
    long pd_not_ref_fix, double pd_not_ref_alpha, double pd_not_ref_beta,
    int pd_not_ref_floored, double pd_not_ref_floor, int pd_not_ref_capped,
    double pd_not_ref_cap,
    /*		calls */
    int *call_type,         /*	0: Call        , 1: KO */
    int ncall, int pay_rec, /*	0: rec pd        , 1: pay pd */
    long *ex_date, long *set_date, double *barrier, /*	KO only */
    int *bar_type, /*	0: up and in        , 1: down and in */
    double *fees,  /*  fees if deal is called in domestic currency */
    double smooth,
    /*	Result */
    CPD_STR cpd) {
  Err err = NULL;

  double *fund_not_s = NULL, *pd_not_s = NULL;

  int i;

  /*	Initialisation */
  cpd->fund_leg = NULL;
  cpd->pd_leg = NULL;
  cpd->call = NULL;

  /*	Funding leg */

  cpd->fund_leg = (PD_FUND_LEG)malloc(sizeof(pd_fund_leg));
  if (!cpd->fund_leg) {
    err = "Memory allocation error (1) in cpd_fill_check_prod_struct";
    goto FREE_RETURN;
  }

  err = cpd_fill_fund_leg(today, eod_fix_flag, fund_not, fund_ccy, fund_ncpn,
                          fund_fix, fund_start, fund_pay, fund_basis, fund_spr,
                          fund_mrg, cpd->fund_leg);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Exotic leg */

  cpd->pd_leg = (PD_EXO_LEG)malloc(sizeof(pd_exo_leg));
  if (!cpd->pd_leg) {
    err = "Memory allocation error (2) in cpd_fill_check_prod_struct";
    goto FREE_RETURN;
  }

  err = cpd_fill_exo_leg(
      today, eod_fix_flag, pd_not, pd_ncpn, pd_fix, pd_start, pd_pay, pd_basis,
      pd_alpha, pd_beta, pd_floored, pd_floor, pd_capped, pd_cap, 0, NULL, NULL,
      NULL, NULL, NULL, pd_not_ref_fix, pd_not_ref_alpha, pd_not_ref_beta,
      pd_not_ref_floored, pd_not_ref_floor, pd_not_ref_capped, pd_not_ref_cap,
      0, 0, 0.0, 0.0, NULL, NULL, cpd->pd_leg);
  if (err) {
    goto FREE_RETURN;
  }

  /*	Calls */

  if (ncall > 0 && ex_date[ncall - 1] >= today + eod_ex_flag) {
    fund_not_s = (double *)calloc(ncall, sizeof(double));
    pd_not_s = (double *)calloc(ncall, sizeof(double));

    if (!fund_not_s || !pd_not_s) {
      goto FREE_RETURN;
    }

    for (i = 0; i < ncall; i++) {
      fund_not_s[i] = fund_not;
      pd_not_s[i] = pd_not;
    }

    err = cpd_fill_calls(today, eod_ex_flag, ncall, call_type, pay_rec, ex_date,
                         set_date, NULL, fund_not_s, pd_not_s, 0, NULL, NULL,
                         NULL, NULL, NULL, barrier, bar_type, fees, 0, smooth,
                         1, 0, 0, cpd);
    if (err) {
      goto FREE_RETURN;
    }
  } else {
    cpd->num_calls = 0;
    cpd->call = NULL;
  }

FREE_RETURN:

  if (err) {
    cpd_free_prod_struct(cpd);
  }

  if (fund_not_s) {
    free(fund_not_s);
  }

  if (pd_not_s) {
    free(pd_not_s);
  }

  return err;
}

/*	Free product structure */
Err cpd_free_prod_struct(CPD_STR cpd) {
  cpd_free_calls(cpd);

  if (cpd->fund_leg) {
    cpd_free_fund_leg(cpd->fund_leg);
    free(cpd->fund_leg);
  }

  if (cpd->pd_leg) {
    cpd_free_exo_leg(cpd->pd_leg);
    free(cpd->pd_leg);
  }

  return NULL;
}

static int set_mask(char *str) {
  int i, n = strlen(str);
  int b, e;
  int dot;

  b = 0;
  while (b < n - 1 && isspace(str[b])) {
    b++;
  }

  e = n - 1;
  while (e > 0 && isspace(str[e])) {
    e--;
  }

  /*	Blank */
  if (e <= b) {
    return GRFNBCELL;
  }

  if (str[b] == '+' || str[b] == '-') {
    b++;
  }

  dot = 0;
  for (i = b; i <= e; i++) {
    if (isdigit(str[i]))
      continue;

    if (str[i] == '.') {
      if (dot) {
        /*	String */
        return GRFNSCELL;
      } else {
        dot = 1;
      }
    } else {
      /*	String */
      return GRFNSCELL;
    }
  }

  /*	Number */
  return GRFNDCELL;
}

/*	Function that produces a GRFN tableau for a CPD deal */
Err cpd_make_grfn_tableau(
    /*	The structure */
    CPD_STR cpd,
    /*	Extra cash-flows coming from past coupons and initial exchange */
    int fund_ncf, long *fund_cfdtes, double *fund_cf, int pd_ncf,
    long *pd_cfdtes, double *pd_cf,
    /*	The underlyings */
    long today, char *dom_name, char *for_name, char *fx_name,
    /*	The tableau and auxiliaries
            Must be allocated with maximum values
                    and initialised with zeros */
    int *num_evt_dtes, int *num_cols, int *aux_rows, int *aux_cols,
    long *evt_dtes, char ***evt, int **mask, double **aux) {
  PD_EXO_LEG exo_leg = cpd->pd_leg;
  PD_FUND_LEG fund_leg = cpd->fund_leg;
  PD_CALL call = cpd->call;
  int i, j, k;
  int ncpn, nexo = exo_leg->num_cpn, nfund = fund_leg->num_cpn,
            ncall = cpd->num_calls, ncols = 5;
  int today_flag, first_ex_flag;
  char *fund_name, tmpstr[MAX_STR_LEN];
  char *dom_yc, *fund_yc;
  double pvcf;
  SrtUndPtr tmpund;
  Err err = NULL;

  /*	Initialise and make event dates = today + exo coupons fixing dates
   */

  if (nexo > 0 && exo_leg->cpn[0].fx_fix_date > today) {
    today_flag = 1;
  } else {
    today_flag = 0;
  }

  if (ncall > 0 && nexo > 0 && call[0].ex_date > today &&
      call[0].ex_date < exo_leg->cpn[0].fx_fix_date - 60) {
    first_ex_flag = 1;
  } else {
    first_ex_flag = 0;
  }

  *num_cols = ncols;
  *aux_cols = 16;

  *num_evt_dtes = ncpn = nexo + today_flag + first_ex_flag;
  evt_dtes[0] = today;
  if (first_ex_flag)
    evt_dtes[1] = call[0].ex_date;
  for (i = 0; i < nexo; i++) {
    evt_dtes[today_flag + first_ex_flag + i] = exo_leg->cpn[i].fx_fix_date;
  }

  /*	Column 0: fx fixings */

  aux_rows[2] = ncpn;
  for (i = 0; i < ncpn; i++) {
    aux[i][2] = add_unit(evt_dtes[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    sprintf(evt[i][0],
            "S(\"%s\")*DF(D[I]        ,A[2        ,I]        ,\"%s\")/DF(D[I]  "
            "      ,A[2        ,I]        ,\"%s\")",
            fx_name, for_name, dom_name);
  }

  /*	Column 1: exo leg */

  if (today_flag) {
    strcpy(evt[0][1], "\0");
  }

  if (first_ex_flag) {
    strcpy(evt[1][1], "\0");
  }

  /*	Coupons */
  aux_rows[3] = aux_rows[4] = aux_rows[5] = aux_rows[6] = aux_rows[7] = ncpn;
  for (i = 0; i < nexo; i++) {
    aux[today_flag + first_ex_flag + i][3] = exo_leg->cpn[i].pay_date;
    aux[today_flag + first_ex_flag + i][4] = exo_leg->cpn[i].alpha;
    aux[today_flag + first_ex_flag + i][5] = exo_leg->cpn[i].beta;
    aux[today_flag + first_ex_flag + i][6] =
        (exo_leg->cpn[i].floored ? exo_leg->cpn[i].floor : NO_FLOOR);
    aux[today_flag + first_ex_flag + i][7] =
        (exo_leg->cpn[i].cxxapped ? exo_leg->cpn[i].cxxap : NO_CAP);

    sprintf(evt[today_flag + first_ex_flag + i][1],
            "MAX(A[6        ,I]        ,MIN(A[7        ,I]        ,A[4        "
            ",I]+A[5        ,I]*C[0  "
            "      ,I]))*DF(D[I]        ,A[3        ,I]        ,\"%s\")",
            dom_name);
  }

  /*	Redemption */
  sprintf(tmpstr,
          "+MAX(A[14        ,0]        ,MIN(A[15        ,0]        ,A[12       "
          " ,0]+A[13        ,0]*C[0  "
          "      ,I]))*DF(D[I]        ,A[3        ,I]        ,\"%s\")",
          dom_name);
  strcat(evt[ncpn - 1][1], tmpstr);

  aux_rows[12] = aux_rows[13] = aux_rows[14] = aux_rows[15] = 1;

  aux[0][12] = exo_leg->not_ref.alpha;
  aux[0][13] = exo_leg->not_ref.beta;
  aux[0][14] = (exo_leg->not_ref.floored ? exo_leg->not_ref.floor : NO_FLOOR);
  aux[0][15] = (exo_leg->not_ref.cxxapped ? exo_leg->not_ref.cxxap : NO_CAP);

  /*	PV of extra cash-flows */
  tmpund = lookup_und(dom_name);
  if (!tmpund) {
    err = "Underlying not found";
    goto FREE_RETURN;
  }
  dom_yc = get_ycname_from_irund(tmpund);
  pvcf = 0.0;
  for (i = 0; i < pd_ncf; i++) {
    pvcf += swp_f_df(today, pd_cfdtes[i], dom_yc) * pd_cf[i];
  }
  if (strlen(evt[0][1]) > 0) {
    if (pvcf < 0) {
      sprintf(tmpstr, "-%.6f", fabs(pvcf));
    } else {
      sprintf(tmpstr, "+%.6f", fabs(pvcf));
    }
  } else {
    sprintf(tmpstr, "%.6f", pvcf);
  }
  strcat(evt[0][1], tmpstr);

  /*	Column 2: funding leg */

  fund_name = (fund_leg->dom_for == 0 ? dom_name : for_name);
  tmpund = lookup_und(fund_name);
  if (!tmpund) {
    err = "Underlying not found";
    goto FREE_RETURN;
  }
  fund_yc = get_ycname_from_irund(tmpund);

  /*	PV of Libor coupons and redemption */

  /*	PV of Libor coupons and final notional */
  pvcf = 0.0;

  if (nfund > 0) {
    pvcf = fund_leg->notional *
           swp_f_df(today, fund_leg->cpn[0].start_date, fund_yc);
  }

  /*	PV of extra cash-flows */
  for (i = 0; i < fund_ncf; i++) {
    pvcf += swp_f_df(today, fund_cfdtes[i], fund_yc) * fund_cf[i];
  }

  /*	Print the total */
  sprintf(evt[0][2], "(%.6f", pvcf);

  /*	PV of funding coupons */
  sprintf(tmpstr, "+PVRNG(0        ,1        ,D[I]        ,\"%s\"))",
          fund_name);
  strcat(evt[0][2], tmpstr);

  /*	Multiply by fx if needed */
  if (fund_leg->dom_for == 1) {
    sprintf(tmpstr, "*S(\"%s\")", fx_name);
    strcat(evt[0][2], tmpstr);
  }

  for (i = 1; i < ncpn; i++) {
    strcpy(evt[i][2], "\0");
  }

  if (fund_leg->num_cpn == 0) {
    aux_rows[0] = aux_rows[1] = 1;
    aux[0][0] = today;
    aux[0][1] = 0.0;
  } else {
    aux_rows[0] = aux_rows[1] = fund_leg->num_cpn;
    for (i = 0; i < fund_leg->num_cpn; i++) {
      aux[i][0] = fund_leg->cpn[i].pay_date;
      aux[i][1] = fund_leg->cpn[i].cxxpn;
    }
  }

  for (i = 0; i < ncpn; i++) {
    strcpy(evt[i][3], "\0");
    strcpy(evt[i][4], "\0");
    aux[i][9] = -999999;
  }

  if (ncall > 0) {
    aux_rows[8] = aux_rows[9] = aux_rows[10] = aux_rows[11] = ncpn;

    /*	Column 3: IV calculation at call date */

    for (i = 0; i < ncall; i++) {
      /* Search the call date to use (closest fx fixing date) */
      j = 0;
      while ((evt_dtes[j] <= call[i].ex_date) && (j < ncpn - 1)) {
        j++;
      }

      /* Correction of
      if (j > 0 && call[i].ex_date - evt_dtes[j-1] < evt_dtes[j] -
      call[i].ex_date)
      {
              j--;
      }
      Pb when
      today<Call1<Fx0<Start1=End0<Call2<Fx1
      assign Call1 to today instead of Fx0

      */
      if (j > 0 && fabs(call[i].ex_date - evt_dtes[j - 1]) <
                       fabs(evt_dtes[j] - call[i].ex_date)) {
        if (((evt_dtes[j - 1] == call[i].ex_date) &&
             (call[i].ex_date < evt_dtes[j] - 60)) ||
            (j > 1) || (first_ex_flag == 1) || (today_flag == 0)) {
          j--;
        }
      }

      /*	First exotic coupon to be called */
      k = call[i].pd_idx + today_flag + first_ex_flag;

      if (k < j)
      /*	This should never happen */
      {
        err = "Serious problem (1) in grfn tableau generation";
        goto FREE_RETURN;
      }

      if (k == j)
      /*	The current coupon is called */
      {
        sprintf(evt[j][3],
                "(PV[1]+C[1        ,I]-A[11        ,I]*DF(D[I]        ,A[10    "
                "    ,I]        ,\"%s\")",
                dom_name);
      }

      if (k == j + 1)
      /*	The next coupon is called */
      {
        sprintf(evt[j][3],
                "(PV[1]-A[11        ,I]*DF(D[I]        ,A[10        ,I]        "
                ",\"%s\")",
                dom_name);
      }

      if (k > j + 1)
      /*	This should never happen */
      {
        err = "Serious problem (2) in grfn tableau generation";
        goto FREE_RETURN;
      }

      /*	First funding coupon to be called */
      k = call[i].fund_idx;

      /*	PV of funding leg */
      sprintf(tmpstr, "-PVRNG(0        ,1        ,%d        ,\"%s\")",
              fund_leg->cpn[k].pay_date - 1, fund_name);
      strcat(evt[j][3], tmpstr);

      /*	Multiply by fx if needed */
      if (fund_leg->dom_for == 1) {
        sprintf(tmpstr, "*S(\"%s\")", fx_name);
        strcat(evt[j][3], tmpstr);
      }

      strcat(evt[j][3], ")");

      /*	Multiply by -1 if payer */
      if (call[i].pay_rec == 1) {
        strcat(evt[j][3], "*(-1)");
      }

      /*	Make call */
      if (cpd->call[i].cxxall_type == 0) {
        sprintf(evt[j][4], "MAX(C[3        ,I]        ,PV[J])");
      } else
      /*	Barrier */
      {
        if (call[i].bar_type == 0) {
          sprintf(evt[j][4],
                  "IF(S(\"%s\")>A[8        ,I]        ,C[3        ,I]        "
                  ",PV[J])",
                  fx_name);
        } else {
          sprintf(evt[j][4],
                  "IF(S(\"%s\")<A[8        ,I]        ,C[3        ,I]        "
                  ",PV[J])",
                  fx_name);
        }

        aux[j][8] = call[i].barrier;
        aux[j][9] = 4;
      }
      aux[j][10] =
          (call[i].set_date >= evt_dtes[j] ? call[i].set_date : evt_dtes[j]);
      aux[j][11] = call[i].pd_not_amt;
    }
  } else {
    aux_rows[8] = aux_rows[9] = aux_rows[10] = aux_rows[11] = 0;
  }

  for (i = 0; i < *num_evt_dtes; i++) {
    for (j = 0; j < *num_cols; j++) {
      mask[i][j] = set_mask(evt[i][j]);
    }
  }

FREE_RETURN:

  return err;
}

/*	Function that produces a GRFN tableau for a CPD deal */
Err cpd_make_grfn_tableau_credit(
    /*	The structure */
    CPD_STR cpd,
    /*	Extra cash-flows coming from past coupons and initial exchange */
    int fund_ncf, long *fund_cfdtes, double *fund_cf, int pd_ncf,
    long *pd_cfdtes, double *pd_cf,
    /*	The underlyings */
    long today, char *dom_name, char *for_name, char *fx_name,

    /*	Credit inputs */
    char *risky_yc, long credit_freq,
    /*	The tableau and auxiliaries
            Must be allocated with maximum values
                    and initialised with zeros */
    int *num_evt_dtes, int *num_cols, int *aux_rows, int *aux_cols,
    long *evt_dtes, char ***evt, int **mask, double **aux) {
  PD_EXO_LEG exo_leg = cpd->pd_leg;
  PD_FUND_LEG fund_leg = cpd->fund_leg;
  PD_CALL call = cpd->call;
  int i, j, k;
  int ncpn, nexo = exo_leg->num_cpn, nfund = fund_leg->num_cpn,
            ncall = cpd->num_calls, ncols = 5;
  int today_flag, first_ex_flag;
  char *fund_name, tmpstr[MAX_STR_LEN];
  char *dom_yc, *fund_yc;
  double pvcf;

  long step;

  SrtUndPtr tmpund;
  Err err = NULL;

  /*	Initialise and make event dates = today + exo coupons fixing dates
   */

  if (exo_leg->cpn[0].fx_fix_date > today) {
    today_flag = 1;
  } else {
    today_flag = 0;
  }

  if (ncall > 0 && call[0].ex_date > today + 15 &&
      call[0].ex_date < exo_leg->cpn[0].fx_fix_date - 60) {
    first_ex_flag = 1;
  } else {
    first_ex_flag = 0;
  }

  *num_cols = ncols;
  *aux_cols = 16;

  /* cpd event dates */
  ncpn = nexo + today_flag + first_ex_flag;

  /* total event dates */
  *num_evt_dtes = nexo + (nexo - 1) * credit_freq + today_flag + first_ex_flag;

  evt_dtes[0] = today;
  if (first_ex_flag)
    evt_dtes[1 + credit_freq] = call[0].ex_date;
  for (i = 0; i < nexo; i++) {
    evt_dtes[today_flag + first_ex_flag + i * (1 + credit_freq)] =
        exo_leg->cpn[i].fx_fix_date;
  }

  /* fill the intermediary dates */
  for (i = 0; i < nexo - 1; i++) {
    step =
        (long)((evt_dtes[today_flag + first_ex_flag +
                         (i + 1) * (1 + credit_freq)] -
                evt_dtes[today_flag + first_ex_flag + i * (1 + credit_freq)]) /
                   (credit_freq + 1) +
               0.5);

    for (j = 1; j <= credit_freq; j++) {
      evt_dtes[today_flag + first_ex_flag + i * (1 + credit_freq) + j] =
          evt_dtes[today_flag + first_ex_flag + i * (1 + credit_freq) + j - 1] +
          step;
    }
  }

  /*	Column 0: fx fixings */

  aux_rows[2] = *num_evt_dtes;

  for (i = 0; i < nexo; i++) {
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][2] =
        add_unit(evt_dtes[today_flag + first_ex_flag + i * (1 + credit_freq)],
                 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    sprintf(evt[today_flag + first_ex_flag + i * (1 + credit_freq)][0],
            "S(\"%s\")*DF(D[I]        ,A[2        ,I]        ,\"%s\")/DF(D[I]  "
            "      ,A[2        ,I]        ,\"%s\")",
            fx_name, for_name, dom_name);
  }

  /* fill the intermediary points */
  for (i = 0; i < nexo - 1; i++) {
    for (j = 1; j <= credit_freq; j++) {
      sprintf(evt[today_flag + first_ex_flag + i * (1 + credit_freq) + j][0],
              "pv[j]");
    }
  }

  /*	Column 1: exo leg */

  if (today_flag) {
    strcpy(evt[0][1], "\0");
  }

  if (first_ex_flag) {
    strcpy(evt[1][1], "\0");
  }

  /*	Coupons */
  aux_rows[3] = aux_rows[4] = aux_rows[5] = aux_rows[6] = aux_rows[7] =
      *num_evt_dtes;

  for (i = 0; i < nexo; i++) {
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][3] =
        exo_leg->cpn[i].pay_date;
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][4] =
        exo_leg->cpn[i].alpha;
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][5] =
        exo_leg->cpn[i].beta;
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][6] =
        (exo_leg->cpn[i].floored ? exo_leg->cpn[i].floor : NO_FLOOR);
    aux[today_flag + first_ex_flag + i * (1 + credit_freq)][7] =
        (exo_leg->cpn[i].cxxapped ? exo_leg->cpn[i].cxxap : NO_CAP);

    sprintf(evt[today_flag + first_ex_flag + i * (1 + credit_freq)][1],
            "MAX(A[6        ,I]        ,MIN(A[7        ,I]        ,A[4        "
            ",I]+A[5        ,I]*C[0  "
            "      ,I]))*DF(D[I]        ,A[3        ,I]        ,\"%s\")",
            dom_name);
  }

  for (i = 0; i < nexo; i++) {
    for (j = 1; j <= credit_freq; j++) {
      aux[today_flag + first_ex_flag + i * (1 + credit_freq) + j][3] = 0.0;
      aux[today_flag + first_ex_flag + i * (1 + credit_freq) + j][4] = 0.0;
      aux[today_flag + first_ex_flag + i * (1 + credit_freq) + j][5] = 0.0;
      aux[today_flag + first_ex_flag + i * (1 + credit_freq) + j][6] = 0.0;
      aux[today_flag + first_ex_flag + i * (1 + credit_freq) + j][7] = 0.0;
      sprintf(evt[today_flag + first_ex_flag + i * (1 + credit_freq) + j][1],
              "pv[j]");
    }
  }

  /*	Redemption */
  sprintf(tmpstr,
          "+MAX(A[14        ,0]        ,MIN(A[15        ,0]        ,A[12       "
          " ,0]+A[13        ,0]*C[0  "
          "      ,I]))*DF(D[I]        ,A[3        ,I]        ,\"%s\")",
          dom_name);
  strcat(evt[*num_evt_dtes - 1][1], tmpstr);

  aux_rows[12] = aux_rows[13] = aux_rows[14] = aux_rows[15] = 1;

  aux[0][12] = exo_leg->not_ref.alpha;
  aux[0][13] = exo_leg->not_ref.beta;
  aux[0][14] = (exo_leg->not_ref.floored ? exo_leg->not_ref.floor : NO_FLOOR);
  aux[0][15] = (exo_leg->not_ref.cxxapped ? exo_leg->not_ref.cxxap : NO_CAP);

  /*	PV of extra cash-flows */
  tmpund = lookup_und(dom_name);
  if (!tmpund) {
    err = "Underlying not found";
    goto FREE_RETURN;
  }
  dom_yc = get_ycname_from_irund(tmpund);
  pvcf = 0.0;
  for (i = 0; i < pd_ncf; i++) {
    pvcf += swp_f_df(today, pd_cfdtes[i], dom_yc) * pd_cf[i];
  }
  if (strlen(evt[0][1]) > 0) {
    if (pvcf < 0) {
      sprintf(tmpstr, "-%.6f", fabs(pvcf));
    } else {
      sprintf(tmpstr, "+%.6f", fabs(pvcf));
    }
  } else {
    sprintf(tmpstr, "%.6f", pvcf);
  }
  strcat(evt[0][1], tmpstr);

  /*	Column 2: funding leg */

  fund_name = (fund_leg->dom_for == 0 ? dom_name : for_name);
  tmpund = lookup_und(fund_name);
  if (!tmpund) {
    err = "Underlying not found";
    goto FREE_RETURN;
  }
  fund_yc = get_ycname_from_irund(tmpund);

  /*	PV of Libor coupons and redemption */

  /*	PV of Libor coupons and final notional */
  pvcf = fund_leg->notional *
         swp_f_df(today, fund_leg->cpn[0].start_date, fund_yc);

  /*	PV of extra cash-flows */
  for (i = 0; i < fund_ncf; i++) {
    pvcf += swp_f_df(today, fund_cfdtes[i], fund_yc) * fund_cf[i];
  }

  /*	Print the total */
  sprintf(evt[0][2], "(%.6f", pvcf);

  /*	PV of funding coupons */
  sprintf(tmpstr, "+PVRNG(0        ,1        ,D[I]        ,\"%s\"))",
          fund_name);
  strcat(evt[0][2], tmpstr);

  /*	Multiply by fx if needed */
  if (fund_leg->dom_for == 1) {
    sprintf(tmpstr, "*S(\"%s\")", fx_name);
    strcat(evt[0][2], tmpstr);
  }

  for (i = 1; i < ncpn; i++) {
    strcpy(evt[i][2], "\0");
  }

  aux_rows[0] = aux_rows[1] = fund_leg->num_cpn;
  for (i = 0; i < fund_leg->num_cpn; i++) {
    aux[i][0] = fund_leg->cpn[i].pay_date;
    aux[i][1] = fund_leg->cpn[i].cxxpn;
  }

  for (i = 0; i < ncpn; i++) {
    strcpy(evt[i][3], "\0");
    strcpy(evt[i][4], "\0");
    aux[i][9] = -999999;
  }

  if (ncall > 0) {
    aux_rows[8] = aux_rows[9] = aux_rows[10] = aux_rows[11] = *num_evt_dtes;

    /*	Column 3: IV calculation at call date */

    for (i = 0; i < ncall; i++) {
      /* Search the call date to use (closest fx fixing date) */
      j = 0;
      while (evt_dtes[j] < call[i].ex_date) {
        j++;
      }

      if (j > 0 &&
          call[i].ex_date - evt_dtes[j - 1] < evt_dtes[j] - call[i].ex_date) {
        j--;
      }

      /*	First exotic coupon to be called */
      k = call[i].pd_idx * (credit_freq + 1) + today_flag + first_ex_flag;

      if (k < j)
      /*	This should never happen */
      {
        err = "Serious problem (1) in grfn tableau generation";
        goto FREE_RETURN;
      }

      if (k == j)
      /*	The current coupon is called */
      {
        sprintf(evt[j][3],
                "(PV[1]+C[1        ,I]-A[11        ,I]*DF(D[I]        ,A[10    "
                "    ,I]        ,\"%s\")",
                dom_name);
      }

      if (k == j + 1 + credit_freq)
      /*	The next coupon is called */
      {
        sprintf(evt[j][3],
                "(PV[1]-A[11        ,I]*DF(D[I]        ,A[10        ,I]        "
                ",\"%s\")",
                dom_name);
      }

      if (k > j + 1 + credit_freq)
      /*	This should never happen */
      {
        err = "Serious problem (2) in grfn tableau generation";
        goto FREE_RETURN;
      }

      /*	First funding coupon to be called */
      k = call[i].fund_idx;

      /*	PV of funding leg */
      sprintf(tmpstr, "-PVRNG(0        ,1        ,%d        ,\"%s\")",
              fund_leg->cpn[k].pay_date - 1, fund_name);
      strcat(evt[j][3], tmpstr);

      /*	Multiply by fx if needed */
      if (fund_leg->dom_for == 1) {
        sprintf(tmpstr, "*S(\"%s\")", fx_name);
        strcat(evt[j][3], tmpstr);
      }

      strcat(evt[j][3], ")");

      /*	Multiply by -1 if payer */
      if (call[i].pay_rec == 1) {
        strcat(evt[j][3], "*(-1)");
      }

      /*	Make call */
      if (cpd->call[i].cxxall_type == 0) {
        sprintf(evt[j][4], "MAX(C[3        ,I]        ,PV[J])");
      } else
      /*	Barrier */
      {
        if (call[i].bar_type == 0) {
          sprintf(evt[j][4],
                  "IF(S(\"%s\")>A[8        ,I]        ,C[3        ,I]        "
                  ",PV[J])",
                  fx_name);
        } else {
          sprintf(evt[j][4],
                  "IF(S(\"%s\")<A[8        ,I]        ,C[3        ,I]        "
                  ",PV[J])",
                  fx_name);
        }

        aux[j][8] = call[i].barrier;
        aux[j][9] = 4;
      }
      aux[j][10] =
          (call[i].set_date >= evt_dtes[j] ? call[i].set_date : evt_dtes[j]);
      aux[j][11] = call[i].pd_not_amt;
    }
  } else {
    aux_rows[8] = aux_rows[9] = aux_rows[10] = aux_rows[11] = 0;
  }

  for (i = 0; i < *num_evt_dtes; i++) {
    for (j = 0; j < *num_cols; j++) {
      mask[i][j] = set_mask(evt[i][j]);
    }
  }

FREE_RETURN:

  return err;
}

/*	Fill underlying structure from a predefined underlying */
Err cpd_fill_smile_und(char *fx3dund, double alpha, double beta,
                       CPD_SMILE_UND und, double dom_vol_shift,
                       double for_vol_shift, double fx_vol_shift) {
  Err err;

  err = cpd_fill_und(fx3dund, &(und->und), NULL, dom_vol_shift, for_vol_shift,
                     fx_vol_shift);
  if (err) {
    return err;
  }

  und->alpha = alpha;
  und->beta = beta;

  return NULL;
}

/*	Fill underlying structure from calibration instruments */
Err cpd_calib_smile_und(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I        , 1: E */
    double fx_spot, long fx_spot_date,
    int dom_calib,        /*	Calibrate domestic underlying */
    char *dom_und,        /*	If no        , domestic underlying to be used */
    char *dom_yc,         /*	Domestic yc */
    char *dom_vc,         /*	Domestic vc (only if calib) */
    char *dom_ref,        /*	Domestic ref rate (only if calib) */
    char *dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char *dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,       /*	Domestic lambda */
    int for_calib,        /*	Same for foreign */
    char *for_und, char *for_yc, char *for_vc, char *for_ref,
    char *for_swap_freq, char *for_swap_basis, double for_lam,
    double *corr_times, double *correl_dom_for, /*	Correlations */
    double *correl_dom_fx, double *correl_for_fx, long corr_n_times,
    double alpha, double beta, CPD_STR cpd, /*	Structure */
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char *vol_curve_name, double start_date,
                           double end_date, double cash_strike, int zero,
                           char *ref_rate_name, double *vol, double *power),
    /*	Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,
    CPD_SMILE_UND und, double dom_vol_shift, double for_vol_shift,
    double fx_vol_shift) {
  int i, nex;
  long *lex = NULL, last;

  double *fx_mkt_vol_time = NULL;

  double *temp_sig_time1 = NULL, *temp_sig_time2 = NULL, *temp_sig1 = NULL,
         *temp_sig2 = NULL;

  double *temp_tau_time = NULL, *temp_tau = NULL;

  int temp_tau_n, temp_sig_n1, temp_sig_n2;
  CPD_UND basic_und;
  int num_cal;
  Err err = NULL;

  /*	Eliminate zero/negative lambdas */

  if (dom_lam < 1.0e-08) {
    dom_lam = 1.0e-08;
  }

  if (for_lam < 1.0e-08) {
    for_lam = 1.0e-08;
  }

  /*	Initialise */

  basic_und = &(und->und);

  basic_und->sigma_date_rates = NULL;
  basic_und->sigma_time_rates = NULL;
  basic_und->sigma_dom = NULL;
  basic_und->sigma_for = NULL;
  basic_und->sigma_date_fx = NULL;
  basic_und->sigma_time_fx = NULL;
  basic_und->sigma_fx = NULL;
  basic_und->model = NULL;
  basic_und->fees = NULL;
  basic_und->fees_dates = NULL;
  basic_und->nb_fees = 0;

  basic_und->today = today;

  basic_und->spot_fx = fx_spot * swp_f_df(today, fx_spot_date, dom_yc) /
                       swp_f_df(today, fx_spot_date, for_yc);

  basic_und->corr_times = corr_times;
  basic_und->correl_dom_for = correl_dom_for;
  basic_und->correl_dom_fx = correl_dom_fx;
  basic_und->correl_for_fx = correl_for_fx;
  basic_und->corr_n_times = corr_n_times;

  und->fx2bdfwd = fx_spot;
  und->alpha = alpha;
  und->beta = beta;

  strcpy(basic_und->dom_yc, dom_yc);
  strcpy(basic_und->for_yc, for_yc);
  strcpy(basic_und->name, "CALIB");

  /*	Find exercise dates for the calibration of interest rates */

  last = cpd->pd_leg->cpn[cpd->pd_leg->num_cpn - 1].pay_date;
  if (dom_calib || for_calib) {
    if (cpd->num_calls > 0 &&
        !(cpd->num_calls == 1 && cpd->call[0].ex_date <= today + eod_flag)) {
      /*	If call dates        , choose call dates as option expiries for
       * calibration */
      nex = cpd->num_calls;
      lex = (long *)calloc(nex, sizeof(long));
      if (!lex) {
        err = "Allocation error (1) in cpd_calib_smile_und";
        goto FREE_RETURN;
      }

      for (i = 0; i < nex; i++) {
        lex[i] = cpd->call[i].ex_date;
      }
    } else {
      /*	If no call dates        , calibrate to last/2 into last/2 */
      nex = 1;
      lex = (long *)calloc(nex, sizeof(long));
      if (!lex) {
        err = "Allocation error (1') in cpd_calib_smile_und";
        goto FREE_RETURN;
      }

      lex[0] = (last + today) / 2;
    }
  }

  /*	Domestic underlying */

  if (dom_calib) {
    basic_und->lda_dom = dom_lam;

    err = cpd_calib_diagonal(
        dom_yc, dom_vc, dom_ref, get_ir_cash_vol, dom_vol_shift, 0, nex, lex,
        last, NULL, NULL, 0, 1.0, 1.0, dom_swap_freq, dom_swap_basis, 1, 0, 1,
        CALPRES, CALPRES, 0, 0, 0, 0, NULL, &dom_lam, 1, 0.0, 0.0, 0.0,
        &temp_sig_n1, &temp_sig_time1, &temp_sig1, NULL);

    if (err) {
      goto FREE_RETURN;
    }
  } else {
    err = Get_LGM_TermStructure(dom_und, &temp_sig_time1, &temp_sig1,
                                &temp_sig_n1, &temp_tau_time, &temp_tau,
                                &temp_tau_n);
    if (err) {
      goto FREE_RETURN;
    }

    err = get_unique_lambda(temp_tau, temp_tau_n, &(basic_und->lda_dom));
    if (err) {
      goto FREE_RETURN;
    }
    free(temp_tau_time);
    temp_tau_time = NULL;
    free(temp_tau);
    temp_tau = NULL;

    for (i = 0; i < temp_sig_n1; i++) {
      temp_sig1[i] += dom_vol_shift;
    }
  }

  /*	Foreign underlying */

  if (for_calib) {
    basic_und->lda_for = for_lam;

    err = cpd_calib_diagonal(
        for_yc, for_vc, for_ref, get_ir_cash_vol, for_vol_shift, 0, nex, lex,
        last, NULL, NULL, 0, 1.0, 1.0, for_swap_freq, for_swap_basis, 1, 0, 1,
        CALPRES, CALPRES, 0, 0, 0, 0, NULL, &for_lam, 1, 0.0, 0.0, 0.0,
        &temp_sig_n2, &temp_sig_time2, &temp_sig2, NULL);

    if (err) {
      goto FREE_RETURN;
    }
  } else {
    err = Get_LGM_TermStructure(for_und, &temp_sig_time2, &temp_sig2,
                                &temp_sig_n2, &temp_tau_time, &temp_tau,
                                &temp_tau_n);
    if (err) {
      goto FREE_RETURN;
    }

    err = get_unique_lambda(temp_tau, temp_tau_n, &(basic_und->lda_for));
    if (err) {
      goto FREE_RETURN;
    }
    free(temp_tau_time);
    temp_tau_time = NULL;
    free(temp_tau);
    temp_tau = NULL;

    for (i = 0; i < temp_sig_n2; i++) {
      temp_sig2[i] += for_vol_shift;
    }
  }

  /*	Merge rates structures */

  err = merge_rates_ts(temp_sig_time1, temp_sig1, temp_sig_n1, temp_sig_time2,
                       temp_sig2, temp_sig_n2, &(basic_und->sigma_time_rates),
                       &(basic_und->sigma_dom), &(basic_und->sigma_for),
                       &(basic_und->sigma_n_rates));

  if (err) {
    goto FREE_RETURN;
  }

  basic_und->sigma_date_rates =
      (double *)calloc(basic_und->sigma_n_rates, sizeof(double));
  if (!basic_und->sigma_date_rates) {
    err = "Allocation error (2) in cpd_calib_smile_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < basic_und->sigma_n_rates; i++) {
    basic_und->sigma_date_rates[i] =
        today + basic_und->sigma_time_rates[i] * DAYS_IN_YEAR + 1.0e-16;
  }

  /*	Calibrate Fx */

  /*	Remove unused fx vol dates */
  i = num_fx_mkt_vol - 1;
  while (i > 0 && fx_mkt_vol_date[i] > last) {
    i--;
  }
  if (i < num_fx_mkt_vol - 1 && fx_mkt_vol_date[i] < last) {
    i++;
  }
  num_fx_mkt_vol = i + 1;

  basic_und->sigma_n_fx = num_fx_mkt_vol;
  fx_mkt_vol_time = (double *)calloc(num_fx_mkt_vol, sizeof(double));
  basic_und->sigma_date_fx = (double *)calloc(num_fx_mkt_vol, sizeof(double));
  basic_und->sigma_time_fx = (double *)calloc(num_fx_mkt_vol, sizeof(double));
  if (!fx_mkt_vol_time || !basic_und->sigma_date_fx ||
      !basic_und->sigma_time_fx) {
    err = "Allocation error (2) in cpd_calib_smile_und";
    goto FREE_RETURN;
  }

  num_cal = 0;
  for (i = 0; i < num_fx_mkt_vol; i++) {
    fx_mkt_vol_time[i] = (fx_mkt_vol_date[i] - today) * YEARS_IN_DAY;
    fx_mkt_vol[i] += fx_vol_shift;
    basic_und->sigma_date_fx[i] =
        add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    basic_und->sigma_time_fx[i] =
        (basic_und->sigma_date_fx[i] - today) * YEARS_IN_DAY;
    if (basic_und->sigma_time_fx[i] >= 10) {
      num_cal++;
    }
  }

  err = Fx3DAlphaBetatsCalibration2_corr(
      basic_und->today, basic_und->sigma_time_fx, fx_mkt_vol_time, fx_mkt_vol,
      basic_und->sigma_n_fx, num_cal, basic_und->sigma_time_rates,
      basic_und->sigma_n_rates, basic_und->sigma_dom, basic_und->lda_dom,
      basic_und->sigma_for, basic_und->lda_for, und->alpha, und->beta,
      und->fx2bdfwd, basic_und->corr_times, basic_und->correl_dom_for,
      basic_und->correl_dom_fx, basic_und->correl_for_fx,
      basic_und->corr_n_times, basic_und->dom_yc, basic_und->for_yc,
      &(basic_und->sigma_fx), 50, 2, 0.1, 0.1, 200);

  if (err) {
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (err)
    cpd_free_und(basic_und);

  if (fx_mkt_vol_time)
    free(fx_mkt_vol_time);

  if (lex)
    free(lex);
  if (temp_sig_time1)
    free(temp_sig_time1);
  if (temp_sig1)
    free(temp_sig1);
  if (temp_sig_time2)
    free(temp_sig_time2);
  if (temp_sig2)
    free(temp_sig2);
  if (temp_tau_time)
    free(temp_tau_time);
  if (temp_tau)
    free(temp_tau);

  return err;
}

void cpd_copy_smile_und(CPD_SMILE_UND src, CPD_SMILE_UND dest) {
  cpd_copy_und(&(src->und), &(dest->und));
  src->alpha = dest->alpha;
  src->beta = dest->beta;
  src->fx2bdfwd = dest->fx2bdfwd;
}

Err cpd_free_smile_und(CPD_SMILE_UND und) { return cpd_free_und(&und->und); }

#if 0

/*	Launch the smile tree */
Err cpd_launch_smile_tree(
CPD_STR		cpd        ,
CPD_SMILE_UND	und        ,
/*	Number of steps */
int			num_stp        ,
/*	GRFN tableau structures */
int			*num_evt_dtes        ,
int			*num_cols        ,
int			*aux_rows        ,
int			*aux_cols        ,
long		*evt_dtes        ,
char		***evt        ,
int			**mask        ,
double		**aux        ,
/*	Results */
double		*fund_leg        ,
double		*pd_leg        ,
double		*call)
{
	Err		err			= NULL;
	int		i        , nsig;
	double	**sig_tab	= NULL        ,
			*corr_date	= NULL        ,
			**correl	= NULL;
	char	***und_tab	= NULL;
	double	tau__        , *tau_ = &tau__        , **tau = &tau_;
	SrtUndPtr	dom_und; 
	SrtUndPtr	for_und;
	char	fx_name1[256]        , fx_name2[256]        , *dom_ccy        , *for_ccy;
	int		num_prod;
	double	*prod_val	= NULL;
	double	*barrier	= NULL;
	int		*bar_col	= NULL;

	/*	1.- Create IR underlyings */
	
	nsig = max (und->und.sigma_n_rates        , und->und.sigma_n_fx);
	sig_tab = dmatrix (0        , nsig-1        , 0        , 1);
	if (!sig_tab)
	{
		err = "Allocation error in cpd_launch_smile_tree";
		goto FREE_RETURN;
	}

	for (i=0; i<und->und.sigma_n_rates; i++)
	{
		sig_tab[i][0] = und->und.sigma_date_rates[i];
		sig_tab[i][1] = und->und.sigma_dom[i];
	}

	tau__ = 1.0 / und->und.lda_dom;
	err = SrtInitIRUnd(
		"DOM"        ,
		und->und.dom_yc        ,
		"LGM"        ,
		und->und.sigma_n_rates        ,
		2        ,
		sig_tab        ,
		1        ,
		1        ,
		tau        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0        ,
		0        ,
		NULL);

	if (err)
	{
		goto FREE_RETURN;
	}

	dom_und = lookup_und ("DOM");
	dom_ccy = get_underlying_ccy (dom_und);

	for (i=0; i<und->und.sigma_n_rates; i++)
	{
		sig_tab[i][1] = und->und.sigma_for[i];
	}

	tau__ = 1.0 / und->und.lda_for;
	err = SrtInitIRUnd(
		"FOR"        ,
		und->und.for_yc        ,
		"LGM"        ,
		und->und.sigma_n_rates        ,
		2        ,
		sig_tab        ,
		1        ,
		1        ,
		tau        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0.0        ,
		0        ,
		0        ,
		NULL);

	if (err)
	{
		goto FREE_RETURN;
	}

	for_und = lookup_und ("FOR");
	for_ccy = get_underlying_ccy (for_und);

	sprintf (fx_name1        ,"%s/%s"        , for_ccy        , dom_ccy);	

	/*	2.- Create corr matrix */

	und_tab = smatrix_size (0        , 2        , 0        , 1        , 256);
	if (!und_tab)
	{
		err = "Allocation error in cpd_launch_smile_tree";
		goto FREE_RETURN;
	}
	
	strcpy (und_tab[0][0]        , "DOM");
	strcpy (und_tab[0][1]        , "FOR");
	strcpy (und_tab[1][0]        , "FOM");
	strcpy (und_tab[1][1]        , fx_name1);
	strcpy (und_tab[2][0]        , "FOR");
	strcpy (und_tab[2][1]        , fx_name1);

	correl = dmatrix (0        , 2        , 0        , und->und.cxxorr_n_times -1);
	corr_date = dvector(0        , und->und.cxxorr_n_times -1);

	if (!correl || !corr_date)
	{
		err = "Allocation error in cpd_launch_smile_tree";
		goto FREE_RETURN;
	}
	
	for (i=0; i<und->und.cxxorr_n_times; i++)
	{
		correl[0][i] = und->und.cxxorrel_dom_for[0];
		correl[1][i] = und->und.cxxorrel_dom_fx[0];
		correl[2][i] = und->und.cxxorrel_for_fx[0];
		corr_date[i] = und->und.today + DAYS_IN_YEAR * und->und.cxxorr_times[i];
	}	

	err = SrtInitCorrelationMatrix(
		und.cxxorr_n_times        ,
		3        ,
		correl        ,
		corr_date        ,
		und_tab);

	if (err)
	{
		goto FREE_RETURN;
	}

	/*	3.- Create Fx underlying */

	for (i=0; i<und->und.sigma_n_fx; i++)
	{
		sig_tab[i][0] = und->und.sigma_date_fx[i];
		sig_tab[i][1] = und->und.sigma_fx[i];
	}
	
	err = SrtInitFXUnd(
		fx_name1        ,
		und->fx2bdfwd        ,
		"STOCH_RATES"        ,
		"DOM"        ,
		"FOR"        ,
		NULL        ,
		und->und.sigma_n_fx        ,
		2        ,
		sig_tab        ,
		fx_name2);

	if (err)
	{
		goto FREE_RETURN;
	}

	/*	4.- Create GRFN tableau */
	
	err = cpd_make_grfn_tableau(
		cpd        ,
		und->und.today        ,
		"DOM"        ,
		"FOR"        ,
		fx_name2        ,
		num_evt_dtes        ,
		num_cols        ,
		aux_rows        ,
		aux_cols        ,
		evt_dtes        ,
		evt        ,
		mask        ,
		aux);

	if (err)
	{
		goto FREE_RETURN;
	}

	barrier = (double*) calloc (*num_evt_dtes        , sizeof (double));
	bar_col = (int*) calloc (*num_evt_dtes        , sizeof (int));

	if (!barrier || !bar_col)
	{
		err = "Allocation error in cpd_launch_smile_tree";
		goto FREE_RETURN;
	}

	for (i=0; i<*num_evt_dtes; i++)
	{
		barrier[i] = aux[i][8];
		bar_col[i] = (int) (aux[i][9] + 1.0e-08);
	}

	/*	5.- Launch the calculation */

	err = SrtGrfn3DFXAlphaBetaTree(
		fx_name2        ,
		und->alpha        ,
		und->beta        ,
		*num_evt_dtes        ,
		evt_dtes        ,
		*num_evt_dtes        ,
		*num_cols        ,
		evt        ,
		mask        ,
		*aux_cols        ,
		aux_rows        ,
		aux        ,
		barrier        ,
		bar_col        ,
		num_stp        ,
		&num_prod        ,
		1        ,
		&prod_val);

	if (err)
	{
		goto FREE_RETURN;
	}

	/*	6.-	Free all & return */
FREE_RETURN:
	
	if (sig_tab)
	{
		free_dmatrix (sig_tab        , 0        , nsig-1        , 0        , 1);
	}

	if (correl)
	{
		free_dmatrix (correl        , 0        , 2        , und->und.cxxorr_n_times -1);
	}

	if (corr_date)
	{
		free_dvector (corr_date        , 0        , und->und.cxxorr_n_times -1);
	}

	if (und_tab)
	{
		free_smatrix_size (und_tab        , 0        , 2        , 0        , 1        , 256);
	}

	if (prod_val)
	{
		free (prod_val);
	}

	if (barrier)
	{
		free (barrier);
	}

	if (bar_col)
	{
		free (bar_col);
	}

	return err;
}

#endif

/* Function  call_GRFN_alphabetaModel

        Price the funding        , power dual leg and the call of a PD in the
   alpha beta model by initialising the dom        , for        , and fx
   underlying if necessary making the grfn tableau from the cpd_str structure of
   the power dual calling the alphabeta grfn tree to price

        Notice that
                the cpd_und underlying and cpd_str structure are filled by the
   cpd_fill_check_all_struct function

                In the cpd_und        , the dom and for instantaneous ts are
   available. the fx ts corresponds either to the fx und passed to the function
   (which in that case is assumed to be a alpha beta calibrated fx ts) or the
   calibrated one by the cpd_fill_check_all_struct that is to say a Lognormal
   Model calibrated fx ts

                In the case where the Fx underlying needs to be calibrated with
   the alpha beta model        , the call_GRFN_alphabetaModel calculates the ts
   and filled the cpd_und with it.

  */
Err call_GRFN_alphabetaModel(
    /* Input */
    CPD_UND cpd_und, /* The CPD underlying */
    CPD_STR cpd_str,

    int use_calib,                /*	0: use fx3dund        , 1: calibrate */
    char *fx_name, int dom_calib, /*	Calibrate domestic underlying */
    char *dom_name, /*	If no        , domestic underlying to be used */
    char *dom_yc,   /*	Domestic yc */
    int for_calib,  /*	Same for foreign */
    char *for_name, char *for_yc,

    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */

    /*		funding */
    double fund_not, /*	Notional */
    int for_fund,    /*	1: if funding has been transform into domestic */
    double
        eq_final_ex, /*	equivalent final exchange use only if for_fund = 1 */
    double
        eq_init_ex, /*	equivalent init exchange use only if for_fund = 1 */
    long fund_start_date, int fund_ncpn, /*	Number of coupons */
    long *fund_fix,                      /*	Fixing dates */
    long *fund_start,                    /*	Start dates */
    long *fund_pay,                      /*	Pay dates */
    char **fund_basis,                   /*	Basis */
    double *fund_mrg,                    /*	Margins */
    double *fund_fix_cpn, /*	Past coupon fixing if relevant        ,
                                                    includes spr        , but
                       not mrg        , cvg and notional */
    /*		pd */
    double pd_not,   /*	Notional */
    int pd_ncpn,     /*	Number of coupons */
    long *pd_fix,    /*	Fx fixing dates */
    long *pd_start,  /*	Start dates */
    long *pd_pay,    /*	Pay dates */
    char **pd_basis, /*	Basis */
    double
        *pd_alpha, /*	Coupon = alpha + beta * fx [capped        , floored] */
    double *pd_beta, int *pd_floored, double *pd_floor, int *pd_capped,
    double *pd_cap, double *pd_fix_fx, /*	Past Fx fixing if relevant */

    /* Call */
    int ncall,       /*	Number of calls */
    long *ex_date,   /*	Call dates */
    double *barrier, /*	in case of a pure KO or a Callable KO */

    /*	EOD Fixing Flag */
    int eod_fix_flag, /*	0: I        , 1: E */
    /*	EOD Payment Flag */
    int eod_pay_flag, /*	0: I        , 1: E */

    /* Model Parameters */
    double alpha, double beta,

    /* Numerical parameter */
    int nb_step, /* Number of step for the tree */

    /* Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol,

    int num_fx_mkt_vol, double fx_spot, /* Spot Fx */

    /* OutPut */
    double *fund_val, /*	Value of the funding leg */
    double *pd_val,   /*	Value of the Power Dual leg */
    double *call_val, /*	Value of the callable feature */
    double
        *call_stdev /*	Standard deviation of the call if applicable */) {

  /* Declaration of local variables */

  Err err = NULL;

  SrtUndPtr fx_und;
  SrtBasisCode bas;

  int iNum;
  int num_cal;

  double **DateSigmaTab = NULL;
  double dom_tau, for_tau;
  double *ptr_dom_tau = &dom_tau;
  double *ptr_for_tau = &for_tau;

  char ***und_names = NULL;
  double *correl_dates = NULL;
  double **correl = NULL;

  int num_evt_dtes, tableauRows, tableauCols, auxWidth;
  int *auxLen = NULL;
  long *evt_dtes = NULL;
  char ***tableauStrings = NULL;
  double **aux = NULL;
  double **aux_transpose = NULL;
  int **mask = NULL;

  int num_prod;
  double *prod_val = NULL;
  double *tab_barrier = NULL;
  int *bar_col = NULL;

  double *maturity_opt = NULL;

  int fund_ncf;
  long *fund_cfdtes;
  double *fund_cf;

  int pd_ncf;
  long *pd_cfdtes;
  double *pd_cf;

  double PdFixingRate;

  int i, j;

  /* Memory allocation */

  DateSigmaTab =
      dmatrix(0, 1, 0, max(cpd_und->sigma_n_fx, cpd_und->sigma_n_rates));
  und_names = smatrix_size(0, 2, 0, 1, 256);
  correl = dmatrix(0, 2, 0, cpd_und->corr_n_times - 1);
  correl_dates = dvector(0, cpd_und->corr_n_times - 1);

  auxLen = ivector(0, MAX_COL - 1);
  evt_dtes = lngvector(0, MAX_NEVT - 1);
  tableauStrings =
      smatrix_size(0, MAX_NEVT - 1, 0, MAX_COL - 1, GRFN_DEF_ARGBUFSZ);

  aux = dmatrix(0, MAX_NEVT - 1, 0, MAX_COL - 1);
  aux_transpose = dmatrix(0, MAX_COL - 1, 0, MAX_NEVT - 1);
  mask = imatrix(0, MAX_NEVT - 1, 0, MAX_COL - 1);

  tab_barrier = dvector(0, MAX_NEVT - 1);
  bar_col = ivector(0, MAX_NEVT - 1);

  maturity_opt = dvector(0, num_fx_mkt_vol - 1);

  fund_cfdtes = lngvector(0, fund_ncpn - 1);
  fund_cf = dvector(0, fund_ncpn - 1);

  pd_cfdtes = lngvector(0, pd_ncpn - 1);
  ;
  pd_cf = dvector(0, pd_ncpn - 1);

  if ((!DateSigmaTab) || (!und_names) || (!correl) || (!correl_dates) ||
      (!auxLen) || (!evt_dtes) || (!tableauStrings) || (!aux) ||
      (!aux_transpose) || (!mask) || (!tab_barrier) || (!bar_col) ||
      (!maturity_opt) || (!fund_cfdtes) || (!fund_cf) || (!pd_cfdtes) ||
      (!pd_cf)) {
    err = "Allocation error in call_GRFN_alphabetaModel ";
    goto FREE_RETURN;
  }

  /* Modify the barrier in the cpd structure
   */
  /* It was transformed into log(barrier/cashFx) by cpd_fill_check_all_struct
   */
  for (iNum = 0; iNum < cpd_str->num_calls; iNum++) {
    if (cpd_str->call[iNum].cxxall_type > 0) {
      /* KO Case */

      /* Find the corresponding barrier and put the barrier back */
      for (i = 0; i < ncall; i++) {
        if (cpd_str->call[iNum].ex_date == ex_date[i])
          cpd_str->call[iNum].barrier = barrier[i];
      }
    }
  }

  /* Underlyings extraction or underlyings initialisation */

  if (!use_calib) {
    /* the fx3dund is passed */
    /* It is supposed to be a alpha beta fxund with the same alpha and beta
     * passed to the function */
    fx_und = lookup_und(fx_name);
    dom_name = get_domname_from_fxund(fx_und);
    for_name = get_forname_from_fxund(fx_und);
  } else {
    /* the fx3dund is not passed */

    if (dom_calib) {
      /* the Domund is not passed */

      /* Fill the DateSigmaTab from the Cpd und */
      for (iNum = 0; iNum < cpd_und->sigma_n_rates; iNum++) {
        DateSigmaTab[0][iNum] = cpd_und->sigma_date_rates[iNum];
        DateSigmaTab[1][iNum] = cpd_und->sigma_dom[iNum];
      }

      /* Init of the Dom Und in memory   */
      /* Name of the underlying is "DOMXXYYZXY" */
      strcpy(dom_name, "DOMXXYYZXY");

      dom_tau = 1.0 / cpd_und->lda_dom;
      err = SrtInitIRUnd(dom_name, dom_yc, "LGM", cpd_und->sigma_n_rates, 2,
                         DateSigmaTab, 1, 1, &ptr_dom_tau, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0, 0, NULL);

      if (err) {
        goto FREE_RETURN;
      }
    }

    if (for_calib) {
      /* the Forund is not passed */

      /* Fill the DateSigmaTab from the Cpd und */
      for (iNum = 0; iNum < cpd_und->sigma_n_rates; iNum++) {
        DateSigmaTab[0][iNum] = cpd_und->sigma_date_rates[iNum];
        DateSigmaTab[1][iNum] = cpd_und->sigma_for[iNum];
      }

      /* Init of the For Und in memory   */
      /* Name of the underlying is "FORXXYYZXY" */
      strcpy(for_name, "FORXXYYZXY");

      for_tau = 1.0 / cpd_und->lda_for;
      err = SrtInitIRUnd(for_name, for_yc, "LGM", cpd_und->sigma_n_rates, 2,
                         DateSigmaTab, 1, 1, &ptr_for_tau, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0, 0, NULL);

      if (err) {
        goto FREE_RETURN;
      }
    }

    /* Name of the Fx Underlying */
    sprintf(fx_name, "%s/%s", for_name, dom_name);

    /* Initialisation of the CorrMatrix */

    /* Fill the underlying names */
    strcpy(und_names[0][0], dom_name);
    strcpy(und_names[0][1], for_name);
    strcpy(und_names[1][0], dom_name);
    strcpy(und_names[1][1], fx_name);
    strcpy(und_names[2][0], for_name);
    strcpy(und_names[2][1], fx_name);

    /* Fill the corresponding correl tableau and the correl dates tableau */
    for (iNum = 0; iNum < cpd_und->corr_n_times; iNum++) {
      correl[0][iNum] = cpd_und->correl_dom_for[iNum];
      correl[1][iNum] = cpd_und->correl_dom_fx[iNum];
      correl[2][iNum] = cpd_und->correl_for_fx[iNum];
      correl_dates[iNum] =
          cpd_und->today + DAYS_IN_YEAR * cpd_und->corr_times[iNum];
    }

    err = SrtInitCorrelationMatrix(cpd_und->corr_n_times, 3, correl,
                                   correl_dates, und_names);

    if (err) {
      goto FREE_RETURN;
    }

    /* Creation of the Fx Underlying */

    /* Calibration of the Fx underlying */
    num_cal = 0; /* Number of opt of maturity > 10 Y */
    for (iNum = 0; iNum < cpd_und->sigma_n_fx; iNum++) {
      maturity_opt[iNum] =
          (fx_mkt_vol_date[iNum] - cpd_und->today) * YEARS_IN_DAY;
      if (maturity_opt[iNum] >= 10) {
        num_cal++;
      }
    }

    /* fill the cpd_und ->sigma_fx by the alpha beta calibrated one */
    err = Fx3DAlphaBetatsCalibration2_corr(
        cpd_und->today, cpd_und->sigma_time_fx, maturity_opt, fx_mkt_vol,
        cpd_und->sigma_n_fx, num_cal, cpd_und->sigma_time_rates,
        cpd_und->sigma_n_rates, cpd_und->sigma_dom, cpd_und->lda_dom,
        cpd_und->sigma_for, cpd_und->lda_for, alpha, beta, fx_spot,
        cpd_und->corr_times, cpd_und->correl_dom_for, cpd_und->correl_dom_fx,
        cpd_und->correl_for_fx, cpd_und->corr_n_times, cpd_und->dom_yc,
        cpd_und->for_yc, &(cpd_und->sigma_fx), 50, 2, 0.1, 0.1, 200);

    if (err) {
      goto FREE_RETURN;
    }

    /* filling the inst term structure of fwd vol fx */
    for (iNum = 0; iNum < cpd_und->sigma_n_fx; iNum++) {
      DateSigmaTab[0][iNum] = cpd_und->sigma_date_fx[iNum];
      DateSigmaTab[1][iNum] = cpd_und->sigma_fx[iNum];
    }

    /* init of the fx underlying */
    err = SrtInitFXUnd(fx_name, fx_spot, "STOCH_RATES", dom_name, for_name,
                       cpd_und->sigma_n_fx, 2, DateSigmaTab, fx_name);

    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Creation of the corresponding GRFN tableau */

  /* Filling the past and not yet paid  coupon */
  /*	Past and other cash-flows */

  /*	Funding */

  /*	Coupon fixed in the past and not yet paid */
  fund_ncf = 0;
  iNum = 0;
  while (iNum < fund_ncpn && fund_fix[iNum] < cpd_und->today + eod_fix_flag) {
    if (fund_pay[iNum] >= cpd_und->today + eod_pay_flag) {
      err = interp_basis(fund_basis[iNum], &bas);
      if (err) {
        goto FREE_RETURN;
      }

      fund_cfdtes[fund_ncf] = fund_pay[iNum];
      fund_cf[fund_ncf] = (fund_fix_cpn[iNum] + fund_mrg[iNum]) *
                          coverage((long)(fund_start[iNum] + 1.0e-08),
                                   (long)(fund_pay[iNum] + 1.0e-08), bas) *
                          fund_not;
      fund_ncf++;
    }

    iNum++;
  }

  /*	Initial exchange */
  if (start_date >= cpd_und->today + eod_pay_flag) {
    fund_cfdtes[fund_ncf] = start_date;
    if (for_fund) {
      fund_cf[fund_ncf] = -eq_init_ex;
    } else {
      fund_cf[fund_ncf] = -fund_not;
    }
    fund_ncf++;
  }

  /*	Final exchange */
  if (for_fund) {
    fund_cfdtes[fund_ncf] = fund_start_date;
    fund_cf[fund_ncf] = eq_final_ex - fund_not;
    fund_ncf++;
  }

  /*	PD */

  /*	Coupon fixed in the past and not yet paid */
  pd_ncf = 0;
  iNum = 0;
  while (iNum < pd_ncpn && pd_fix[iNum] < cpd_und->today + eod_fix_flag) {
    if (pd_pay[iNum] >= cpd_und->today + eod_pay_flag) {
      err = interp_basis(pd_basis[iNum], &bas);
      if (err) {
        goto FREE_RETURN;
      }

      pd_cfdtes[pd_ncf] = (long)(pd_pay[iNum] + 1.0e-08);

      /* Now pd_past_fix contains the fx fixing for the coupon */

      /* Old Method

        pd_cf[pd_ncf] = pd_past_fix[i]
              * coverage ((long) (pd_start[i] + 1.0e-08)        , (long)
        (pd_pay[i]
        + 1.0e-08)        , bas)
              * pd_not;
      */

      /* New Method */
      PdFixingRate = pd_alpha[iNum] + pd_beta[iNum] * pd_fix_fx[iNum];

      if (pd_floored[iNum])
        PdFixingRate = max(PdFixingRate, pd_floor[iNum]);

      if (pd_capped[iNum])
        PdFixingRate = min(PdFixingRate, pd_cap[iNum]);

      pd_cf[pd_ncf] = PdFixingRate *
                      coverage((long)(pd_start[iNum] + 1.0e-08),
                               (long)(pd_pay[iNum] + 1.0e-08), bas) *
                      pd_not;
      /* */
      pd_ncf++;
    }

    iNum++;
  }

  /*	Initial exchange */
  if (start_date >= cpd_und->today + eod_pay_flag) {
    pd_cfdtes[pd_ncf] = start_date;
    pd_cf[pd_ncf] = -pd_not;
    pd_ncf++;
  }

  /* Making the GRFN Tableau */
  err = cpd_make_grfn_tableau(
      /*	The structure */
      cpd_str,
      /*	Extra cash-flows coming from past coupons and initial exchange
       */
      fund_ncf, fund_cfdtes, fund_cf, pd_ncf, pd_cfdtes, pd_cf,
      /*	The underlyings */
      cpd_und->today, dom_name, for_name, fx_name,
      /*	The tableau and auxiliaries
              Must be allocated with maximum values
                      and initialised with zeros */
      &num_evt_dtes, &tableauCols, auxLen, &auxWidth, evt_dtes, tableauStrings,
      mask, aux);

  /* Number of rows in the GRFN tableau = number of event dates */
  tableauRows = num_evt_dtes;

  if (err) {
    goto FREE_RETURN;
  }

  /* Pricing with the alpha beta model */

  /* Fill the tab_barrier and bar_col array */
  for (iNum = 0; iNum < num_evt_dtes; iNum++) {
    tab_barrier[iNum] = aux[iNum][8];
    bar_col[iNum] = (int)(aux[iNum][9] + 1.0e-08);
  }

  /* Tranpose the auxilliaries */
  for (i = 0; i < MAX_COL; i++)
    for (j = 0; j < MAX_NEVT; j++)
      aux_transpose[i][j] = aux[j][i];

  /* launch the alpha beta tree on the grfn tableau */
  err = SrtGrfn3DFXAlphaBetaTree(
      fx_name, alpha, beta, num_evt_dtes, evt_dtes, tableauRows, tableauCols,
      tableauStrings, mask, auxWidth, auxLen, aux_transpose, 0, 0, tab_barrier,
      bar_col, nb_step, &num_prod, 1, &prod_val);

  if (err) {
    goto FREE_RETURN;
  }

  /* Filling the output */
  *fund_val = prod_val[2];
  *pd_val = prod_val[1];
  *call_val = prod_val[4];
  *call_stdev = 0;

FREE_RETURN:

  if (DateSigmaTab) {
    free_dmatrix(DateSigmaTab, 0, 1, 0,
                 max(cpd_und->sigma_n_fx, cpd_und->sigma_n_rates));
  }

  if (und_names) {
    free_smatrix_size(und_names, 0, 2, 0, 1, 256);
  }

  if (correl) {
    free_dmatrix(correl, 0, 2, 0, cpd_und->corr_n_times - 1);
  }

  if (correl_dates) {
    free_dvector(correl_dates, 0, cpd_und->corr_n_times - 1);
  }

  if (auxLen) {
    free_ivector(auxLen, 0, MAX_COL - 1);
  }

  if (evt_dtes) {
    free_lngvector(evt_dtes, 0, MAX_NEVT - 1);
  }

  if (tableauStrings) {
    free_smatrix_size(tableauStrings, 0, MAX_NEVT - 1, 0, MAX_COL - 1,
                      GRFN_DEF_ARGBUFSZ);
  }

  if (aux) {
    free_dmatrix(aux, 0, MAX_NEVT - 1, 0, MAX_COL - 1);
  }

  if (aux_transpose) {
    free_dmatrix(aux_transpose, 0, MAX_COL - 1, 0, MAX_NEVT - 1);
  }

  if (mask) {
    free_imatrix(mask, 0, MAX_NEVT - 1, 0, MAX_COL - 1);
  }

  if (tab_barrier) {
    free_dvector(tab_barrier, 0, MAX_NEVT - 1);
  }

  if (bar_col) {
    free_ivector(bar_col, 0, MAX_NEVT - 1);
  }

  if (maturity_opt) {
    free_dvector(maturity_opt, 0, num_fx_mkt_vol - 1);
  }

  if (fund_cfdtes) {
    free_lngvector(fund_cfdtes, 0, fund_ncpn - 1);
  }

  if (fund_cf) {
    free_dvector(fund_cf, 0, fund_ncpn - 1);
  }

  if (pd_cfdtes) {
    free_lngvector(pd_cfdtes, 0, pd_ncpn - 1);
  }

  if (pd_cf) {
    free_dvector(pd_cf, 0, pd_ncpn - 1);
  }

  return err;
}