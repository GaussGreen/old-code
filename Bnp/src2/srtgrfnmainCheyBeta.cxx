/**********************************************************************
 *      Name: SrtGrfnMainCheyBeta.cxx * Function: Entry point to GRFN with raw
 *data                       * Copyright: (C) Paribas Capital Markets Ltd. *
 *--------------------------------------------------------------------*
 *    Author: J.L. * Date: 27/06/01 *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "CheyBetaPricing.h"
#include "CheyBetagrfn.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allfx3f.h"

char *SrtGrfnCHEYBETApde(char *underlying, int numeventdates, long *eventdates,
                         long tableauRows, long tableauCols,
                         char ***tableauStrings, int **tableauMask,
                         long auxWidth, long *auxLen, double **aux, int nstept,
                         int nstepx, int nstepphi, int *nb_prod,
                         double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  GRFNPARMCHEYBETA grfn_prm;
  int forback;
  int flag = 0;
  long nstp;

  double next_d;

  double *evt_tms = NULL, *time = NULL, *ifr = NULL, *date = NULL;

  long *evt_dts = NULL;

  int *is_event = NULL;
  void **void_prm = NULL;

  double *dff = NULL, *gam = NULL, *gam_sqr = NULL;

  long today, spot_date;

  int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

  char *domestic_name, *yc;

  SrtUndPtr *und_ptr = NULL, und = NULL;

  int i, j;

  Err err = NULL;

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */
  err = srt_f_set_default_GrfnParams(&defParm);
  defParm.force_mc = 0;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           underlying, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  free_str = 1;

  /*	Now        , lookup underlyings involved */
  err = FIRSTGetUndFromDeal(&xStr, &num_und, &und_ptr);

  if (err) {
    goto FREE_RETURN;
  }

  if (num_und != 1) {
    err = "Product should involve only one underlying";
    goto FREE_RETURN;
  }

  und = und_ptr[0];

  /* look for the underlying name */
  und = lookup_und(underlying);
  if (!und) {
    err = "cannot find the underlying";
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) == CHEY_BETA) {
    domestic_name = und->underl_name;
  } else {
    err = "Model must be CHEYBETA";
    goto FREE_RETURN;
  }

  if (strcmp(domestic_name, und->underl_name)) {
    err = "Tableau uses different underlying";
    goto FREE_RETURN;
  }

  /* look for the today date */
  today = get_today_from_underlying(und);
  spot_date = get_spotdate_from_underlying(und);

  yc = (char *)get_ycname_from_irund(und);

  /* Get number of columns */
  err = FIRSTGetNumColFromDeal(&xStr, &num_col);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Get the maximum number of dfs required	*/
  err = FIRSTGetMaxNumDfFromDeal(&xStr, &max_num_df);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Next        , get the time steps */
  err = FIRSTGetEvtDatesFromDeal(&xStr, &num_evt, &evt_dts, &evt_tms);
  if (err) {
    goto FREE_RETURN;
  }

  /* discretise in time			*/

  nstp = num_evt;

  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  /*	Fill the time vector */
  err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, nstept);
  if (err) {
    goto FREE_RETURN;
  }

  void_prm = (void **)calloc(nstp, sizeof(void *));
  is_event = (int *)calloc(nstp, sizeof(int));
  date = (double *)calloc(nstp, sizeof(double));
  ifr = (double *)calloc(nstp, sizeof(double));

  if (max_num_df > 0) {
    dff = dvector(0, max_num_df - 1);
    gam = dvector(0, max_num_df - 1);
    gam_sqr = dvector(0, max_num_df - 1);
  }

  if (!void_prm || !is_event || !ifr || !date ||
      ((!dff || !gam || !gam_sqr) && (max_num_df > 0))) {
    err = "Memory allocation failure";
    goto FREE_RETURN;
  }

  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  j = xStr.num_evt - 1;
  next_d = evt_dts[j] + 1;

  for (i = nstp - 1; i >= 0; i--) {

    ifr[i] = swp_f_zr(date[i], next_d, yc);

    if (j >= 0 && fabs(time[i] - evt_tms[j]) < 1.0E-08) {
      grfn_prm = malloc(sizeof(grfn_parm_cheybeta));

      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;

      grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
      grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
      grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

      grfn_prm->dff = dff;
      grfn_prm->gam = gam;
      grfn_prm->gam_sqr = gam_sqr;

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;

      j--;
      while (j >= 0 && xStr.evt[j].evt == NULL) {
        j--;
      }
    } else if (j > 0 && xStr.am[j - 1]) {
      grfn_prm = malloc(sizeof(grfn_parm_cheybeta));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j - 1;

      grfn_prm->num_df = xStr.evt[j - 1].evt->dflen[0];
      grfn_prm->df_tms = xStr.evt[j - 1].evt->dft[0];
      grfn_prm->df_dts = xStr.evt[j - 1].evt->dfd[0];

      grfn_prm->dff = dff;
      grfn_prm->gam = gam;
      grfn_prm->gam_sqr = gam_sqr;

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;
    } else {
      is_event[i] = 0;
      void_prm[i] = NULL;
    }
    next_d = date[i];
  }

  /*	Eventually! call to function */

  *prod_val = dvector(0, num_col - 1);

  if (!*prod_val) {
    err = "Memory allocation failure";
    goto FREE_RETURN;
  }

  err = cheybeta_pricing_pde(
      /*	Time data		*/
      nstp, time, date,
      /*	Discretisation	*/
      nstepx, nstepphi,
      /*	Model data		*/
      und,
      /*	Product data */
      void_prm, is_event,
      /*	Market data */
      ifr, yc,
      /*	Payoff function */
      payoff_cheybeta_pde,
      /*	Result */
      num_col, *prod_val);

  if (err) {
    goto FREE_RETURN;
  }

  *nb_prod = num_col;

  /*	Add PV of Past */
  (*prod_val)[num_col - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (free_str) {
    FIRSTFreeUndFromDeal(num_und, &und_ptr);

    FIRSTFreeEvtDatesFromDeal(nstp, &evt_dts, &evt_tms);

    FIRSTFreeMktStruct(&xStr);
  }

  if (max_num_df > 0) {
    if (dff)
      free_dvector(dff, 0, max_num_df - 1);
    if (gam)
      free_dvector(gam, 0, max_num_df - 1);
    if (gam_sqr)
      free_dvector(gam_sqr, 0, max_num_df - 1);

    dff = NULL;
    gam = NULL;
    gam_sqr = NULL;
  }

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMCHEYBETA)void_prm[i];
        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (is_event)
    free(is_event);
  if (date)
    free(date);
  if (ifr)
    free(ifr);
  if (time)
    free(time);

  is_event = NULL;
  date = NULL;
  ifr = NULL;
  time = NULL;

  return err;
}

long Get_Up_Index(double T, double *Maturity, long nbrMat) {
  static int i;

  for (i = 0; (i < nbrMat) && (Maturity[i] < T); i++)
    ;
  return i;
}

char *SrtGrfnCHEYBETA(char *underlying, int numeventdates, long *eventdates,
                      long tableauRows, long tableauCols,
                      char ***tableauStrings, int **tableauMask, long auxWidth,
                      long *auxLen, double **aux,
                      /* param */
                      int method, /* 0 MC        , 1 MC adj        , 2 PDE */
                      int nstept,
                      /* for MC */
                      long numpaths, SrtMCSamType gen_method,
                      /* for PDE */
                      int nstepx, int nstepphi, int *nb_prod,
                      double **prod_val) {
  int free_str = 0;
  FIRSTAllMkts xStr;
  SrtGrfnParam defParm;
  GRFNPARMCHEYBETA grfn_prm;
  int forback;
  int flag = 0;
  long nstp;

  double next_d;

  double *evt_tms = NULL, *time = NULL, *ifr = NULL, *date = NULL;

  long *evt_dts = NULL;

  int *is_event = NULL;
  void **void_prm = NULL;

  int *stop_vol = NULL, *stop_lambda = NULL;
  long index;
  double *dff = NULL, *gam = NULL, *gam_sqr = NULL;

  long today, spot_date;

  int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

  char *domestic_name, *yc;

  SrtUndPtr *und_ptr = NULL, und = NULL;

  /* used in MC case */
  CHEYBETA_MDL mdl;

  long tempnt;
  int i, j;
  Err err = NULL;

  /*	Initialise the GRFN tableau */

  /*	First        , initialise the param struct */
  err = srt_f_set_default_GrfnParams(&defParm);
  if (method == 2)
    defParm.force_mc = 0;
  else
    defParm.force_mc = 1;

  err = FIRSTInitMktStruct(numeventdates, eventdates, tableauRows, tableauCols,
                           tableauStrings, tableauMask, auxWidth, auxLen, aux,
                           underlying, &defParm, &forback, &xStr);

  if (err) {
    goto FREE_RETURN;
  }

  free_str = 1;

  /*	Now        , lookup underlyings involved */
  err = FIRSTGetUndFromDeal(&xStr, &num_und, &und_ptr);

  if (err) {
    goto FREE_RETURN;
  }

  if (num_und != 1) {
    err = "Product should involve only one underlying";
    goto FREE_RETURN;
  }

  und = und_ptr[0];

  /* look for the underlying name */
  und = lookup_und(underlying);
  if (!und) {
    err = "cannot find the underlying";
    goto FREE_RETURN;
  }

  if (get_mdltype_from_irund(und) == CHEY_BETA) {
    domestic_name = und->underl_name;
  } else {
    err = "Model must be CHEYBETA";
    goto FREE_RETURN;
  }

  if (strcmp(domestic_name, und->underl_name)) {
    err = "Tableau uses different underlying";
    goto FREE_RETURN;
  }

  /* look for the today date */
  today = get_today_from_underlying(und);
  spot_date = get_spotdate_from_underlying(und);

  yc = (char *)get_ycname_from_irund(und);

  /* Get number of columns */
  err = FIRSTGetNumColFromDeal(&xStr, &num_col);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Get the maximum number of dfs required	*/
  err = FIRSTGetMaxNumDfFromDeal(&xStr, &max_num_df);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Next        , get the time steps */
  err = FIRSTGetEvtDatesFromDeal(&xStr, &num_evt, &evt_dts, &evt_tms);
  if (err) {
    goto FREE_RETURN;
  }

  /* discretise in time			*/

  nstp = num_evt;

  time = (double *)calloc(nstp, sizeof(double));
  if (!time) {
    err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
    goto FREE_RETURN;
  }
  memcpy(time, xStr.tms, nstp * sizeof(double));

  /* nt is the number of time step per year */
  tempnt = (long)(nstept * (long)(time[nstp - 1]));
  if (tempnt >= nstp + 1)
    nstept = tempnt;
  else
    nstept = nstp + 1;

  /*	Fill the time vector */
  err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, nstept);
  if (err) {
    goto FREE_RETURN;
  }

  /* if MC method merge Vol - Tau dates with Times */
  if (method < 2) {
    chey_beta_mdl_init(&mdl);
    chey_beta_mdl_build_from_und(&mdl, und);

    num_f_concat_vector(&nstp, &time, mdl.num_sigma, mdl.sigma_times);
    num_f_concat_vector(&nstp, &time, mdl.num_lambda, mdl.lambda_times);

    num_f_sort_vector(nstp, time);
    num_f_unique_vector(&nstp, time);

    stop_vol = (int *)calloc(nstp, sizeof(int));
    stop_lambda = (int *)calloc(nstp, sizeof(int));

    if (!(stop_vol && stop_lambda)) {
      err = "Memory allocation failure";
      goto FREE_RETURN;
    }

    memset(stop_vol, 0, nstp * sizeof(int));
    memset(stop_lambda, 0, nstp * sizeof(int));
  }

  void_prm = (void **)calloc(nstp, sizeof(void *));
  is_event = (int *)calloc(nstp, sizeof(int));
  date = (double *)calloc(nstp, sizeof(double));
  ifr = (double *)calloc(nstp, sizeof(double));

  if (!(is_event && date && ifr)) {
    err = "Memory allocation failure";
    goto FREE_RETURN;
  }

  /* modify the dates and time */
  for (i = 0; i < nstp; i++) {
    date[i] = today + DAYS_IN_YEAR * time[i];

    if (i > 0 && date[i] - date[i - 1] >= 1) {
      date[i] = (long)(date[i] + 1.0e-08);
      time[i] = YEARS_IN_DAY * (date[i] - today);
    }
  }

  /* populate the stop_vol and stop_lambda */
  if (method < 2) {
    /* initialise the stop_vol */
    for (i = 0; i < mdl.num_sigma - 1; i++) {
      index = Get_Up_Index(mdl.sigma_times[i], time, nstp);
      if (index == nstp)
        break;
      stop_vol[index] = i + 1;
    }

    /* initialise the stop_lambda */
    for (i = 0; i < mdl.num_lambda - 1; i++) {
      index = Get_Up_Index(mdl.lambda_times[i], time, nstp);
      if (index == nstp)
        break;
      stop_lambda[index] = i + 1;
    }
  }

  if (max_num_df > 0) {
    dff = dvector(0, max_num_df - 1);
    gam = dvector(0, max_num_df - 1);
    gam_sqr = dvector(0, max_num_df - 1);
  }

  if (!void_prm || !is_event ||
      ((!dff || !gam || !gam_sqr) && (max_num_df > 0))) {
    err = "Memory allocation failure";
    goto FREE_RETURN;
  }

  j = xStr.num_evt - 1;
  next_d = evt_dts[j] + 1;

  for (i = nstp - 1; i >= 0; i--) {
    ifr[i] = swp_f_zr(date[i], next_d, yc);

    if (j >= 0 && fabs(time[i] - evt_tms[j]) < 1.0E-08) {
      grfn_prm = malloc(sizeof(grfn_parm_cheybeta));

      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j;

      grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
      grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
      grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

      grfn_prm->dff = dff;
      grfn_prm->gam = gam;
      grfn_prm->gam_sqr = gam_sqr;

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;

      j--;
      while (j >= 0 && xStr.evt[j].evt == NULL) {
        j--;
      }
    } else if (j > 0 && xStr.am[j - 1]) {
      grfn_prm = malloc(sizeof(grfn_parm_cheybeta));
      grfn_prm->global = &xStr;
      grfn_prm->local = xStr.evt + j - 1;

      grfn_prm->num_df = xStr.evt[j - 1].evt->dflen[0];
      grfn_prm->df_tms = xStr.evt[j - 1].evt->dft[0];
      grfn_prm->df_dts = xStr.evt[j - 1].evt->dfd[0];

      grfn_prm->dff = dff;
      grfn_prm->gam = gam;
      grfn_prm->gam_sqr = gam_sqr;

      is_event[i] = 1;
      void_prm[i] = (void *)grfn_prm;
    } else {
      is_event[i] = 0;
      void_prm[i] = NULL;
    }
    next_d = date[i];
  }

  /*	Eventually! call to function */

  if (method == 2) {

    err = cheybeta_pricing_pde(
        /*	Time data		*/
        nstp, time, date,
        /*	Discretisation	*/
        nstepx, nstepphi,
        /*	Model data		*/
        und,
        /*	Product data */
        void_prm, is_event,
        /*	Market data */
        ifr, yc,
        /*	Payoff function */
        payoff_cheybeta_pde,
        /*	Result */
        num_col, prod_val[0]);

    if (err) {
      goto FREE_RETURN;
    }

  } else if ((method == 0) || (method == 1)) {

    err = cheybeta_pricing_mc(
        /*	Time data		*/
        nstp, time, date,
        /*	Discretisation	*/
        numpaths, method, /* 0 balantisam        , 1 balantisam adjusted */
        gen_method,
        /*	Model data		*/
        und,
        /*	Product data */
        void_prm, is_event,
        /*	Market data */
        ifr, yc, stop_vol, stop_lambda,
        /*	Payoff function */
        payoff_cheybeta_mc,
        /*	Result */
        num_col, prod_val);
    if (err) {
      goto FREE_RETURN;
    }
  } else {
    err = "Unknown method";
    goto FREE_RETURN;
  }

  *nb_prod = num_col;

  /*	Add PV of Past */

  prod_val[0][num_col - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

  if (free_str) {
    FIRSTFreeUndFromDeal(num_und, &und_ptr);

    FIRSTFreeEvtDatesFromDeal(nstp, &evt_dts, &evt_tms);

    FIRSTFreeMktStruct(&xStr);
  }

  if (void_prm) {
    for (i = 0; i < nstp; i++) {
      if (void_prm[i]) {
        grfn_prm = (GRFNPARMCHEYBETA)void_prm[i];
        free(grfn_prm);
      }
    }

    free(void_prm);
  }

  if (max_num_df > 0) {
    if (dff)
      free_dvector(dff, 0, max_num_df - 1);
    if (gam)
      free_dvector(gam, 0, max_num_df - 1);
    if (gam_sqr)
      free_dvector(gam_sqr, 0, max_num_df - 1);
  }

  if (is_event)
    free(is_event);
  if (date)
    free(date);
  if (ifr)
    free(ifr);
  if (time)
    free(time);
  if (stop_vol)
    free(stop_vol);
  if (stop_lambda)
    free(stop_lambda);

  return err;
}
