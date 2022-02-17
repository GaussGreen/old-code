/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_fx_init.c

   PURPOSE:     Provide a few utilities when dealing with a TermStruct:
                    - initialise a new TS:      srt_f_init_FX_TermStruct
                    - output an existing TS:    srt_f_display_FX_TermStruct
   ------------------------------------------------------------------------ */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts_fx.h"
#include "srt_h_ts_irm.h"

static Err make_fx_cumvol_vector(TermStruct *l);

static Err make_M_fx_vector(TermStruct *l);
static Err make_N_fx_vector(TermStruct *l);
static Err make_V_dx_vector(TermStruct *l);
static Err make_W_dx_vector(TermStruct *l);
static Err make_V_fx_vector(TermStruct *l);
static Err make_W_fx_vector(TermStruct *l);
static Err make_H_fd_vector(TermStruct *l);
static Err make_N_fx_vector(TermStruct *l);
static Err make_O_fd_vector(TermStruct *l);
static Err make_P_fd_vector(TermStruct *l);
static Err make_R_fd_vector(TermStruct *l);
static Err make_M_fd_vector(TermStruct *l);
static Err make_Q_fd_vector(TermStruct *l);
static Err make_S_fd_vector(TermStruct *l);
static Err make_P_fd_vector(TermStruct *l);
static Err make_Phi_fd_vector(TermStruct *l);
static Err make_T_fd_vector(TermStruct *l);
static Err make_U_fd_vector(TermStruct *l);
static Err make_V_fd_vector(TermStruct *l);
static Err make_W_fd_vector(TermStruct *l);
static Err make_X_dx_vector(TermStruct *l);
static Err make_Y_dx_vector(TermStruct *l);

/* -------------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------------
   The MAIN function to call when a TS has to be initialised for an FX
   underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_init_FX_TermStruct(Date today, double **sig, /* sig[col][row] */
                             int sig_cols, int nsig,

                             SrtMdlType mdl_type,

                             /* SMILE */
                             double beta,

                             /* STOCH VOL */
                             double vovol,

                             /* For FX_STOCH_RATE */
                             char *undName, char *domirundName,
                             char *forirundName,

                             /* OUTPUT */
                             TermStruct **ts) {
  Err err;
  SrtUndPtr dom_irund, for_irund;

  if (mdl_type == FX_STOCH_RATES) {
    /* The Domestic Name SHOULD be an underlying */
    dom_irund = lookup_und(domirundName);
    if (!dom_irund)
      return serror("Underlying %s not found for fx currencies", domirundName);

    /* The Foreign Name SHOULD be an underlying */
    for_irund = lookup_und(forirundName);
    if (!for_irund)
      return serror("Underlying %s not found for fx currencies", forirundName);

    /* Gets the Model Type for Both Underlyings */
    err = srt_f_init_FX_STOCH_RATES_TermStruct(undName, dom_irund, for_irund,
                                               ts, today, sig, sig_cols, nsig);

    if (err)
      return err;
  } else if ((mdl_type == BLACK_SCHOLES) || (mdl_type == NORMAL_BS)) {

    err = srt_f_init_FX_NORMLOG_TermStruct(ts, today, sig, sig_cols, nsig);
  }

  return err;
}

/* ------------------------------------------------------------------------- */
/* A function to output an existing TermStructure (Memory Allocation done
 * inside) */

Err srt_f_display_FX_TermStruct(char *szFxUndName,

                                long *sigma_n, double **sigma_date,
                                double **sigma,

                                long *corr_date_n, double **corr_date,
                                double **corr)

{

  SrtLst *l;
  TermStruct *ts;
  Err err = NULL;
  SrtUndPtr sFxUnd;
  String szDomUndName, szForUndName;
  SrtCorrLstPtr sCorrlist;
  SrtCorrLstVal *psCorrval;
  SrtListAtom *psCorr_lc;

  long corr_start_index;

  double today, time;

  /*-------------- Get the local volatility term
   * structure---------------------------*/

  sFxUnd = lookup_und(szFxUndName);
  if (!sFxUnd)
    return serror("Can not get the FX underlying name ");

  today = get_today_from_underlying(sFxUnd);
  err = get_underlying_ts(sFxUnd, &ts);
  if (err)
    serror("Can not get the underlying term struct ");

  l = ts->head;
  *sigma_n = 0;

  /* Allocate memory for the initial pointers */
  *sigma_date = (double *)malloc(sizeof(double));
  *sigma = (double *)malloc(sizeof(double));

  /* Loop on all the elements of the Term Sturcture */
  while (l != NULL) {

    (*sigma_n)++;
    *sigma_date = realloc(*sigma_date, (*sigma_n) * sizeof(double));
    *sigma = realloc(*sigma, (*sigma_n) * sizeof(double));
    (*sigma_date)[*sigma_n - 1] =
        (double)(((FxTermStructVal *)l->element->val.pval)->date);
    (*sigma)[*sigma_n - 1] = ((FxTermStructVal *)l->element->val.pval)->sigx;

    l = l->next;
  }

  /*------------- Get the correlation term
   * structure-----------------------------*/

  szForUndName = get_forname_from_fxund(sFxUnd);
  szDomUndName = get_domname_from_fxund(sFxUnd);

  sCorrlist = srt_f_GetTheCorrelationList();
  if (!sCorrlist->head->element)
    return serror("correlation list improperly initialised");

  psCorr_lc = sCorrlist->head;

  /* Allocate memory for the initial pointers */

  *corr_date = (double *)malloc(sizeof(double));
  *corr = (double *)malloc(sizeof(double));

  (*corr_date_n) = 0;
  while (psCorr_lc != NULL) {
    psCorrval = (SrtCorrLstVal *)psCorr_lc->element->val.pval;
    if (psCorrval->ncorrel < 3)
      return serror("Fatal error in srt_f_display_FX_TermStruct: Wrong number "
                    "of correlation");

    (*corr_date_n)++;
    *corr_date = realloc(*corr_date, (*corr_date_n) * sizeof(double));

    /* old version wrong when corr matrix contains more than 3 underlyings... */
    /*
     *corr = realloc(*corr  ,(*corr_date_n)*psCorrval->ncorrel*sizeof(double));
     */

    /* obviously we only ask for three correlations per date */
    *corr = realloc(*corr, (*corr_date_n) * 3 * sizeof(double));

    (*corr_date)[*corr_date_n - 1] = psCorrval->date;
    time = (psCorrval->date - today) * YEARS_IN_DAY;

    /* old version wrong for the same raisons */
    /*
    corr_start_index = (*corr_date_n-1)*psCorrval->ncorrel;
    */

    corr_start_index = (*corr_date_n - 1) * 3;

    err = srt_f_get_corr_from_CorrList(sCorrlist, szDomUndName, szForUndName,
                                       time, &(*corr)[corr_start_index]);

    if (err)
      return serror("Fatal Error in srt_f_get_corr_from_CorrList");

    err = srt_f_get_corr_from_CorrList(sCorrlist, szDomUndName, szFxUndName,
                                       time, &(*corr)[corr_start_index + 1]);

    if (err)
      return serror("Fatal Error in srt_f_get_corr_from_CorrList");

    err = srt_f_get_corr_from_CorrList(sCorrlist, szForUndName, szFxUndName,
                                       time, &(*corr)[corr_start_index + 2]);

    if (err)
      return serror("Fatal Error in srt_f_get_corr_from_CorrList");

    psCorr_lc = psCorr_lc->next;
  }

  return NULL;
}

/* ------------------------------------------------------------------------- */
/* A function to output an existing TermStructure without the correlation part
   and without the check. This is for Fx TS defined with BS type */

Err srt_f_display_FXBS_TermStruct(char *szFxUndName,

                                  long *sigma_n, double **sigma_date,
                                  double **sigma)

{

  SrtLst *l;
  TermStruct *ts;
  Err err = NULL;
  SrtUndPtr sFxUnd;

  double today;

  /*-------------- Get the local volatility term
   * structure---------------------------*/

  sFxUnd = lookup_und(szFxUndName);
  if (!sFxUnd)
    return serror("Can not get the FX underlying name ");

  today = get_today_from_underlying(sFxUnd);
  err = get_underlying_ts(sFxUnd, &ts);
  if (err)
    serror("Can not get the underlying term struct ");

  l = ts->head;
  *sigma_n = 0;

  /* Allocate memory for the initial pointers */
  *sigma_date = (double *)malloc(sizeof(double));
  *sigma = (double *)malloc(sizeof(double));

  /* Loop on all the elements of the Term Sturcture */
  while (l != NULL) {

    (*sigma_n)++;
    *sigma_date = realloc(*sigma_date, (*sigma_n) * sizeof(double));
    *sigma = realloc(*sigma, (*sigma_n) * sizeof(double));
    (*sigma_date)[*sigma_n - 1] =
        (double)(((FxTermStructVal *)l->element->val.pval)->date);
    (*sigma)[*sigma_n - 1] = ((FxTermStructVal *)l->element->val.pval)->sigx;

    l = l->next;
  }

  return NULL;
}

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------
   The MAIN function to free a Term Structure of Volatility initialised
   for a FOREIGN EXCHANGE underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_free_FX_TermStruct(TermStruct **l) {
  Err err;

  if (!l)
    return NULL;

  err = srt_f_lstfree(*l, SRT_YES);
  if (err)
    return err;

  srt_free(*l);

  return NULL;
}

/* ----------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------
                              NON SMILE SPECIFIC FUNCTIONS
   -------------------------------------------------------------------------------
 */

/* Function required to free a TermStructObject once attached into a linked list
 */

Err srt_f_FxTermStructValFree(void *tsvalptr) {
  FxTermStructVal *tsval = (FxTermStructVal *)tsvalptr;
  Err err = NULL;

  srt_free(tsval);

  return err;
}

/* ------------------------------------------------------------------------ */

/* Initialises the Volatility Term Structure for Non smile FX Models */

Err srt_f_init_FX_NORMLOG_TermStruct(TermStruct **ts, Date today,
                                     double **sig_data, int sig_cols,
                                     int num_sig) {
  int cur_sig;
  FxTermStructVal *tsval;
  Err err = 0;
  long i;
  long ticker;

  /* Check on the number of sigma */
  if (sig_cols != 2)
    return serror("Need TWO columns for the Sigma Curve");

  /* Create a blank double linked ts that will become the TermStruct */
  srt_f_lstcreate(ts, "TermStruct");

  cur_sig = 0;

  /* If the data exists but is NULL: consider garbage */
  for (i = 0; i < num_sig; i++) {
    if (!sig_data[0][i] || !sig_data[1][i])
      break;
  }
  num_sig = i;

  /* Current index for tau and sigma start at a date >= today */
  while ((cur_sig < num_sig) && ((sig_data[0][cur_sig] - today) <= EPSILON))
    cur_sig++;

  /* Return error if there is no acceptable sigma or tau value */
  if (cur_sig == num_sig)
    return serror("Init_TermStruct err: no current data !!!");

  for (; cur_sig < num_sig; cur_sig++) {

    /* Creation of a new TermStruct2 */
    tsval = (FxTermStructVal *)srt_calloc(1, sizeof(FxTermStructVal));

    /* transfer of information */
    tsval->date = DTOL(sig_data[0][cur_sig]);
    tsval->time = (sig_data[0][cur_sig] - today) * YEARS_IN_DAY;
    tsval->sigx = sig_data[1][cur_sig];

    /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
    srt_f_lstins(*ts, "BsTsAtom", tsval->date, OBJ_PTR_FX_TermStruct,
                 (void *)tsval, &srt_f_FxTermStructValFree, &ticker);
  }

  /* Computes cum_vol values */
  err = make_fx_cumvol_vector(*ts);
  if (err) {
    srt_f_lstfree(*ts, SRT_YES);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_fx_cumvol_vector(TermStruct *l) {
  FxTermStructVal *c, *p;
  SrtLst *lc;

  Err err = 0;

  /* Starts at the beginning of the Term Structure */
  lc = l->head;

  /* Gets the TermStructVal attached to this element */
  c = (FxTermStructVal *)lc->element->val.pval;

  /* Computes the value of the cumulative volatility */
  c->int_sig2_dt = c->sigx * c->sigx * c->time;
  c->int_sig_dt = c->sigx * c->time;

  /* Moves on to the next point */
  lc = lc->next;

  /* Loops on all the elements of the TermStructure */
  while (lc != NULL) {
    p = c;
    c = (FxTermStructVal *)lc->element->val.pval;

    /* The Cumulative Volatility as the integral of Sigma^2*t */
    c->int_sig2_dt = p->int_sig2_dt + c->sigx * c->sigx * (c->time - p->time);

    /* The Integral of Sigma * t */
    c->int_sig_dt = p->int_sig_dt + c->sigx * (c->time - p->time);

    /* Moves on to the next element in the Term Structure */
    lc = lc->next;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

/*-----------------For Jumping Numeraire
 * ------------------------------------------*/

Err srt_f_init_FX_STOCH_RATES_TermStruct(char *undName, SrtUndPtr dom_und,
                                         SrtUndPtr for_und, TermStruct **fx_ts,
                                         Date today, double **sig_data,
                                         int sig_cols, int num_sig)

{
  FxTermStructVal *fx_tsval, *fx_tsval_p;
  Err err = 0;
  long i;
  long ticker;
  TermStruct *subls;
  SrtCorrLstPtr corrlist;
  TermStruct *domts, *forts;
  String dom_und_name, for_und_name;
  SrtMdlType dom_mdl_type;
  SrtMdlType for_mdl_type;
  SrtListAtom *fx_atom, *dom_atom, *for_atom;
  IrmTermStructVal *dom_tsval1F, *for_tsval1F;
  TwoFacIrmTermStructVal *dom_tsval2F, *for_tsval2F;
  double time;

  /* First Creates a Term Structure that is like a normal FX one */
  err = srt_f_init_FX_NORMLOG_TermStruct(fx_ts, today, sig_data, sig_cols,
                                         num_sig);
  if (err) {
    return err;
  }
  /* Gets the Model Type for Both Underlyings */
  err = get_underlying_mdltype(dom_und, &dom_mdl_type);
  if (err) {

    return err;
  }
  err = get_underlying_mdltype(for_und, &for_mdl_type);
  if (err) {

    return err;
  }

  /* If not lgm or vasicek / lgm or vasicek : nothing else to be done */
  if ((dom_mdl_type != LGM && dom_mdl_type != VASICEK) ||
      (for_mdl_type != LGM && for_mdl_type != VASICEK))
    return NULL;

  /* ------------------- Workswith the Domestic Term Structure -------------- */

  /* Gets the Domestic Underlying */
  err = get_underlying_ts(dom_und, &domts);
  if (err) {
    return serror("Can not get the domestic term struct in "
                  "srt_f_init_FX_STOCH_RATES_TermStruct.");
  }

  /* Starts at the top of the FX Term Struct */
  fx_atom = (*fx_ts)->head;

  /* Starts at the top of the DOM LGM Term Struct */
  dom_atom = domts->head;

  /* Loops on the Domestic Term Struct dates to add them to the FX one */

  if (get_mdldim_from_irund(dom_und) != TWO_FAC) {
    while (dom_atom != NULL) {
      /* Gets the Domestic LGm Term Struct Val attached */
      dom_tsval1F = (IrmTermStructVal *)dom_atom->element->val.pval;

      /* Moves on to the FX date just after (or on) this date */
      while ((fx_atom != NULL) &&
             (((FxTermStructVal *)fx_atom->element->val.pval)->time <
              dom_tsval1F->time))
        fx_atom = fx_atom->next;
      if (fx_atom == NULL)
        fx_atom = (*fx_ts)->tail;

      fx_tsval = ((FxTermStructVal *)fx_atom->element->val.pval);

      /* If not on the same date  , we add the date onto the FX ts */
      if (fx_tsval->time != dom_tsval1F->time) {
        /* Creation of a new IrmTermStructVal (simple element of the Term
         * Struct) */
        fx_tsval = (FxTermStructVal *)srt_calloc(1, sizeof(FxTermStructVal));

        /* Adds the FX specific information */
        fx_tsval->date = dom_tsval1F->date;
        fx_tsval->time = dom_tsval1F->time;
        fx_tsval->sigx = find_fx_sig(fx_tsval->time, *fx_ts);

        /* Insert the TermStructAtom in the TS (==linked ts sorted by key =
         * date)  */
        srt_f_lstins(*fx_ts, "FXTsAtom", fx_tsval->date, OBJ_PTR_FX_TermStruct,
                     (void *)fx_tsval, &srt_f_FxTermStructValFree, &ticker);
      }

      /* Moves on to the next element in the Domestic Term Structure */
      dom_atom = dom_atom->next;

    } /* END of loop on Domestic Term Struct Dates */
  } else {
    while (dom_atom != NULL) {
      /* Gets the Domestic LGm Term Struct Val attached */
      dom_tsval2F = (TwoFacIrmTermStructVal *)dom_atom->element->val.pval;

      /* Moves on to the FX date just after (or on) this date */
      while ((fx_atom != NULL) &&
             (((FxTermStructVal *)fx_atom->element->val.pval)->time <
              dom_tsval2F->time))
        fx_atom = fx_atom->next;
      if (fx_atom == NULL)
        fx_atom = (*fx_ts)->tail;

      fx_tsval = ((FxTermStructVal *)fx_atom->element->val.pval);

      /* If not on the same date  , we add the date onto the FX ts */
      if (fx_tsval->time != dom_tsval2F->time) {
        /* Creation of a new IrmTermStructVal (simple element of the Term
         * Struct) */
        fx_tsval = (FxTermStructVal *)srt_calloc(1, sizeof(FxTermStructVal));

        /* Adds the FX specific information */
        fx_tsval->date = dom_tsval2F->date;
        fx_tsval->time = dom_tsval2F->time;
        fx_tsval->sigx = find_fx_sig(fx_tsval->time, *fx_ts);

        /* Insert the TermStructAtom in the TS (==linked ts sorted by key =
         * date)  */
        srt_f_lstins(*fx_ts, "FXTsAtom", fx_tsval->date, OBJ_PTR_FX_TermStruct,
                     (void *)fx_tsval, &srt_f_FxTermStructValFree, &ticker);
      }

      /* Moves on to the next element in the Domestic Term Structure */
      dom_atom = dom_atom->next;

    } /* END of loop on Domestic Term Struct Dates */
  }

  /* ------------------- Workswith the Foreign Term Structure -------------- */

  /* Gets the Foreign Underlying */
  err = get_underlying_ts(for_und, &forts);
  if (err)

    /* Starts at the top of the FX Term Struct */
    fx_atom = (*fx_ts)->head;

  /* Starts at the top of the FOR LGM Term Struct */
  for_atom = forts->head;

  /* Loops on the Foreign Term Struct dates to add them to the FX one */

  if (get_mdldim_from_irund(for_und) != TWO_FAC) {
    while (for_atom != NULL) {
      /* Gets the Foreign LGm Term Struct Val attached */
      for_tsval1F = (IrmTermStructVal *)for_atom->element->val.pval;

      /* Moves on to the FX date just after (or on) this date */
      while ((fx_atom != NULL) &&
             (((FxTermStructVal *)fx_atom->element->val.pval)->time <
              for_tsval1F->time))
        fx_atom = fx_atom->next;
      if (fx_atom == NULL)
        fx_atom = (*fx_ts)->tail;

      fx_tsval = ((FxTermStructVal *)fx_atom->element->val.pval);

      /* If not on the same date  , we add the date onto the FX ts */
      if (fx_tsval->time != for_tsval1F->time) {
        /* Creation of a new IrmTermStructVal (simple element of the Term
         * Struct) */
        fx_tsval = (FxTermStructVal *)srt_calloc(1, sizeof(FxTermStructVal));

        /* Adds the FX specific information */
        fx_tsval->date = for_tsval1F->date;
        fx_tsval->time = for_tsval1F->time;
        fx_tsval->sigx = find_fx_sig(fx_tsval->time, *fx_ts);

        /* Insert the TermStructAtom in the TS (==linked ts sorted by key =
         * date)  */
        srt_f_lstins(*fx_ts, "FXTsAtom", fx_tsval->date, OBJ_PTR_FX_TermStruct,
                     (void *)fx_tsval, &srt_f_FxTermStructValFree, &ticker);
      }

      /* Moves on to the next element in the Foreign Term Structure */
      for_atom = for_atom->next;

    } /* END of loop on Foreign Term Struct Dates */
  } else {
    while (for_atom != NULL) {
      /* Gets the Foreign LGm Term Struct Val attached */
      for_tsval2F = (TwoFacIrmTermStructVal *)for_atom->element->val.pval;

      /* Moves on to the FX date just after (or on) this date */
      while ((fx_atom != NULL) &&
             (((FxTermStructVal *)fx_atom->element->val.pval)->time <
              for_tsval2F->time))
        fx_atom = fx_atom->next;
      if (fx_atom == NULL)
        fx_atom = (*fx_ts)->tail;

      fx_tsval = ((FxTermStructVal *)fx_atom->element->val.pval);

      /* If not on the same date  , we add the date onto the FX ts */
      if (fx_tsval->time != for_tsval2F->time) {
        /* Creation of a new IrmTermStructVal (simple element of the Term
         * Struct) */
        fx_tsval = (FxTermStructVal *)srt_calloc(1, sizeof(FxTermStructVal));

        /* Adds the FX specific information */
        fx_tsval->date = for_tsval2F->date;
        fx_tsval->time = for_tsval2F->time;
        fx_tsval->sigx = find_fx_sig(fx_tsval->time, *fx_ts);

        /* Insert the TermStructAtom in the TS (==linked ts sorted by key =
         * date)  */
        srt_f_lstins(*fx_ts, "FXTsAtom", fx_tsval->date, OBJ_PTR_FX_TermStruct,
                     (void *)fx_tsval, &srt_f_FxTermStructValFree, &ticker);
      }

      /* Moves on to the next element in the Foreign Term Structure */
      for_atom = for_atom->next;

    } /* END of loop on Foreign Term Struct Dates */
  }

  /* ---------  Add Domestic and Foreign specific information in the FX Term
   * Struct ------- */

  /* Gets the Underlying Names (for their correlation) */
  dom_und_name = get_underlying_name(dom_und);
  for_und_name = get_underlying_name(for_und);

  /* Gets the Global Correlation List */
  corrlist = srt_f_GetTheCorrelationList();
  if (!corrlist->head->element)
    return serror("Correlation list improperly initialised...");

  /* Starts at the beginning of the FX Term Structure */
  fx_atom = (*fx_ts)->head;

  /* Loops on all the Dates of the Fx Term Struct to add ... */
  while (fx_atom != NULL) {
    /* Gets the TermStructVal attached to the atom of the TermStruct */
    fx_tsval = (FxTermStructVal *)fx_atom->element->val.pval;

    time = 0;
    if (fx_atom != (*fx_ts)->head) {
      fx_tsval_p = (FxTermStructVal *)fx_atom->previous->element->val.pval;
      time = fx_tsval_p->time;
    }

    /* Gets all the Correlations for this date */
    err = srt_f_get_corr_from_CorrList(corrlist, dom_und_name, for_und_name,
                                       fx_tsval->time, &fx_tsval->rhofd);
    err = srt_f_get_corr_from_CorrList(corrlist, undName, for_und_name,
                                       fx_tsval->time, &fx_tsval->rhofx);
    err = srt_f_get_corr_from_CorrList(corrlist, dom_und_name, undName,
                                       fx_tsval->time, &fx_tsval->rhodx);
    if (err) {
      srt_f_lstfree(*fx_ts, SRT_YES);
      return err;
    }

    /* Stores the information in the TermStructVal */
    for (i = 0; i < 2; i++) {
      if (i == 0)
        subls = domts;
      else if (i == 1)
        subls = forts;

      fx_tsval->StochRates_Ts[i].sig = find_sig(fx_tsval->time, subls);
      fx_tsval->StochRates_Ts[i].tau = find_tau(fx_tsval->time, subls);
      fx_tsval->StochRates_Ts[i].lambda =
          (1.0 / fx_tsval->StochRates_Ts[i].tau);
      fx_tsval->StochRates_Ts[i].F = F_func(time, subls);
      G_H_func(time, subls, &fx_tsval->StochRates_Ts[i].G,
               &fx_tsval->StochRates_Ts[i].H);
      fx_tsval->StochRates_Ts[i].Psi = Psi_func(time, subls);
      fx_tsval->StochRates_Ts[i].I = I_func(time, subls);
      fx_tsval->StochRates_Ts[i].J = J_func(time, subls);
      fx_tsval->StochRates_Ts[i].K = K_func(time, subls);
      fx_tsval->StochRates_Ts[i].L = L_func(time, subls);
      fx_tsval->StochRates_Ts[i].O = O_func(time, subls);
      fx_tsval->StochRates_Ts[i].Q = Q_func(time, subls);
      fx_tsval->StochRates_Ts[i].Phi = Phi_func(time, subls);

    } /* END of loop on dom and und */

    /* Moves on to the next element in the Forex Term Structure */
    fx_atom = fx_atom->next;

  } /* END of loop on FX TermStruct dates */

  /* Builds the Fx Cumulative Volatilities (not used in this case  , but why
   * not) */
  err = make_fx_cumvol_vector(*fx_ts);
  if (err) {
    srt_f_lstfree(*fx_ts, SRT_YES);
    return err;
  }

  /* Builds the Nice Integrals we need for the Fx Stoch Rates with Jumping
   * (LGM/LGM) */
  err = make_M_fx_vector(*fx_ts);
  err = make_N_fx_vector(*fx_ts);
  err = make_O_fd_vector(*fx_ts);
  err = make_P_fd_vector(*fx_ts);
  err = make_R_fd_vector(*fx_ts);
  err = make_M_fd_vector(*fx_ts);
  err = make_H_fd_vector(*fx_ts);
  err = make_Q_fd_vector(*fx_ts);
  err = make_S_fd_vector(*fx_ts);
  err = make_V_dx_vector(*fx_ts);
  err = make_W_dx_vector(*fx_ts);
  err = make_V_fx_vector(*fx_ts);
  err = make_W_fx_vector(*fx_ts);
  err = make_Phi_fd_vector(*fx_ts);
  err = make_S_fd_vector(*fx_ts);
  err = make_T_fd_vector(*fx_ts);
  err = make_U_fd_vector(*fx_ts);
  err = make_V_fd_vector(*fx_ts);
  err = make_W_fd_vector(*fx_ts);
  err = make_X_dx_vector(*fx_ts);
  err = make_Y_dx_vector(*fx_ts);

  if (err) {
    srt_f_lstfree(*fx_ts, SRT_YES);
  }

  /* Return an error if there is any */
  return err;

} /* END  Err srt_f_init_FX_STOCH_RATES_TermStruct(... ) */

/* --------------------------------------------------------------------------------------
 */

/* ------------------------------ STATIC FUNCTIONS
 * ------------------------------------- */

static Err make_M_fx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp;
  double sig, tau, F;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->M_fx = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    sig = p->StochRates_Ts[1].sig;
    tau = p->StochRates_Ts[1].tau;
    F = p->StochRates_Ts[1].F;

    temp = (exp(time / tau) - 1) * tau;

    c->M_fx = p->M_fx + p->rhofx * p->sigx * (sig / F) * temp;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_N_fx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp;
  double fortau, forsig, forF;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->N_fx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    fortau = p->StochRates_Ts[1].tau;
    forsig = p->StochRates_Ts[1].sig;
    forF = p->StochRates_Ts[1].F;
    temp = (1 - exp(-time / fortau)) * fortau;

    c->N_fx = p->N_fx + p->M_fx * forF * temp +
              p->rhofx * p->sigx * forsig * fortau * (time - temp);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_M_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double domsig, forsig, domtau, fortau, lambda, domPsi, forPsi, domF, forF;
  double temp, temp1, temp2, temp3;
  double time;

  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->P_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;
    domPsi = p->StochRates_Ts[0].Psi;
    forPsi = p->StochRates_Ts[1].Psi;

    temp = domsig * forsig;
    temp1 = (exp(time / domtau) - 1) * domtau;
    temp2 = (exp(time / fortau) - 1) * fortau;
    lambda = (1.0 / domtau) + (1.0 / fortau);
    temp3 = (exp(lambda * time) - 1) / lambda;

    c->M_fd =
        p->M_fd + p->rhofd * (temp / (domF * forF)) * domPsi * forPsi * temp3 +
        p->rhofd * (temp / domF) * fortau * domPsi * (temp3 - temp1) +
        p->rhofd * (temp / forF) * domtau * forPsi * (temp3 - temp2) +
        p->rhofd * temp * domtau * fortau * (temp3 + time - temp1 - temp2);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

/* --------------------------------------------------------------------------------------
 */

static Err make_O_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1;
  double dom_sig, dom_tau, dom_F, for_sig, for_tau, for_F;
  double lambda;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->O_fd = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    dom_tau = p->StochRates_Ts[0].tau;
    dom_sig = p->StochRates_Ts[0].sig;
    dom_F = p->StochRates_Ts[0].F;

    for_tau = p->StochRates_Ts[1].tau;
    for_sig = p->StochRates_Ts[1].sig;
    for_F = p->StochRates_Ts[1].F;

    lambda = 1.0 / dom_tau + 1.0 / for_tau;

    temp = (dom_sig / dom_F) * (for_sig / for_F);
    temp1 = (exp(time * lambda) - 1) / lambda;

    c->O_fd = p->O_fd + p->rhofd * temp * temp1;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_P_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1, temp2;
  double dom_sig, dom_tau, dom_F, for_Psi;
  double for_sig, for_tau, for_F;
  double lambda;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->P_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    dom_sig = p->StochRates_Ts[0].sig;
    for_sig = p->StochRates_Ts[1].sig;
    dom_tau = p->StochRates_Ts[0].tau;
    for_tau = p->StochRates_Ts[1].tau;
    dom_F = p->StochRates_Ts[0].F;
    for_F = p->StochRates_Ts[1].F;
    for_Psi = p->StochRates_Ts[1].Psi;

    temp = (dom_sig * for_sig) / (dom_F * for_F);
    lambda = (1.0 / dom_tau) + (1.0 / for_tau);
    temp1 = (exp(lambda * time) - 1) / lambda;
    temp2 = (exp(time / dom_tau) - 1) * dom_tau;

    c->P_fd = p->P_fd + p->rhofd * for_Psi * temp * temp1 +
              p->rhofd * temp * for_F * (temp1 - temp2) * for_tau;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_R_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1, temp2;
  double domsig, domtau, domF, domPsi;
  double forsig, fortau, forF;
  double lambda;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->P_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;
    domPsi = p->StochRates_Ts[0].Psi;

    temp = (domsig * forsig) / (domF * forF);
    lambda = (1.0 / domtau) + (1.0 / fortau);
    temp1 = (exp(lambda * time) - 1) / lambda;
    temp2 = (exp(time / fortau) - 1) * fortau;

    c->R_fd = p->R_fd + p->rhofd * domPsi * temp * temp1 +
              p->rhofd * temp * domF * (temp1 - temp2) * domtau;
  }
  return err;
}

static Err make_S_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1;
  double domsig, forsig, domtau, fortau, domF, forF, lambda;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->S_fd = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;

    temp = (domsig * forsig) / (domF * forF);
    lambda = 1 / domtau + 1 / fortau;
    temp1 = (exp(time * lambda) - 1) / lambda;

    c->S_fd = p->S_fd + p->rhofd * temp * temp1;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_Phi_fd_vector(TermStruct *l) {

  FxTermStructVal *c;
  SrtLst *lc;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->Phi_fd = 0.0;

  while (lc->next != NULL) {
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;
    c->Phi_fd = c->StochRates_Ts[0].F * c->StochRates_Ts[1].F * c->S_fd;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_H_fd_vector(TermStruct *l) {

  FxTermStructVal *c;
  SrtLst *lc;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->H_fd = 0.0;

  while (lc->next != NULL) {
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;
    c->H_fd =
        c->StochRates_Ts[1].F * (c->R_fd - c->StochRates_Ts[0].Psi * c->O_fd);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_Q_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time, temp, temp1, temp2;
  double domsig, forsig, domtau, fortau, domF, forF, domPsi, forPsi, lambda;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->Q_fd = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;
    domPsi = p->StochRates_Ts[0].Psi;
    forPsi = p->StochRates_Ts[1].Psi;

    temp = (1.0 - exp(-time / fortau)) * fortau;
    lambda = 1.0 / domtau + 1.0 / fortau;
    temp1 = (1.0 - exp(-lambda * time)) / lambda;
    temp2 = (exp(time / domtau) - 1.0) * domtau;

    c->Q_fd =
        p->Q_fd - forF * p->O_fd * domPsi * temp -
        domF * forF * p->O_fd * domtau * (temp - temp1) -
        p->rhofd * domsig * forsig * domPsi / (domF * lambda) * (temp2 - temp) -
        p->rhofd * domsig * forsig * domtau / lambda *
            (temp2 - temp - time + temp1) +
        forF * p->P_fd * temp +
        p->rhofd * domPsi * domsig * forsig / (domF * lambda) * (temp2 - temp) +
        p->rhofd * domsig * forsig *
            (temp2 / lambda - temp / lambda - time * fortau + temp * fortau) *
            domtau;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_T_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time, temp, temp1, temp2;
  double domtau, fortau, domsig, forsig, domF, forF, forPsi, lambda;
  Err err = 0;
  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->T_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;
    forPsi = p->StochRates_Ts[1].Psi;
    lambda = 1 / domtau + 1 / fortau;

    temp = (1 - exp(-time / fortau)) * fortau;
    temp1 = (domsig * forsig) / (domF * forF);
    temp2 = (1 - exp(time / domtau)) * domtau;

    c->T_fd = p->T_fd + forF * p->S_fd * temp +
              p->rhofd * forF * temp1 * (-temp2 / lambda - temp / lambda);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_U_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double domtau, fortau, lambda, domsig, forsig, domF, forF, domPsi, temp,
      temp1, temp2;
  Err err = 0;
  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->U_fd = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domtau = p->StochRates_Ts[0].tau;
    fortau = p->StochRates_Ts[1].tau;
    domsig = p->StochRates_Ts[0].sig;
    forsig = p->StochRates_Ts[1].sig;
    domF = p->StochRates_Ts[0].F;
    forF = p->StochRates_Ts[1].F;
    domPsi = p->StochRates_Ts[0].Psi;
    temp = (1 - exp(-time / fortau)) * fortau;
    temp1 = (exp(time / domtau) - 1) * domtau;
    lambda = 1 / domtau + 1 / fortau;
    temp2 = (1 - exp(-lambda * time)) / lambda;

    c->U_fd = p->U_fd + forF * domPsi * p->O_fd * temp +
              p->O_fd * domF * forF * domtau * (temp - temp2) +
              domPsi * (p->rhofd / lambda) * (domsig * forsig / (domF)) *
                  (temp1 - temp) +
              p->rhofd * (domsig * forsig * domtau / lambda) *
                  (temp1 - temp - time + temp2);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_V_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double dDomSig, dForSig, dDomTau, dDomF, dTemp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->V_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    dDomSig = p->StochRates_Ts[0].sig;
    dForSig = p->StochRates_Ts[1].sig;
    dDomTau = p->StochRates_Ts[0].tau;
    dDomF = p->StochRates_Ts[0].F;
    dTemp = (exp(time / dDomTau) - 1) * dDomTau;

    c->V_fd = p->V_fd + p->rhofd * (dDomSig * dForSig / dDomF) * dTemp;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_W_fd_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double dDomSig, dForSig, dDomTau, dDomF, dDomPsi, dTemp;

  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->W_fd = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    dDomSig = p->StochRates_Ts[0].sig;
    dForSig = p->StochRates_Ts[1].sig;
    dDomTau = p->StochRates_Ts[0].tau;
    dDomF = p->StochRates_Ts[0].F;
    dDomPsi = p->StochRates_Ts[0].Psi;
    dTemp = (exp(time / dDomTau) - 1) * dDomTau;

    c->W_fd = p->W_fd +
              p->rhofd * (dDomSig * dForSig / dDomF) * dDomPsi * dTemp +
              p->rhofd * dDomSig * dForSig * dDomTau * (dTemp - time);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_X_dx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double domsig, domtau, domF, sqdomF, temp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;
  c->X_dx = 0.0;
  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    domtau = p->StochRates_Ts[0].tau;
    domF = p->StochRates_Ts[0].F;
    sqdomF = domF * domF;
    temp = 0.5 * (exp(2 * time / domtau) - 1) * domtau;

    c->X_dx = p->X_dx + p->rhodx * p->sigx * (domsig / sqdomF) * temp;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

/* --------------------------------------------------------------------------------------
 */

static Err make_Y_dx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double domsig, domtau, domF, sqdomF, domPsi;
  double temp, temp1;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;

  c->Y_dx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    domtau = p->StochRates_Ts[0].tau;
    domF = p->StochRates_Ts[0].F;
    sqdomF = domF * domF;
    domPsi = p->StochRates_Ts[0].Psi;
    temp = (exp(time / domtau) - 1) * domtau;
    temp1 = 0.5 * (exp(2 * time / domtau) - 1) * domtau;

    c->Y_dx = p->Y_dx +
              p->rhodx * p->sigx * domsig * (domPsi / sqdomF) * temp1 +
              p->rhodx * p->sigx * domsig * domtau * (temp1 - temp);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

/* --------------------------------------------------------------------------------------
 */

static Err make_V_dx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double domsig, domF, domtau;
  double temp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;

  c->V_dx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    domtau = p->StochRates_Ts[0].tau;
    domF = p->StochRates_Ts[0].F;
    temp = (exp(time / domtau) - 1) * domtau;

    c->V_dx = p->V_dx + p->rhodx * p->sigx * (domsig / domF) * temp;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_W_dx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double domsig, domtau, domF, domPsi;
  double temp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;

  c->W_dx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    domsig = p->StochRates_Ts[0].sig;
    domtau = p->StochRates_Ts[0].tau;
    domF = p->StochRates_Ts[0].F;
    domPsi = p->StochRates_Ts[0].Psi;
    temp = (exp(time / domtau) - 1) * domtau;

    c->W_dx = p->W_dx + p->rhodx * p->sigx * domsig * (domPsi / domF) * temp +
              p->rhodx * p->sigx * domsig * domtau * (temp - time);
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_V_fx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double forsig, forF, fortau;
  double temp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;

  c->V_fx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    forsig = p->StochRates_Ts[1].sig;
    fortau = p->StochRates_Ts[1].tau;
    forF = p->StochRates_Ts[1].F;
    temp = (exp(time / fortau) - 1) * fortau;

    c->V_fx = p->V_fx + p->rhofx * p->sigx * (forsig / forF) * temp;
  }
  return err;
}

/* --------------------------------------------------------------------------------------
 */

static Err make_W_fx_vector(TermStruct *l) {

  FxTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double forsig, fortau, forF, forPsi;
  double temp;
  Err err = 0;

  lc = l->head;
  c = (FxTermStructVal *)lc->element->val.pval;

  c->W_fx = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (FxTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (FxTermStructVal *)lc->element->val.pval;

    forsig = p->StochRates_Ts[1].sig;
    fortau = p->StochRates_Ts[1].tau;
    forF = p->StochRates_Ts[1].F;
    forPsi = p->StochRates_Ts[1].Psi;
    temp = (exp(time / fortau) - 1) * fortau;

    c->W_fx = p->W_fx + p->rhofx * p->sigx * forsig * (forPsi / forF) * temp +
              p->rhofx * p->sigx * forsig * fortau * (temp - time);
  }
  return err;
}
