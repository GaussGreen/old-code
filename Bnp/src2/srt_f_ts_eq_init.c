/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_eq_init.c

   PURPOSE:     Function to initialise an EQ TermStruct:
   ------------------------------------------------------------------------ */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts_eq.h"

static Err make_eq_cumvol_vector(TermStruct *l);
static Err make_V_ir_e_vector(TermStruct *l);
static Err make_W_ir_e_vector(TermStruct *l);

/* ----------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------
   The MAIN function to call when a TS has to be initialised for an EQD
   underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_init_EQ_TermStruct(Date today, double **sig, /* sig[col][row] */
                             int sig_cols, int nsig,

                             char *undname, char *ir_or_yc_name,
                             SrtMdlType mdl_type,

                             /* SRVGS Model Parameter */
                             double omega, double beta, double gamma,

                             double voldrift, double vovol, double rho,

                             /* Output */
                             TermStruct **ts) {
  Err err;

  if ((mdl_type == BLACK_SCHOLES) || (mdl_type == NORMAL_BS) ||
      (mdl_type == FX_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES) ||
      (mdl_type == EQ_STOCH_RATES_SRVGS)) {
    err = srt_f_init_EQ_NORMLOG_TermStruct(undname, ir_or_yc_name, mdl_type, ts,
                                           today, sig, sig_cols, nsig,

                                           /* SRVGS Model Parameter */
                                           omega, beta, gamma,

                                           voldrift, vovol, rho);
  }
  return err;
}

/* ------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   A function to output an existing TermStructure (Memory Allocation done
   inside)
   ----------------------------------------------------------------------- */

Err srt_f_display_EQ_TermStruct(TermStruct *ts, double **sigma_date,
                                double **sigma, long *sigma_n)

{
  SrtLst *l = ts->head;

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
        (double)(((EquityTermStructVal *)l->element->val.pval)->date);
    (*sigma)[*sigma_n - 1] = ((EquityTermStructVal *)l->element->val.pval)->sig;

    l = l->next;
  }

  return NULL;
}

/* ----------------------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------------
   The MAIN function to free a Term Structure of Volatility initialised
   for an EQD underlying
   -------------------------------------------------------------------------------
 */

Err srt_f_free_EQ_TermStruct(TermStruct **l) {
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
Err srt_f_EquityTermStructValFree(void *tsvalptr) {
  EquityTermStructVal *tsval = (EquityTermStructVal *)tsvalptr;
  Err err = NULL;

  srt_free(tsval);

  return err;
}

/* ------------------------------------------------------------------------- */

/* Initialises the Volatility Term Structure for Eq Models */

Err srt_f_init_EQ_NORMLOG_TermStruct(char *undname, char *ir_und_or_yc_name,
                                     SrtMdlType mdl_type, TermStruct **ts,
                                     Date today, double **sig_data,
                                     int sig_cols, int num_sig,

                                     double omega, double beta, double gamma,

                                     double voldrift, double vovol,
                                     double rho) {
  int cur_sig;
  EquityTermStructVal *tsval, *tsval_p;
  IrmTermStructVal *ir_tsval;

  Err err = 0;
  long i;
  long ticker;
  SrtUndPtr ir_und;
  TermStruct *ir_ts;
  SrtListAtom *atom, *ir_atom;
  SrtCorrLstPtr corrlist;
  SrtMdlType ir_mdl_type;
  double time;

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
    return serror("init_TermStruct err: no current data !!!");

  for (; cur_sig < num_sig; cur_sig++) {

    /* Creation of a new TermStruct2 */
    tsval = (EquityTermStructVal *)srt_calloc(1, sizeof(EquityTermStructVal));

    /* transfer of information */
    tsval->date = DTOL(sig_data[0][cur_sig]);
    tsval->time = (sig_data[0][cur_sig] - today) * YEARS_IN_DAY;
    tsval->sig = sig_data[1][cur_sig];

    /* SRVGS model parameters */
    tsval->omega = omega;
    tsval->beta = beta;
    tsval->gamma = gamma;

    tsval->basevol = sig_data[1][cur_sig];
    tsval->voldrift = voldrift;
    tsval->vovol = vovol;

    tsval->rho_spot_vol = rho;

    /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
    srt_f_lstins(*ts, "BsTsAtom", tsval->date, OBJ_PTR_EQ_TermStruct,
                 (void *)tsval, &srt_f_EquityTermStructValFree, &ticker);
  }

  if ((mdl_type == EQ_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES_SRVGS)) {
    /* gets the interest rate underlying name && its term structure */
    ir_und = lookup_und(ir_und_or_yc_name);
    if (!ir_und)
      return serror("can not get interest rate underlying");

    err = get_underlying_mdltype(ir_und, &ir_mdl_type);
    if (err)
      return err;

    if (ir_mdl_type == LGM) {
      err = get_underlying_ts(ir_und, &ir_ts);
      if (err)
        return serror("can not get the interest rate term structure");

      /* starts at the top of the lgm term struct */
      ir_atom = ir_ts->head;
      /* starts at the top of the equity term sruct */
      atom = (*ts)->head;

      /* loops on the interest rate term structure dates to add them to the
       * equity one */
      while (ir_atom != NULL) {
        /* get the lgm term struct val attached */
        ir_tsval = (IrmTermStructVal *)ir_atom->element->val.pval;

        /* moves on the equity date just after or on this date */
        while ((atom != NULL) &&
               (((EquityTermStructVal *)atom->element->val.pval)->time <
                ir_tsval->time))
          atom = atom->next;
        if (atom == NULL)
          atom = (*ts)->tail;

        tsval = ((EquityTermStructVal *)atom->element->val.pval);
        /* if not on the same date   , add the date on the equity term struct */
        if (tsval->time != ir_tsval->time) {
          tsval =
              (EquityTermStructVal *)srt_calloc(1, sizeof(EquityTermStructVal));
          /* adds the equity information */
          tsval->date = ir_tsval->date;
          tsval->time = ir_tsval->time;
          tsval->sig = find_eq_sig(tsval->time, *ts);

          tsval->beta = find_eq_beta(tsval->time, *ts);
          tsval->gamma = find_eq_gamma(tsval->time, *ts);
          tsval->omega = find_eq_omega(tsval->time, *ts);
          tsval->basevol = find_eq_basevol(tsval->time, *ts);
          tsval->voldrift = find_eq_voldrift(tsval->time, *ts);
          tsval->vovol = find_eq_vovol(tsval->time, *ts);
          tsval->rho_spot_vol = find_eq_rho(tsval->time, *ts);

          /* Insert the TermStructAtom in the TS (==linked ts sorted by key =
           * date)  */
          srt_f_lstins(*ts, "EQTsAtom", tsval->date, OBJ_PTR_EQ_TermStruct,
                       (void *)tsval, &srt_f_free_EQ_TermStruct, &ticker);
        }

        /* moves on the next element of the ir ts */
        ir_atom = ir_atom->next;
      } /* end of loop on the ir ts */

      /* gets the global correlation list */
      corrlist = srt_f_GetTheCorrelationList();
      if (!corrlist->head->element)
        return serror("Correlation list improperly initialised...");

      /* loops on all the dates of the EQ term stucture */
      atom = (*ts)->head;
      while (atom != NULL) {
        tsval = (EquityTermStructVal *)atom->element->val.pval;
        time = 0;
        if (atom != (*ts)->head) {
          tsval_p = (EquityTermStructVal *)atom->previous->element->val.pval;
          time = tsval_p->time;
        }

        err = srt_f_get_corr_from_CorrList(corrlist, ir_und_or_yc_name, undname,
                                           tsval->time, &tsval->rho_e_ir);
        if (err) {
          srt_f_lstfree(*ts, SRT_YES);
          return err;
        }

        tsval->StochRate_Ts.sig = find_sig(tsval->time, ir_ts);
        tsval->StochRate_Ts.tau = find_tau(tsval->time, ir_ts);
        tsval->StochRate_Ts.F = F_func(time, ir_ts);
        tsval->StochRate_Ts.Psi = Psi_func(time, ir_ts);

        atom = atom->next;
      } /* end of loop */
    }
  } /* end of if(ir_mdl_type == LGM) */

  /* Computes cum_vol values */
  err = make_eq_cumvol_vector(*ts);
  if (err) {
    srt_f_lstfree(*ts, SRT_YES);
    return err;
  }

  err = make_V_ir_e_vector(*ts);
  if (err) {
    srt_f_lstfree(*ts, SRT_YES);
    return err;
  }

  err = make_W_ir_e_vector(*ts);
  if (err) {
    srt_f_lstfree(*ts, SRT_YES);
    return err;
  }

  return NULL;
}

/* ------------------------------------------------------------------------- */

/* Computes both integral of sigma and sigm^2 */

static Err make_eq_cumvol_vector(TermStruct *l) {
  EquityTermStructVal *c, *p;
  SrtLst *lc;

  Err err = 0;

  /* Starts at the beginning of the Term Structure */
  lc = l->head;

  /* Gets the TermStructVal attached to this element */
  c = (EquityTermStructVal *)lc->element->val.pval;

  /* Computes the value of the cumulative volatility */
  c->int_sig2_dt = c->sig * c->sig * c->time;
  c->int_sig_dt = c->sig * c->time;

  /* Moves on to the next point */
  lc = lc->next;

  /* Loops on all the elements of the TermStructure */
  while (lc != NULL) {
    p = c;
    c = (EquityTermStructVal *)lc->element->val.pval;

    /* The Cumulative Volatility as the integral of Sigma^2*t */
    c->int_sig2_dt = p->int_sig2_dt + c->sig * c->sig * (c->time - p->time);

    /* The Integral of Sigma * t */
    c->int_sig_dt = p->int_sig_dt + c->sig * (c->time - p->time);

    /* Moves on to the next element in the Term Structure */
    lc = lc->next;
  }
  return err;
}

/* ------------------------------------------------------------------------- */

static Err make_V_ir_e_vector(TermStruct *l) {
  EquityTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double tau, ir_sig, F;

  lc = l->head;
  c = (EquityTermStructVal *)lc->element->val.pval;
  c->V_ir_e = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (EquityTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (EquityTermStructVal *)lc->element->val.pval;

    tau = p->StochRate_Ts.tau;
    ir_sig = p->StochRate_Ts.sig;
    F = p->StochRate_Ts.F;

    c->V_ir_e = p->V_ir_e + p->rho_e_ir * p->sig * (ir_sig / F) *
                                (exp(time / tau) - 1) * tau;
  }

  return NULL;
}

static Err make_W_ir_e_vector(TermStruct *l) {
  EquityTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double tau, ir_sig, F, Psi;

  lc = l->head;
  c = (EquityTermStructVal *)lc->element->val.pval;
  c->W_ir_e = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (EquityTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (EquityTermStructVal *)lc->element->val.pval;

    ir_sig = p->StochRate_Ts.sig;
    tau = p->StochRate_Ts.tau;
    F = p->StochRate_Ts.F;
    Psi = p->StochRate_Ts.Psi;

    c->W_ir_e = p->W_ir_e +
                p->rho_e_ir * p->sig * ir_sig * (Psi / F) *
                    (exp(time / tau) - 1) * tau +
                p->rho_e_ir * p->sig * ir_sig * tau *
                    ((exp(time / tau) - 1) * tau - time);
  }

  return NULL;
}
