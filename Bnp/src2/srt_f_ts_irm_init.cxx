/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_init.cxx

   PURPOSE:     Provide a few utilities when dealing with a TermStruct:
                    - initialise a new TS:      srt_f_init_IRM_TermStruct
                    - output an existing TS:    srt_f_display_ts
   ------------------------------------------------------------------------ */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts_irm.h"

/* ----------------------------------------------------------------------- */

static Err make_new_G_H_vector(TermStruct *l);
static Err make_I_vector(TermStruct *l);
static Err make_J_vector(TermStruct *l);
static Err make_K_vector(TermStruct *l);
static Err make_L_vector(TermStruct *l);
static Err make_O_vector(TermStruct *l);
static Err make_Q_vector(TermStruct *l);
static Err make_R_vector(TermStruct *l);
static Err make_S_vector(TermStruct *l);
static Err make_T_vector(TermStruct *l);
static Err make_Phi_vector(TermStruct *l);
static Err make_Zeta_vector(TermStruct *l);
static Err make_vasicek_mean_sr_vector(TermStruct *ts);
static Err make_vasicek_mean_int_sr_vector(TermStruct *ts);
static Err make_int_phi_vector(TermStruct *ts);

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   The MAIN function to call when a Interest Rate Model TermStructure has
   to be initialised
   ------------------------------------------------------------------------- */

Err srt_f_init_IRM_TermStruct(Date today, double **sig, /* sig[col][row] */
                              int sig_cols, int nsig,
                              double **tau, /* tau[col][row] */
                              int tau_cols, int ntau,

                              SrtMdlType mdl_type, SrtMdlDim mdl_dim,

                              /* SMILE */
                              double beta,

                              /* TWO FACTOR */
                              double alpha, double gamma, double rho,

                              /* STOCH VOL */
                              double vovol,

                              /* ETABETA  or MIXED */
                              double eta_or_omega,

                              /* VASICEK MODEL */
                              double vasicek_init_cond, int num_mean_rev_level,
                              int num_mean_rev_level_cols,
                              double **vasicek_mean_rev_level_data,

                              /* OUTPUT */
                              TermStruct **ts) {
  Err flag;

  if (mdl_dim == TWO_FAC) {
    flag = srt_f_init_IRM_TwoFac_TermStruct(ts, today, sig, sig_cols, nsig, tau,
                                            tau_cols, ntau, mdl_type, beta,
                                            alpha, gamma, rho, eta_or_omega);
  } else if (mdl_dim == ONE_FAC) {
    flag = srt_f_init_IRM_OneFac_TermStruct(
        ts, today, sig, sig_cols, nsig, tau, tau_cols, ntau, mdl_type, alpha,
        rho, gamma, beta, vasicek_init_cond, num_mean_rev_level,
        num_mean_rev_level_cols, vasicek_mean_rev_level_data);
  }

  return flag;

} /* END Err srt_f_init_IRM_TermStruct(...) */

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   The MAIN function to call to FREE an Interest Rate Model Term Structure
   ------------------------------------------------------------------------- */

Err srt_f_free_IRM_TermStruct(TermStruct **l) {
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

/* -------------------------------------------------------------------------
                  ONE FACTOR SPECIFIC FUNCTIONS
   ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */

/* Function required to free a TermStructObject once attached into a linked list
 */
Err srt_f_irmtsvalfree(void *tsvalptr) {
  Err err = NULL;

  IrmTermStructVal *tsval;

  tsval = (IrmTermStructVal *)tsvalptr;

  if (tsval->is_M_allocated_here == SRT_YES) {
    if (tsval->M_beta_eta)
      free_dmatrix(tsval->M_beta_eta, 0, SMAX - 1, 0, THETAMAX - 1);
  }
  tsval->M_beta_eta = NULL;

  srt_free(tsval);

  return err;
}

/* ------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------
   The MAIN function to call when a Interest Rate Model TermStructure has
   to be initialised with a ONE FACTOR MODEL
   ------------------------------------------------------------------------- */

Err srt_f_init_IRM_OneFac_TermStruct(
    TermStruct **ts, Date today, double **sig_data, int sig_cols, int num_sig,
    double **tau_data, int tau_cols, int num_tau, SrtMdlType mdl_type,
    double vovol, double rho, double meanvol, double beta,
    double vasicek_init_cond, int num_mean_rev_level,
    int num_mean_rev_level_cols, double **vasicek_mean_rev_level_data) {
  int cur_sig, cur_tau, cur_mean_rev_level;
  IrmTermStructVal *tsval;
  Err err = 0;
  int over_sig, over_tau;
  long i, j;
  double **ppdSigmaValues;
  double **ppdTauValues;
  long ticker;
  SrtLst *ls;

  /* Checks on the number of columns in Sigma : one value        , two columns
   * or three columns (with Beta) */
  if ((sig_cols != 2) && (sig_cols != 3) && (sig_cols != 6)) {
    if ((sig_cols != 1) || (num_sig > 1))
      return serror(
          "Need TWO/THREE/SIX columns or ONE SINGLE VALUE for the Sigma Curve");
  }

  /* Checks on the number of columns in Tau : one value or two columns */
  if (tau_cols != 2) {
    if ((tau_cols != 1) || (num_tau > 1))
      return serror("Need TWO columns or ONE SINGLE VALUE for the Tau Curve");
  }

  /* Allocate memory for the real size matrix */
  ppdSigmaValues = dmatrix(0, 5, 0, num_sig - 1);
  ppdTauValues = dmatrix(0, 1, 0, num_tau - 1);

  /* If only one single value for sigma        , create a fake Term Struct with
   * constants */
  if ((num_sig == 1) && (sig_cols == 1)) {
    ppdSigmaValues[0][0] = today + 365.0;
    ppdSigmaValues[1][0] = sig_data[0][0];
  } else
  /* Transfer the sigma data from the input array */
  {
    for (i = 0; i < num_sig; i++) {
      ppdSigmaValues[0][i] = sig_data[0][i];
      ppdSigmaValues[1][i] = sig_data[1][i];
    }
  }

  /* If the third column (the Beta) is missing        , fill it in with the
   * input beta or according to model type */
  if (sig_cols < 3) {
    for (i = 0; i < num_sig; i++) {
      ppdSigmaValues[2][i] = beta;
      ppdSigmaValues[3][i] = vovol;
      ppdSigmaValues[4][i] = meanvol;
      ppdSigmaValues[5][i] = rho;
    }
  } else if (sig_cols == 3) {
    /* Transfer the input data */
    for (i = 0; i < num_sig; i++) {
      ppdSigmaValues[2][i] = sig_data[2][i];
      ppdSigmaValues[3][i] = vovol;
      ppdSigmaValues[4][i] = meanvol;
      ppdSigmaValues[5][i] = rho;
    }
  } else {
    /* Transfer the input data SV */
    for (i = 0; i < num_sig; i++) {
      ppdSigmaValues[2][i] = sig_data[2][i];
      ppdSigmaValues[3][i] = sig_data[3][i];
      ppdSigmaValues[4][i] = sig_data[4][i];
      ppdSigmaValues[5][i] = sig_data[5][i];
    }
  }

  /* Overwrites the Values of Beta for explicit models */
  if ((mdl_type == LGM) || (mdl_type == NEWLGM) ||
      (mdl_type == LGM_STOCH_VOL)) {
    for (i = 0; i < num_sig; i++)
      ppdSigmaValues[2][i] = 0.0;
  } else if ((mdl_type == CHEY) || (mdl_type == CHEY_STOCH_VOL)) {
    for (i = 0; i < num_sig; i++)
      ppdSigmaValues[2][i] = 1.0;
  }

  /* If only one single value for tau        , create a fake Term Struct with
   * constants */
  if ((num_tau == 1) && (tau_cols == 1)) {
    ppdTauValues[0][0] = today + 365;
    ppdTauValues[1][0] = tau_data[0][0];
  } else
  /* Transfer the tau data from the input array */
  {
    for (i = 0; i < num_tau; i++) {
      ppdTauValues[0][i] = tau_data[0][i];
      ppdTauValues[1][i] = tau_data[1][i];
    }
  }

  /* Starts the extraction of the data from the inputs */
  cur_tau = 0;
  cur_sig = 0;
  cur_mean_rev_level = 0;

  over_sig = 0;
  over_tau = 0;

  /* If the data exists but is NULL: consider what follows as garbage */
  for (i = 0; i < num_sig; i++) {
    if (!ppdSigmaValues[0][i] || !ppdSigmaValues[1][i])
      break;
  }
  num_sig = i;
  for (i = 0; i < num_tau; i++) {
    if (!ppdTauValues[0][i] || !ppdTauValues[1][i])
      break;
  }
  num_tau = i;
  for (i = 0; i < num_mean_rev_level; i++) {
    if (!vasicek_mean_rev_level_data[0][i] ||
        !vasicek_mean_rev_level_data[1][i])
      break;
  }
  num_mean_rev_level = i;

  /* Current index for tau and sigma start at a date >= today */
  while ((cur_sig < num_sig) &&
         ((ppdSigmaValues[0][cur_sig] - today) <= EPSILON))
    cur_sig++;
  while ((cur_tau < num_tau) && ((ppdTauValues[0][cur_tau] - today) <= EPSILON))
    cur_tau++;
  while (
      (cur_mean_rev_level < num_mean_rev_level) &&
      ((vasicek_mean_rev_level_data[0][cur_mean_rev_level] - today) <= EPSILON))
    cur_mean_rev_level++;

  /* Return error if there is no acceptable sigma or tau value */
  if (cur_sig == num_sig) {
    free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
    free_dmatrix(ppdTauValues, 0, 1, 0, num_tau - 1);
    return serror("No Current Sigma Values in TS");
  }
  if (cur_tau == num_tau) {
    free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
    free_dmatrix(ppdTauValues, 0, 1, 0, num_tau - 1);
    return serror("No Current Tau Values in TS");
  }

  if ((cur_mean_rev_level == num_mean_rev_level) && (num_mean_rev_level > 0)) {
    free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
    free_dmatrix(ppdTauValues, 0, 1, 0, num_tau - 1);
    return serror("no current mean reversion level values in term struct");
  }

  /* Make sure there is no negative values of ppdSigmaValues */
  for (i = 0; i < num_sig; i++) {
    if (ppdSigmaValues[1][i] < 0.0) {
      free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
      free_dmatrix(ppdTauValues, 0, 1, 0, num_tau - 1);
      return serror("Cannot input negative volatility in an IRM Term Struct");
    }
  }

  /* Create a blank double linked ts that will become the TermStruct */
  srt_f_lstcreate(ts, "TermStruct");

  /* Loops to fill in the TS as long as there is something to put in */
  while ((over_tau == 0) || (over_sig == 0)) {
    /* Creation of a new IrmTermStructVal (simple element of the Term Struct) */
    tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

    /* Transfer of initial information */
    tsval->sig = ppdSigmaValues[1][cur_sig];
    tsval->tau = ppdTauValues[1][cur_tau];
    tsval->beta = ppdSigmaValues[2][cur_sig];
    tsval->vovol = ppdSigmaValues[3][cur_sig];
    tsval->meanvol = ppdSigmaValues[4][cur_sig];
    tsval->rho = ppdSigmaValues[5][cur_sig];

    /* 1st case: date_sigma = date_tau */
    if (((over_sig == 0) && (over_tau == 0)) &&
        (ppdSigmaValues[0][cur_sig] == ppdTauValues[0][cur_tau])) {

      /* Transfer of specific information */
      tsval->time = (ppdSigmaValues[0][cur_sig] - today) * YEARS_IN_DAY;
      tsval->date = DTOL(ppdSigmaValues[0][cur_sig]);
      tsval->val_origin = BOTH_DATE;

      /* Update the indexes */
      if (cur_sig < num_sig - 1)
        cur_sig++;
      else
        over_sig = 1;

      if (cur_tau < num_tau - 1)
        cur_tau++;
      else
        over_tau = 1;

    } else
        /* 2nd case: date_sigma < date_tau */
        if ((over_tau == 1) ||
            (((over_sig == 0) && (over_tau == 0)) &&
             (ppdSigmaValues[0][cur_sig] < ppdTauValues[0][cur_tau]))) {

      /* Transfer of specific information */
      tsval->date = DTOL(ppdSigmaValues[0][cur_sig]);
      tsval->time = (ppdSigmaValues[0][cur_sig] - today) * YEARS_IN_DAY;
      tsval->val_origin = SIGMA_DATE;

      /* Update the indexes */
      if (cur_sig < num_sig - 1)
        cur_sig++;
      else
        over_sig = 1;

    } else
        /* 3rd case: date_tau < date_sigma */
        if ((over_sig == 1) ||
            (((over_sig == 0) && (over_tau == 0)) &&
             (ppdSigmaValues[0][cur_sig] > ppdTauValues[0][cur_tau]))) {
      /* Transfer of specific information */
      tsval->date = DTOL(ppdTauValues[0][cur_tau]);
      tsval->time = (ppdTauValues[0][cur_tau] - today) * YEARS_IN_DAY;
      tsval->val_origin = TAU_DATE;

      /* Update the indexes */
      if (cur_tau < num_tau - 1)
        cur_tau++;
      else
        over_tau = 1;
    }

    /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
    srt_f_lstins(*ts, "OneFacTsAtom", tsval->date, OBJ_PTR_IRM_TermStruct,
                 (void *)tsval, &srt_f_irmtsvalfree, &ticker);
  }

  if (mdl_type == VASICEK) {
    /* loop on the mean reversion level dates to add them to the interest
     * underlying one */

    while (cur_mean_rev_level < num_mean_rev_level) {
      /* moves on to the i.r dates just after or on this date */
      ls = (*ts)->head;
      tsval = (IrmTermStructVal *)ls->element->val.pval;
      while (
          (ls != NULL) &&
          (tsval->date < vasicek_mean_rev_level_data[0][cur_mean_rev_level])) {
        ls = ls->next;
        if (ls)
          tsval = (IrmTermStructVal *)ls->element->val.pval;
      }
      if (ls == NULL)
        ls = (*ts)->tail;

      tsval = (IrmTermStructVal *)ls->element->val.pval;
      /* if not on the same date        , we add the date onto the i.r term
       * structure
       */
      if (tsval->date != vasicek_mean_rev_level_data[0][cur_mean_rev_level]) {
        /* creation of a new IrmTermStructVal (simple element of the Term
         * Struct) */
        tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

        tsval->val_origin = VASICEK_MEAN_REV_DATE;

        /* adds the specific information */
        tsval->date = vasicek_mean_rev_level_data[0][cur_mean_rev_level];
        tsval->time =
            (vasicek_mean_rev_level_data[0][cur_mean_rev_level] - today) *
            YEARS_IN_DAY;

        tsval->sig = find_sig(tsval->time, *ts);
        tsval->tau = find_tau(tsval->time, *ts);

        tsval->vasicek_init_cond = vasicek_init_cond;
        tsval->mean_rev_level =
            vasicek_mean_rev_level_data[1][cur_mean_rev_level];

        /* insert the term struct atom in the ts */
        srt_f_lstins(*ts, "OneFacTsAtom", tsval->date, OBJ_PTR_IRM_TermStruct,
                     (void *)tsval, &srt_f_irmtsvalfree, &ticker);

      } else /* If it tsval->date =
                vasicek_mean_rev_level_data[0][cur_mean_rev_level] */
      {
        /* adds the specific information */
        tsval->vasicek_init_cond = vasicek_init_cond;
        tsval->mean_rev_level =
            vasicek_mean_rev_level_data[1][cur_mean_rev_level];
      }

      cur_mean_rev_level += 1;
    }

    /* loops on all the dates of the i.r term struct to add ..*/
    ls = (*ts)->head;
    while (ls != NULL) {
      tsval = (IrmTermStructVal *)ls->element->val.pval;

      tsval->sig = find_sig(tsval->time, *ts);
      tsval->tau = find_tau(tsval->time, *ts);

      tsval->vasicek_init_cond = vasicek_init_cond;

      j = 0;
      while ((j < num_mean_rev_level) &&
             (tsval->date > vasicek_mean_rev_level_data[0][j]))
        j++;
      if (j == num_mean_rev_level)
        tsval->mean_rev_level =
            vasicek_mean_rev_level_data[1][num_mean_rev_level - 1];
      else
        tsval->mean_rev_level = vasicek_mean_rev_level_data[1][j];

      ls = ls->next;
    }
  }

  /* Free the allocated memory */
  free_dmatrix(ppdSigmaValues, 0, 5, 0, num_sig - 1);
  free_dmatrix(ppdTauValues, 0, 1, 0, num_tau - 1);

  /* computes F and Psi values for LGM        , Cheyette...*/
  err = make_F_Psi_vector(*ts);
  if (err)
    return err;

  /* computes G and H values for NEWLGM */
  if (is_model_New_type(mdl_type)) {
    err = make_new_G_H_vector(*ts);
    if (err)
      return err;
  }

  /* computes G and H values used in LGM*/
  if (mdl_type == LGM || mdl_type == LGM_STOCH_VOL) {
    err = make_G_H_vector(*ts);
    if (err)
      return err;
    err = make_I_vector(*ts);
    if (err)
      return err;
    err = make_J_vector(*ts);
    if (err)
      return err;
    err = make_K_vector(*ts);
    if (err)
      return err;
    err = make_L_vector(*ts);
    if (err)
      return err;
    err = make_O_vector(*ts);
    if (err)
      return err;
    err = make_Q_vector(*ts);
    if (err)
      return err;
    err = make_R_vector(*ts);
    if (err)
      return err;
    err = make_S_vector(*ts);
    if (err)
      return err;

    err = make_T_vector(*ts);
    if (err)
      return err;

    err = make_Phi_vector(*ts);
    if (err)
      return err;

    err = make_int_phi_vector(*ts);
    if (err)
      return err;

  } else if (mdl_type == ETABETA) {
    err = make_G_H_vector(*ts);
    if (err)
      return err;
    err = make_Zeta_vector(*ts);
    if (err)
      return err;
    err = make_M_matrix(*ts);
    if (err)
      return err;
    err = make_Lambda_vector(*ts);
    if (err)
      return err;
  } else if (mdl_type == VASICEK) {
    err = make_G_H_vector(*ts);
    if (err)
      return err;

    err = make_L_vector(*ts);
    if (err)
      return err;

    err = make_O_vector(*ts);
    if (err)
      return err;

    err = make_S_vector(*ts);
    if (err)
      return err;

    err = make_T_vector(*ts);
    if (err)
      return err;

    err = make_vasicek_mean_sr_vector(*ts);
    if (err)
      return err;

    err = make_vasicek_mean_int_sr_vector(*ts);
    if (err)
      return err;
  }

  return err;

} /* END Err srt_f_init_IRM_OneFac_TermStruct(...) */

/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
   A function to output an existing TermStructure (Memory Allocation done
   inside)
   ----------------------------------------------------------------------- */
Err srt_f_display_IRM_OneFac_TermStruct(TermStruct *ts, double **sigma_date,
                                        double **sigma, double **beta,
                                        double **vovol, double **rho,
                                        double **meanvol, long *plNumSigmas,
                                        double **tau_date, double **tau,
                                        long *plNumTaus)

{
  SrtLst *l = ts->head;

  *plNumSigmas = *plNumTaus = 0;

  /* Allocate memory for the initial pointers */
  *sigma_date = (double *)srt_malloc(sizeof(double));
  *tau_date = (double *)srt_malloc(sizeof(double));

  *sigma = (double *)srt_malloc(sizeof(double));
  *beta = (double *)srt_malloc(sizeof(double));
  *vovol = (double *)srt_malloc(sizeof(double));
  *rho = (double *)srt_malloc(sizeof(double));
  *meanvol = (double *)srt_malloc(sizeof(double));
  *tau = (double *)srt_malloc(sizeof(double));

  /* Loop on all the elements of the Term Sturcture */
  while (l != NULL) {
    if ((((IrmTermStructVal *)l->element->val.pval)->val_origin == TAU_DATE) ||
        (((IrmTermStructVal *)l->element->val.pval)->val_origin == BOTH_DATE)) {

      (*plNumTaus)++;
      *tau_date = realloc(*tau_date, (*plNumTaus) * sizeof(double));
      *tau = realloc(*tau, (*plNumTaus) * sizeof(double));
      (*tau_date)[*plNumTaus - 1] =
          (double)(((IrmTermStructVal *)l->element->val.pval)->date);
      (*tau)[*plNumTaus - 1] = ((IrmTermStructVal *)l->element->val.pval)->tau;
    }

    if ((((IrmTermStructVal *)l->element->val.pval)->val_origin ==
         SIGMA_DATE) ||
        (((IrmTermStructVal *)l->element->val.pval)->val_origin == BOTH_DATE)) {

      (*plNumSigmas)++;
      *sigma_date = realloc(*sigma_date, (*plNumSigmas) * sizeof(double));
      *sigma = realloc(*sigma, (*plNumSigmas) * sizeof(double));
      *beta = realloc(*beta, (*plNumSigmas) * sizeof(double));
      *vovol = realloc(*vovol, (*plNumSigmas) * sizeof(double));
      *rho = realloc(*rho, (*plNumSigmas) * sizeof(double));
      *meanvol = realloc(*meanvol, (*plNumSigmas) * sizeof(double));
      (*sigma_date)[*plNumSigmas - 1] =
          (double)(((IrmTermStructVal *)l->element->val.pval)->date);
      (*sigma)[*plNumSigmas - 1] =
          ((IrmTermStructVal *)l->element->val.pval)->sig;
      (*beta)[*plNumSigmas - 1] =
          ((IrmTermStructVal *)l->element->val.pval)->beta;
      (*vovol)[*plNumSigmas - 1] =
          ((IrmTermStructVal *)l->element->val.pval)->vovol;
      (*rho)[*plNumSigmas - 1] =
          ((IrmTermStructVal *)l->element->val.pval)->rho;
      (*meanvol)[*plNumSigmas - 1] =
          ((IrmTermStructVal *)l->element->val.pval)->meanvol;
    }

    l = l->next;
  }

  return NULL;
}

/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------

                FUNCTIONS MORE SPECIFIC TO THE HANDLING OF THE TERM STRUCT

   ------------------------------------------------------------------------------
 */

Err srt_f_tsupdate(TermStruct **l) {
  Err err = 0;
  double cur_tau, cur_sig;
  int done_sig, done_tau;
  SrtLst *ls;
  IrmTermStructVal *tsval;

  ls = (*l)->tail;
  done_sig = 0;
  done_tau = 0;

  /* We search for the 2 last original values of sigma and tau */
  /* Search completed from the end of the ts */
  while (((done_sig != 1) || (done_tau != 1)) && (ls != NULL)) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;
    switch (tsval->val_origin) {
    case SIGMA_DATE:
      cur_sig = tsval->sig;
      done_sig = 1;
      break;
    case TAU_DATE:
      cur_tau = tsval->tau;
      done_tau = 1;
      break;
    case BOTH_DATE:
      cur_tau = tsval->tau;
      done_tau = 1;
      cur_sig = tsval->sig;
      done_sig = 1;
      break;
    default:
      break;
    }
    ls = ls->previous;
  }

  ls = (*l)->tail;
  while (ls != NULL) {
    tsval = (IrmTermStructVal *)ls->element->val.pval;

    switch (tsval->val_origin) {
    case SIGMA_DATE:
      cur_sig = tsval->sig;
      tsval->tau = cur_tau;
      break;
    case TAU_DATE:
      cur_tau = tsval->tau;
      tsval->sig = cur_sig;
      break;
    case BOTH_DATE:
      cur_tau = tsval->tau;
      cur_sig = tsval->sig;
      break;
    default:
      break;
    }

    ls = ls->previous;
  }

  /* recomputes F and Psi values */
  make_F_Psi_vector(*l);

  /* recomputes G and H values */
  err = make_G_H_vector(*l);
  make_J_vector(*l);

  return err;
}

/* ------------------------------------------------------------------------- */

Err srt_f_tsaddtau(TermStruct *ts, Ddate date, double today, double tau) {
  IrmTermStructVal *tsval;
  SrtObject *objPtr;
  Err err = 0;
  long ticker;

  if (date <= today)
    return "Tau date < today";

  if (srt_f_lstgetobj(*ts, "", date, &objPtr) == NULL) {
    tsval = (IrmTermStructVal *)objPtr->val.pval;
    tsval->tau = tau;

    if (tsval->val_origin == SIGMA_DATE)
      tsval->val_origin = BOTH_DATE;
    else
      tsval->val_origin = TAU_DATE;

    return err;
  }

  tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

  /* transfer of information */
  tsval->date = DTOL(date);
  tsval->time = (date - today) * YEARS_IN_DAY;
  tsval->tau = tau;
  tsval->val_origin = TAU_DATE;

  /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
  srt_f_lstins(ts, "OneFacTsAtom", tsval->date, OBJ_PTR_IRM_TermStruct,
               (void *)tsval, &srt_f_irmtsvalfree, &ticker);

  /* !!! The ts is not accurate anymore !!! */

  return err;
}

/* ------------------------------------------------------------------------- */

Err srt_f_tsaddsig(TermStruct *ts, Ddate date, double today, double sig) {
  IrmTermStructVal *tsval;
  SrtObject *objPtr;
  Err err = 0;
  long ticker;

  if (date <= today)
    return "Sigma date < today";

  if (srt_f_lstgetobj(*ts, "", date, &objPtr) == NULL) {
    tsval = (IrmTermStructVal *)objPtr->val.pval;
    tsval->sig = sig;
    if (tsval->val_origin == TAU_DATE)
      tsval->val_origin = BOTH_DATE;
    else
      tsval->val_origin = SIGMA_DATE;

    return err;
  }

  tsval = (IrmTermStructVal *)srt_calloc(1, sizeof(IrmTermStructVal));

  /* transfer of information */
  tsval->date = DTOL(date);
  tsval->time = (date - today) * YEARS_IN_DAY;
  tsval->sig = sig;
  tsval->val_origin = SIGMA_DATE;

  /* Insert the TermStructAtom in the TS (==linked ts sorted by key = date) */
  srt_f_lstins(ts, "OneFacTsAtom", tsval->date, OBJ_PTR_IRM_TermStruct,
               (void *)tsval, &srt_f_irmtsvalfree, &ticker);

  /* !!! The ts is not accurate anymore !!! */

  return err;
}

Err make_vasicek_mean_sr_vector(TermStruct *ts) {
  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  Err err = NULL;
  double time;

  lc = ts->head;
  c = (IrmTermStructVal *)lc->element->val.pval;

  c->vasicek_mean_sr = 0; /* initialisation */

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;

    c->vasicek_mean_sr = p->vasicek_mean_sr +
                         p->mean_rev_level * (exp(time / p->tau) - 1) / p->F;
  }

  return err;
}

Err make_vasicek_mean_int_sr_vector(TermStruct *ts) {
  SrtLst *lp, *lc;
  IrmTermStructVal *p, *c;
  Err err = NULL;
  double time, exp_term;

  lc = ts->head;
  c = (IrmTermStructVal *)lc->element->val.pval;

  c->vasicek_mean_int_sr = 0; /* initialisation */

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;

    exp_term = (1 - exp(-time / p->tau)) * p->tau;

    c->vasicek_mean_int_sr =
        p->vasicek_mean_int_sr +
        (p->vasicek_mean_sr * p->F - p->mean_rev_level) * exp_term +
        p->mean_rev_level * time;
  }

  return err;
}

/* ----------------------------------------------------------------------- */

/*
Description: in order to make interpolation easier        , the value of F
                (resp. Psi) running up to date Ti-1 is stored on the ith
                element of the ts.
*/

Err make_F_Psi_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc : lc = current ts element */
  double time, temp, temp1;
  Err err = NULL;

  lc = l->head; /* first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->F = 1.0;
  c->Psi = 0;

  while (lc->next != NULL) /* At least 2 elements in the ts */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    temp = exp(-time / p->tau);
    if (fabs(1 / p->tau) > EPS)
      temp1 = p->tau * (1 - temp);
    else
      temp1 = time;

    c->F = p->F * temp;
    c->Psi = p->Psi + p->F * temp1;
  }
  return err;
}
/* ------------------------------------------------------------------------- */
/*
Description: in order to make interpolation easier        , the value of G
                (resp. H) running up to date Ti-1 is stored on the ith
                element of the ts.
*/
Err make_G_H_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c  : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc  with lc current ts element */
                           /* p: previous        ,  c: current */
  double temp, temp2, temp3;
  double time;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->G = 0;
  c->H = 0;

  while (lc->next != NULL) /* At least 2 elements in the ts */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    temp = exp(-time / p->tau);
    if (fabs(1 / p->tau) > EPS) {
      temp2 = (1 - temp) * p->tau;
      temp3 = (1.0 / (temp * temp) - 1.0) * p->tau;
    } else {
      temp2 = time;
      temp3 = 2.0 * time;
    }
    c->G = p->G + 0.5 * temp3 * (p->sig * p->sig) / (p->F * p->F);
    c->H = temp * p->H + temp2 * temp * p->F * p->F * p->G +
           0.5 * p->sig * p->sig * temp2 * temp2;
  }
  return err;
}
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/*
Description: in order to make interpolation easier        , the value of G
                (resp. H) running up to date Ti-1 is stored on the ith
                element of the ts for the new LGM model they are in fact equal
to sig*sig*t and tau !!!.
*/
static Err make_new_G_H_vector(TermStruct *l) {
  IrmTermStructVal *c, *p, *n; /* p < c  : c = pval of current ts element */
  SrtLst *lc, *lp, *ln;        /* lp < lc  with lc current ts element */
                               /* p: previous        ,  c: current */
  double diff = 0.0;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->G = 0;
  c->H = 0;

  while (lc != NULL) /* At least 2 elements in the ts */
  {
    if (c->val_origin == BOTH_DATE) {
      c->G = c->sig * c->sig * c->time;
      c->H = c->tau;
    } else if (c->val_origin == SIGMA_DATE) {
      c->G = c->sig * c->sig * c->time;

      ln = lc;
      while ((ln != NULL) &&
             (((IrmTermStructVal *)ln->element->val.pval)->val_origin !=
              TAU_DATE) &&
             (((IrmTermStructVal *)ln->element->val.pval)->val_origin !=
              BOTH_DATE)) {
        ln = ln->next;
      }

      lp = lc;
      while ((lp != NULL) &&
             (((IrmTermStructVal *)lp->element->val.pval)->val_origin !=
              TAU_DATE) &&
             (((IrmTermStructVal *)lp->element->val.pval)->val_origin !=
              BOTH_DATE)) {
        lp = lp->previous;
      }

      if ((lp != NULL) && (ln != NULL)) {
        p = (IrmTermStructVal *)lp->element->val.pval;
        n = (IrmTermStructVal *)ln->element->val.pval;

        c->H = p->tau +
               (n->tau - p->tau) * (c->time - p->time) / (n->time - p->time);
      }
      if (ln == NULL)
        ln = l->tail;
      if (lp == NULL)
        lp = l->head;

      if ((((IrmTermStructVal *)ln->element->val.pval)->val_origin ==
           TAU_DATE) ||
          (((IrmTermStructVal *)ln->element->val.pval)->val_origin ==
           BOTH_DATE)) {
        n = (IrmTermStructVal *)ln->element->val.pval;
        c->H = n->tau * c->time / n->time;
      }

      if ((((IrmTermStructVal *)lp->element->val.pval)->val_origin ==
           TAU_DATE) ||
          (((IrmTermStructVal *)lp->element->val.pval)->val_origin ==
           BOTH_DATE)) {
        p = (IrmTermStructVal *)ln->element->val.pval;
        c->H = p->tau * c->time / p->time;
      }

    } else if (c->val_origin == TAU_DATE) {
      c->H = c->tau;

      ln = lc;
      while ((ln != NULL) &&
             (((IrmTermStructVal *)ln->element->val.pval)->val_origin !=
              TAU_DATE) &&
             (((IrmTermStructVal *)ln->element->val.pval)->val_origin !=
              BOTH_DATE)) {
        ln = ln->next;
      }

      lp = lc;
      while ((lp != NULL) &&
             (((IrmTermStructVal *)lp->element->val.pval)->val_origin !=
              SIGMA_DATE) &&
             (((IrmTermStructVal *)lp->element->val.pval)->val_origin !=
              BOTH_DATE)) {
        lp = lp->previous;
      }

      if ((lp != NULL) && (ln != NULL)) {
        p = (IrmTermStructVal *)lp->element->val.pval;
        n = (IrmTermStructVal *)ln->element->val.pval;

        c->G = p->sig * p->sig * p->time +
               (n->sig * n->sig * n->time - p->sig * p->sig * p->time) *
                   (c->time - p->time) / (n->time - p->time);
      }

      if (ln == NULL)
        ln = l->tail;
      if (lp == NULL)
        lp = l->head;

      if ((((IrmTermStructVal *)ln->element->val.pval)->val_origin ==
           TAU_DATE) ||
          (((IrmTermStructVal *)ln->element->val.pval)->val_origin ==
           BOTH_DATE)) {
        n = (IrmTermStructVal *)ln->element->val.pval;
        c->G = n->sig * n->sig * c->time;
      }

      if ((((IrmTermStructVal *)lp->element->val.pval)->val_origin ==
           TAU_DATE) ||
          (((IrmTermStructVal *)lp->element->val.pval)->val_origin ==
           BOTH_DATE)) {
        p = (IrmTermStructVal *)ln->element->val.pval;
        c->G = p->sig * p->sig * c->time;
      }
    }

    lc = lc->next;
    if (lc != NULL) {
      c = (IrmTermStructVal *)lc->element->val.pval;
    }
  }

  /* CHECK THAT at least G is an non negative increasing function */
  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;

  if (c->G <= 0.0) {
    err = "term Struct for New Model has negative value";
    return err;
  }

  while (lc->next != NULL) /* At least 2 elements in the ts */
  {
    diff = -c->G;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;
    diff += c->G;
    /* G has to be possitive */
    if (c->G <= 0.0) {
      err = "term Struct for New Model has negative value";
      return err;
    }
    /* G has to be increasing */
    if (diff < 0.0) {
      err = "Vol term Struct for New Model has to be cumulative vol :: check "
            "the input";
      return err;
    }
  }

  return err;
}
/* ------------------------------------------------------------------------- */

static Err make_J_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c  : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc  with lc current ts element */
                           /* p: previous        , c: current */
  double time;
  double temp, temp2;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->J = 0;

  while (lc->next != NULL)
  /* lc ranges from second element to next to the last */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;                             /* <=> [i-1] */
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    temp = exp(-time / p->tau);
    if (fabs(1 / p->tau) > EPS)
      temp2 = (temp - 1) * p->tau;
    else
      temp2 = time;
    c->J = p->J + temp2 * p->sig / p->F;
  }
  return err;
}

/* ------------------------------------------------------------------------- */
static Err make_I_vector(TermStruct *l) {

  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1;
  double tau, sig, F, G;
  Err err = 0;

  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->I = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    tau = p->tau;
    sig = p->sig;
    F = p->F;
    G = p->G;

    temp = (1 - exp(-time / tau)) * tau;
    temp1 = (exp(time / tau) + exp(-time / tau) - 2);

    c->I = p->I + F * G * temp + (sig * tau) * (sig * tau) * (temp1) / (2 * F);
  }
  return err;
}

/* ------------------------------------------------------------------------- */

static Err make_K_vector(TermStruct *l) {

  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp;
  double tau, sig, F, J;

  Err err = 0;
  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->K = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    tau = p->tau;
    sig = p->sig;
    F = p->F;
    J = p->J;
    temp = (1 - exp(-time / tau)) * tau;

    c->K = p->K + F * J * temp + sig * tau * (time - temp);
  }
  return err;
}
/* ------------------------------------------------------------------------- */

static Err make_L_vector(TermStruct *l) {

  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;

  Err err = 0;
  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->L = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    c->L = p->L + p->F * p->G * (1 - exp(-time / p->tau)) * p->tau +
           0.5 * pow(p->sig * p->tau, 2) *
               (exp(-time / p->tau) + exp(time / p->tau) - 2) / p->F;
  }
  return err;
}
/* ------------------------------------------------------------------------- */

static Err make_O_vector(TermStruct *l) {

  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  Err err = NULL;

  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->O = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;
    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;

    c->O = p->O + p->F * p->G * p->Psi * (1 - exp(-time / p->tau)) * p->tau +
           p->F * p->F * p->G * p->tau *
               ((1 - exp(-time / p->tau)) * p->tau -
                0.5 * (1 - exp(-2 * time / p->tau)) * p->tau) +
           0.5 * (p->Psi / p->F) * p->tau * p->sig * p->sig *
               ((exp(time / p->tau) - 1) * p->tau -
                (1 - exp(-time / p->tau)) * p->tau) +
           0.5 * pow(p->sig * p->tau, 2) *
               ((exp(time / p->tau) - 1) * p->tau -
                (1 - exp(-time / p->tau)) * p->tau) +
           0.5 * pow(p->sig * p->tau, 2) *
               (0.5 * (1 - exp(-2 * time / p->tau)) * p->tau - time);
  }
  return err;
}
/* ------------------------------------------------------------------------- */

static Err make_Q_vector(TermStruct *l) {

  IrmTermStructVal *p, *c;
  SrtLst *lp, *lc;
  double time;
  double temp, temp1;
  Err err = 0;

  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;

  c->Q = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;

    temp = (1.0 - exp(-time / p->tau)) * p->tau;
    temp1 = 0.5 * (1.0 - exp(-2 * time / p->tau)) * p->tau;

    c->Q = p->Q + p->H * temp + p->F * p->F * p->G * p->tau * (temp - temp1) +
           0.5 * p->sig * p->sig * p->tau * p->tau * (time - 2 * temp + temp1);
  }
  return err;
}

/*------------------------------------------------------------------*/

/* ------------------------------------------------------------------------- */

static Err make_R_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c  : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc  with lc current ts element */
                           /* p: previous        , c: current */
  double time;
  double temp;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->R = 0;

  while (lc->next != NULL)
  /* lc ranges from second element to next to the last */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;                             /* <=> [i-1] */
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    if (fabs(1 / p->tau) > EPS)
      temp = 0.5 * (exp(2 * time / p->tau) - 1) * p->tau;
    else
      temp = time;

    c->R = p->R + (p->sig * p->sig / (p->F * p->F)) * temp;
  }
  return err;
}

/*--------------------------------------------------------------------*/

static Err make_S_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c  : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc  with lc current ts element */
                           /* p: previous        , c: current */
  double time;
  double temp, temp1;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->S = 0;

  while (lc->next != NULL)
  /* lc ranges from second element to next to the last */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;                             /* <=> [i-1] */
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    if (fabs(1 / p->tau) > EPS)
      temp = (exp(time / p->tau) - 1) * p->tau;
    else
      temp = time;

    if (fabs(2 / p->tau) > EPS)
      temp1 = 0.5 * (exp(2 * time / p->tau) - 1) * p->tau;
    else
      temp1 = time;

    c->S = p->S + (p->sig * p->sig / (p->F * p->F)) * p->Psi * temp1 +
           (p->sig * p->sig * p->tau) * (temp1 - temp) / (p->F);
  }
  return err;
}

static Err make_T_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c  : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc  with lc current ts element */
                           /* p: previous        , c: current */
  double time;
  double temp, temp1, temp2;
  Err err = 0;

  lc = l->head; /* lc as the first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->T = 0;

  while (lc->next != NULL)
  /* lc ranges from second element to next to the last */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;                             /* <=> [i-1] */
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    if (fabs(1 / p->tau) > EPS)
      temp = (1 - exp(-time / p->tau)) * p->tau;
    else
      temp = time;

    if (fabs(1 / p->tau) > EPS)
      temp1 = (exp(time / p->tau) - 1) * p->tau;
    else
      temp1 = time;

    if (fabs(2 / p->tau) > EPS)
      temp2 = 0.5 * (exp(2 * time / p->tau) - 1) * p->tau;
    else
      temp2 = time;

    c->T = p->T + p->F * p->S * temp +
           0.5 * p->sig * p->sig * p->Psi / p->F * p->tau * (temp1 - temp) +
           0.5 * p->sig * p->sig * p->tau * p->tau * (temp1 - temp) -
           p->sig * p->sig * p->tau * p->tau * (time - temp);
  }

  return err;
}

/* ------------------------------------------------------------------------- */
static Err make_Phi_vector(TermStruct *l) {

  IrmTermStructVal *c;
  SrtLst *lc;

  Err err = 0;

  lc = l->head; /* first element of the list */
  c = (IrmTermStructVal *)lc->element->val.pval;

  c->Phi = 0.0;

  while (lc->next != NULL) /* At least 2 elements in the list */
  {
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    c->Phi = c->F * c->F * c->G;
  }
  return err;
}

static Err make_int_phi_vector(TermStruct *l) {
  IrmTermStructVal *c, *p;
  SrtLst *lc, *lp;
  Err err = NULL;
  double time, temp;

  lc = l->head;
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->int_phi = 0.0;

  while (lc->next != NULL) {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval;
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval;

    temp = 0.5 * (1 - exp(-2 * time / p->tau)) * p->tau;

    c->int_phi = p->int_phi + p->F * p->F * p->G * temp +
                 0.5 * p->sig * p->sig * p->tau * (time - temp);
  }

  return err;
}

/* ------------------------------------------------------------------------- */

static Err make_Zeta_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc : lc = current ts element */
  double time;
  Err err = NULL;

  lc = l->head; /* first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->Zeta = 0.0;

  while (lc->next != NULL) /* At least 2 elements in the ts */
  {
    if (lc->previous != NULL) {
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */
      time = c->time - p->time;
    } else
      time = c->time;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    c->Zeta = p->Zeta + c->sig * c->sig * time;
  }
  return err;
}

/* -------------------------------------------------------------------------- */

Err make_Lambda_vector(TermStruct *l) {
  IrmTermStructVal *p, *c; /* p < c : c = pval of current ts element */
  SrtLst *lp, *lc;         /* lp < lc : lc = current ts element */
  double time_prev;
  Err err = NULL;

  lc = l->head; /* first element of the ts */
  c = (IrmTermStructVal *)lc->element->val.pval;
  c->Lambda = 0.0;

  while (lc->next != NULL) /* At least 2 elements in the ts */
  {
    if (lc->previous != NULL) {
      time_prev = p->time;
      lp = lc->previous;
      p = (IrmTermStructVal *)lp->element->val.pval; /* <=> [i-1] */

    } else
      time_prev = 0.;

    p = c;
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */

    c->Lambda = p->Lambda +
                p->tau * (exp(-time_prev / p->tau) - exp(-p->time / p->tau));
  }
  return err;
}

/* ------------------------------------------------------------------------- */
