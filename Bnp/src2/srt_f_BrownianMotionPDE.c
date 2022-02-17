
#include "SRT_H_ALL.H>
#include "math.h"
#include "srt_h_BrownianMotion.h"

#define ONE_HALF 0.5
#define PDE_THETA 0.5
#define TREE_MESH_SPACING 1.0
#define TREE_MESH_SPACING_SQUARE 1.0
#define LIM_IN_STDEV 6
#define MAX_TIME 0.0416666666666667
#define MIN_NODE 30
#define MAX_NODE 1000

static srt_black_scholes_parms static_black_scholes_parms;

static Err set_pde_bound(long pde_dim, double *mesh_bound,
                         Err (*pde_boundary_func)(), double **pde_bound) {
  Err err = NULL;
  long i;

  for (i = 0; i < pde_dim; i++) {
    err = pde_boundary_func(mesh_bound[i], &(*pde_bound)[i]);
    if (err)
      return err;
  }

  return err;
}

static Err srt_f_eur_option_payoff(double dwt, double *payoff);
static Err srt_f_amer_option_payoff(Date eval_date, Date ExerDate, double dwt,
                                    double *payoff);

Err srt_f_brownain_motion_eur_pde(
    Date eval_date,            /* EVALUATION DATE */
    Date exp_date,             /* EXPIRY DATE */
    Err (*term_payoff_func)(), /* FINAL MATURITY payoff */
    long min_node, long num_mesh, SrtPdeSolvingScheme pde_solving_scheme,
    SrtPdeBoundaryCond pde_boundary_cond, double *pde_bound,

    double **div_mat,

    /*OUTPUTS */
    double *option_pv) {

  Err err = NULL;
  double *cur_premium_ptr = NULL, *prev_premium_ptr = NULL, *payoff = NULL,
         *X = NULL, prob[3], yr_to_exp_date, time_step = MAX_TIME, dX;
  long max_index, num_time_step, i, j;
  double *diag_A = NULL, *diag_B = NULL, *diag_C = NULL, *vect_R = NULL,
         *vect_U = NULL, pde_params[2];

  /* FORCE THE BOUNDARY CONDITION TYPE */
  pde_boundary_cond = DIRICHLET;

  yr_to_exp_date = (double)(exp_date - eval_date) * YEARS_IN_DAY;

  time_step = min(time_step, yr_to_exp_date / min_node);
  num_time_step = ((long)(yr_to_exp_date / time_step)) + 1;

  /* IF THE SCHEME IS EXPLICIT (TREE) */
  if (pde_solving_scheme == EXPLICIT) {
    dX = TREE_MESH_SPACING * sqrt(time_step);
    max_index = (long)DTOL(LIM_IN_STDEV * sqrt(yr_to_exp_date) / dX) + 1;

    prev_premium_ptr = dvector(-max_index, max_index);
    cur_premium_ptr = dvector(-max_index, max_index);
    payoff = dvector(-max_index, max_index);
    X = dvector(-max_index, max_index);

    /* SET THE TRANSITION PROBABILITIES */
    prob[0] = ONE_HALF * TREE_MESH_SPACING;
    prob[2] = ONE_HALF * TREE_MESH_SPACING;
    prob[1] = 1.0 - prob[0] - prob[2];

    /* FILL THE STATE VARIABLES */
    for (j = -max_index; j <= max_index; j++)
      X[j] = (j * dX);

    /*payoff FUNCTION AT THE LAST STEP*/
    for (j = -max_index; j <= max_index; j++) {
      err = (*term_payoff_func)(X[j], &payoff[j]);
      if (err)
        return err;

      prev_premium_ptr[j] = payoff[j];
    }

    /* start_date OF THE forward ALGORITHM */
    for (i = 1; i < num_time_step; i++) /* LOOP ON THE TIME STEP */
    {
      if (pde_boundary_cond ==
          DIRICHLET) /*CASE OF OF DIRICHLET BOUNDARIES CONDITIONS */
      {
        for (j = (-max_index + 1); j <= (max_index - 1);
             j++) /* LOOP ON THE SPACE POINTS */
        {
          cur_premium_ptr[j] = prob[0] * prev_premium_ptr[j - 1] +
                               prob[1] * prev_premium_ptr[j] +
                               prob[2] * prev_premium_ptr[j + 1];
        }

        /* APPLY THE DIRICHLET CONDITIONS */
        cur_premium_ptr[-max_index] = pde_bound[0];
        cur_premium_ptr[max_index] = pde_bound[1];

        for (j = -max_index; j <= max_index; j++)
          prev_premium_ptr[j] = cur_premium_ptr[j];

      } /* end_date OF pde_boundary_cond == DIRICHLET */

      else if (pde_boundary_cond == NEUMANN) /*CASE OF OF NEUMANN CONDITIONS */
      {

      } /* end_date OF pde_boundary_cond == NEUMANN */

    } /*end_date OF LOOP ON THE SPACE POINTS */

    (*option_pv) = prev_premium_ptr[0];

  } /* end_date OF THE EXPLICIT SCHEME */

  else if (pde_solving_scheme == THETASCHEME_CENTRED) {
    dX = LIM_IN_STDEV * sqrt(yr_to_exp_date) / num_mesh;
    max_index = (long)DTOL(LIM_IN_STDEV * sqrt(yr_to_exp_date) / dX);

    /* ALLOCATIONS OF THE TRIDAG STRUCTURES */
    diag_A = dvector(1, max_index + max_index - 1);
    diag_B = dvector(1, max_index + max_index - 1);
    diag_C = dvector(1, max_index + max_index - 1);
    vect_R = dvector(1, max_index + max_index - 1);
    vect_U = dvector(1, max_index + max_index - 1);

    prev_premium_ptr = dvector(-max_index, max_index);
    payoff = dvector(-max_index, max_index);
    X = dvector(-max_index, max_index);

    /* FILL THE STATE VARIABLES */
    for (j = -max_index; j <= max_index; j++)
      X[j] = (j * dX);

    pde_params[0] = ONE_HALF * PDE_THETA * time_step / (dX * dX);
    pde_params[1] = ONE_HALF * (PDE_THETA - 1.0) * time_step / (dX * dX);

    /* FILL THE DIAG MATRIX FOR THE TRIDAG */
    for (j = 1; j < (max_index + max_index); j++) {
      diag_A[j] = pde_params[0];
      diag_B[j] = -1.0 - 2 * pde_params[0];
      diag_C[j] = pde_params[0];
    }

    for (j = -max_index; j <= max_index; j++) {
      err = (*term_payoff_func)(X[j], &payoff[j]);
      if (err)
        return err;

      prev_premium_ptr[j] = payoff[j];
    }

    /* start_date THE CONVOLUTION - LOOP ON TIME STEP */
    for (i = 1; i < num_time_step; i++) /* LOOP ON THE TIME STEP */
    {
      for (j = 1; j < (max_index + max_index);
           j++) /* LOOP ON THE SPACE POINTS */
      {
        if (j == 1) /* BOUNDARY CONDITION TO APPLY */
        {
          vect_R[j] =
              pde_params[1] * prev_premium_ptr[j + 1 - max_index] -
              (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
              pde_params[1] * pde_bound[0] - pde_params[0] * pde_bound[0]

              + (1 - PDE_THETA) * div_mat[i - 1][j] + PDE_THETA * div_mat[i][j];
        } else if (j == (max_index + max_index -
                         1)) /* BOUNDARY CONDITION TO APPLY */
        {
          vect_R[j] =
              pde_params[1] * pde_bound[1] - pde_params[0] * pde_bound[1] -
              (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
              pde_params[1] * prev_premium_ptr[j - 1 - max_index]

              + (1 - PDE_THETA) * div_mat[i - 1][j] + PDE_THETA * div_mat[i][j];

        } else /* NO BOUNDARY CONDITION TO APPLY */
        {
          vect_R[j] =
              pde_params[1] * prev_premium_ptr[j + 1 - max_index] -
              (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
              pde_params[1] * prev_premium_ptr[j - 1 - max_index]

              + (1 - PDE_THETA) * div_mat[i - 1][j] + PDE_THETA * div_mat[i][j];
        }

      } /* end_date LOOP ON THE SPACE POINTS */

      err = tridag(diag_A, diag_B, diag_C, vect_R, vect_U,
                   (max_index + max_index - 1));
      if (err)
        return err;

      for (j = (-max_index + 1); j <= (max_index - 1); j++)
        prev_premium_ptr[j] = vect_U[j + max_index];

    } /*end_date OF LOOP ON TIME STEP */

    (*option_pv) = prev_premium_ptr[0];

    if (diag_A)
      free_dvector(diag_A, 1, (max_index + max_index - 1));
    diag_A = NULL;
    if (diag_B)
      free_dvector(diag_B, 1, (max_index + max_index - 1));
    diag_B = NULL;
    if (diag_C)
      free_dvector(diag_C, 1, (max_index + max_index - 1));
    diag_C = NULL;
    if (vect_R)
      free_dvector(vect_R, 1, (max_index + max_index - 1));
    vect_R = NULL;
    if (vect_U)
      free_dvector(vect_U, 1, (max_index + max_index - 1));
    vect_U = NULL;
    if (prev_premium_ptr)
      free_dvector(prev_premium_ptr, -max_index, max_index);
    prev_premium_ptr = NULL;
    if (payoff)
      free_dvector(payoff, -max_index, max_index);
    payoff = NULL;
  }

  return (err);
}

Err srt_f_brownain_motion_amer_pde(
    Date eval_date, /* EVALUATION DATE */
    Date *tEx,      /* FROM 0 TO nEx-1 */
    long nEx, long min_node, long num_mesh,
    Err (*amer_payoff_func)(), /* AMERICAN payoff */
    SrtPdeSolvingScheme pde_solving_scheme,
    SrtPdeBoundaryCond pde_boundary_cond, double **pde_bound,

    /*OUTPUTS */
    double *option_pv) {

  Err err = NULL;
  double *cur_premium_ptr = NULL, *prev_premium_ptr = NULL, *payoff = NULL,
         *X = NULL, prob[3], yr_to_exp_date, time_step = MAX_TIME, dX, Weight;
  long max_index, num_time_step, i, j, n, pde_dim = 1, k, index, N;
  long max_node = MAX_NODE, index_exercise = 0, num_interval;
  double *diag_A = NULL, *diag_B = NULL, *diag_C = NULL, *vect_R = NULL,
         *vect_U = NULL, pde_params[2];
  Date *tDeal;
  double start_date, end_date, dt;
  double **pde_boundary_deal, Xmin, Xmax;

  /** FIRST STEP INITIALISATION OF THE DATES **/
  yr_to_exp_date = (double)(tEx[nEx - 1] - eval_date) * YEARS_IN_DAY;
  time_step = min(time_step, yr_to_exp_date / min_node);

  tDeal = (Date *)malloc((max_node) * sizeof(Date));

  n = 0;
  tDeal[0] = eval_date;
  start_date = (double)eval_date;
  for (i = 0; i < nEx; i++) {
    end_date = (double)tEx[i];
    num_interval =
        ((long)((end_date - start_date) * YEARS_IN_DAY / time_step)) + 1;
    num_interval = min(num_interval, (long)(end_date - start_date));
    dt = (end_date - start_date) * YEARS_IN_DAY / num_interval;
    time_step = min(time_step, dt);
    for (j = 1; j < num_interval; j++) {
      n++;
      tDeal[n] = (Date)((double)tDeal[n - 1] + dt / YEARS_IN_DAY);
    }
    n++;
    tDeal[n] = tEx[i];
    start_date = end_date;
  }

  num_time_step = ((long)(yr_to_exp_date / time_step)) +
                  1; /** end_date OF THE FIRST STEP */

  /* GET THE pde_bound AT EACH EXERCISE DATE */

  if (pde_solving_scheme == EXPLICIT)
    dX = TREE_MESH_SPACING * sqrt(time_step);
  else if (pde_solving_scheme == THETASCHEME_CENTRED)
    dX = LIM_IN_STDEV * sqrt(yr_to_exp_date) / num_mesh;

  max_index = (long)DTOL(LIM_IN_STDEV * sqrt(yr_to_exp_date) / dX) + 1;

  if (pde_boundary_cond == DIRICHLET) {
    Xmin = -max_index * dX;
    Xmax = max_index * dX;

    for (i = 0; i < nEx; i++) {
      err = (*amer_payoff_func)(eval_date, tEx[i], Xmin, &pde_bound[0][i]);
      if (err)
        return err;

      err = (*amer_payoff_func)(eval_date, tEx[i], Xmax, &pde_bound[1][i]);
      if (err)
        return err;

    } /* end_date OF GETTING THE pde_bound AT EACH EXERCISE DATE */
  } else if (pde_boundary_cond == NEUMANN) {
    for (i = 0; i < nEx; i++) {
      pde_bound[0][i] = 0.;
      pde_bound[1][i] = 1.;
    }
  }

  /* GET THE pde_bound AT EACH DEAL DATE BY INTERPOLATION */
  pde_boundary_deal = dmatrix(0, pde_dim, 0, num_time_step);
  k = 1;
  while ((k < num_time_step) && (tDeal[k] < tEx[0])) {
    Weight = ((double)(tDeal[k] - eval_date)) / (tEx[0] - eval_date);
    pde_boundary_deal[0][k] = Weight * pde_bound[0][0];
    pde_boundary_deal[1][k] = Weight * pde_bound[1][0];

    k++;
  }

  pde_boundary_deal[0][k] = pde_bound[0][0];
  pde_boundary_deal[1][k] = pde_bound[1][0];

  index = 1;
  while ((k < n)) {
    if (tDeal[k] == tEx[index]) {
      pde_boundary_deal[0][k] = pde_bound[0][index];
      pde_boundary_deal[1][k] = pde_bound[1][index];
    } else {
      while (tDeal[k] > tEx[index])
        index++;
      Weight =
          ((double)(tDeal[k] - tEx[index - 1])) / (tEx[index] - tEx[index - 1]);
      pde_boundary_deal[0][k] =
          pde_bound[0][index - 1] +
          Weight * (pde_bound[0][index] - pde_bound[0][index - 1]);
      pde_boundary_deal[1][k] =
          pde_bound[1][index - 1] +
          Weight * (pde_bound[1][index] - pde_bound[1][index - 1]);
    }

    k++;
  }
  /* end_date OF THE INTERPOLATION OF THE BOUNDARY CONDITIONS */

  /* IF THE SCHEME IS EXPLICIT (TREE) */
  if (pde_solving_scheme == EXPLICIT) {

    prev_premium_ptr = dvector(-max_index, max_index);
    cur_premium_ptr = dvector(-max_index, max_index);
    payoff = dvector(-max_index, max_index);
    X = dvector(-max_index, max_index);

    /* SET THE TRANSITION PROBABILITIES */
    prob[0] = ONE_HALF * TREE_MESH_SPACING;
    prob[2] = ONE_HALF * TREE_MESH_SPACING;
    prob[1] = 1.0 - prob[0] - prob[2];

    /* FILL THE STATE VARIABLES */
    for (j = -max_index; j <= max_index; j++)
      X[j] = (j * dX);

    /*payoff FUNCTION AT THE LAST EXERCISE DATE: tEx[nEx-1] */
    for (j = -max_index; j <= max_index; j++) {
      err = (*amer_payoff_func)(eval_date, tEx[nEx - 1], X[j], &payoff[j]);
      if (err)
        return err;

      prev_premium_ptr[j] = payoff[j];
    }

    /* start_date OF THE forward ALGORITHM */
    index_exercise = nEx - 2;
    for (i = (num_time_step - 1); i > 0; i--) /* LOOP ON THE TIME STEP */
    {
      if (pde_boundary_cond ==
          DIRICHLET) /*CASE OF OF DIRICHLET BOUNDARIES CONDITIONS */
      {
        for (j = (-max_index + 1); j <= (max_index - 1);
             j++) /* LOOP ON THE SPACE POINTS */
        {
          cur_premium_ptr[j] = prob[0] * prev_premium_ptr[j - 1] +
                               prob[1] * prev_premium_ptr[j] +
                               prob[2] * prev_premium_ptr[j + 1];
        }

        /* APPLY THE DIRICHLET CONDITIONS */
        cur_premium_ptr[-max_index] = pde_boundary_deal[0][i];
        cur_premium_ptr[max_index] = pde_boundary_deal[1][i];

        if (tDeal[i] == tEx[index_exercise]) {
          for (j = -max_index; j <= max_index; j++) {
            err = (*amer_payoff_func)(eval_date, tEx[index_exercise], X[j],
                                      &payoff[j]);
            if (err)
              return err;

            prev_premium_ptr[j] = max(cur_premium_ptr[j], payoff[j]);
          }

          index_exercise--;

        } else {
          for (j = -max_index; j <= max_index; j++)
            prev_premium_ptr[j] = cur_premium_ptr[j];
        }

      } /* end_date OF pde_boundary_cond == DIRICHLET */

      else if (pde_boundary_cond == NEUMANN) /*CASE OF OF NEUMANN CONDITIONS */
      {

      } /* end_date OF pde_boundary_cond == NEUMANN */

    } /*end_date OF LOOP ON THE SPACE POINTS */

    (*option_pv) = prev_premium_ptr[0];

  } /* end_date OF THE EXPLICIT SCHEME */

  else if (pde_solving_scheme == THETASCHEME_CENTRED) {
    if (pde_boundary_cond == DIRICHLET)
      N = (max_index + max_index - 1);
    else if (pde_boundary_cond == NEUMANN)
      N = (max_index + max_index + 1);

    /* ALLOCATIONS OF THE TRIDAG STRUCTURES */
    diag_A = dvector(1, N);
    diag_B = dvector(1, N);
    diag_C = dvector(1, N);
    vect_R = dvector(1, N);
    vect_U = dvector(1, N);

    prev_premium_ptr = dvector(-max_index, max_index);
    payoff = dvector(-max_index, max_index);
    X = dvector(-max_index, max_index);

    /* FILL THE STATE VARIABLES */
    for (j = -max_index; j <= max_index; j++)
      X[j] = (j * dX);

    pde_params[0] = ONE_HALF * PDE_THETA * time_step / (dX * dX);
    pde_params[1] = ONE_HALF * (PDE_THETA - 1.0) * time_step / (dX * dX);

    /* FILL THE DIAG MATRIX FOR THE TRIDAG */
    if (pde_boundary_cond == DIRICHLET) {
      for (j = 1; j <= N; j++) {
        diag_A[j] = pde_params[0];
        diag_B[j] = -1.0 - 2 * pde_params[0];
        diag_C[j] = pde_params[0];
      }

    } else if (pde_boundary_cond == NEUMANN) {
      diag_A[1] = pde_params[0];
      diag_B[1] = -1.0 - pde_params[0]; /* DUE TO NEUMANN CONDITION */
      diag_C[1] = pde_params[0];

      for (j = 2; j < N; j++) {
        diag_A[j] = pde_params[0];
        diag_B[j] = -1.0 - 2 * pde_params[0];
        diag_C[j] = pde_params[0];
      }

      diag_A[N] = pde_params[0];
      diag_B[N] = -1.0 - pde_params[0]; /* DUE TO NEUMANN CONDITION */
      diag_C[N] = pde_params[0];
    }

    for (j = -max_index; j <= max_index; j++) {
      err = (*amer_payoff_func)(eval_date, tEx[nEx - 1], X[j], &payoff[j]);
      if (err)
        return err;

      prev_premium_ptr[j] = payoff[j];
    }

    /* start_date THE CONVOLUTION - LOOP ON TIME STEP */
    index_exercise = nEx - 2;
    for (i = (num_time_step - 1); i > 0; i--) /* LOOP ON THE TIME STEP */
    {
      for (j = 1; j <= N; j++) /* LOOP ON THE SPACE POINTS */
      {
        if (j == 1) /* BOUNDARY CONDITION TO APPLY */
        {
          if (pde_boundary_cond == DIRICHLET) {
            vect_R[j] =
                pde_params[1] * prev_premium_ptr[j + 1 - max_index] -
                (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
                pde_params[1] * pde_boundary_deal[0][i] -
                pde_params[0] * pde_boundary_deal[0][i];
          } else if (pde_boundary_cond == NEUMANN) {
            vect_R[j] =
                pde_params[1] * prev_premium_ptr[j + 1 - max_index] -
                (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index]

                + pde_params[1] * (prev_premium_ptr[j - max_index] -
                                   pde_boundary_deal[0][i] * dX) +
                pde_params[0] * pde_boundary_deal[0][i] * dX;
          }
        } else if (j == N) /* BOUNDARY CONDITION TO APPLY */
        {
          if (pde_boundary_cond == DIRICHLET) {
            vect_R[j] =
                pde_params[1] * pde_boundary_deal[1][i] -
                pde_params[0] * pde_boundary_deal[1][i] -
                (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
                pde_params[1] * prev_premium_ptr[j - 1 - max_index];
          } else if (pde_boundary_cond == NEUMANN) {
            vect_R[j] =
                pde_params[1] * (pde_boundary_deal[1][i] * dX +
                                 prev_premium_ptr[j - max_index]) -
                pde_params[0] * pde_boundary_deal[1][i] * dX -
                (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
                pde_params[1] * prev_premium_ptr[j - 1 - max_index];
          }

        } else /* NO BOUNDARY CONDITION TO APPLY */
        {
          vect_R[j] =
              pde_params[1] * prev_premium_ptr[j + 1 - max_index] -
              (1.0 + 2.0 * pde_params[1]) * prev_premium_ptr[j - max_index] +
              pde_params[1] * prev_premium_ptr[j - 1 - max_index];
        }

      } /* end_date LOOP ON THE SPACE POINTS */

      err = tridag(diag_A, diag_B, diag_C, vect_R, vect_U, N);
      if (err)
        return err;

      if (tDeal[i] == tEx[index_exercise]) {
        if (pde_boundary_cond == DIRICHLET) {
          for (j = (-max_index + 1); j <= (max_index - 1); j++) {
            err = (*amer_payoff_func)(eval_date, tEx[index_exercise], X[j],
                                      &payoff[j]);
            if (err)
              return err;

            prev_premium_ptr[j] = max(vect_U[j + max_index], payoff[j]);
          }
        } else if (pde_boundary_cond == NEUMANN) {
          for (j = -max_index; j <= max_index; j++) {
            err = (*amer_payoff_func)(eval_date, tEx[index_exercise], X[j],
                                      &payoff[j]);
            if (err)
              return err;

            prev_premium_ptr[j] = max(vect_U[j + max_index], payoff[j]);
          }
        }

        index_exercise--;
      } else {
        if (pde_boundary_cond == DIRICHLET) {
          for (j = (-max_index + 1); j <= (max_index - 1); j++)
            prev_premium_ptr[j] = vect_U[j + max_index];
        } else if (pde_boundary_cond == NEUMANN) {
          for (j = -max_index; j <= max_index; j++)
            prev_premium_ptr[j] = vect_U[j + max_index];
        }
      }

    } /*end_date OF LOOP ON TIME STEP */

    (*option_pv) = prev_premium_ptr[0];

    if (diag_A)
      free_dvector(diag_A, 1, N);
    diag_A = NULL;
    if (diag_B)
      free_dvector(diag_B, 1, N);
    diag_B = NULL;
    if (diag_C)
      free_dvector(diag_C, 1, N);
    diag_C = NULL;
    if (vect_R)
      free_dvector(vect_R, 1, N);
    vect_R = NULL;
    if (vect_U)
      free_dvector(vect_U, 1, N);
    vect_U = NULL;
    if (prev_premium_ptr)
      free_dvector(prev_premium_ptr, -max_index, max_index);
    prev_premium_ptr = NULL;
    if (payoff)
      free_dvector(payoff, -max_index, max_index);
    payoff = NULL;
    if (pde_boundary_deal)
      free_dmatrix(pde_boundary_deal, 0, pde_dim, 0, num_time_step);
    pde_boundary_deal = NULL;
  }

  return (err);
}

Err srt_f_pde_eur_option(Date eval_date, Date exp_date, double forward,
                         double equity_strike, double equity_vol,
                         char *rec_pay_str, char *pde_solving_scheme_str,
                         long min_node, long num_mesh,

                         /*OUTPUTS */
                         double *option_pv)

{
  Err err = NULL;
  double yr_to_exp_date, time_step, dX;
  long num_time_step, max_index;
  long pde_dim = 1;
  double yr_to_exp, *pde_bound, *mesh_bound, **div_mat;
  SrtReceiverType pay_rec;
  SrtPdeSolvingScheme pde_solving_scheme;

  /* FILL THE EUROPEAN PARAMS STRUCTURE */
  static_black_scholes_parms.equity_vol = equity_vol;
  static_black_scholes_parms.forward = forward;
  static_black_scholes_parms.equity_strike = equity_strike;
  yr_to_exp = (double)(exp_date - eval_date) * YEARS_IN_DAY;
  static_black_scholes_parms.yr_to_exp = yr_to_exp;

  err = interp_rec_pay(rec_pay_str, &pay_rec);
  if (err)
    return err;

  /* GET THE pde_solving_scheme*/
  strupper(pde_solving_scheme_str);
  strip_white_space(pde_solving_scheme_str);

  if (!strcmp(pde_solving_scheme_str, "EXPLICIT"))
    pde_solving_scheme = EXPLICIT;
  else if (!strcmp(pde_solving_scheme_str, "THETASCHEME_CENTRED"))
    pde_solving_scheme = THETASCHEME_CENTRED;

  static_black_scholes_parms.pay_rec = pay_rec;

  /* GET THE DIRICHLET CONDITION */
  pde_bound = dvector(0, pde_dim);
  mesh_bound = dvector(0, pde_dim);

  mesh_bound[0] = -LIM_IN_STDEV * sqrt(yr_to_exp);
  mesh_bound[1] = LIM_IN_STDEV * sqrt(yr_to_exp);

  err =
      set_pde_bound(pde_dim, mesh_bound, &srt_f_eur_option_payoff, &pde_bound);
  if (err)
    return err;

  /* ALLOCATE THE DIVIDEND MATRIX */
  yr_to_exp_date = (double)(exp_date - eval_date) * YEARS_IN_DAY;
  time_step = min(MAX_TIME, yr_to_exp_date / min_node);
  num_time_step = ((long)(yr_to_exp_date / time_step)) + 1;

  dX = LIM_IN_STDEV * sqrt(yr_to_exp_date) / num_mesh;
  max_index = (long)DTOL(LIM_IN_STDEV * sqrt(yr_to_exp_date) / dX);

  div_mat = dmatrix(0, num_time_step - 1, 0, max_index + max_index - 1);

  err = srt_f_brownain_motion_eur_pde(
      eval_date, exp_date, &srt_f_eur_option_payoff, min_node, num_mesh,
      pde_solving_scheme, DIRICHLET, /*BOUNDARY CONDITION TYPE */
      pde_bound,

      div_mat,

      option_pv);

  if (err)
    return err;

  if (pde_bound)
    free_dvector(pde_bound, 0, pde_dim);
  pde_bound = NULL;
  if (mesh_bound)
    free_dvector(mesh_bound, 0, pde_dim);
  mesh_bound = NULL;
  if (div_mat)
    free_dmatrix(div_mat, 0, num_time_step - 1, 0, max_index + max_index - 1);
  div_mat = NULL;

  return err;
}

Err srt_f_pde_amer_option(Date eval_date, Date *tEx, long nEx, double forward,
                          double equity_strike, double equity_vol,
                          char *rec_pay_str, char *pde_solving_scheme_str,
                          char *pde_boundary_cond_str, long min_node,
                          long num_mesh,

                          /*OUTPUTS */
                          double *option_pv)

{
  Err err = NULL;
  long pde_dim = 1;
  double **pde_bound;
  SrtReceiverType pay_rec;
  SrtPdeSolvingScheme pde_solving_scheme;
  SrtPdeBoundaryCond pde_boundary_cond;

  /* FILL THE EUROPEAN PARAMS STRUCTURE */
  static_black_scholes_parms.equity_vol = equity_vol;
  static_black_scholes_parms.forward = forward;
  static_black_scholes_parms.equity_strike = equity_strike;

  err = interp_rec_pay(rec_pay_str, &pay_rec);
  if (err)
    return err;

  /* GET THE pde_solving_scheme*/
  strupper(pde_solving_scheme_str);
  strip_white_space(pde_solving_scheme_str);
  strip_white_space(pde_boundary_cond_str);

  if (!strcmp(pde_solving_scheme_str, "EXPLICIT"))
    pde_solving_scheme = EXPLICIT;
  else if (!strcmp(pde_solving_scheme_str, "THETASCHEME_CENTRED"))
    pde_solving_scheme = THETASCHEME_CENTRED;

  if (!strcmp(pde_boundary_cond_str, "NEUMANN"))
    pde_boundary_cond = NEUMANN;
  else
    pde_boundary_cond = DIRICHLET;

  static_black_scholes_parms.pay_rec = pay_rec;

  /* GET THE DIRICHLET CONDITION FOR EACH EXERCISE DATE */
  pde_bound = dmatrix(0, pde_dim, 0, nEx - 1);

  err = srt_f_brownain_motion_amer_pde(
      eval_date, tEx, nEx, min_node, num_mesh, &srt_f_amer_option_payoff,
      pde_solving_scheme, pde_boundary_cond, pde_bound,

      option_pv);

  if (err)
    return err;

  if (pde_bound)
    free_dmatrix(pde_bound, 0, pde_dim, 0, nEx - 1);
  pde_bound = NULL;

  return err;
}

static Err srt_f_eur_option_payoff(double dwt, double *payoff) {
  double sq_equity_vol, equity_vol, forward, equity_strike, yr_to_exp;
  SrtReceiverType pay_rec;

  /* GET STATIC INFORMATIONS */
  equity_vol = static_black_scholes_parms.equity_vol;
  sq_equity_vol = equity_vol * equity_vol;
  forward = static_black_scholes_parms.forward;
  yr_to_exp = static_black_scholes_parms.yr_to_exp;
  pay_rec = static_black_scholes_parms.pay_rec;
  equity_strike = static_black_scholes_parms.equity_strike;

  if (pay_rec == SRT_RECEIVER)
    (*payoff) =
        max(0, equity_strike - forward * exp(equity_vol * dwt -
                                             0.5 * sq_equity_vol * yr_to_exp));
  else if (pay_rec == SRT_RECEIVER)
    (*payoff) = max(
        0, forward * exp(equity_vol * dwt - 0.5 * sq_equity_vol * yr_to_exp) -
               equity_strike);

  return NULL;
}

static Err srt_f_amer_option_payoff(Date eval_date, Date ExerDate, double dwt,
                                    double *payoff) {
  double sq_equity_vol, equity_vol, forward, equity_strike, time, Eps;
  SrtReceiverType pay_rec;

  /* GET STATIC INFORMATIONS */
  equity_vol = static_black_scholes_parms.equity_vol;
  sq_equity_vol = equity_vol * equity_vol;
  forward = static_black_scholes_parms.forward;
  time = (double)(ExerDate - eval_date) * YEARS_IN_DAY;
  pay_rec = static_black_scholes_parms.pay_rec;
  equity_strike = static_black_scholes_parms.equity_strike;

  if (pay_rec == SRT_RECEIVER)
    Eps = -1.0;
  else
    Eps = +1.0;

  (*payoff) = max(
      0, Eps * (forward * exp(equity_vol * dwt - 0.5 * sq_equity_vol * time) -
                equity_strike));

  return NULL;
}

#undef ONE_HALF
#undef PDE_THETA
#undef TREE_MESH_SPACING
#undef TREE_MESH_SPACING_SQUARE
#undef LIM_IN_STDEV
#undef MAX_TIME
#undef MIN_NODE
#undef MAX_NODE
